####################################################################################################
# - File for the MCMC sampler of the Accelerated Failure Time model where there is a fixed intercept
####################################################################################################

library(CASdatasets)
library(survival)
library(segmented)
library(mvtnorm)
library(invgamma)
library(truncnorm)
library(LaplacesDemon)
library(dplyr)
library(MASS)
#library(TruncatedDistributions)


#===================== - Final dataset for inference - =====================
R1 <- 6
R2 <- 4
R3 <- 3
# - Covariates storage

X1_cov <- matrix(NA, sample_size, R1)
X2_cov <- matrix(NA, sample_size, R2)
X3_cov <- matrix(NA, sample_size, R3)

# - Surrending
X1_cov[,1] <- uslapseagent$annual.premium
X1_cov[,2] <- uslapseagent$termination.cause
X1_cov[,3] <- ifelse(uslapseagent$acc.death.rider=="Rider", 1, 0)
X1_cov[,4] <- ifelse(uslapseagent$gender=="Female", 1, 0)
X1_cov[,5] <- ifelse(uslapseagent$premium.frequency %in% c("Annual", 'Other'), 1, 0)
X1_cov[,6] <- ifelse(uslapseagent$underwriting.age %in% c('Young',"Middle"), 1, 0)

# - Death
X2_cov[,1] <- ifelse(uslapseagent$gender=="Female", 1, 0)
X2_cov[,2] <- ifelse(uslapseagent$underwriting.age=="Old", 1, 0)
X2_cov[,3] <- ifelse(uslapseagent$living.place=="Other", 1, 0)
X2_cov[,4] <- ifelse(uslapseagent$risk.state=="Smoker", 1, 0)

# - Other cause
X3_cov[,1] <- uslapseagent$annual.premium
X3_cov[,2] <- ifelse(uslapseagent$premium.frequency %in% c("Annual", 'Other'), 1, 0)
X3_cov[,3] <- ifelse(uslapseagent$acc.death.rider=="Rider", 1, 0)

t_i <- uslapseagent$duration
d_ci <- uslapseagent[,c("surrender", "death", "other")]

#================================ - Parameters storage - =============================================

n_iter <- 50000 #number of iterations, set by the researcher

K = 25 # Dunson (2010) claims a number of pieces between 20 and 50
M <- 3
R <- 11


# - beta regression parameters
beta1 <- matrix(0, n_iter, R1)
colnames(beta1) <- c(sprintf("beta_1%d", c(1:R1)))
beta2 <- matrix(0, n_iter, R2)
colnames(beta2) <- c(sprintf("beta_2%d", c(1:R2)))
beta3 <- matrix(0, n_iter, R3)
colnames(beta3) <- c(sprintf("beta_3%d", c(1:R3)))

sigma2_c <- matrix(3, n_iter, M)

# - theta frailty parameters
theta_star <- matrix(0, n_iter, M*K)
colnames(theta_star) <- rep(c(sprintf("theta_s_%d", c(1:M))), K)

# - pi mixture weight
pi_mw <- matrix(NA, n_iter, K)
colnames(pi_mw) <- c(sprintf("pi_%d", c(1:K)))

# - mu_theta (parameter of the base distribution)
mu_theta <- matrix(NA, n_iter, M)
colnames(mu_theta) <- c(sprintf("mu_theta_%d", c(1:M)))

# - Covariance_theta
Sigma_theta <- matrix(NA, n_iter, (M^2))

# - psi from stick-breaking procedure
psi_sbp <- matrix(NA, n_iter, K)
colnames(psi_sbp) <- c(sprintf("psi_%d", c(1:K)))

# phi concentration parameter of the Dirichlet Process
phi <- matrix(NA, n_iter, 1)

# - Number of classes
N_k_tb <- matrix(NA, n_iter, K)
N_c_tb <- rep(0, n_iter)
lpd <- rep(NA, n_iter)

############################### - Set prior distributions - #####################################

M=3
K=25

# - alpha
omega1 <- 1
omega2 <- 1

# - beta
mu_beta1_prior <- rep(0, R1)
mu_beta2_prior <- rep(0, R2)
mu_beta3_prior <- rep(0, R3)

Sigma_beta1_prior <- rep(9, R1) #rep(1, 2) # assuming independently distributed beta_cl
Sigma_beta2_prior <- rep(9, R2) #rep(1, 2)
Sigma_beta3_prior <- rep(9, R3) #rep(1, 2)

# - theta
## - mu_theta
lambda1 <- rep(0, M) # - indicated in the paper as m
lambda2 <-  9 * matrix(diag(1, M), M, M)# - 0.0001 * matrix(diag(1, M), M, M) # - indicated in the paper as B

## - Covariance theta
lambda3 <- M + 5
lambda4 <- 0.001 * diag(1, M) %*% matrix(c(1,0.5,0.2,0.5,1,0.3,0.2,0.3,1), M, M, byrow=TRUE) %*% diag(1, M)

# - phi
lambda5 <- 1
lambda6 <- 1


####################################### - Starting values - ################################################

# beta
beta1[1,] <- rnorm(R1, 0, 1)# c(0, 0.14, 0.8, -0.2, -0.07, -0.3, 0.01) #st_val
beta2[1,] <- rnorm(R2, 0, 1)# c(0, -0.0784, 0.2, 0.121, 0.2) #st_val
beta3[1,] <- rnorm(R3, 0, 1)# c(0,-0.0784, 0.2, 0.121) #st_val


# - theta frailty parameters
theta_star[1,] <- rnorm(M*K, 0, 0.5)# c(0.2, 0.3, 0.1) #st_val

# - mu_theta (parameter of the base distribution)
mu_theta[1,] <- rnorm(M, 0, 0.5)# c(0.2, 0.3, 0.1) #st_val

# - Covariance_theta
Sigma_theta[1,] <- c(1,0.5,0.2,0.5,1,0.3,0.2,0.3,1)

# - psi from stick-breaking procedure
# - Initial value for the stick-breaking weights V (random)
phi[1,1] <- 1

psi_sbp[1,K] <- 1
for(k in 1:(K-1)){
  psi_sbp[1,k] <- rbeta(1, shape1 = 1, shape2 = phi[1,1], ncp = 0)
}

# - pi mixture weight
pi_mw[1, 1] <- psi_sbp[1]
for(k in 2:K){
  pi_mw[1,k] <- psi_sbp[1,k] * prod(1-psi_sbp[1,1:(k-1)])
}

######################################## - MCMC sampling - ###################################
#  tryCatch({

# - Preparatory
N_k <- rep(0, K)
prob_alloc_matrix <- matrix(NA, 1, K) # storage of the posterior weights from which to sample K

S_unit <- rep(0, sample_size)
S_unit <- sample(c(1:K), sample_size, prob=rep(1/K, K), replace=TRUE)

y_i <- matrix(NA, sample_size, M)
y_i[,1] <- log(t_i)
y_i[,2] <- log(t_i)
y_i[,3] <- log(t_i)


S_unit_fc <- function(unit){
  
  # - use log sum exp trick to deal with very small numbers which R inappropriately set to zero
  pf_vector <- rep(0, K)
  th_k_pos <- c(1:M) # position vector for theta
  
  piece1 <- sum(beta1[i-1,] * X1_cov[unit,])
  piece2 <- sum(beta2[i-1,] * X2_cov[unit,]) # t_i[unit] * exp(alpha2[i-1,1] + sum(beta2[i-1,] * X2_cov[unit,]))
  piece3 <- sum(beta3[i-1,] * X3_cov[unit,])
  
  pf_vector <- log(pi_mw[i-1,]) + dnorm(y_i[unit,1], mean = (piece1 + theta_s1_lp), sd = sqrt(sigma2_c[i-1,1]), log=TRUE) + 
                                  dnorm(y_i[unit,2], mean = (piece2 + theta_s2_lp), sd = sqrt(sigma2_c[i-1,2]), log=TRUE) + 
                                  dnorm(y_i[unit,3], mean = (piece3 + theta_s3_lp), sd = sqrt(sigma2_c[i-1,3]), log=TRUE) 
  
  max_pf <- max(pf_vector)
  
  denom_calc_lse <- max_pf + log(sum(exp(pf_vector-max_pf)))
  
  prob_alloc <- exp(pf_vector - denom_calc_lse) # - in vectorized form
  
  # - Sample cluster allocation
  S_unit_output <- sample(c(1:K), 1, prob=prob_alloc)
  
  return(S_unit_output)
}

pnorm1_fc <- function(unit){
  
  y_output <- ifelse(d_ci[unit,1]==1, log(t_i[unit]), rtruncnorm(1, a=log(t_i[unit]), b=6, mean = (sum(beta1[i-1,] * X1_cov[unit,]) + theta_star[i-1,(1 + M * (S_unit[unit] - 1))]), sd = sqrt(sigma2_c[i-1,1])))
  
  return(y_output)
}

pnorm2_fc <- function(unit){

  y_output <- ifelse(d_ci[unit,2]==1, log(t_i[unit]), rtruncnorm(1, a=log(t_i[unit]), b=6, mean = (sum(beta2[i-1,] * X2_cov[unit,]) + theta_star[i-1,(2 + M * (S_unit[unit] - 1))]), sd = sqrt(sigma2_c[i-1,2])))
  
  return(y_output)
}

pnorm3_fc <- function(unit){

  y_output <- ifelse(d_ci[unit,3]==1, log(t_i[unit]), rtruncnorm(1, a=log(t_i[unit]), b=6, mean = (sum(beta3[i-1,] * X3_cov[unit,]) + theta_star[i-1,(3 + M * (S_unit[unit] - 1))]), sd = sqrt(sigma2_c[i-1,3])))
  
  return(y_output)
}


i <- 2
for(i in i:n_iter){
  #  start_time <- Sys.time()
  
  # - Step 1: Sample new mixture component
  
  theta_s1_lp <- theta_star[i-1,seq(1, 73, 3)]
  theta_s2_lp <- theta_star[i-1,seq(2, 74, 3)]
  theta_s3_lp <- theta_star[i-1,seq(3, 75, 3)]
  
  # - Step 0: Data Augment y_i
  y_i[,1] <- sapply(c(1:sample_size), FUN = pnorm1_fc)
  y_i[,2] <- sapply(c(1:sample_size), FUN = pnorm2_fc)
  y_i[,3] <- sapply(c(1:sample_size), FUN = pnorm3_fc)
  
  S_unit <- sapply(c(1:sample_size), FUN = S_unit_fc)
  
  #    profvis({
  
  # - Step 2.1: Sample stick-breaking weights
  for(k in 1:(K-1)){
    psi_sbp[i,k] <- rbeta(1, shape1 = 1 + sum(ifelse(S_unit==k,1,0)), shape2 = phi[i-1,1] + sum(ifelse(S_unit > k,1,0)), ncp = 0)
  }
  psi_sbp[,K] <- 1
  # - Step 2.2: Update mixture distribution pi
  ## - Computation of pi
  pi_mw[i, 1] <- psi_sbp[i,1]
  for(k in 2:K){
    pi_mw[i,k] <- psi_sbp[i,k] * prod(1-psi_sbp[i,1:(k-1)])
  }
  
  # - Step 3: Sample Beta_c
  Sigma_beta1_post <- solve(t(X1_cov / sigma2_c[i-1,1]) %*% X1_cov + diag(1/Sigma_beta1_prior))
  
  Sigma_beta2_post <- solve(t(X2_cov / sigma2_c[i-1,2]) %*% X2_cov + diag(1/Sigma_beta2_prior))
  
  Sigma_beta3_post <- solve(t(X3_cov / sigma2_c[i-1,3]) %*% X3_cov + diag(1/Sigma_beta3_prior))
  
  mu_beta1_post <- Sigma_beta1_post %*% (t(X1_cov / sigma2_c[i-1,1]) %*% (y_i[,1] - theta_star[i-1,(1 + M * (S_unit - 1))]))
  
  mu_beta2_post <- Sigma_beta2_post %*% (t(X2_cov / sigma2_c[i-1,2]) %*% (y_i[,2] - theta_star[i-1,(2 + M * (S_unit - 1))]))
  
  mu_beta3_post <- Sigma_beta3_post %*% (t(X3_cov / sigma2_c[i-1,3]) %*% (y_i[,3] - theta_star[i-1,(3 + M * (S_unit - 1))]))
  
  beta1[i,] <- mvrnorm(n = 1, mu_beta1_post, Sigma_beta1_post)
  beta2[i,] <- mvrnorm(n = 1, mu_beta2_post, Sigma_beta2_post)
  beta3[i,] <- mvrnorm(n = 1, mu_beta3_post, Sigma_beta3_post)
  
  # - Step 4: Sample theta*
  th_k_pos <- c(1:M)
  
  solve_Sigma_theta <- solve(matrix(Sigma_theta[i-1,], M,M))
  for(k in 1:K){
    N_k[k] <- sum(ifelse(S_unit==k,1,0)) # number of people with the same k - or within the same cluster
    
    list1 <- ifelse(S_unit==k,1,0) * (y_i[,1] - (X1_cov %*% beta1[i,])) 
    list2 <- ifelse(S_unit==k,1,0) * (y_i[,2] - (X2_cov %*% beta2[i,])) 
    list3 <- ifelse(S_unit==k,1,0) * (y_i[,3] - (X3_cov %*% beta3[i,])) 
    
    diff_vector <- c(sum(list1), sum(list2), sum(list3))
    
    var_k <- solve(solve_Sigma_theta + N_k[k] * diag(1 / sigma2_c[i-1,]))
    
    #mean_k <- var_k %*% (solve_Sigma_theta %*% mu_theta[i-1,] + diag(1 / sigma2_c[i-1,]) %*% diff_vector)
    mean_k <- var_k %*% (solve_Sigma_theta %*% mu_theta[i-1,] + diff_vector / sigma2_c[i-1,])
    
    theta_star[i,th_k_pos] <- mvrnorm(n = 1, mean_k, var_k)
    
    th_k_pos <- th_k_pos + M
    
  }
  N_k_tb[i,] <- N_k
  
  # - Step 4 and 5: sample mu_theta and Sigma_theta
  
  N_c <- sum(ifelse(N_k[1:K]>0,1,0)) # Number of clusters with at least one element
  
  # - Posterior variance of mu_0
  Lambda2 <- solve(solve(lambda2) + N_c * solve(matrix(Sigma_theta[i-1,], M,M)))
  
  # - Posterior mean of mu_0 and preparatory calculation for the posterior parameter of the Wishart matrix
  mu_theta_sum_over_k <- rep(0, M)
  th_k_pos <- c(1:M)
  
  for(k_count in 1:K){ ## - check if it performs matrix operations (unlikely)
    mu_theta_sum_over_k[1] <- mu_theta_sum_over_k[1] + ifelse(N_k[k_count] > 0, theta_star[i,th_k_pos[1]], 0)
    mu_theta_sum_over_k[2] <- mu_theta_sum_over_k[2] + ifelse(N_k[k_count] > 0, theta_star[i,th_k_pos[2]], 0)
    mu_theta_sum_over_k[3] <- mu_theta_sum_over_k[3] + ifelse(N_k[k_count] > 0, theta_star[i,th_k_pos[3]], 0)
    th_k_pos <- th_k_pos + M
  }
  
  Lambda1 <- Lambda2 %*% (solve(lambda2) %*% lambda1 + solve(matrix(Sigma_theta[i-1,], M,M)) %*% mu_theta_sum_over_k)
  mu_theta[i,] <- rmvnorm(1, mean = Lambda1, sigma = Lambda2) # - update of mu_0
  
  # - Update of Sigma_theta
  means_cross_prod <- matrix(0, M, M)
  th_k_pos <- c(1:M)
  
  for(k in 1:K){
    if(N_k[k] > 0){
      means_cross_prod <- means_cross_prod + (theta_star[i,th_k_pos] - mu_theta[i,]) %*% t(theta_star[i,th_k_pos] - mu_theta[i,])
    }
    th_k_pos <- th_k_pos + M
  }
  
  Lambda3 <- N_c + lambda3
  Lambda4 <- lambda3 * lambda4 + means_cross_prod
  
  Sigma_theta[i,] <- matrix(rinvwishart(nu = Lambda3, S = Lambda4), nrow=1, ncol=ncol(Sigma_theta), byrow = TRUE)
  
  # - Step 6: Update sigma_c from a conjugate Inverse-Gamma distribution
  suff_stat1 <- 0.5 * sum((y_i[,1] - beta1[i,1] * X1_cov[,1] - beta1[i,2] * X1_cov[,2] - beta1[i,3] * X1_cov[,3] - beta1[i,4] * X1_cov[,4] - beta1[i,5] * X1_cov[,5] - beta1[i,6] * X1_cov[,6] - theta_star[i,(1 + M * (S_unit - 1))])^2)
  suff_stat2 <- 0.5 * sum((y_i[,2] - beta2[i,1] * X2_cov[,1] - beta2[i,2] * X2_cov[,2] - beta2[i,3] * X2_cov[,3] - beta2[i,4] * X2_cov[,4]                                                     - theta_star[i,(2 + M * (S_unit - 1))])^2)
  suff_stat3 <- 0.5 * sum((y_i[,3] - beta3[i,1] * X3_cov[,1] - beta3[i,2] * X3_cov[,2] - beta3[i,3] * X3_cov[,3]                                                                               - theta_star[i,(3 + M * (S_unit - 1))])^2)
  
  sigma2_c[i,1] <- 1 / rgamma(1, shape = (1 + 0.5 * sample_size), rate = (1 + suff_stat1) )
  sigma2_c[i,2] <- 1 / rgamma(1, shape = (1 + 0.5 * sample_size), rate = (1 + suff_stat2) )
  sigma2_c[i,3] <- 1 / rgamma(1, shape = (1 + 0.5 * sample_size), rate = (1 + suff_stat3) )
  
  # - Step 9: Update concentration parameter phi
  N_c <- sum(ifelse(N_k>0,1,0)) # Number of clusters with at least one element
  N_c_tb[i] <- N_c
  
  zeta <- rbeta(1, shape1 = phi[i-1,1] + 1, shape2 = sample_size, ncp=0)
  pi_zeta <- (lambda5 + N_c - 1) / (lambda5 + N_c - 1 + sample_size * (lambda6 - log(zeta) ))
  
  mix_subs <- sample(c(1, 0), 1, prob = c(pi_zeta, 1 - pi_zeta))
  phi[i,1] <- rgamma(1, shape = (lambda5 + N_c - ifelse(mix_subs==0,0,1)), rate = lambda6 - log(zeta) )
  
  
  #     })
  #  end_time <- Sys.time()
  #  end_time - start_time
  print(paste(i, "% iteration"))
  
  if(i %in% seq(5000, n_iter, 5000)){
    save.image("~/Desktop/TUM_Research/AFT_MCMC_sv1_hp1.RData")
  }
  
}

plot(c(1:(i-1)), beta1[1:(i-1),1], type="l")
plot(c(1:(i-1)), beta1[1:(i-1),2], type="l")
plot(c(1:(i-1)), beta1[1:(i-1),3], type="l")
plot(c(1:(i-1)), beta1[1:(i-1),4], type="l")
plot(c(1:(i-1)), beta1[1:(i-1),5], type="l")
plot(c(1:(i-1)), beta1[1:(i-1),6], type="l")

plot(c(1:(i-1)), beta2[1:(i-1),1], type="l")
plot(c(1:(i-1)), beta2[1:(i-1),2], type="l")
plot(c(1:(i-1)), beta2[1:(i-1),3], type="l")
plot(c(1:(i-1)), beta2[1:(i-1),4], type="l")

plot(c(1:(i-1)), beta3[1:(i-1),1], type="l")
plot(c(1:(i-1)), beta3[1:(i-1),2], type="l")
plot(c(1:(i-1)), beta3[1:(i-1),3], type="l")

plot(c(1:(i-1)), mu_theta[1:(i-1),1], type="l")
plot(c(1:(i-1)), mu_theta[1:(i-1),2], type="l")
plot(c(1:(i-1)), mu_theta[1:(i-1),3], type="l")

plot(c(1:(i-1)), pi_mw[1:(i-1),1], type="l")
plot(c(1:(i-1)), theta_star[1:(i-1),1], type="l")
plot(c(1:(i-1)), theta_star[1:(i-1),2], type="l")
plot(c(1:(i-1)), theta_star[1:(i-1),3], type="l")

plot(c(1:(i-1)), pi_mw[1:(i-1),2], type="l")
plot(c(1:(i-1)), theta_star[1:(i-1),4], type="l")
plot(c(1:(i-1)), theta_star[1:(i-1),5], type="l")
plot(c(1:(i-1)), theta_star[1:(i-1),6], type="l")

plot(c(1:(i-1)), pi_mw[1:(i-1),3], type="l")
plot(c(1:(i-1)), theta_star[1:(i-1),7], type="l")
plot(c(1:(i-1)), theta_star[1:(i-1),8], type="l")
plot(c(1:(i-1)), theta_star[1:(i-1),9], type="l")

plot(c(1:(i-1)), pi_mw[1:(i-1),2], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),3], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),4], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),5], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),6], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),7], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),8], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),9], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),10], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),11], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),12], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),13], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),14], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),15], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),16], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),17], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),18], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),19], type="l")
plot(c(1:(i-1)), pi_mw[1:(i-1),20], type="l")

