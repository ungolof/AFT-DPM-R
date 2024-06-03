####################################################################################################
# - File for the AFT model where there is a fixed intercept and an univariate DP distributed theta
# Set m_theta=0, and alpha fixed intercept over mixture components. Advantage: more parameters updated
# in a single sweep
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
R1 <- 7
R2 <- 5
R3 <- 4
# - Covariates storage

X1_cov <- matrix(NA, sample_size, R1)
X2_cov <- matrix(NA, sample_size, R2)
X3_cov <- matrix(NA, sample_size, R3)

# - Surrending
X1_cov[,1] <- 1
X1_cov[,2] <- uslapseagent$annual.premium
X1_cov[,3] <- uslapseagent$termination.cause
X1_cov[,4] <- ifelse(uslapseagent$acc.death.rider=="Rider", 1, 0)
X1_cov[,5] <- ifelse(uslapseagent$gender=="Female", 1, 0)
X1_cov[,6] <- ifelse(uslapseagent$premium.frequency %in% c("Annual", 'Other'), 1, 0)
X1_cov[,7] <- ifelse(uslapseagent$underwriting.age %in% c('Young',"Middle"), 1, 0)

# - Death
X2_cov[,1] <- 1
X2_cov[,2] <- ifelse(uslapseagent$gender=="Female", 1, 0)
X2_cov[,3] <- ifelse(uslapseagent$underwriting.age=="Old", 1, 0)
X2_cov[,4] <- ifelse(uslapseagent$living.place=="Other", 1, 0)
X2_cov[,5] <- ifelse(uslapseagent$risk.state=="Smoker", 1, 0)

# - Other cause
X3_cov[,1] <- 1
X3_cov[,2] <- uslapseagent$annual.premium
X3_cov[,3] <- ifelse(uslapseagent$premium.frequency %in% c("Annual", 'Other'), 1, 0)
X3_cov[,4] <- ifelse(uslapseagent$acc.death.rider=="Rider", 1, 0)

t_i <- uslapseagent$duration
d_ci <- uslapseagent[,c("surrender", "death", "other")]

#================================ - Parameters storage - =============================================

n_iter <- 50000 #number of iterations, set by the researcher

K = 25 # Dunson (2010) claims a number of pieces between 20 and 50
M <- 3

# - beta regression parameters
beta1 <- matrix(0, n_iter, R1)
colnames(beta1) <- c(sprintf("beta_1%d", c(1:R1)))
beta2 <- matrix(0, n_iter, R2)
colnames(beta2) <- c(sprintf("beta_2%d", c(1:R2)))
beta3 <- matrix(0, n_iter, R3)
colnames(beta3) <- c(sprintf("beta_3%d", c(1:R3)))

sigma2_c <- matrix(3, n_iter, M)

# - theta frailty parameters
theta_star <- matrix(0, n_iter, K)
colnames(theta_star) <- sprintf("theta_s_%d", c(1:K))

# - gamma 2 and gamma 3
gamma <- matrix(NA, n_iter, ncol = M-1)

# - pi mixture weight
pi_mw <- matrix(NA, n_iter, K)
colnames(pi_mw) <- c(sprintf("pi_%d", c(1:K)))

# - Covariance_theta
sigma2_theta <- matrix(NA, n_iter, 1)

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

# - beta (including the intercept alpha)
mu_beta1_prior <- rep(0, R1)
mu_beta2_prior <- rep(0, R2)
mu_beta3_prior <- rep(0, R3)

Sigma_beta1_prior <- c(16, rep(9, R1-1)) #rep(1, 2) # assuming independently distributed beta_cl
Sigma_beta2_prior <- c(16, rep(9, R2-1)) #rep(1, 2)
Sigma_beta3_prior <- c(16, rep(9, R3-1)) #rep(1, 2)

# - theta
## - sigma^2 theta
lambda3 <- 1
lambda4 <- 1

# - gamma (set a normal prior with variance equal to 4)

# - phi
lambda5 <- 1
lambda6 <- 1


####################################### - Starting values - ################################################

# beta
beta1[1,] <- rnorm(R1, 0, 1)# c(0, 0.14, 0.8, -0.2, -0.07, -0.3, 0.01) #st_val
beta2[1,] <- rnorm(R2, 0, 1)# c(0, -0.0784, 0.2, 0.121, 0.2) #st_val
beta3[1,] <- rnorm(R3, 0, 1)# c(0,-0.0784, 0.2, 0.121) #st_val

# - gamma
gamma[1,] <- rnorm(M-1, 0, 1)

# - theta frailty parameters
theta_star[1,] <- rnorm(K, 0, 0.5)# c(0.2, 0.3, 0.1) #st_val

# - Covariance_theta
sigma2_theta[1,] <- 1

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
  
  #pf_vector <- sapply(c(1:K), FUN = pf_vector_fc)
  # - pf_vector in vectorized form
  #  pf_vector <- log(pi_mw[i-1,]) + d_ci[unit,1] * dnorm(y_i[unit,1], mean = (piece1 + theta_s1_lp), sd = sqrt(sigma2_c[i-1,1]), log=TRUE) + (1 - d_ci[unit,1]) * pnorm(y_i[unit], mean = (piece1 + theta_s1_lp), sd = sqrt(sigma2_c[i-1,1]), log=TRUE, lower.tail = FALSE) +
  #                                  d_ci[unit,2] * dnorm(y_i[unit,2], mean = (piece2 + theta_s2_lp), sd = sqrt(sigma2_c[i-1,2]), log=TRUE) + (1 - d_ci[unit,2]) * pnorm(y_i[unit], mean = (piece2 + theta_s2_lp), sd = sqrt(sigma2_c[i-1,2]), log=TRUE, lower.tail = FALSE) +
  #                                  d_ci[unit,3] * dnorm(y_i[unit,3], mean = (piece3 + theta_s3_lp), sd = sqrt(sigma2_c[i-1,3]), log=TRUE) + (1 - d_ci[unit,3]) * pnorm(y_i[unit], mean = (piece3 + theta_s3_lp), sd = sqrt(sigma2_c[i-1,3]), log=TRUE, lower.tail = FALSE)
  
  pf_vector <- log(pi_mw[i-1,]) + dnorm(y_i[unit,1], mean = (piece1 +                theta_lp), sd = sqrt(sigma2_c[i-1,1]), log=TRUE) + 
                                  dnorm(y_i[unit,2], mean = (piece2 + gamma[i-1,1] * theta_lp), sd = sqrt(sigma2_c[i-1,2]), log=TRUE) + 
                                  dnorm(y_i[unit,3], mean = (piece3 + gamma[i-1,2] * theta_lp), sd = sqrt(sigma2_c[i-1,3]), log=TRUE) 
  
  max_pf <- max(pf_vector)
  
  denom_calc_lse <- max_pf + log(sum(exp(pf_vector-max_pf)))
  
  prob_alloc <- exp(pf_vector - denom_calc_lse) # - in vectorized form
  
  # - Sample cluster allocation
  S_unit_output <- sample(c(1:K), 1, prob=prob_alloc)
  
  return(S_unit_output)
}

pnorm1_fc <- function(unit){
  
  #  prob1 <- pnorm(log(t_i[unit]), mean = (sum(beta1[i-1,] * X1_cov[unit,]) + theta_star[i-1,(1 + M * (S_unit[unit] - 1))]), sd = sqrt(sigma2_c[i-1,1]), lower.tail = TRUE)
  
  #  unif1 <- runif(1, prob1, 1)
  
  #  y_output <- ifelse(d_ci[unit,1]==1, log(t_i[unit]), qnorm(unif1, mean = (sum(beta1[i-1,] * X1_cov[unit,]) + theta_star[i-1,(1 + M * (S_unit[unit] - 1))]), sd = sqrt(sigma2_c[i-1,1]), lower.tail = TRUE))
  y_output <- ifelse(d_ci[unit,1]==1, log(t_i[unit]), rtruncnorm(1, a=log(t_i[unit]), b=6, mean = (sum(beta1[i-1,] * X1_cov[unit,]) + theta_star[i-1,S_unit[unit]]), sd = sqrt(sigma2_c[i-1,1])))
  
  return(y_output)
}

pnorm2_fc <- function(unit){

  y_output <- ifelse(d_ci[unit,2]==1, log(t_i[unit]), rtruncnorm(1, a=log(t_i[unit]), b=6, mean = (sum(beta2[i-1,] * X2_cov[unit,]) + gamma[i-1,1] * theta_star[i-1,S_unit[unit]]), sd = sqrt(sigma2_c[i-1,2])))
  
  return(y_output)
}

pnorm3_fc <- function(unit){

  y_output <- ifelse(d_ci[unit,3]==1, log(t_i[unit]), rtruncnorm(1, a=log(t_i[unit]), b=6, mean = (sum(beta3[i-1,] * X3_cov[unit,]) + gamma[i-1,2] * theta_star[i-1,S_unit[unit]]), sd = sqrt(sigma2_c[i-1,3])))
  
  return(y_output)
}


i <- 2
for(i in i:n_iter){
  #  start_time <- Sys.time()
  
  # - Step 1: Sample new mixture component
  theta_lp <- theta_star[i-1,]

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
  
  # - Step 3: Sample Beta_c and gamma

  Sigma_beta1_post <- solve(t(X1_cov / sigma2_c[i-1,1]) %*% X1_cov + diag(1/Sigma_beta1_prior))

  X2_cov_aug <- cbind(X2_cov, theta_star[i-1,S_unit])
  Sigma_beta2_post <- solve(t(X2_cov_aug / sigma2_c[i-1,2]) %*% X2_cov_aug + diag(1/c(Sigma_beta2_prior, 16)))
  
  X3_cov_aug <- cbind(X3_cov, theta_star[i-1,S_unit])
  Sigma_beta3_post <- solve(t(X3_cov_aug / sigma2_c[i-1,3]) %*% X3_cov_aug + diag(1/c(Sigma_beta3_prior, 16)))
  
  mu_beta1_post <- Sigma_beta1_post %*% (t(X1_cov / sigma2_c[i-1,1]) %*% (y_i[,1] - theta_star[i-1,S_unit]))
  
  mu_beta2_post <- Sigma_beta2_post %*% (t(X2_cov_aug / sigma2_c[i-1,2]) %*% y_i[,2])
  
  mu_beta3_post <- Sigma_beta3_post %*% (t(X3_cov_aug / sigma2_c[i-1,3]) %*% y_i[,3])
  
  beta1[i,] <- mvrnorm(n = 1, mu_beta1_post, Sigma_beta1_post)
  
  vec2 <- mvrnorm(n = 1, mu_beta2_post, Sigma_beta2_post)
  beta2[i,] <- vec2[1:R2]

  vec3 <- mvrnorm(n = 1, mu_beta3_post, Sigma_beta3_post)
  beta3[i,] <- mvrnorm(n = 1, mu_beta3_post, Sigma_beta3_post)[1:R3]
  
  gamma[i,1] <- vec2[R2+1]
  gamma[i,2] <- vec3[R3+1]
#  Sigma_alpha1_post <- 1 / (sample_size / sigma2_c[i-1,1] + 1 / 16)
#  Sigma_alpha2_post <- 1 / (sample_size / sigma2_c[i-1,2] + 1 / 16)
#  Sigma_alpha3_post <- 1 / (sample_size / sigma2_c[i-1,3] + 1 / 16)
  
#  mu_alpha1_post <- Sigma_alpha1_post * (sum(y_i[,1] - beta1[i-1,2:R1] %*% t(X1_cov[,2:R1]) - theta_star[i-1,S_unit]) / sigma2_c[i-1,1])
#  mu_alpha2_post <- Sigma_alpha2_post * (sum(y_i[,2] - beta2[i-1,2:R2] %*% t(X2_cov[,2:R2]) - gamma[i-1,1] * theta_star[i-1,S_unit]) / sigma2_c[i-1,2])
#  mu_alpha3_post <- Sigma_alpha3_post * (sum(y_i[,3] - beta3[i-1,2:R3] %*% t(X3_cov[,2:R3]) - gamma[i-1,2] * theta_star[i-1,S_unit]) / sigma2_c[i-1,3])
  
#  beta1[i,1] <- rnorm(n = 1, mu_alpha1_post, sd=sqrt(Sigma_alpha1_post))
#  beta2[i,1] <- rnorm(n = 1, mu_alpha2_post, sd=sqrt(Sigma_alpha2_post))
#  beta3[i,1] <- rnorm(n = 1, mu_alpha3_post, sd=sqrt(Sigma_alpha3_post))
  
#  Sigma_beta1_post <- solve(t(X1_cov[,2:R1] / sigma2_c[i-1,1]) %*% X1_cov[,2:R1] + diag(1/Sigma_beta1_prior[2:R1]))
  
#  Sigma_beta2_post <- solve(t(X2_cov[,2:R2] / sigma2_c[i-1,2]) %*% X2_cov[,2:R2] + diag(1/Sigma_beta2_prior[2:R2]))
  
#  Sigma_beta3_post <- solve(t(X3_cov[,2:R3] / sigma2_c[i-1,3]) %*% X3_cov[,2:R3] + diag(1/Sigma_beta3_prior[2:R3]))
  
#  mu_beta1_post <- Sigma_beta1_post %*% (t(X1_cov[,2:R1] / sigma2_c[i-1,1]) %*% (y_i[,1] - theta_star[i-1,S_unit]))
  
#  mu_beta2_post <- Sigma_beta2_post %*% (t(X2_cov[,2:R2] / sigma2_c[i-1,2]) %*% (y_i[,2] - gamma[i-1,1] * theta_star[i-1,S_unit]))
  
#  mu_beta3_post <- Sigma_beta3_post %*% (t(X3_cov[,2:R3] / sigma2_c[i-1,3]) %*% (y_i[,3] - gamma[i-1,2] * theta_star[i-1,S_unit]))
  
#  beta1[i,c(2:R1)] <- mvrnorm(n = 1, mu_beta1_post, Sigma_beta1_post)
  
#  beta2[i,c(2:R2)] <- mvrnorm(n = 1, mu_beta2_post, Sigma_beta2_post)
  
#  beta3[i,c(2:R3)] <- mvrnorm(n = 1, mu_beta3_post, Sigma_beta3_post)
  
  
#  Sigma_gamma2_post <- 1/(sum(theta_star[i-1,S_unit]^2) / sigma2_c[i-1,2] + 1/16 )
#  mu_gamma2_post <- Sigma_gamma2_post * (t(theta_star[i-1,S_unit] / sigma2_c[i-1,2]) %*% (y_i[,2] - X2_cov %*% beta2[i,]))

#  Sigma_gamma3_post <- 1/(sum(theta_star[i-1,S_unit]^2) / sigma2_c[i-1,3] + 1/16 )
#  mu_gamma3_post <- Sigma_gamma3_post * (t(theta_star[i-1,S_unit] / sigma2_c[i-1,3]) %*% (y_i[,3] - X3_cov %*% beta3[i,]))
  
#  gamma[i,1] <- rnorm(n = 1, mu_gamma2_post, sd=sqrt(Sigma_gamma2_post))
#  gamma[i,2] <- rnorm(n = 1, mu_gamma3_post, sd=sqrt(Sigma_gamma3_post))

  # - Step 4: Sample theta*

  for(k in 1:K){
    N_k[k] <- sum(ifelse(S_unit==k,1,0)) # number of people with the same k - or within the same cluster

    list1 <- ifelse(S_unit==k,1,0) * (y_i[,1] - (X1_cov %*% beta1[i,])) #  - beta1[i,1] * X1_cov[,1] - beta1[i,2] * X1_cov[,2] - beta1[i,3] * X1_cov[,3] - beta1[i,4] * X1_cov[,4] - beta1[i,5] * X1_cov[,5] - beta1[i,6] * X1_cov[,6] - beta1[i,7] * X1_cov[,7])
    list2 <- ifelse(S_unit==k,1,0) * (y_i[,2] - (X2_cov %*% beta2[i,])) #  - beta2[i,1] * X2_cov[,1] - beta2[i,2] * X2_cov[,2] - beta2[i,3] * X2_cov[,3] - beta2[i,4] * X2_cov[,4] - beta2[i,5] * X2_cov[,5])
    list3 <- ifelse(S_unit==k,1,0) * (y_i[,3] - (X3_cov %*% beta3[i,])) # beta3[i,1] * X3_cov[,1] - beta3[i,2] * X3_cov[,2] - beta3[i,3] * X3_cov[,3] - beta3[i,4] * X3_cov[,4])
    
    diff_vector <- c(sum(list1), sum(list2), sum(list3))
    
    var_k <- 1/(1 / sigma2_theta[i-1,1] + N_k[k] * (1 / sigma2_c[i-1,1] + (gamma[i,1]^2) / sigma2_c[i-1,2] + (gamma[i,2]^2) / sigma2_c[i-1,3]))
    
    mean_k <- var_k * (diff_vector[1] / sigma2_c[i-1,1] + gamma[i,1] * diff_vector[2] / sigma2_c[i-1,2] + gamma[i,2] * diff_vector[3] / sigma2_c[i-1,3])    
    
    theta_star[i,k] <- mvrnorm(n = 1, mean_k, var_k)


  }
  N_k_tb[i,] <- N_k
  
  # - Step 4 and 5: sample mu_theta and Sigma_theta
  
  N_c <- sum(ifelse(N_k[1:K]>0,1,0)) # Number of clusters with at least one element
  
  # - Posterior variance of mu_0
  
  # - Update of sigma2_theta
  shape_theta <- 2 + N_c
  rate_theta <- 2 + sum(ifelse(N_k>0,1,0) * (theta_star[i,]^2))
  sigma2_theta[i,1] <- 1 / rgamma(1, shape = shape_theta, rate = rate_theta)  
  
  # - Step 6: Update sigma_c from a conjugate Inverse-Gamma distribution
  suff_stat1 <- 0.5 * sum((y_i[,1] - beta1[i,1] * X1_cov[,1] - beta1[i,2] * X1_cov[,2] - beta1[i,3] * X1_cov[,3] - beta1[i,4] * X1_cov[,4] - beta1[i,5] * X1_cov[,5] - beta1[i,6] * X1_cov[,6] - beta1[i,7] * X1_cov[,7] -              theta_star[i,S_unit])^2)
  suff_stat2 <- 0.5 * sum((y_i[,2] - beta2[i,1] * X2_cov[,1] - beta2[i,2] * X2_cov[,2] - beta2[i,3] * X2_cov[,3] - beta2[i,4] * X2_cov[,4] - beta2[i,5] * X2_cov[,5]                                                     - gamma[i,1] * theta_star[i,S_unit])^2)
  suff_stat3 <- 0.5 * sum((y_i[,3] - beta3[i,1] * X3_cov[,1] - beta3[i,2] * X3_cov[,2] - beta3[i,3] * X3_cov[,3] - beta3[i,4] * X3_cov[,4]                                                                               - gamma[i,2] * theta_star[i,S_unit])^2)
  
  sigma2_c[i,1] <- 1 / rgamma(1, shape = (1 + 0.5 * sample_size), rate = (1 + suff_stat1) )
  sigma2_c[i,2] <- 1 / rgamma(1, shape = (1 + 0.5 * sample_size), rate = (1 + suff_stat2) )
  sigma2_c[i,3] <- 1 / rgamma(1, shape = (1 + 0.5 * sample_size), rate = (1 + suff_stat3) )
  
  # - Step 9: Update concentration parameter phi
  N_c <- sum(ifelse(N_k>0,1,0)) # Number of clusters with at least one element
  N_c_tb[i] <- N_c
  
  zeta <- rbeta(1, shape1 = phi[i-1,1] + 1, shape2 = sample_size, ncp=0)
  pi_zeta <- (lambda5 + N_c - 1) / (lambda5 + N_c - 1 + sample_size * (lambda6 - log(zeta) ))
  
  mix_subs <- sample(c(1, 0), 1, prob = c(pi_zeta, 1 - pi_zeta))
  phi[i,1] <- rgamma(1, shape = (lambda5 + N_c - ifelse(mix_subs==0,0,1)), scale = lambda6 - log(zeta) )
  
  
  #     })
  #  end_time <- Sys.time()
  #  end_time - start_time
  print(paste(i, "% iteration"))
  
  if(i %in% seq(5000, n_iter, 5000)){
    save.image("~/Desktop/LFCR/Revision/AFT_HW_DPM.RData")
  }
  
}

#==========================================================

## - 1) Discard first 40,000 iterations (burn-in) and thinning every 50
burn_in <- 25000
thinning <- 20
n_iter <- 30000
seq_bi_thn <- seq((burn_in+thinning), n_iter, thinning)

# - beta regression parameters
beta1_thn <- beta1[seq_bi_thn,]
beta2_thn <- beta2[seq_bi_thn,]
beta3_thn <- beta3[seq_bi_thn,]

# - theta frailty parameters
theta_star_thn <- theta_star[seq_bi_thn,]

# - gamma
gamma_thn <- gamma[seq_bi_thn,]

# - pi mixture weight
pi_mw_thn <- pi_mw[seq_bi_thn,]

# - Covariance_theta
sigma2_theta_thn <- sigma2_theta[seq_bi_thn,]

# - psi from stick-breaking procedure
psi_sbp_thn <- psi_sbp[seq_bi_thn,]

# phi concentration parameter of the Dirichlet Process
phi_thn <- phi[seq_bi_thn,]

sigma2_c_thn <- sigma2_c[seq_bi_thn,]




##### - Build WAIC

## - Speeding up computations

WAIC_pieces <- function(unit){ # returns first log(mean(f)) and then mean(log(f))
  
  M_theta <- length(seq_bi_thn)
  sum_over_m <- 0
  log_sum_over_m <- 0
  for(m in 1:M_theta){
    pf_vector_num <- log(pi_mw_thn[m,]) + d_ci[unit,1] * dlnorm(t_i[unit], meanlog = sum(beta1_thn[m,] * X1_cov[unit,]) +                  theta_star_thn[m,], sdlog = sqrt(sigma2_c_thn[m,1]), log = TRUE) + (1 - d_ci[unit,1]) * plnorm(t_i[unit], meanlog = sum(beta1_thn[m,] * X1_cov[unit,]) +                  theta_star_thn[m,], sdlog = sqrt(sigma2_c_thn[m,1]), lower.tail = FALSE, log.p = TRUE) + 
                                          d_ci[unit,2] * dlnorm(t_i[unit], meanlog = sum(beta2_thn[m,] * X2_cov[unit,]) + gamma_thn[m,1] * theta_star_thn[m,], sdlog = sqrt(sigma2_c_thn[m,2]), log = TRUE) + (1 - d_ci[unit,2]) * plnorm(t_i[unit], meanlog = sum(beta2_thn[m,] * X2_cov[unit,]) + gamma_thn[m,1] * theta_star_thn[m,], sdlog = sqrt(sigma2_c_thn[m,2]), lower.tail = FALSE, log.p = TRUE) + 
                                          d_ci[unit,3] * dlnorm(t_i[unit], meanlog = sum(beta3_thn[m,] * X3_cov[unit,]) + gamma_thn[m,2] * theta_star_thn[m,], sdlog = sqrt(sigma2_c_thn[m,3]), log = TRUE) + (1 - d_ci[unit,3]) * plnorm(t_i[unit], meanlog = sum(beta3_thn[m,] * X3_cov[unit,]) + gamma_thn[m,2] * theta_star_thn[m,], sdlog = sqrt(sigma2_c_thn[m,3]), lower.tail = FALSE, log.p = TRUE)  
    
    max_pf_n <- max(pf_vector_num)
    
    pf_vector_den <- log(pi_mw_thn[m,]) + d_ci[unit,2] * dlnorm(t_i[unit], meanlog = sum(beta2_thn[m,] * X2_cov[unit,]) + gamma_thn[m,1] * theta_star_thn[m,], sdlog = sqrt(sigma2_c_thn[m,2]), log = TRUE) + (1 - d_ci[unit,2]) * plnorm(t_i[unit], meanlog = sum(beta2_thn[m,] * X2_cov[unit,]) + gamma_thn[m,1] * theta_star_thn[m,], sdlog = sqrt(sigma2_c_thn[m,2]), lower.tail = FALSE, log.p = TRUE) + 
                                          d_ci[unit,3] * dlnorm(t_i[unit], meanlog = sum(beta3_thn[m,] * X3_cov[unit,]) + gamma_thn[m,2] * theta_star_thn[m,], sdlog = sqrt(sigma2_c_thn[m,3]), log = TRUE) + (1 - d_ci[unit,3]) * plnorm(t_i[unit], meanlog = sum(beta3_thn[m,] * X3_cov[unit,]) + gamma_thn[m,2] * theta_star_thn[m,], sdlog = sqrt(sigma2_c_thn[m,3]), lower.tail = FALSE, log.p = TRUE)  
    
    max_pf_d <- max(pf_vector_den)
    
    # - Average of log(f)
    log_f_lapse <- max_pf_n + log(sum(exp(pf_vector_num - max_pf_n))) - max_pf_d - log(sum(exp(pf_vector_den - max_pf_d)))
    log_sum_over_m <- log_sum_over_m + log_f_lapse
    
    # - Average of f
    f_lapse <- exp(log_f_lapse)
    sum_over_m <- sum_over_m + f_lapse
    
  }
  return(c(-log(M_theta) + log(sum_over_m), log_sum_over_m/M_theta))
  
}


WAIC_insample <- matrix(NA, sample_size, 2)

for(a in 1:sample_size){
  start_time <- Sys.time()
  WAIC_insample[a,] <- WAIC_pieces(a)
  end_time <- Sys.time()
  
  print(end_time - start_time)
  print(a)
  
}

p_WAIC_i <- 2 * (colSums(WAIC_insample)[1] - colSums(WAIC_insample)[2])
WAIC_i_AFT_HW_DPM <- -2 * colSums(WAIC_insample)[1] + 2 * p_WAIC_i


sample_size_test <- nrow(uslapseagent_test)
X1_cov <- matrix(NA, sample_size_test, R1)
X2_cov <- matrix(NA, sample_size_test, R2)
X3_cov <- matrix(NA, sample_size_test, R3)

# - Surrending
X1_cov[,1] <- 1
X1_cov[,2] <- uslapseagent_test$annual.premium
X1_cov[,3] <- uslapseagent_test$termination.cause
X1_cov[,4] <- ifelse(uslapseagent_test$acc.death.rider=="Rider", 1, 0)
X1_cov[,5] <- ifelse(uslapseagent_test$gender=="Female", 1, 0)
X1_cov[,6] <- ifelse(uslapseagent_test$premium.frequency %in% c("Annual", 'Other'), 1, 0)
X1_cov[,7] <- ifelse(uslapseagent_test$underwriting.age %in% c('Young',"Middle"), 1, 0)

# - Death
X2_cov[,1] <- 1
X2_cov[,2] <- ifelse(uslapseagent_test$gender=="Female", 1, 0)
X2_cov[,3] <- ifelse(uslapseagent_test$underwriting.age=="Old", 1, 0)
X2_cov[,4] <- ifelse(uslapseagent_test$living.place=="Other", 1, 0)
X2_cov[,5] <- ifelse(uslapseagent_test$risk.state=="Smoker", 1, 0)

# - Other cause
X3_cov[,1] <- 1
X3_cov[,2] <- uslapseagent_test$annual.premium
X3_cov[,3] <- ifelse(uslapseagent_test$premium.frequency %in% c("Annual", 'Other'), 1, 0)
X3_cov[,4] <- ifelse(uslapseagent_test$acc.death.rider=="Rider", 1, 0)


#=================================== - Prediction - =============================

# - Parameters for analysis - Take the posterior mean
beta1_fin <- colMeans(beta1_thn)
beta2_fin <- colMeans(beta2_thn)
beta3_fin <- colMeans(beta3_thn)
gamma_fin <- colMeans(gamma)
sigma2_1 <- mean(sigma2_c_thn[,1])
sigma2_2 <- mean(sigma2_c_thn[,2])
sigma2_3 <- mean(sigma2_c_thn[,3])
theta_star_fin <- colMeans(theta_star_thn)
pi_mw_m <- colMeans(pi_mw_thn)


# - Support function

integrand <- function(t){
    pi_mw_m[1]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[1] ), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[1]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[1] ), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[2]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[2] ), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[2]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[2] ), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[3]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[3] ), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[3]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[3] ), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[4]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[4]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[4]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[4]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[5]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[5]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[5]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[5]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[6]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[6]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[6]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[6]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[7]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[7]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[7]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[7]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[8]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[8]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[8]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[8]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[9]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[9]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[9]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[9]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[10] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[10]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[10]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[10]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[11] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[11]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[11]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[11]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[12] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[12]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[12]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[12]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[13] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[13]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[13]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[13]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[14] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[14]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[14]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[14]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[15] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[15]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[15]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[15]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[16] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[16]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[16]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[16]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[17] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[17]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[17]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[17]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[18] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[18]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[18]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[18]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[19] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[19]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[19]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[19]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[20] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[20]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[20]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[20]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[21] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[21]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[21]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[21]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[22] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[22]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[22]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[22]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[23] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[23]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[23]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[23]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[24] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[24]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[24]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[24]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[25] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[25]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + gamma_fin[1] * theta_star_fin[25]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + gamma_fin[2] * theta_star_fin[25]), sd=sqrt(sigma2_3), lower.tail = FALSE)
  
}




#========================= - Prediction survival rates - formula (10) - =======================

y_i <- log(t_i)

# - Determination of the risk set R_tq
n_quarters <- round(max(t_i), 0) + 1
R_tq <- rep(0, n_quarters)

R_tq[1] <- sample_size_test
for(q in 2:n_quarters){
  R_tq[q] <- sum(ifelse(t_i>q-1,1,0)) ## - changed
}

F_t_table_AFT_HW_DPM <- matrix(NA, sample_size_test, n_quarters) # - takes F(d_q; X)

rates_tb_qt_AFT_HW_DPM <- function(j){
  integrate(integrand, lower=-Inf, upper=log(j))$value
}

for(i in 1:sample_size_test){ # - Calculate the F in formula (9)
  F_t_table_AFT_HW_DPM[i,] <- sapply(1:n_quarters, rates_tb_qt_AFT_HW_DPM)
}

# \hat{r}_{t,i} - individually predicted rates for computation of eq. (10) in Milhaud-Dutang
r_pred_table_AFT_HW_DPM <- matrix(NA, sample_size_test, n_quarters)
r_pred_table_AFT_HW_DPM[,1] <- F_t_table_AFT_HW_DPM[,1]
r_pred_table_AFT_HW_DPM[,2:n_quarters] <- (F_t_table_AFT_HW_DPM[,2:n_quarters] - F_t_table_AFT_HW_DPM[,1:(n_quarters-1)]) / (1 - F_t_table_AFT_HW_DPM[,1:(n_quarters-1)])



