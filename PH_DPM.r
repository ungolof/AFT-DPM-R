#========================================================================================================
# - Proportional hazard model; Consider other implementations, such as the case where we data-augment T
# - This is the implemetation where we consider alpha and a mean (m_theta) equal to zero
#========================================================================================================


library(CASdatasets)
library(survival)
library(segmented)
library(mvtnorm)
library(invgamma)
library(truncnorm)
library(LaplacesDemon)
library(dplyr)
library(MASS)
library(readxl)

#==================================== - LFCR - ==================

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

# 2) - Implementation
sample_size <- nrow(uslapseagent)
n_iter <- 50000 #number of iterations, set by the researcher

K = 20 # Dunson (2010) claims a number of pieces between 20 and 50
M <- 3

# - alpha piecewise constant hazard function
alpha1 <- matrix(NA, n_iter, 1)
alpha2 <- matrix(NA, n_iter, 1)
alpha3 <- matrix(NA, n_iter, 1)

# - beta regression parameters
beta1 <- matrix(0, n_iter, R1)
colnames(beta1) <- c(sprintf("beta_1%d", c(1:R1)))
beta2 <- matrix(0, n_iter, R2)
colnames(beta2) <- c(sprintf("beta_2%d", c(1:R2)))
beta3 <- matrix(0, n_iter, R3)
colnames(beta3) <- c(sprintf("beta_3%d", c(1:R3)))

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

# 3) - Set prior distributions

## - alpha
omega11 <- 1
omega12 <- 1
omega21 <- 1
omega22 <- 1
omega31 <- 1
omega32 <- 1

## - gamma conjugate prior distribution for exp(beta)
omega1_beta1 <- rep(1, R1)
omega1_beta2 <- rep(1, R2)
omega1_beta3 <- rep(1, R3)

omega2_beta1 <- rep(1, R1)
omega2_beta2 <- rep(1, R2)
omega2_beta3 <- rep(1, R3)

# - beta
mu_beta1 <- rep(0, 4)
mu_beta2 <- rep(0, 4)
mu_beta3 <- rep(0, 2)

sigma_beta1 <- rep(1, 4) # assuming independently distributed beta_cl
sigma_beta2 <- rep(1, 4) 
sigma_beta3 <- rep(1, 2) 


## - mu_theta
lambda1 <- rep(0, M) # - indicated in the paper as m
lambda2 <- matrix(diag(0.5, M), M, M) # - indicated in the paper as B

mu_theta[1,] <- 0

N_k_tb <- matrix(NA, n_iter, K)

S_unit_fc <- function(unit){
  
  # - use log sum exp trick to deal with very small numbers which R inappropriately set to zero
  pf_vector <- rep(0, K)
  th_k_pos <- c(1:M) # position vector for theta
  
  piece1 <- t_i[unit] * exp(sum(beta1[i-1,] * X1_cov[unit,])) #sum(beta1[i-1,] * X_cov[unit,]))
  piece2 <- t_i[unit] * exp(sum(beta2[i-1,] * X2_cov[unit,])) #exp(sum(beta2[i-1,] * X_cov[unit,]))
  piece3 <- t_i[unit] * exp(sum(beta3[i-1,] * X3_cov[unit,])) #exp(sum(beta3[i-1,] * X_cov[unit,]))
  
  #pf_vector <- sapply(c(1:K), FUN = pf_vector_fc)
  # - pf_vector in vectorized form
  pf_vector <- log(pi_mw[i-1,]) - piece1 * exp(theta_s1_lp) - piece2 * exp(theta_s2_lp) - piece3 * exp(theta_s3_lp) + d_ci[unit,1] * theta_s1_lp + d_ci[unit,2] * theta_s2_lp + d_ci[unit,3] * theta_s3_lp
  
  max_pf <- max(pf_vector)
  
  denom_calc_lse <- max_pf + log(sum(exp(pf_vector-max_pf)))
  
  prob_alloc <- exp(pf_vector - denom_calc_lse) # - in vectorized form
  
  # - Sample cluster allocation
  S_unit_output <- sample(c(1:K), 1, prob=prob_alloc)
  
  return(S_unit_output)
}


i <- 2
for(i in i:n_iter){
  #  start_time <- Sys.time()
  
  #     profvis({
  
  # - Step 1: Sample new mixture component
  theta_s1_lp <- theta_star[i-1,seq(1, K*M-2, M)]
  theta_s2_lp <- theta_star[i-1,seq(2, K*M-1, M)]
  theta_s3_lp <- theta_star[i-1,seq(3, K*M  , M)]
  
  S_unit <- sapply(c(1:sample_size), FUN = S_unit_fc)
  
  # - Step 2.1: Sample stick-breaking weights
  for(k in 1:K){
    psi_sbp[i,k] <- rbeta(1, shape1 = 1 + sum(ifelse(S_unit==k,1,0)), shape2 = phi[i-1,1] + sum(ifelse(S_unit > k,1,0)), ncp = 0)
  }
  
  # - Step 2.2: Update mixture distribution pi
  ## - Computation of pi
  pi_mw[i, 1] <- psi_sbp[i,1]
  for(k in 2:K){
    pi_mw[i,k] <- psi_sbp[i,k] * prod(1-psi_sbp[i,1:(k-1)])
  }
  
  # - Step 3: Sample theta*
  th_k_pos <- c(1:M)
  log_sum_Si <- rep(0, K)
  
  for(k in 1:K){
    log_sum_Si <- 0
    N_k[k] <- sum(ifelse(S_unit==k,1,0)) # number of people with the same k - or within the same cluster
    
    # - Sample theta* using the Adaptive Metropolis Hastings of Vihola (2012)
#    r_i <- mvrnorm(n = 1, rep(0, M), diag(1, M))
#    theta_proposed <- matrix(theta_star[i-1,th_k_pos],M,1) + Sigma_ch_th[, ((k-1) * M + 1):(k * M)] %*% matrix(r_i, M,1)

    theta_proposed <- rtmvnorm(1, mean=theta_star[i-1,th_k_pos], sigma = Sigma_th[, ((k-1) * M + 1):(k * M)], lower=rep(-20, M), upper=rep(20,M))
      
    # - Get likelihood of the data conditional to the clustering variable S_i
    log_sum_Si <- log_sum_Si - dmvnorm(theta_proposed, mean = theta_star[i-1,th_k_pos], sigma = Sigma_th[, ((k-1) * M + 1):(k * M)], log = TRUE) +
                               dmvnorm(theta_proposed, mean = mu_theta[i-1,], sigma = matrix(Sigma_theta[i-1,], M,M), log = TRUE) + # this line corresponds to the prior
                               dmvnorm(theta_star[i-1,th_k_pos], mean = t(theta_proposed), sigma = Sigma_th[, ((k-1) * M + 1):(k * M)], log = TRUE) -
                               dmvnorm(theta_star[i-1,th_k_pos], mean = mu_theta[i-1,], sigma = matrix(Sigma_theta[i-1,], M,M), log = TRUE) # this line corresponds to the prior
    
    log_sum_Si <- log_sum_Si + sum(ifelse(S_unit==k, - t_i * (exp(beta1[i-1,1] * X1_cov[,1] + beta1[i-1,2] * X1_cov[,2] + beta1[i-1,3] * X1_cov[,3] + beta1[i-1,4] * X1_cov[,4] + beta1[i-1,5] * X1_cov[,5] + beta1[i-1,6] * X1_cov[,6]) * (exp(theta_proposed[1]) - exp(theta_star[i-1,th_k_pos[1]])) + 
                                                              exp(beta2[i-1,1] * X2_cov[,1] + beta2[i-1,2] * X2_cov[,2] + beta2[i-1,3] * X2_cov[,3] + beta2[i-1,4] * X2_cov[,4]                                                        ) * (exp(theta_proposed[2]) - exp(theta_star[i-1,th_k_pos[2]])) + 
                                                              exp(beta3[i-1,1] * X3_cov[,1] + beta3[i-1,2] * X3_cov[,2] + beta3[i-1,3] * X3_cov[,3]                                                                                    ) * (exp(theta_proposed[3]) - exp(theta_star[i-1,th_k_pos[3]]))) + d_ci[,1] * (theta_proposed[1] - theta_star[i-1,th_k_pos[1]]) + d_ci[,2] * (theta_proposed[2] - theta_star[i-1,th_k_pos[2]]) + d_ci[,3] * (theta_proposed[3] - theta_star[i-1,th_k_pos[3]]), 0))
    
    
    # Acceptance step
    accept_i <- min(1, exp(log_sum_Si))
    
    u_unif <- runif(1, 0, 1)
    for(j in 1:M){
      theta_star[i,th_k_pos[j]] <- ifelse(u_unif > accept_i, theta_star[i-1,th_k_pos[j]], theta_proposed[j])
    }
    
    # Adjsut Sigma_delta_chol
    Sigma_th[, ((k-1) * M + 1):(k * M)] <- Sigma_ch_th[, ((k-1) * M + 1):(k * M)] %*% (diag(1, M) + (i^(-gamma_Vihola)) * (accept_i - opt_accept) * (r_i %*% t(r_i)) / (sum(r_i^2))) %*% t(Sigma_ch_th[, ((k-1) * M + 1):(k * M)])
    Sigma_ch_th[, ((k-1) * M + 1):(k * M)] <- t(chol(Sigma_th[, ((k-1) * M + 1):(k * M)]))
    
    th_k_pos <- th_k_pos + M
  }
  N_k_tb[i,] <- N_k
  # - Step 4 and 5: sample mu_theta and Sigma_theta
  
  N_c <- sum(ifelse(N_k>0,1,0)) # Number of clusters with at least one element
  
  # - Posterior variance of mu_0
  Lambda2 <- solve(solve(lambda2) + N_c * solve(matrix(Sigma_theta[i-1,], M,M)))
  
  # - Posterior mean of mu_0 and preparatory calculation for the posterior parameter of the Wishart matrix
  mu_theta_sum_over_k <- rep(0, M)
  th_k_pos <- c(1:M)
  
  for(k_count in 1:K){ 
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
      #      means_cross_prod <- means_cross_prod + (theta_star[i,th_k_pos] - mu_theta[i,]) %*% t(theta_star[i,th_k_pos] - mu_theta[i,])
      means_cross_prod <- means_cross_prod + (theta_star[i,th_k_pos]) %*% t(theta_star[i,th_k_pos])
    }
    th_k_pos <- th_k_pos + M
  }
  
  Lambda3 <- N_c + lambda3
  Lambda4 <- lambda3 * lambda4 + means_cross_prod
  
  Sigma_theta[i,] <- matrix(rinvwishart(nu = Lambda3, S = Lambda4), nrow=1, ncol=ncol(Sigma_theta), byrow = TRUE)
  
  # - Step 7: Sample beta_c
  
  ## - Acceptance - Rejection sampling
  
  ## - Beta1
  ### - r=1
  # beta*
  beta_star1 <- rtruncnorm(1, a=-10, b=10, mean=beta1[i-1,1], sd=sigma_beta1_pr[1]) # rnorm(n = 1, beta1[i-1,1], sigma_beta1_pr[1])
  # Acceptance step
  log_ef1 <- dnorm(beta_star1, mean=mu_beta1[1], sd=sigma_beta1[1], log=TRUE) - dnorm(beta_star1, mean=beta1[i-1,1], sd=sigma_beta1_pr[1], log=TRUE) - dnorm(beta1[i-1,1], mean=mu_beta1[1], sd=sigma_beta1[1], log=TRUE) + dnorm(beta1[i-1,1], mean=beta_star1, sd=sigma_beta1_pr[1], log=TRUE)
  log_ef1 <- log_ef1 + sum( - t_i * exp(beta1[i-1,2] * X1_cov[,2] + beta1[i-1,3] * X1_cov[,3] + beta1[i-1,4] * X1_cov[,4] + beta1[i-1,5] * X1_cov[,5] + beta1[i-1,6] * X1_cov[,6] + theta_star[i,(1 + M * (S_unit - 1))]) * (exp(beta_star1 * X1_cov[,1]) - exp(beta1[i-1,1] * X1_cov[,1])) + d_ci[,1] * X1_cov[,1] * (beta_star1 - beta1[i-1,1]) )
  
  accept_beta1 <- min(1, exp(log_ef1))
  beta_unif1 <- runif(1, 0, 1)
  beta1[i,1] <- ifelse(beta_unif1 > accept_beta1, beta1[i-1,1], beta_star1)
  # - Adjust sigma_beta
  sigma_beta1_pr[1] <- sigma_beta1_pr[1] * (1 + (i^(-gamma_Vihola)) * (accept_beta1 - beta_tg_accept))^0.5
  
  ### - r=2
  # beta*
  beta_star1 <- rtruncnorm(1, a=-10, b=10, mean=beta1[i-1,2], sd=sigma_beta1_pr[2]) #rnorm(n = 1, beta1[i-1,2], sigma_beta1_pr[2])
  # Acceptance step
  log_ef1 <- dnorm(beta_star1, mean=mu_beta1[2], sd=sigma_beta1[2], log=TRUE) - dnorm(beta_star1, mean=beta1[i-1,2], sd=sigma_beta1_pr[2], log=TRUE) - dnorm(beta1[i-1,2], mean=mu_beta1[2], sd=sigma_beta1[2], log=TRUE) + dnorm(beta1[i-1,2], mean=beta_star1, sd=sigma_beta1_pr[2], log=TRUE)
  log_ef1 <- log_ef1 + sum( - t_i * exp(beta1[i,1] * X1_cov[,1] + beta1[i-1,3] * X1_cov[,3] + beta1[i-1,4] * X1_cov[,4] + beta1[i-1,5] * X1_cov[,5] + beta1[i-1,6] * X1_cov[,6] + theta_star[i,(1 + M * (S_unit - 1))]) * (exp(beta_star1 * X1_cov[,2]) - exp(beta1[i-1,2] * X1_cov[,2])) + d_ci[,1] * X1_cov[,2] * (beta_star1 - beta1[i-1,2]) )
  
  accept_beta1 <- min(1, exp(log_ef1))
  beta_unif1 <- runif(1, 0, 1)
  beta1[i,2] <- ifelse(beta_unif1 > accept_beta1, beta1[i-1,2], beta_star1)
  # - Adjust sigma_beta
  sigma_beta1_pr[2] <- sigma_beta1_pr[2] * (1 + (i^(-gamma_Vihola)) * (accept_beta1 - beta_tg_accept))^0.5
  
  
  ## - Beta3
  ### - r=1
  # beta*
  beta_star3 <- rtruncnorm(1, a=-10, b=10, mean=beta3[i-1,1], sd=sigma_beta3_pr[1]) #rnorm(n = 1, beta3[i-1,1], sigma_beta3_pr[1])
  # Acceptance step
  log_ef3 <- dnorm(beta_star3, mean=mu_beta3[1], sd=sigma_beta3[1], log=TRUE) - dnorm(beta_star3, mean=beta3[i-1,1], sd=sigma_beta3_pr[1], log=TRUE) - dnorm(beta3[i-1,1], mean=mu_beta3[1], sd=sigma_beta3[1], log=TRUE) + dnorm(beta3[i-1,1], mean=beta_star3, sd=sigma_beta3_pr[1], log=TRUE)
  log_ef3 <- log_ef3 + sum( - t_i * exp(beta3[i-1,2] * X3_cov[,2] + beta3[i-1,3] * X3_cov[,3] + theta_star[i,(3 + M * (S_unit - 1))]) * (exp(beta_star3 * X3_cov[,1]) - exp(beta3[i-1,1] * X3_cov[,1])) + d_ci[,3] * X3_cov[,1] * (beta_star3 - beta3[i-1,1]) )
  
  accept_beta3 <- min(1, exp(log_ef3))
  beta_unif3 <- runif(1, 0, 1)
  beta3[i,1] <- ifelse(beta_unif3 > accept_beta3, beta3[i-1,1], beta_star3)
  # - Adjust sigma_beta
  sigma_beta3_pr[1] <- sigma_beta3_pr[1] * (1 + (i^(-gamma_Vihola)) * (accept_beta3 - beta_tg_accept))^0.5
  
  
  ## - beta 1 (conjugate modified gamma)
  ### - r=3
  Omega1_beta1 <- omega1_beta1[3] + sum(d_ci[,1] * X1_cov[,3])
  lc1 <- sum(t_i * X1_cov[,3] * exp(beta1[i,1] * X1_cov[,1] + beta1[i,2] * X1_cov[,2] + beta1[i-1,4] * X1_cov[,4] + beta1[i-1,5] * X1_cov[,5] + beta1[i-1,6] * X1_cov[,6] + theta_star[i,(1 + M * (S_unit - 1))]))
  Omega2_beta1 <- omega2_beta1[3] + lc1
  beta1[i,3] <- log(rgamma(1, shape = Omega1_beta1, rate = Omega2_beta1))
  
  ### - r=4
  Omega1_beta1 <- omega1_beta1[4] + sum(d_ci[,1] * X1_cov[,4])
  lc1 <- sum(t_i * X1_cov[,4] * exp(beta1[i,1] * X1_cov[,1] + beta1[i,2] * X1_cov[,2] + beta1[i,3] * X1_cov[,3] + beta1[i-1,5] * X1_cov[,5] + beta1[i-1,6] * X1_cov[,6] + theta_star[i,(1 + M * (S_unit - 1))]))
  Omega2_beta1 <- omega2_beta1[4] + lc1
  beta1[i,4] <- log(rgamma(1, shape = Omega1_beta1, rate = Omega2_beta1))
  
  ### - r=5
  Omega1_beta1 <- omega1_beta1[5] + sum(d_ci[,1] * X1_cov[,5])
  lc1 <- sum(t_i * X1_cov[,5] * exp(beta1[i,1] * X1_cov[,1] + beta1[i,2] * X1_cov[,2] + beta1[i,3] * X1_cov[,3] + beta1[i,4] * X1_cov[,4] + beta1[i-1,6] * X1_cov[,6] + theta_star[i,(1 + M * (S_unit - 1))]))
  Omega2_beta1 <- omega2_beta1[5] + lc1
  beta1[i,5] <- log(rgamma(1, shape = Omega1_beta1, rate = Omega2_beta1))
  
  ### - r=6
  Omega1_beta1 <- omega1_beta1[6] + sum(d_ci[,1] * X1_cov[,6])
  lc1 <- sum(t_i * X1_cov[,6] * exp(beta1[i,1] * X1_cov[,1] + beta1[i,2] * X1_cov[,2] + beta1[i,3] * X1_cov[,3] + beta1[i,4] * X1_cov[,4] + beta1[i,5] * X1_cov[,5] + theta_star[i,(1 + M * (S_unit - 1))]))
  Omega2_beta1 <- omega2_beta1[6] + lc1
  beta1[i,6] <- log(rgamma(1, shape = Omega1_beta1, rate = Omega2_beta1))
  
  
  ## - beta 2 (conjugate modified gamma)
  ### - r=1
  Omega1_beta2 <- omega1_beta2[1] + sum(d_ci[,2] * X2_cov[,1])
  lc2 <- sum(t_i * X2_cov[,1] * exp(beta2[i-1,2] * X2_cov[,2] + beta2[i-1,3] * X2_cov[,3] + beta2[i-1,4] * X2_cov[,4] + theta_star[i,(2 + M * (S_unit - 1))]))
  Omega2_beta2 <- omega2_beta2[1] + lc2
  beta2[i,1] <- log(rgamma(1, shape = Omega1_beta2, rate = Omega2_beta2))
  
  ### - r=2
  Omega1_beta2 <- omega1_beta2[2] + sum(d_ci[,2] * X2_cov[,2])
  lc2 <- sum(t_i * X2_cov[,2] * exp(beta2[i,1] * X2_cov[,1] + beta2[i-1,3] * X2_cov[,3] + beta2[i-1,4] * X2_cov[,4] + theta_star[i,(2 + M * (S_unit - 1))]))
  Omega2_beta2 <- omega2_beta2[2] + lc2
  beta2[i,2] <- log(rgamma(1, shape = Omega1_beta2, rate = Omega2_beta2))
  
  ### - r=3
  Omega1_beta2 <- omega1_beta2[3] + sum(d_ci[,2] * X2_cov[,3])
  lc2 <- sum(t_i * X2_cov[,3] * exp(beta2[i,1] * X2_cov[,1] + beta2[i,2] * X2_cov[,2] + beta2[i-1,4] * X2_cov[,4] + theta_star[i,(2 + M * (S_unit - 1))]))
  Omega2_beta2 <- omega2_beta2[3] + lc2
  beta2[i,3] <- log(rgamma(1, shape = Omega1_beta2, rate = Omega2_beta2))
  
  ### - r=4
  Omega1_beta2 <- omega1_beta2[4] + sum(d_ci[,2] * X2_cov[,4])
  lc2 <- sum(t_i * X2_cov[,4] * exp(beta2[i,1] * X2_cov[,1] + beta2[i,2] * X2_cov[,2] + beta2[i,3] * X2_cov[,3] + theta_star[i,(2 + M * (S_unit - 1))]))
  Omega2_beta2 <- omega2_beta2[4] + lc2
  beta2[i,4] <- log(rgamma(1, shape = Omega1_beta2, rate = Omega2_beta2))
  
  ## - beta3
  ### - r=2
  Omega1_beta3 <- omega1_beta3[2] + sum(d_ci[,3] * X3_cov[,2])
  lc2 <- sum(t_i * X3_cov[,2] * exp(beta3[i,1] * X3_cov[,1] + beta3[i-1,3] * X3_cov[,3] + theta_star[i,(3 + M * (S_unit - 1))]))
  Omega2_beta3 <- omega2_beta3[2] + lc2
  beta3[i,2] <- log(rgamma(1, shape = Omega1_beta3, rate = Omega2_beta3))
  
  ### - r=3
  Omega1_beta3 <- omega1_beta3[3] + sum(d_ci[,3] * X3_cov[,3])
  lc2 <- sum(t_i * X3_cov[,3] * exp(beta3[i,1] * X3_cov[,1] + beta3[i,2] * X3_cov[,2] + theta_star[i,(3 + M * (S_unit - 1))]))
  Omega2_beta3 <- omega2_beta3[3] + lc2
  beta3[i,3] <- log(rgamma(1, shape = Omega1_beta3, rate = Omega2_beta3))
  
  
  # - Step 9: Update concentration parameter phi
  zeta <- rbeta(1, shape1 = phi[i-1,1] + 1, shape2 = sample_size, ncp=0)
  pi_zeta <- (lambda5 + N_c - 1) / (lambda5 + N_c - 1 + sample_size * (lambda6 - log(zeta) ))
  
  mix_subs <- sample(c(1, 0), 1, prob = c(pi_zeta, 1 - pi_zeta))
  phi[i,1] <- rgamma(1, shape = (lambda5 + N_c - ifelse(mix_subs==0,0,1)), rate = lambda6 - log(zeta) )
  
  print(paste(i, "% iteration"))
  
  if(i %in% seq(5000, n_iter, 5000)){
    save.image("~/Desktop/LFCR/Revision/PH_DPM_noalpha.RData")
  }
  
  #     })
  #  end_time <- Sys.time()
  #  end_time - start_time
}

#=====================================

# - Convergence check
plot(c(1:(i-1)), mu_theta[1:(i-1),1], type="l")
plot(c(1:(i-1)), mu_theta[1:(i-1),2], type="l")
plot(c(1:(i-1)), mu_theta[1:(i-1),3], type="l")

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

plot(c(1:(i-1)), theta_star[1:(i-1),1], type="l")
plot(c(1:(i-1)), theta_star[1:(i-1),2], type="l")
plot(c(1:(i-1)), theta_star[1:(i-1),3], type="l")
plot(c(1:(i-1)), theta_star[1:(i-1),4], type="l")
plot(c(1:(i-1)), theta_star[1:(i-1),5], type="l")
plot(c(1:(i-1)), theta_star[1:(i-1),6], type="l")


## - 1) Discard first 40,000 iterations (burn-in) and thinning every 50
burn_in <- 40000
thinning <- 50
n_iter <- 50000
seq_bi_thn <- seq((burn_in+thinning), n_iter, thinning)

# - beta regression parameters
beta1_thn <- beta1[seq_bi_thn,]
beta2_thn <- beta2[seq_bi_thn,]
beta3_thn <- beta3[seq_bi_thn,]

# - theta frailty parameters
theta_star_thn <- theta_star[seq_bi_thn,]

# - pi mixture weight
pi_mw_thn <- pi_mw[seq_bi_thn,]

# - mu_theta (parameter of the base distribution)
mu_theta_thn <- mu_theta[seq_bi_thn,]


##### - Build WAIC

# - computing the average over the posterior sample

## - Speeding up computations

WAIC_pieces <- function(unit){ # returns first log(mean(f)) and then mean(log(f))
  
  M_theta <- length(seq_bi_thn)
  sum_over_m <- 0
  log_sum_over_m <- 0
  
  for(m in 1:M_theta){
    pf_vector_num <- log(pi_mw_thn[m,]) - t_i[unit] * exp(sum(beta1_thn[m,] * X1_cov[unit,]) + theta_star_thn[m,seq(1,M*K-2,M)]) + d_ci[unit,1] * (sum(beta1_thn[m,] * X1_cov[unit,]) + theta_star_thn[m,seq(1,M*K-2,M)]) - 
                                          t_i[unit] * exp(sum(beta2_thn[m,] * X2_cov[unit,]) + theta_star_thn[m,seq(2,M*K-1,M)]) + d_ci[unit,2] * (sum(beta2_thn[m,] * X2_cov[unit,]) + theta_star_thn[m,seq(2,M*K-1,M)]) - 
                                          t_i[unit] * exp(sum(beta3_thn[m,] * X3_cov[unit,]) + theta_star_thn[m,seq(3,M*K  ,M)]) + d_ci[unit,3] * (sum(beta3_thn[m,] * X3_cov[unit,]) + theta_star_thn[m,seq(3,M*K  ,M)])
    
    max_pf_n <- max(pf_vector_num)
    
    pf_vector_den <- log(pi_mw_thn[m,]) - t_i[unit] * exp(sum(beta2_thn[m,] * X2_cov[unit,]) + theta_star_thn[m,seq(2,M*K-1,M)]) + d_ci[unit,2] * (sum(beta2_thn[m,] * X2_cov[unit,]) + theta_star_thn[m,seq(2,M*K-1,M)]) - 
                                          t_i[unit] * exp(sum(beta3_thn[m,] * X3_cov[unit,]) + theta_star_thn[m,seq(3,M*K  ,M)]) + d_ci[unit,3] * (sum(beta3_thn[m,] * X3_cov[unit,]) + theta_star_thn[m,seq(3,M*K  ,M)])
    
    max_pf_d <- max(pf_vector_den)
    
    # - Average of f
    f_lapse <- exp(max_pf_n + log(sum(exp(pf_vector_num - max_pf_n))) - max_pf_d - log(sum(exp(pf_vector_den - max_pf_d))))
    
    sum_over_m <- sum_over_m + f_lapse
    
    # - Average of log(f)
    log_f_lapse <- max_pf_n + log(sum(exp(pf_vector_num - max_pf_n))) - max_pf_d - log(sum(exp(pf_vector_den - max_pf_d)))
    
    log_sum_over_m <- log_sum_over_m + log_f_lapse
    
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
WAIC_i_PH_DPM <- -2 * colSums(WAIC_insample)[1] + 2 * p_WAIC_i


sample_size_test <- nrow(uslapseagent_test)
X1_cov <- matrix(NA, sample_size_test, R1)
X2_cov <- matrix(NA, sample_size_test, R2)
X3_cov <- matrix(NA, sample_size_test, R3)

# - Surrending
X1_cov[,1] <- uslapseagent_test$annual.premium
X1_cov[,2] <- uslapseagent_test$termination.cause
X1_cov[,3] <- ifelse(uslapseagent_test$acc.death.rider=="Rider", 1, 0)
X1_cov[,4] <- ifelse(uslapseagent_test$gender=="Female", 1, 0)
X1_cov[,5] <- ifelse(uslapseagent_test$premium.frequency %in% c("Annual", 'Other'), 1, 0)
X1_cov[,6] <- ifelse(uslapseagent_test$underwriting.age %in% c('Young',"Middle"), 1, 0)

# - Death
X2_cov[,1] <- ifelse(uslapseagent_test$gender=="Female", 1, 0)
X2_cov[,2] <- ifelse(uslapseagent_test$underwriting.age=="Old", 1, 0)
X2_cov[,3] <- ifelse(uslapseagent_test$living.place=="Other", 1, 0)
X2_cov[,4] <- ifelse(uslapseagent_test$risk.state=="Smoker", 1, 0)

# - Other cause
X3_cov[,1] <- uslapseagent_test$annual.premium
X3_cov[,2] <- ifelse(uslapseagent_test$premium.frequency %in% c("Annual", 'Other'), 1, 0)
X3_cov[,3] <- ifelse(uslapseagent_test$acc.death.rider=="Rider", 1, 0)

t_i <- uslapseagent_test$duration
d_ci <- uslapseagent_test[,c("surrender", "death", "other")]


# - Predicted lapse rates
beta1_pm <- colMeans(beta1_thn)
beta2_pm <- colMeans(beta2_thn)
beta3_pm <- colMeans(beta3_thn)

pi_mw_m <- colMeans(pi_mw_thn)

theta_star_fin <- colMeans(theta_star_thn)

F_t_table_PH_DPM <- F_t_table

for(i in 1:sample_size_test){
  for(j in 1:n_quarters){
    F_t_table_PH_DPM[i,j] <-  (pi_mw_m[1]  * (1 - exp(-j * (exp(sum(beta1_pm * X1_cov[i,]) + theta_star_fin[1]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[2]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[3])))) * exp(theta_star_fin[1]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[1]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[2]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[3])) + 
                                 pi_mw_m[2]  * (1 - exp(-j * (exp(sum(beta1_pm * X1_cov[i,]) + theta_star_fin[4]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[5]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[6])))) * exp(theta_star_fin[4]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[4]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[5]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[6])) + 
                                 pi_mw_m[3]  * (1 - exp(-j * (exp(sum(beta1_pm * X1_cov[i,]) + theta_star_fin[7]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[8]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[9])))) * exp(theta_star_fin[7]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[7]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[8]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[9])) + 
                                 pi_mw_m[4]  * (1 - exp(-j * (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[10]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[11]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[12])))) * exp(theta_star_fin[10]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[10]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[11]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[12])) + 
                                 pi_mw_m[5]  * (1 - exp(-j * (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[13]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[14]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[15])))) * exp(theta_star_fin[13]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[13]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[14]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[15])) + 
                                 pi_mw_m[6]  * (1 - exp(-j * (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[16]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[17]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[18])))) * exp(theta_star_fin[16]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[16]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[17]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[18])) + 
                                 pi_mw_m[7]  * (1 - exp(-j * (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[19]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[20]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[21])))) * exp(theta_star_fin[19]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[19]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[20]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[21])) + 
                                 pi_mw_m[8]  * (1 - exp(-j * (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[22]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[23]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[24])))) * exp(theta_star_fin[22]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[22]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[23]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[24])) + 
                                 pi_mw_m[9]  * (1 - exp(-j * (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[25]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[26]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[27])))) * exp(theta_star_fin[25]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[25]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[26]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[27])) + 
                                 pi_mw_m[10] * (1 - exp(-j * (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[28]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[29]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[30])))) * exp(theta_star_fin[28]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[28]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[29]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[30])) + 
                                 pi_mw_m[11] * (1 - exp(-j * (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[31]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[32]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[33])))) * exp(theta_star_fin[31]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[31]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[32]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[33])) + 
                                 pi_mw_m[12] * (1 - exp(-j * (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[34]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[35]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[36])))) * exp(theta_star_fin[34]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[34]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[35]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[36])) + 
                                 pi_mw_m[13] * (1 - exp(-j * (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[37]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[38]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[39])))) * exp(theta_star_fin[37]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[37]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[38]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[39])) + 
                                 pi_mw_m[14] * (1 - exp(-j * (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[40]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[41]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[42])))) * exp(theta_star_fin[40]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[40]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[41]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[42])) + 
                                 pi_mw_m[15] * (1 - exp(-j * (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[43]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[44]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[45])))) * exp(theta_star_fin[43]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[43]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[44]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[45])) + 
                                 pi_mw_m[16] * (1 - exp(-j * (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[46]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[47]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[48])))) * exp(theta_star_fin[46]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[46]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[47]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[48])) + 
                                 pi_mw_m[17] * (1 - exp(-j * (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[49]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[50]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[51])))) * exp(theta_star_fin[49]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[49]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[50]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[51])) + 
                                 pi_mw_m[18] * (1 - exp(-j * (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[52]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[53]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[54])))) * exp(theta_star_fin[52]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[52]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[53]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[54])) + 
                                 pi_mw_m[19] * (1 - exp(-j * (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[55]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[56]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[57])))) * exp(theta_star_fin[55]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[55]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[56]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[57])) + 
                                 pi_mw_m[20] * (1 - exp(-j * (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[58]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[59]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[60])))) * exp(theta_star_fin[58]) / (exp( sum(beta1_pm * X1_cov[i,]) + theta_star_fin[58]) + exp( sum(beta2_pm * X2_cov[i,]) + theta_star_fin[59]) + exp( sum(beta3_pm * X3_cov[i,]) + theta_star_fin[60]))) * exp( sum(beta1_pm * X1_cov[i,]))
  }
}

# \hat{r}_{t,i} - individually predicted rates for computation of eq. (10) in Milhaud-Dutang
r_pred_table_PH_DPM <- matrix(NA, sample_size_test, n_quarters)
r_pred_table_PH_DPM[,1] <- F_t_table_PH_DPM[,1]
r_pred_table_PH_DPM[,2:n_quarters] <- (F_t_table_PH_DPM[,2:n_quarters] - F_t_table_PH_DPM[,1:(n_quarters-1)]) / (1 - F_t_table_PH_DPM[,1:(n_quarters-1)])

r_t_hat_PH_DPM <- colSums(r_pred_table_PH_DPM * risk_set_ind) / R_tq



