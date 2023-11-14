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

# - gamma 2 and gamma 3
gamma <- matrix(NA, n_iter, ncol = M-1)


# - mu_theta (parameter of the base distribution)
#mu_theta <- matrix(NA, n_iter, M)
#colnames(mu_theta) <- c(sprintf("mu_theta_%d", c(1:M)))

# - Covariance_theta
sigma2_theta <- matrix(NA, n_iter, 1)

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

# - 4) Starting values

# - alpha
alpha1[1,1] <- rnorm(1, 0, 1)
alpha2[1,1] <- rnorm(1, 0, 1)
alpha3[1,1] <- rnorm(1, 0, 1)

# beta
beta1[1,] <- rnorm(R1, 0, 1)# c(0, 0.14, 0.8, -0.2, -0.07, -0.3, 0.01) #st_val
beta2[1,] <- rnorm(R2, 0, 1)# c(0, -0.0784, 0.2, 0.121, 0.2) #st_val
beta3[1,] <- rnorm(R3, 0, 1)# c(0,-0.0784, 0.2, 0.121) #st_val

# - gamma
gamma[1,] <- rnorm(M-1, 0, 1)

# - theta frailty parameters
theta_i <- rep(0, sample_size)

# - Covariance_theta
sigma2_theta[1,] <- 1



#================================ - MCMC sampling - ===============================
#  tryCatch({ 

# - Preparatory

N_k <- rep(0, K)
S_unit <- sample(c(1:K), sample_size, replace=TRUE) #rep(0, sample_size)

## - Adaptive Metropolis-H. Vihola (2012)
sigma_theta_pr <- rep(1, K)

sigma_beta1_pr <- rep(1, 4) * 0.001
sigma_beta2_pr <- rep(1, 4) * 0.001
sigma_beta3_pr <- rep(1, 2) * 0.001
gamma_Vihola <- 0.6 # - should be between 0.5 and 1
opt_accept <- 0.234 # Efficient acceptance rate from Roberts and Rosenthal (2001)
beta_tg_accept <- 0.234

sigma_gamma_pr<- rep(1, 2)

theta_unit_f <- function(unit){
  target_d <- function(theta){
    #      mean_i <- c(sum(beta1[i,] * X1_cov[unit,]) + theta, sum(beta2[i,] * X2_cov[unit,]) + theta * gamma[i,1], sum(beta3[i,] * X3_cov[unit,]) + theta * gamma[i,1])
    inv_cov_i <- diag(1/sigma2_c[i-1,])
    
    #f_theta <- dmvnorm(y_i[unit,], mean = c(sum(beta1[i,] * X1_cov[unit,]) + theta, sum(beta2[i,] * X2_cov[unit,]) + theta * gamma[i,1], sum(beta3[i,] * X3_cov[unit,]) + theta * gamma[i,1]), sigma = cov_i, log=TRUE) - 0.5 * (theta^2) / sigma2_theta[i-1,1]
    
    f_theta <- - t_i[unit] * exp(alpha1[i-1,1] + beta1[i-1,] * X1_cov[unit,] +                theta) + d_ci[unit,1] * (alpha1[i-1,1] + beta1[i-1,] * X1_cov[unit,] +                theta) - 
                 t_i[unit] * exp(alpha2[i-1,1] + beta2[i-1,] * X2_cov[unit,] + gamma[i-1,1] * theta) + d_ci[unit,2] * (alpha2[i-1,1] + beta2[i-1,] * X2_cov[unit,] + gamma[i-1,1] * theta) - 
                 t_i[unit] * exp(alpha3[i-1,1] + beta3[i-1,] * X3_cov[unit,] + gamma[i-1,2] * theta) + d_ci[unit,3] * (alpha3[i-1,1] + beta3[i-1,] * X3_cov[unit,] + gamma[i-1,2] * theta) - 
                 0.5 * (theta^2) / sigma2_theta[i-1,1]
    
    return(f_theta)
  }
  
  theta_unit <- arms(1, target_d, -10, 10,  metropolis = FALSE, include_n_evaluations = TRUE)$samples   # adaptive_MH(i, n_samples, warm_up, thin)
  return(theta_unit)
}


i <- 2
for(i in i:n_iter){
  #  start_time <- Sys.time()
  
  #     profvis({
  
  # - Step 3: Sample theta_i
  theta_i <- sapply(c(1:sample_size), FUN = theta_unit_f)  
  
  # - Update of sigma2_theta
  shape_theta <- 2 + 0.5 * sample_size
  rate_theta <- 2 + sum(theta_i^2)
  sigma2_theta[i,1] <- 1 / rgamma(1, shape = shape_theta, rate = rate_theta)  
  
  
  # - Step 6: Sample alpha_c
  ## - alpha1
  Omega11 <- omega11 + sum(d_ci[,1])
  Omega12 <- omega12
  
  Omega12 <- Omega12 + sum(t_i * exp(beta1[i-1,1] * X1_cov[,1] + beta1[i-1,2] * X1_cov[,2] + beta1[i-1,3] * X1_cov[,3] + beta1[i-1,4] * X1_cov[,4] + beta1[i-1,5] * X1_cov[,5] + beta1[i-1,6] * X1_cov[,6] + theta_i))
  
  vu1 <- rgamma(1, shape = Omega11, rate = Omega12)
  alpha1[i,1] <- log(vu1)
  
  ## - alpha2
  Omega21 <- omega21 + sum(d_ci[,2])
  Omega22 <- omega22
  
  Omega22 <- Omega22 + sum(t_i * exp(beta2[i-1,1] * X2_cov[,1] + beta2[i-1,2] * X2_cov[,2] + beta2[i-1,3] * X2_cov[,3] + beta2[i-1,4] * X2_cov[,4] + gamma[i-1,1] * theta_i))
  
  vu2 <- rgamma(1, shape = Omega21, rate = Omega22)
  alpha2[i,1] <- log(vu2)
  
  ## - alpha3
  Omega31 <- omega31 + sum(d_ci[,3])
  Omega32 <- omega32
  
  Omega32 <- Omega32 + sum(t_i * exp(beta3[i-1,1] * X3_cov[,1] + beta3[i-1,2] * X3_cov[,2] + beta3[i-1,3] * X3_cov[,3] + gamma[i-1,2] * theta_i))
  
  vu3 <- rgamma(1, shape = Omega31, rate = Omega32)
  alpha3[i,1] <- log(vu3)
  
  
  # - Step 7: Sample beta_c
  
  ## - Acceptance - Rejection sampling
  
  ## - Beta1
  ### - r=1
  # beta*
  beta_star1 <- rnorm(n = 1, beta1[i-1,1], sigma_beta1_pr[1])
  # Acceptance step
  log_ef1 <- dnorm(beta_star1, mean=mu_beta1[1], sd=sigma_beta1[1], log=TRUE) - dnorm(beta_star1, mean=beta1[i-1,1], sd=sigma_beta1_pr[1], log=TRUE) - dnorm(beta1[i-1,1], mean=mu_beta1[1], sd=sigma_beta1[1], log=TRUE) + dnorm(beta1[i-1,1], mean=beta_star1, sd=sigma_beta1_pr[1], log=TRUE)
  log_ef1 <- log_ef1 + sum( - t_i * exp(alpha1[i,1] + beta1[i-1,2] * X1_cov[,2] + beta1[i-1,3] * X1_cov[,3] + beta1[i-1,4] * X1_cov[,4] + beta1[i-1,5] * X1_cov[,5] + beta1[i-1,6] * X1_cov[,6] + theta_i) * (exp(beta_star1 * X1_cov[,1]) - exp(beta1[i-1,1] * X1_cov[,1])) + d_ci[,1] * X1_cov[,1] * (beta_star1 - beta1[i-1,1]) )
  
  accept_beta1 <- min(1, exp(log_ef1))
  beta_unif1 <- runif(1, 0, 1)
  beta1[i,1] <- ifelse(beta_unif1 > accept_beta1, beta1[i-1,1], beta_star1)
  # - Adjust sigma_beta
  sigma_beta1_pr[1] <- sigma_beta1_pr[1] * (1 + (i^(-gamma_Vihola)) * (accept_beta1 - beta_tg_accept))^0.5
  
  ### - r=2
  # beta*
  beta_star1 <- rnorm(n = 1, beta1[i-1,2], sigma_beta1_pr[2])
  # Acceptance step
  log_ef1 <- dnorm(beta_star1, mean=mu_beta1[2], sd=sigma_beta1[2], log=TRUE) - dnorm(beta_star1, mean=beta1[i-1,2], sd=sigma_beta1_pr[2], log=TRUE) - dnorm(beta1[i-1,2], mean=mu_beta1[2], sd=sigma_beta1[2], log=TRUE) + dnorm(beta1[i-1,2], mean=beta_star1, sd=sigma_beta1_pr[2], log=TRUE)
  log_ef1 <- log_ef1 + sum( - t_i * exp(alpha1[i,1] + beta1[i,1] * X1_cov[,1] + beta1[i-1,3] * X1_cov[,3] + beta1[i-1,4] * X1_cov[,4] + beta1[i-1,5] * X1_cov[,5] + beta1[i-1,6] * X1_cov[,6] + theta_i) * (exp(beta_star1 * X1_cov[,2]) - exp(beta1[i-1,2] * X1_cov[,2])) + d_ci[,1] * X1_cov[,2] * (beta_star1 - beta1[i-1,2]) )
  
  accept_beta1 <- min(1, exp(log_ef1))
  beta_unif1 <- runif(1, 0, 1)
  beta1[i,2] <- ifelse(beta_unif1 > accept_beta1, beta1[i-1,2], beta_star1)
  # - Adjust sigma_beta
  sigma_beta1_pr[2] <- sigma_beta1_pr[2] * (1 + (i^(-gamma_Vihola)) * (accept_beta1 - beta_tg_accept))^0.5
  
  
  ## - Beta3
  ### - r=1
  # beta*
  beta_star3 <- rnorm(n = 1, beta3[i-1,1], sigma_beta3_pr[1])
  # Acceptance step
  log_ef3 <- dnorm(beta_star3, mean=mu_beta3[1], sd=sigma_beta3[1], log=TRUE) - dnorm(beta_star3, mean=beta3[i-1,1], sd=sigma_beta3_pr[1], log=TRUE) - dnorm(beta3[i-1,1], mean=mu_beta3[1], sd=sigma_beta3[1], log=TRUE) + dnorm(beta3[i-1,1], mean=beta_star3, sd=sigma_beta3_pr[1], log=TRUE)
  log_ef3 <- log_ef3 + sum( - t_i * exp(alpha3[i,1] + beta3[i-1,2] * X3_cov[,2] + beta3[i-1,3] * X3_cov[,3] + gamma[i-1,2] * theta_i) * (exp(beta_star3 * X3_cov[,1]) - exp(beta3[i-1,1] * X3_cov[,1])) + d_ci[,3] * X3_cov[,1] * (beta_star3 - beta3[i-1,1]) )
  
  accept_beta3 <- min(1, exp(log_ef3))
  beta_unif3 <- runif(1, 0, 1)
  beta3[i,1] <- ifelse(beta_unif3 > accept_beta3, beta3[i-1,1], beta_star3)
  # - Adjust sigma_beta
  sigma_beta3_pr[1] <- sigma_beta3_pr[1] * (1 + (i^(-gamma_Vihola)) * (accept_beta3 - beta_tg_accept))^0.5
  
  
  ## - beta 1 (conjugate modified gamma)
  ### - r=3
  Omega1_beta1 <- omega1_beta1[3] + sum(d_ci[,1] * X1_cov[,3])
  lc1 <- sum(t_i * X1_cov[,3] * exp(alpha1[i,1] + beta1[i,1] * X1_cov[,1] + beta1[i,2] * X1_cov[,2] + beta1[i-1,4] * X1_cov[,4] + beta1[i-1,5] * X1_cov[,5] + beta1[i-1,6] * X1_cov[,6] + theta_i))
  Omega2_beta1 <- omega2_beta1[3] + lc1
  beta1[i,3] <- log(rgamma(1, shape = Omega1_beta1, rate = Omega2_beta1))
  
  ### - r=4
  Omega1_beta1 <- omega1_beta1[4] + sum(d_ci[,1] * X1_cov[,4])
  lc1 <- sum(t_i * X1_cov[,4] * exp(alpha1[i,1] + beta1[i,1] * X1_cov[,1] + beta1[i,2] * X1_cov[,2] + beta1[i,3] * X1_cov[,3] + beta1[i-1,5] * X1_cov[,5] + beta1[i-1,6] * X1_cov[,6] + theta_i))
  Omega2_beta1 <- omega2_beta1[4] + lc1
  beta1[i,4] <- log(rgamma(1, shape = Omega1_beta1, rate = Omega2_beta1))
  
  ### - r=5
  Omega1_beta1 <- omega1_beta1[5] + sum(d_ci[,1] * X1_cov[,5])
  lc1 <- sum(t_i * X1_cov[,5] * exp(alpha1[i,1] + beta1[i,1] * X1_cov[,1] + beta1[i,2] * X1_cov[,2] + beta1[i,3] * X1_cov[,3] + beta1[i,4] * X1_cov[,4] + beta1[i-1,6] * X1_cov[,6] + theta_i))
  Omega2_beta1 <- omega2_beta1[5] + lc1
  beta1[i,5] <- log(rgamma(1, shape = Omega1_beta1, rate = Omega2_beta1))
  
  ### - r=6
  Omega1_beta1 <- omega1_beta1[6] + sum(d_ci[,1] * X1_cov[,6])
  lc1 <- sum(t_i * X1_cov[,6] * exp(alpha1[i,1] + beta1[i,1] * X1_cov[,1] + beta1[i,2] * X1_cov[,2] + beta1[i,3] * X1_cov[,3] + beta1[i,4] * X1_cov[,4] + beta1[i,5] * X1_cov[,5] + theta_i))
  Omega2_beta1 <- omega2_beta1[6] + lc1
  beta1[i,6] <- log(rgamma(1, shape = Omega1_beta1, rate = Omega2_beta1))
  
  
  ## - beta 2 (conjugate modified gamma)
  ### - r=1
  Omega1_beta2 <- omega1_beta2[1] + sum(d_ci[,2] * X2_cov[,1])
  lc2 <- sum(t_i * X2_cov[,1] * exp(alpha2[i,1] + beta2[i-1,2] * X2_cov[,2] + beta2[i-1,3] * X2_cov[,3] + beta2[i-1,4] * X2_cov[,4] + gamma[i-1,1] * theta_i))
  Omega2_beta2 <- omega2_beta2[1] + lc2
  beta2[i,1] <- log(rgamma(1, shape = Omega1_beta2, rate = Omega2_beta2))
  
  ### - r=2
  Omega1_beta2 <- omega1_beta2[2] + sum(d_ci[,2] * X2_cov[,2])
  lc2 <- sum(t_i * X2_cov[,2] * exp(alpha2[i,1] + beta2[i,1] * X2_cov[,1] + beta2[i-1,3] * X2_cov[,3] + beta2[i-1,4] * X2_cov[,4] + gamma[i-1,1] * theta_i))
  Omega2_beta2 <- omega2_beta2[2] + lc2
  beta2[i,2] <- log(rgamma(1, shape = Omega1_beta2, rate = Omega2_beta2))
  
  ### - r=3
  Omega1_beta2 <- omega1_beta2[3] + sum(d_ci[,2] * X2_cov[,3])
  lc2 <- sum(t_i * X2_cov[,3] * exp(alpha2[i,1] + beta2[i,1] * X2_cov[,1] + beta2[i,2] * X2_cov[,2] + beta2[i-1,4] * X2_cov[,4] + gamma[i-1,1] * theta_i))
  Omega2_beta2 <- omega2_beta2[3] + lc2
  beta2[i,3] <- log(rgamma(1, shape = Omega1_beta2, rate = Omega2_beta2))
  
  ### - r=4
  Omega1_beta2 <- omega1_beta2[4] + sum(d_ci[,2] * X2_cov[,4])
  lc2 <- sum(t_i * X2_cov[,4] * exp(alpha2[i,1] + beta2[i,1] * X2_cov[,1] + beta2[i,2] * X2_cov[,2] + beta2[i,3] * X2_cov[,3] + gamma[i-1,1] * theta_i))
  Omega2_beta2 <- omega2_beta2[4] + lc2
  beta2[i,4] <- log(rgamma(1, shape = Omega1_beta2, rate = Omega2_beta2))
  
  ## - beta3
  ### - r=2
  Omega1_beta3 <- omega1_beta3[2] + sum(d_ci[,3] * X3_cov[,2])
  lc2 <- sum(t_i * X3_cov[,2] * exp(alpha3[i,1] + beta3[i,1] * X3_cov[,1] + beta3[i-1,3] * X3_cov[,3] + gamma[i-1,2] * theta_i))
  Omega2_beta3 <- omega2_beta3[2] + lc2
  beta3[i,2] <- log(rgamma(1, shape = Omega1_beta3, rate = Omega2_beta3))
  
  ### - r=3
  Omega1_beta3 <- omega1_beta3[3] + sum(d_ci[,3] * X3_cov[,3])
  lc2 <- sum(t_i * X3_cov[,3] * exp(alpha3[i,1] + beta3[i,1] * X3_cov[,1] + beta3[i,2] * X3_cov[,2] + gamma[i-1,2] * theta_i))
  Omega2_beta3 <- omega2_beta3[3] + lc2
  beta3[i,3] <- log(rgamma(1, shape = Omega1_beta3, rate = Omega2_beta3))
  
  # - Step 8: Sample gamma
  # gamma2*
  gamma_star <- rnorm(n = 1, gamma[i-1,1], sigma_gamma_pr[1])
  # Acceptance step
  log_ef <- dnorm(gamma_star, mean=0, sd=3, log=TRUE) - dnorm(gamma_star, mean=gamma[i-1,1], sd=sigma_gamma_pr[1], log=TRUE) - dnorm(gamma[i-1,1], mean=0, sd=3, log=TRUE) + dnorm(gamma[i-1,1], mean=gamma_star, sd=sigma_gamma_pr[1], log=TRUE)
  log_ef <- log_ef + sum( - t_i * exp(alpha2[i,1] + beta2[i,1] * X2_cov[,1] + beta2[i,2] * X2_cov[,2] + beta2[i,3] * X2_cov[,3] + beta2[i,4] * X2_cov[,4]) * (exp(gamma_star * theta_i) - exp(gamma[i-1,1] * theta_i)) + d_ci[,2] * theta_i * (gamma_star - gamma[i-1,1]) )
  
  accept_gamma <- min(1, exp(log_ef))
  gamma_unif <- runif(1, 0, 1)
  gamma[i,1] <- ifelse(gamma_unif > accept_gamma, gamma[i-1,1], gamma_star)
  # - Adjust sigma_beta
  sigma_gamma_pr[1] <- sigma_gamma_pr[1] * (1 + (i^(-gamma_Vihola)) * (accept_gamma - opt_accept))^0.5
  
  # gamma3*
  gamma_star <- rnorm(n = 1, gamma[i-1,2], sigma_gamma_pr[2])
  # Acceptance step
  log_ef <- dnorm(gamma_star, mean=0, sd=3, log=TRUE) - dnorm(gamma_star, mean=gamma[i-1,2], sd=sigma_gamma_pr[2], log=TRUE) - dnorm(gamma[i-1,2], mean=0, sd=3, log=TRUE) + dnorm(gamma[i-1,2], mean=gamma_star, sd=sigma_gamma_pr[2], log=TRUE)
  log_ef <- log_ef + sum( - t_i * exp(alpha3[i,1] + beta3[i,1] * X3_cov[,1] + beta3[i,2] * X3_cov[,2] + beta3[i,3] * X3_cov[,3]) * (exp(gamma_star * theta_i) - exp(gamma[i-1,2] * theta_i)) + d_ci[,3] * theta_i * (gamma_star - gamma[i-1,2]) )
  
  accept_gamma <- min(1, exp(log_ef))
  gamma_unif <- runif(1, 0, 1)
  gamma[i,2] <- ifelse(gamma_unif > accept_gamma, gamma[i-1,2], gamma_star)
  # - Adjust sigma_beta
  sigma_gamma_pr[2] <- sigma_gamma_pr[2] * (1 + (i^(-gamma_Vihola)) * (accept_gamma - opt_accept))^0.5
  
  print(paste(i, "% iteration"))
  
  if(i %in% seq(5000, n_iter, 5000)){
    save.image("~/Desktop/LFCR/Revision/PH_HW.RData")
  }
  
  #     })
  #  end_time <- Sys.time()
  #  end_time - start_time
}

# - Convergence check
plot(c(1:(i-1)), alpha1[1:(i-1),1], type="l")
plot(c(1:(i-1)), alpha2[1:(i-1),1], type="l")
plot(c(1:(i-1)), alpha3[1:(i-1),1], type="l")

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




##### - Build WAIC

# - computing the average over the posterior sample
burn_in <- 4000
thinning <- 20
n_iter <- 6000
seq_bi_thn <- seq((burn_in+thinning), n_iter, thinning)

alpha1_thn <- alpha1[seq_bi_thn,1]
alpha2_thn <- alpha2[seq_bi_thn,1]
alpha3_thn <- alpha3[seq_bi_thn,1]

# - beta regression parameters
beta1_thn <- beta1[seq_bi_thn,]
beta2_thn <- beta2[seq_bi_thn,]
beta3_thn <- beta3[seq_bi_thn,]

# - theta frailty parameters
sigma2_theta_thn <- sigma2_theta[seq_bi_thn]

# - gamma
gamma_thn <- gamma[seq_bi_thn,]

## - log(mean(f))


WAIC_pieces <- function(unit){ # returns first log(mean(f)) and then mean(log(f))
  
  M_theta <- length(seq_bi_thn)
  sum_over_m <- 0
  log_sum_over_m <- 0
  
  for(m in 1:M_theta){
    integrand_WAIC_PH_HW_num <- function(theta){
      dnorm(theta, mean=0, sd=sqrt(sigma2_theta_thn[m])) * exp(- t_i[unit] * (exp(alpha1_thn[m] + sum(beta1_thn[m,] * X1_cov[unit,]) +                  theta)) + d_ci[unit,1] * (alpha1_thn[m] + sum(beta1_thn[m,] * X1_cov[unit,]) +                  theta) -
                                                                 t_i[unit] * (exp(alpha2_thn[m] + sum(beta2_thn[m,] * X2_cov[unit,]) + gamma_thn[m,1] * theta)) + d_ci[unit,2] * (alpha2_thn[m] + sum(beta2_thn[m,] * X2_cov[unit,]) + gamma_thn[m,1] * theta) -
                                                                 t_i[unit] * (exp(alpha3_thn[m] + sum(beta3_thn[m,] * X3_cov[unit,]) + gamma_thn[m,2] * theta)) + d_ci[unit,3] * (alpha3_thn[m] + sum(beta3_thn[m,] * X3_cov[unit,]) + gamma_thn[m,2] * theta))
    }
    
    integrand_WAIC_PH_HW_den <- function(theta){
      dnorm(theta, mean=0, sd=sqrt(sigma2_theta_thn[m])) * exp(- t_i[unit] * (exp(alpha2_thn[m] + sum(beta2_thn[m,] * X2_cov[unit,]) + gamma_thn[m,1] * theta)) + d_ci[unit,2] * (alpha2_thn[m] + sum(beta2_thn[m,] * X2_cov[unit,]) + gamma_thn[m,1] * theta) -
                                                                 t_i[unit] * (exp(alpha3_thn[m] + sum(beta3_thn[m,] * X3_cov[unit,]) + gamma_thn[m,2] * theta)) + d_ci[unit,3] * (alpha3_thn[m] + sum(beta3_thn[m,] * X3_cov[unit,]) + gamma_thn[m,2] * theta))
    }
    
    pf_vector_num <- integrate(integrand_WAIC_PH_HW_num, lower=-30, upper=30)$value    
    pf_vector_den <- integrate(integrand_WAIC_PH_HW_den, lower=-30, upper=30)$value
    
    # - Average of f
    f_lapse <- pf_vector_num / pf_vector_den
    
    sum_over_m <- sum_over_m + f_lapse
    
    # - Average of log(f)
    log_f_lapse <- log(f_lapse)
    
    log_sum_over_m <- log_sum_over_m + log_f_lapse
    
  }
  return(c(-log(M_theta) + log(sum_over_m), log_sum_over_m/M_theta))
  
}

WAIC_insample <- matrix(NA, sample_size, 2)
a=1

for(a in 1:sample_size){
  start_time <- Sys.time()
  #  tryCatch({
  WAIC_insample[a,] <- WAIC_pieces(a)
  #  },
  #  error=function(error_message) {
  #    message("This is my custom message.")
  #  })
  end_time <- Sys.time()
  
  print(end_time - start_time)
  print(a)
  print(WAIC_insample[a,])
  
  #print(WAIC_insample[a,])
  
  
}

p_WAIC_i <- 2 * (colSums(WAIC_insample)[1] - colSums(WAIC_insample)[2])
WAIC_i_AFT_DPM <- -2 * colSums(WAIC_insample)[1] + 2 * p_WAIC_i


# - Simplified version of the HW calculation

WAIC_pieces <- function(unit){ # returns first log(mean(f)) and then mean(log(f))
  
  M_theta <- length(seq_bi_thn)
  sum_over_m <- 0
  log_sum_over_m <- 0
  
  for(m in 1:M_theta){
    theta_sample <- rnorm(10, mean=0, sd=sqrt(sigma2_theta_thn[m]))
    f_lapse <- mean(exp(- exp(alpha1_thn[m] + sum(beta1_thn[m,] * X1_cov[unit,]) +                  theta_sample) + d_ci[unit,1] * (alpha1_thn[m] + sum(beta1_thn[m,] * X1_cov[unit,]) +                  theta_sample)))
    
    # - Average of f
    sum_over_m <- sum_over_m + f_lapse
    
    # - Average of log(f)
    log_f_lapse <- log(f_lapse)
    
    log_sum_over_m <- log_sum_over_m + log_f_lapse
    
  }
  return(c(-log(M_theta) + log(sum_over_m), log_sum_over_m/M_theta))
  
}

# - Conservative estimation of the number of parameters (8: alpha, beta for the six covariates and sigma_theta)
## - It eventually makes more sense to estimate a DIC-like information criteria

# - Compute prediction rates

# - Covariates storage
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
## - 1) Discard first 40,000 iterations (burn-in) and thinning every 50
burn_in <- 4000
thinning <- 20
n_iter <- 6000
seq_bi_thn <- seq((burn_in+thinning), n_iter, thinning)

alpha1_thn <- alpha1[seq_bi_thn,1]
alpha2_thn <- alpha2[seq_bi_thn,1]
alpha3_thn <- alpha3[seq_bi_thn,1]

# - beta regression parameters
beta1_thn <- beta1[seq_bi_thn,]
beta2_thn <- beta2[seq_bi_thn,]
beta3_thn <- beta3[seq_bi_thn,]

# - theta frailty parameters
sigma2_theta_thn <- sigma2_theta[seq_bi_thn]

# - gamma
gamma_thn <- gamma[seq_bi_thn,]


alpha1_pm <- mean(alpha1_thn)
alpha2_pm <- mean(alpha2_thn)
alpha3_pm <- mean(alpha3_thn)

beta1_pm <- colMeans(beta1_thn)
beta2_pm <- colMeans(beta2_thn)
beta3_pm <- colMeans(beta3_thn)

gamma_pm <- colMeans(gamma_thn)

sigma2_theta_pm <- mean(sigma2_theta_thn)

integrand_PH_HW <- function(theta){
  dnorm(theta, mean=0, sd=sqrt(sigma2_theta_pm)) * ((1 - exp(-j * (exp(alpha1_pm + sum(beta1_pm * X1_cov[i,]) + theta) + exp(alpha2_pm + sum(beta2_pm * X2_cov[i,]) + gamma_pm[1] * theta) + exp(alpha3_pm + sum(beta3_pm * X3_cov[i,]) + gamma_pm[2] * theta)))) * exp(theta_star_fin[1]) / (exp(alpha1_pm + sum(beta1_pm * X1_cov[i,]) + theta) + exp(alpha2_pm + sum(beta2_pm * X2_cov[i,]) + gamma_pm[1] * theta) + exp(alpha3_pm + sum(beta3_pm * X3_cov[i,]) + gamma_pm[2] * theta))  ) * exp(alpha1_pm + sum(beta1_pm * X1_cov[i,]) + theta)
}



F_t_table_PH_HW <- F_t_table

for(i in 1:sample_size_test){
  for(j in 1:n_quarters){
    integrand_PH_HW <- function(theta){
      dnorm(theta, mean=0, sd=sqrt(sigma2_theta_pm)) * ((1 - exp(-j * (exp(alpha1_pm + sum(beta1_pm * X1_cov[i,]) + theta) + exp(alpha2_pm + sum(beta2_pm * X2_cov[i,]) + gamma_pm[1] * theta) + exp(alpha3_pm + sum(beta3_pm * X3_cov[i,]) + gamma_pm[2] * theta)))) * exp(theta) / (exp(alpha1_pm + sum(beta1_pm * X1_cov[i,]) + theta) + exp(alpha2_pm + sum(beta2_pm * X2_cov[i,]) + gamma_pm[1] * theta) + exp(alpha3_pm + sum(beta3_pm * X3_cov[i,]) + gamma_pm[2] * theta))  ) * exp(alpha1_pm + sum(beta1_pm * X1_cov[i,]))
    }
    
    F_t_table_PH_HW[i,j] <- integrate(integrand_PH_HW, lower=-10, upper=10)$value
  }
}

# \hat{r}_{t,i} - individually predicted rates for computation of eq. (10) in Milhaud-Dutang
r_pred_table_PH_HW <- matrix(NA, sample_size_test, n_quarters)
r_pred_table_PH_HW[,1] <- F_t_table_PH_HW[,1]
r_pred_table_PH_HW[,2:n_quarters] <- (F_t_table_PH_HW[,2:n_quarters] - F_t_table_PH_HW[,1:(n_quarters-1)]) / (1 - F_t_table_PH_HW[,1:(n_quarters-1)])

r_t_hat_PH_HW <- colSums(r_pred_table_PH_HW * risk_set_ind) / R_tq







