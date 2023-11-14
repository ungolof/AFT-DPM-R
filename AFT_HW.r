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
library(armspp)
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

K = 1 # Dunson (2010) claims a number of pieces between 20 and 50
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
theta_i <- rep(0, sample_size)

# - gamma 2 and gamma 3
gamma <- matrix(NA, n_iter, ncol = M-1)

# - Covariance_theta
sigma2_theta <- matrix(NA, n_iter, 1)

lpd <- rep(NA, n_iter)

############################### - Set prior distributions - #####################################

M=3
K=1

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



####################################### - Starting values - ################################################

# beta
beta1[1,] <- rnorm(R1, 0, 1)# c(0, 0.14, 0.8, -0.2, -0.07, -0.3, 0.01) #st_val
beta2[1,] <- rnorm(R2, 0, 1)# c(0, -0.0784, 0.2, 0.121, 0.2) #st_val
beta3[1,] <- rnorm(R3, 0, 1)# c(0,-0.0784, 0.2, 0.121) #st_val

# - gamma
gamma[1,] <- rnorm(M-1, 0, 1)

# - Covariance_theta
sigma2_theta[1:n_iter] <- 1


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




pnorm1_fc <- function(unit){
  
  y_output <- ifelse(d_ci[unit,1]==1, log(t_i[unit]), rtruncnorm(1, a=log(t_i[unit]), b=6, mean = (sum(beta1[i-1,] * X1_cov[unit,]) + theta_i[unit]), sd = sqrt(sigma2_c[i-1,1])))
  
  return(y_output)
}

pnorm2_fc <- function(unit){
  
  y_output <- ifelse(d_ci[unit,2]==1, log(t_i[unit]), rtruncnorm(1, a=log(t_i[unit]), b=6, mean = (sum(beta2[i-1,] * X2_cov[unit,]) + gamma[i-1,1] * theta_i[unit]), sd = sqrt(sigma2_c[i-1,2])))
  
  return(y_output)
}

pnorm3_fc <- function(unit){
  
  y_output <- ifelse(d_ci[unit,3]==1, log(t_i[unit]), rtruncnorm(1, a=log(t_i[unit]), b=6, mean = (sum(beta3[i-1,] * X3_cov[unit,]) + gamma[i-1,2] * theta_i[unit]), sd = sqrt(sigma2_c[i-1,3])))
  
  return(y_output)
}

theta_unit_f <- function(unit){
  target_d <- function(theta){
    #      mean_i <- c(sum(beta1[i,] * X1_cov[unit,]) + theta, sum(beta2[i,] * X2_cov[unit,]) + theta * gamma[i,1], sum(beta3[i,] * X3_cov[unit,]) + theta * gamma[i,1])
    inv_cov_i <- diag(1/sigma2_c[i-1,])
    
    #f_theta <- dmvnorm(y_i[unit,], mean = c(sum(beta1[i,] * X1_cov[unit,]) + theta, sum(beta2[i,] * X2_cov[unit,]) + theta * gamma[i,1], sum(beta3[i,] * X3_cov[unit,]) + theta * gamma[i,1]), sigma = cov_i, log=TRUE) - 0.5 * (theta^2) / sigma2_theta[i-1,1]
    
    f_theta <- - 0.5 * (y_i[unit,] - c(sum(beta1[i,] * X1_cov[unit,]) + theta, sum(beta2[i,] * X2_cov[unit,]) + theta * gamma[i,1], sum(beta3[i,] * X3_cov[unit,]) + theta * gamma[i,2])) %*% inv_cov_i %*% (y_i[unit,] - c(sum(beta1[i,] * X1_cov[unit,]) + theta, sum(beta2[i,] * X2_cov[unit,]) + theta * gamma[i,1], sum(beta3[i,] * X3_cov[unit,]) + theta * gamma[i,2])) - 0.5 * (theta^2) / sigma2_theta[i-1,1]
    
    return(f_theta)
  }
  
  theta_unit <- arms(1, target_d, -10, 10,  metropolis = FALSE, include_n_evaluations = TRUE)$samples   # adaptive_MH(i, n_samples, warm_up, thin)
  return(theta_unit)
}

theta_unit_sample <- function(unit){
  betax <- c(sum(beta1[i,] * X1_cov[unit,]), sum(beta2[i,] * X2_cov[unit,]), sum(beta3[i,] * X3_cov[unit,]))
  
  mean_i <- variance_i * (gamma_vector %*% Sigma_y_1 %*% y_i[unit,] - betax %*% Sigma_y_1 %*% gamma_vector)
  
  theta_unit <- rnorm(1, mean=mean_i, sd=sqrt(variance_i))
  return(theta_unit)
}



i <- 2
for(i in i:n_iter){
  #  start_time <- Sys.time()
  
  # - Step 0: Data Augment y_i
  y_i[,1] <- sapply(c(1:sample_size), FUN = pnorm1_fc)
  y_i[,2] <- sapply(c(1:sample_size), FUN = pnorm2_fc)
  y_i[,3] <- sapply(c(1:sample_size), FUN = pnorm3_fc)
  
  #    profvis({
  
  # - Step 3: Sample Beta_c and gamma
  
  Sigma_beta1_post <- solve(t(X1_cov / sigma2_c[i-1,1]) %*% X1_cov + diag(1/Sigma_beta1_prior))
  
  X2_cov_aug <- cbind(X2_cov, theta_i)
  Sigma_beta2_post <- solve(t(X2_cov_aug / sigma2_c[i-1,2]) %*% X2_cov_aug + diag(1/c(Sigma_beta2_prior, 16)))
  
  X3_cov_aug <- cbind(X3_cov, theta_i)
  Sigma_beta3_post <- solve(t(X3_cov_aug / sigma2_c[i-1,3]) %*% X3_cov_aug + diag(1/c(Sigma_beta3_prior, 16)))
  
  mu_beta1_post <- Sigma_beta1_post %*% (t(X1_cov / sigma2_c[i-1,1]) %*% (y_i[,1] - theta_i))
  
  mu_beta2_post <- Sigma_beta2_post %*% (t(X2_cov_aug / sigma2_c[i-1,2]) %*% y_i[,2])
  
  mu_beta3_post <- Sigma_beta3_post %*% (t(X3_cov_aug / sigma2_c[i-1,3]) %*% y_i[,3])
  
  beta1[i,] <- mvrnorm(n = 1, mu_beta1_post, Sigma_beta1_post)
  
  vec2 <- mvrnorm(n = 1, mu_beta2_post, Sigma_beta2_post)
  beta2[i,] <- vec2[1:R2]
  
  vec3 <- mvrnorm(n = 1, mu_beta3_post, Sigma_beta3_post)
  beta3[i,] <- mvrnorm(n = 1, mu_beta3_post, Sigma_beta3_post)[1:R3]
  
  gamma[i,1] <- vec2[R2+1]
  gamma[i,2] <- vec3[R3+1]
  
  # - Step 4: Sample theta_i
 
  gamma_vector <- c(1, gamma[i,1], gamma[i,2])
  Sigma_y_1 <- diag(1/sigma2_c[i-1,])
  variance_i <- 1 / (gamma_vector %*% Sigma_y_1 %*% gamma_vector + 1 / sigma2_theta[i,1])
  theta_i <- sapply(c(1:sample_size), FUN = theta_unit_sample)  
  
  # - Update of sigma2_theta
  shape_theta <- 2 + 0.5 * sample_size
  rate_theta <- 2 + sum(theta_i^2)
  sigma2_theta[i,1] <- 1 / rgamma(1, shape = shape_theta, rate = rate_theta)  
  
  # - Step 6: Update sigma_c from a conjugate Inverse-Gamma distribution
  suff_stat1 <- 0.5 * sum((y_i[,1] - beta1[i,1] * X1_cov[,1] - beta1[i,2] * X1_cov[,2] - beta1[i,3] * X1_cov[,3] - beta1[i,4] * X1_cov[,4] - beta1[i,5] * X1_cov[,5] - beta1[i,6] * X1_cov[,6] - beta1[i,7] * X1_cov[,7] -              theta_i)^2)
  suff_stat2 <- 0.5 * sum((y_i[,2] - beta2[i,1] * X2_cov[,1] - beta2[i,2] * X2_cov[,2] - beta2[i,3] * X2_cov[,3] - beta2[i,4] * X2_cov[,4] - beta2[i,5] * X2_cov[,5]                                                     - gamma[i,1] * theta_i)^2)
  suff_stat3 <- 0.5 * sum((y_i[,3] - beta3[i,1] * X3_cov[,1] - beta3[i,2] * X3_cov[,2] - beta3[i,3] * X3_cov[,3] - beta3[i,4] * X3_cov[,4]                                                                               - gamma[i,2] * theta_i)^2)
  
  sigma2_c[i,1] <- 1 / rgamma(1, shape = (1 + 0.5 * sample_size), rate = (1 + suff_stat1) )
  sigma2_c[i,2] <- 1 / rgamma(1, shape = (1 + 0.5 * sample_size), rate = (1 + suff_stat2) )
  sigma2_c[i,3] <- 1 / rgamma(1, shape = (1 + 0.5 * sample_size), rate = (1 + suff_stat3) )
  

  #     })
  #  end_time <- Sys.time()
  #  end_time - start_time
  print(paste(i, "% iteration"))
  
  if(i %in% seq(5000, n_iter, 5000)){
    save.image("~/Desktop/LFCR/Revision/AFT_HW.RData")
  }
  
}






plot(c(1:(i-1)), beta1[1:(i-1),1], type="l")
plot(c(1:(i-1)), beta1[1:(i-1),2], type="l")
plot(c(1:(i-1)), beta1[1:(i-1),3], type="l")
plot(c(1:(i-1)), beta1[1:(i-1),4], type="l")
plot(c(1:(i-1)), beta1[1:(i-1),5], type="l")
plot(c(1:(i-1)), beta1[1:(i-1),6], type="l")
plot(c(1:(i-1)), beta1[1:(i-1),7], type="l")

plot(c(1:(i-1)), beta2[1:(i-1),1], type="l")
plot(c(1:(i-1)), beta2[1:(i-1),2], type="l")
plot(c(1:(i-1)), beta2[1:(i-1),3], type="l")
plot(c(1:(i-1)), beta2[1:(i-1),4], type="l")
plot(c(1:(i-1)), beta2[1:(i-1),5], type="l")

plot(c(1:(i-1)), beta3[1:(i-1),1], type="l")
plot(c(1:(i-1)), beta3[1:(i-1),2], type="l")
plot(c(1:(i-1)), beta3[1:(i-1),3], type="l")
plot(c(1:(i-1)), beta3[1:(i-1),4], type="l")


## - 1) Discard first 40,000 iterations (burn-in) and thinning every 50
n_iter <- 10000
burn_in <- 5000
thinning <- 20

seq_bi_thn <- seq((burn_in+thinning), n_iter, thinning)

# - beta regression parameters
beta1_thn <- beta1[seq_bi_thn,]
beta2_thn <- beta2[seq_bi_thn,]
beta3_thn <- beta3[seq_bi_thn,]

# - gamma
gamma_thn <- gamma[seq_bi_thn,]

# - Covariance_theta
sigma2_theta_thn <- sigma2_theta[seq_bi_thn,]

sigma2_c_thn <- sigma2_c[seq_bi_thn,]

# - Simplified WAIC implementation (just considering the extra-variation in the lapse distribution)
log_lik_m_simple_ln <- function(unit){
  log_lik <- d_ci[unit,1] * dlnorm(t_i[unit], meanlog=sum(beta1_thn[m,] * X1_cov[unit,]), sdlog=sqrt(sigma2_c_thn[m,1] + sigma2_theta_thn[m]), log=TRUE) + (1 - d_ci[unit,1]) * plnorm(t_i[unit], meanlog=sum(beta1_thn[m,] * X1_cov[unit,]), sdlog=sqrt(sigma2_c_thn[m,1] + sigma2_theta_thn[m]), log.p=TRUE, lower.tail = FALSE)
  return(log_lik)
}

WAIC_pieces <- function(unit){ # returns first log(mean(f)) and then mean(log(f))
  
  M_theta <- length(seq_bi_thn)
  sum_over_m <- 0
  log_sum_over_m <- 0
  
  for(m in 1:M_theta){
    log_f_lapse <- log_lik_m_simple_ln(unit)
    
    log_sum_over_m <- log_sum_over_m + log_f_lapse
    
    # - Average of f
    f_lapse <- exp(log_f_lapse)      
    sum_over_m <- sum_over_m + f_lapse
    
  }
  return(c(-log(M_theta) + log(sum_over_m), log_sum_over_m/M_theta))
  
}

a=1
for(a in a:sample_size){
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
WAIC_i_AFT_HW <- -2 * colSums(WAIC_insample)[1] + 2 * p_WAIC_i
