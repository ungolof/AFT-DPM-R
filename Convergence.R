library(matrixStats)


# - Traceplots analysis

## - beta surrending
plot(c(1:(i-1)), beta1[c(1:(i-1)),1], type='l')
plot(c(1:(i-1)), beta1[c(1:(i-1)),2], type='l')
plot(c(1:(i-1)), beta1[c(1:(i-1)),3], type='l')
plot(c(1:(i-1)), beta1[c(1:(i-1)),4], type='l')
plot(c(1:(i-1)), beta1[c(1:(i-1)),5], type='l')
plot(c(1:(i-1)), beta1[c(1:(i-1)),6], type='l')
plot(c(1:(i-1)), beta1[c(1:(i-1)),7], type='l')

## - beta death
plot(c(1:(i-1)), beta2[c(1:(i-1)),1], type='l')
plot(c(1:(i-1)), beta2[c(1:(i-1)),2], type='l')
plot(c(1:(i-1)), beta2[c(1:(i-1)),3], type='l')
plot(c(1:(i-1)), beta2[c(1:(i-1)),4], type='l')
plot(c(1:(i-1)), beta2[c(1:(i-1)),5], type='l')

## - beta others
plot(c(1:(i-1)), beta3[c(1:(i-1)),1], type='l')
plot(c(1:(i-1)), beta3[c(1:(i-1)),2], type='l')
plot(c(1:(i-1)), beta3[c(1:(i-1)),3], type='l')
plot(c(1:(i-1)), beta3[c(1:(i-1)),4], type='l')

## - phi
plot(c(1:(i-1)), phi[c(1:(i-1)),1], type='l')


## - mu_theta
plot(c(1:(i-1)), mu_theta[c(1:(i-1)),1], type='l')
plot(c(1:(i-1)), mu_theta[c(1:(i-1)),2], type='l')
plot(c(1:(i-1)), mu_theta[c(1:(i-1)),3], type='l')


## - theta* as drawn from the Dirichlet Process
plot(c(1:(i-1)), theta_star[c(1:(i-1)),4], type='l')
plot(c(1:(i-1)), theta_star[c(1:(i-1)),5], type='l')
plot(c(1:(i-1)), theta_star[c(1:(i-1)),6], type='l')

plot(c(1:(i-1)), theta_star[c(1:(i-1)),7], type='l')
plot(c(1:(i-1)), theta_star[c(1:(i-1)),8], type='l')
plot(c(1:(i-1)), theta_star[c(1:(i-1)),9], type='l')

plot(c(1:(i-1)), theta_star[c(1:(i-1)),10], type='l')
plot(c(1:(i-1)), theta_star[c(1:(i-1)),11], type='l')
plot(c(1:(i-1)), theta_star[c(1:(i-1)),12], type='l')

# - pi (mixture weights from the Dirichlet Process)
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),1], type='l')
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),2], type='l')
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),3], type='l')
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),4], type='l')
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),5], type='l')
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),6], type='l')
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),7], type='l')
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),8], type='l')
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),9], type='l')
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),10], type='l')
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),11], type='l')

## - Number of clusters and cluster occupancy
plot(c(1:(i-1)), N_c_tb[c(1:(i-1))], type='l')
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),1], type='l')
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),2], type='l')
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),3], type='l')
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),4], type='l')
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),5], type='l')
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),6], type='l')


## - Log-posterior density

lpd_fnc <- function(unit){
  lpd_i <- 0
  pf_vector <- rep(0, K)
  th_k_pos <- c(1:M) # position vector for theta

  piece1 <- sum(beta1[i,] * X1_cov[unit,])
  piece2 <- sum(beta2[i,] * X2_cov[unit,]) # t_i[unit] * exp(alpha2[i-1,1] + sum(beta2[i-1,] * X2_cov[unit,]))
  piece3 <- sum(beta3[i,] * X3_cov[unit,])

  pf_vector <- log(pi_mw[i,]) + ifelse(d_ci[unit,1]==1, dnorm(log(t_i[unit]), mean = (piece1 + theta_star[i,seq(1, 73, 3)]), sd = sqrt(sigma2_c[i,1]), log=TRUE), pnorm(log(t_i[unit]), mean = (piece1 + theta_star[i,seq(1, 73, 3)]), sd = sqrt(sigma2_c[i,1]), log=TRUE, lower.tail = FALSE)) +
                                ifelse(d_ci[unit,2]==1, dnorm(log(t_i[unit]), mean = (piece2 + theta_star[i,seq(2, 74, 3)]), sd = sqrt(sigma2_c[i,2]), log=TRUE), pnorm(log(t_i[unit]), mean = (piece2 + theta_star[i,seq(2, 74, 3)]), sd = sqrt(sigma2_c[i,2]), log=TRUE, lower.tail = FALSE)) +
                                ifelse(d_ci[unit,3]==1, dnorm(log(t_i[unit]), mean = (piece3 + theta_star[i,seq(3, 75, 3)]), sd = sqrt(sigma2_c[i,3]), log=TRUE), pnorm(log(t_i[unit]), mean = (piece3 + theta_star[i,seq(3, 75, 3)]), sd = sqrt(sigma2_c[i,3]), log=TRUE, lower.tail = FALSE))

  max_pf <- max(pf_vector)

  lpd_i <- max_pf + log(sum(exp(pf_vector-max_pf)))
  return(lpd_i)
}

lpd <- rep(0, n_iter)
for(i in 1:n_iter){
  #  profvis({
  # - log-likelihood
  lpd[i] <- 0
  lpd_i_ct <- sapply(c(1:sample_size), FUN = lpd_fnc)

  lpd[i] <- lpd[i] + sum(lpd_i_ct)

  # - variance
  lpd[i] <- lpd[i] + sum(dgamma(1/sigma2_c[i,], shape = 1, rate = 1, log=TRUE))

  # - theta*
  th_k_pos <- c(1:M)
  log_sum <- 0
  for(k in 1:K){
    log_sum <- log_sum + dmvnorm(theta_star[i,th_k_pos], mu_theta[i,], sigma = matrix(Sigma_theta[i,], M,M), log = TRUE)
    th_k_pos <- th_k_pos + M
  }

  lpd[i] <- lpd[i] + log_sum

  # - psi from the stick-breaking process
  lpd[i] <- lpd[i] + sum(dbeta(psi_sbp[i,], 1, phi[i,1], ncp = 0, log = TRUE))

  # beta regression parameters
  lpd[i] <- lpd[i] + sum(dnorm(beta1[i,], mean=0, sd=1, log=TRUE))
  lpd[i] <- lpd[i] + sum(dnorm(beta2[i,], mean=0, sd=1, log=TRUE))
  lpd[i] <- lpd[i] + sum(dnorm(beta3[i,], mean=0, sd=1, log=TRUE))

  lpd[i] <- lpd[i] + dmvnorm(mu_theta[i,], lambda1, sigma = lambda2, log = TRUE)
  lpd[i] <- lpd[i] + dinvwishart(matrix(Sigma_theta[i,], M,M), nu = Lambda3, S = Lambda4, log=TRUE)
  lpd[i] <- lpd[i] + dgamma(phi[i,1], shape=lambda5, rate = lambda6, log = TRUE)
  #  })
  print(paste(i, "% iteration"))

  if(i %in% seq(5000, n_iter, 5000)){
      save.image("C:/Users/z3532149/Desktop/LFCR/AFT-DP/AFT_MCMC_v0_Conv.RData")
  }

}

plot(c(1:(i-1)), lpd[c(1:(i-1))], type='l')

par(mfrow = c(1, 2))
plot(c(1:50000), lpd, type='l', ylab='log-posterior density value', xlab="Iteration", main="Log-posterior density traceplot")
plot(density(lpd[seq(40050,50000,50)]), type='l', xlab='log-posterior d. value', main='Lpd d. (burn-in 40000, thin. e. 50')







cov2cor(cov(theta_star[,1:3]))
cov2cor(cov(theta_star[,4:6]))
cov2cor(cov(theta_star[,7:9]))
cov2cor(cov(theta_star[,10:12]))
cov2cor(cov(theta_star[,13:15]))
cov2cor(cov(theta_star[,16:18]))


#============================= - Get final draws - ========================
## - 1) Discard first 40,000 iterations (burn-in) and thinning every 50
burn_in <- 40000
thinning <- 20

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

# - Covariance_theta
Sigma_theta_thn <- Sigma_theta[seq_bi_thn,]

# - psi from stick-breaking procedure
psi_sbp_thn <- psi_sbp[seq_bi_thn,]

# phi concentration parameter of the Dirichlet Process
phi_thn <- phi[seq_bi_thn,]

sigma2_c_thn <- sigma2_c[seq_bi_thn,]

#============= - Posterior density of the parameters

# - mu_theta
plot(density(mu_theta_thn[,1]), xlab = "beta", main=" ")
abline(v=mean(mu_theta_thn[,1]), lty=1) # abline(v=log(phi0_p), lty=2)
abline(v=quantile(mu_theta_thn[,1], prob=0.975), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=quantile(mu_theta_thn[,1], prob=0.025), lty=2) # abline(v=log(phi0_p), lty=2)

plot(density(mu_theta_thn[,2]), xlab = "beta", main=" ") #, ylim=c(0, 10.3),
abline(v=mean(mu_theta_thn[,2]), lty=1) # abline(v=log(phi0_p), lty=2)
abline(v=quantile(mu_theta_thn[,2], prob=0.975), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=quantile(mu_theta_thn[,2], prob=0.025), lty=2) # abline(v=log(phi0_p), lty=2)

plot(density(mu_theta_thn[,3]), xlab = "beta", main=" ") #, ylim=c(0, 10.3),
abline(v=mean(mu_theta_thn[,3]), lty=1) # abline(v=log(phi0_p), lty=2)
abline(v=quantile(mu_theta_thn[,3], prob=0.975), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=quantile(mu_theta_thn[,3], prob=0.025), lty=2) # abline(v=log(phi0_p), lty=2)

plot(density(sigma2_c_thn[,1]))
plot(density(sigma2_c_thn[,2]))
plot(density(sigma2_c_thn[,3]))


plot(density(theta_star_thn[,1]))
plot(density(theta_star_thn[,2]))
plot(density(theta_star_thn[,3]))
plot(density(theta_star_thn[,4]))
plot(density(theta_star_thn[,5]))
plot(density(theta_star_thn[,6]))
plot(density(theta_star_thn[,7]))
plot(density(theta_star_thn[,8]))
plot(density(theta_star_thn[,9]))
plot(density(theta_star_thn[,10]))
plot(density(theta_star_thn[,11]))
plot(density(theta_star_thn[,12]))
plot(density(theta_star_thn[,13]))
plot(density(theta_star_thn[,14]))
plot(density(theta_star_thn[,15]))
plot(density(theta_star_thn[,16]))
plot(density(theta_star_thn[,17]))
plot(density(theta_star_thn[,18]))
plot(density(theta_star_thn[,19]))
plot(density(theta_star_thn[,20]))
plot(density(theta_star_thn[,21]))
plot(density(theta_star_thn[,22]))
plot(density(theta_star_thn[,23]))
plot(density(theta_star_thn[,24]))
plot(density(theta_star_thn[,25]))
plot(density(theta_star_thn[,26]))
plot(density(theta_star_thn[,27]))
plot(density(theta_star_thn[,28]))
plot(density(theta_star_thn[,29]))
plot(density(theta_star_thn[,30]))
plot(density(theta_star_thn[,31]))
plot(density(theta_star_thn[,32]))
plot(density(theta_star_thn[,33]))
plot(density(theta_star_thn[,34]))
plot(density(theta_star_thn[,35]))
plot(density(theta_star_thn[,36]))
plot(density(theta_star_thn[,37]))
plot(density(theta_star_thn[,38]))
plot(density(theta_star_thn[,39]))
plot(density(theta_star_thn[,40]))

# - theta_star

plot(density(theta_star_thn[,4]), xlab = "beta", main=" ")
abline(v=mean(theta_star_thn[,4]), lty=1) # abline(v=log(phi0_p), lty=2)
abline(v=quantile(theta_star_thn[,4], prob=0.975), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=quantile(theta_star_thn[,4], prob=0.025), lty=2) # abline(v=log(phi0_p), lty=2)

plot(density(theta_star_thn[,5]), xlab = "beta", main=" ") #, ylim=c(0, 10.3),
abline(v=mean(theta_star_thn[,5]), lty=1) # abline(v=log(phi0_p), lty=2)
abline(v=quantile(theta_star_thn[,5], prob=0.975), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=quantile(theta_star_thn[,5], prob=0.025), lty=2) # abline(v=log(phi0_p), lty=2)

plot(density(theta_star_thn[,6]), xlab = "beta", main=" ") #, ylim=c(0, 10.3),
abline(v=mean(theta_star_thn[,6]), lty=1) # abline(v=log(phi0_p), lty=2)
abline(v=quantile(theta_star_thn[,6], prob=0.975), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=quantile(theta_star_thn[,6], prob=0.025), lty=2) # abline(v=log(phi0_p), lty=2)

plot(density(phi_thn), xlab = "beta", main=" ") #, ylim=c(0, 10.3),
mean(phi_thn)

# - Analysis of covariates
beta1_lq <- round(colQuantiles(beta1_thn, rows = NULL, cols = NULL, probs = 0.025),4)
beta1_uq <- round(colQuantiles(beta1_thn, rows = NULL, cols = NULL, probs = 0.975),4)

beta2_lq <- round(colQuantiles(beta2_thn, rows = NULL, cols = NULL, probs = 0.025),4)
beta2_uq <- round(colQuantiles(beta2_thn, rows = NULL, cols = NULL, probs = 0.975),4)

beta3_lq <- round(colQuantiles(beta3_thn, rows = NULL, cols = NULL, probs = 0.025),4)
beta3_uq <- round(colQuantiles(beta3_thn, rows = NULL, cols = NULL, probs = 0.975),4)


round(colMeans(beta1_thn),4)
round(colMeans(beta2_thn),4)
round(colMeans(beta3_thn),4)

round(beta1_lq,4)
round(beta1_uq,4)

theta_s_lq <- round(colQuantiles(theta_star_thn, rows = NULL, cols = NULL, probs = 0.025),4)
theta_s_uq <- round(colQuantiles(theta_star_thn, rows = NULL, cols = NULL, probs = 0.975),4)
round(colMeans(theta_star_thn),4)[c(1,4,7,10,43)]
pi_mw_lq <- round(colQuantiles(pi_mw_thn, rows = NULL, cols = NULL, probs = 0.025),4)
pi_mw_uq <- round(colQuantiles(pi_mw_thn, rows = NULL, cols = NULL, probs = 0.975),4)


round(colMeans(sigma2_c_thn),4)
round(colQuantiles(sigma2_c_thn, rows = NULL, cols = NULL, probs = 0.025),4)
round(colQuantiles(sigma2_c_thn, rows = NULL, cols = NULL, probs = 0.975),4)

cbind(round(colMeans(beta1_thn),4), beta1_lq, beta1_uq, round(colMeans(beta2_thn),4), beta2_lq, beta2_uq, round(colMeans(beta3_thn),4), beta3_lq, beta3_uq)

cbind(round(colMeans(sigma2_c_thn),4), theta_s_lq, theta_s_uq, round(colMeans(beta2_thn),4), beta2_lq, beta2_uq, round(colMeans(beta3_thn),4), beta3_lq, beta3_uq)

cbind(round(colMeans(theta_star_thn[,seq(1,73,3)]),4), theta_s_lq[seq(1,73,3)], theta_s_uq[seq(1,73,3)], 
      round(colMeans(theta_star_thn[,seq(2,74,3)]),4), theta_s_lq[seq(2,74,3)], theta_s_uq[seq(2,74,3)], 
      round(colMeans(theta_star_thn[,seq(3,75,3)]),4), theta_s_lq[seq(3,75,3)], theta_s_uq[seq(3,75,3)])

cbind(round(colMeans(pi_mw_thn[,seq(3,24,3)]),4), pi_mw_lq[seq(3,24,3)], pi_mw_uq[seq(3,24,3)], 
      round(colMeans(pi_mw_thn[,seq(4,25,3)]),4), pi_mw_lq[seq(4,25,3)], pi_mw_uq[seq(4,25,3)],
      round(colMeans(pi_mw_thn[,seq(5,23,3)]),4), pi_mw_lq[seq(5,25,3)], pi_mw_uq[seq(5,25,3)])

cbind(round(mean(phi_thn),4), round(quantile(phi_thn, probs = 0.025),4), round(quantile(phi_thn, probs = 0.975),4))
cbind(round(mean(pi_mw_thn[,1]),4), round(quantile(pi_mw_thn[,1], probs = 0.025),4), round(quantile(pi_mw_thn[,1], probs = 0.975),4))
cbind(round(mean(pi_mw_thn[,2]),4), round(quantile(pi_mw_thn[,2], probs = 0.025),4), round(quantile(pi_mw_thn[,2], probs = 0.975),4))


#===================================== - Ergodic average - ================================
mu_theta_erg <- matrix(NA, nrow=n_iter, ncol=ncol(mu_theta))
for (i in 25000:n_iter){
  mu_theta_erg[i,] <- colSums(mu_theta[c(1:i), , drop = FALSE])/i
  # mu_theta_erg[i,] <- colSums(mu_theta[c(1:i), , drop = FALSE])/i
}

plot(c(1:n_iter), mu_theta_erg[,1], type='l')
lines(c(1:n_iter), mu_theta_erg[,2], col='blue')
lines(c(1:n_iter), mu_theta_erg[,3], col='red')


# - Plots for presentation
colMeans(beta1_thn)

par(mfrow = c(3, 2))
plot(density(beta1_thn[,1]), main="Annual premium (mean = -0.0410)", xlab="beta annual premium")
abline(v=mean(beta1_thn[,1]))
abline(v=beta1_lq[1], lty=3)
abline(v=beta1_uq[1], lty=3)
abline(v=0, lwd=2, col='red')
legend("topright", legend=c("Mean", "95% Cred. Int.", "0"), lty=c(1,3,1), lwd=c(1,1,2), col=c('black', 'black','red'), bty='n')


plot(density(beta1_thn[,2]), main="DJIA (mean = -0.6450)", xlab="beta DJIA")
abline(v=mean(beta1_thn[,2]))
abline(v=beta1_lq[2], lty=3)
abline(v=beta1_uq[2], lty=3)
abline(v=0, lwd=2, col='red')

plot(density(beta1_thn[,3]), main="Accidental death rider (mean = 0.0534)", xlab="beta ADR")
abline(v=mean(beta1_thn[,3]))
abline(v=beta1_lq[3], lty=3)
abline(v=beta1_uq[3], lty=3)
abline(v=0, lwd=2, col='red')

plot(density(beta1_thn[,4]), main="Gender (mean = 0.0310)", xlab="beta Gender")
abline(v=mean(beta1_thn[,4]))
abline(v=beta1_lq[4], lty=3)
abline(v=beta1_uq[4], lty=3)
abline(v=0, lwd=2, col='red')

plot(density(beta1_thn[,5]), main="Premium frequency (mean = 0.0744)", xlab="beta premium frequency")
abline(v=mean(beta1_thn[,5]))
abline(v=beta1_lq[5], lty=3)
abline(v=beta1_uq[5], lty=3)
abline(v=0, lwd=2, col='red')

plot(density(beta1_thn[,6]), main="Underwriting age (mean = -0.1036)", xlab="beta UW age")
abline(v=mean(beta1_thn[,6]))
abline(v=beta1_lq[6], lty=3)
abline(v=beta1_uq[6], lty=3)
abline(v=0, lwd=2, col='red')
legend("topright", legend=c("Mean", "95% Cred. Int.", "0"), lty=c(1,3,1), lwd=c(1,1,2), col=c('black', 'black','red'), bty='n')





#===================== - Group analysis - ==========================

N_c_tb_thn <- mean(N_c_tb[seq_bi_thn])
N_k_tb_thn <- colMeans(N_k_tb[seq_bi_thn,])
100 * N_k_tb_thn / sample_size

sum((100 * N_k_tb_thn / sample_size)[c(1:4,8,10,11,15)])

par(mfrow=c(1,1), mar = c(4.5,4,1,0.5))
barplot(N_k_tb_thn, main="", horiz=FALSE, xlab = "Mixture component (k)", ylab="Post. mean n_k", cex.names=0.7, names.arg=c(1:K))
abline(h=550)
# - Covariance in log(T)

n_post_draws <- 1000
## 1 - Sample the row d
row_draw <- sample(c(1:length(seq_bi_thn)), size=n_post_draws, replace=TRUE)

alpha_c_thn <- cbind(beta1_thn[,1], beta2_thn[,1], beta3_thn[,1])

n_t_draws <- 1000

t_draws <- matrix(NA, n_post_draws * n_t_draws, M)
count <- 1
for(i in 1:n_post_draws){
  for(c in 1:M){
    t_draws[count:(count + n_t_draws - 1),c] <- rnorm(n_t_draws, mean=(alpha_c_thn[row_draw[i],c] + theta_star_thn[row_draw[i],c]), sd=sqrt(sigma2_c_thn[row_draw[i],c]))
  }
  count <- count + n_t_draws
}
  
cov(t_draws)
cov2cor(cov(t_draws))

# - 4th class

count <- 1
for(i in 1:n_post_draws){
  for(c in 1:M){
    t_draws[count:(count + n_t_draws - 1),c] <- rnorm(n_t_draws, mean=(alpha_c_thn[row_draw[i],c] + theta_star_thn[row_draw[i],M*3 + c]), sd=sqrt(sigma2_c_thn[row_draw[i],c]))
  }
  count <- count + n_t_draws
}

cor(t_draws)
cov2cor(cov(t_draws))


cov2cor(matrix(colMeans(Sigma_theta_thn),3,3))

cor(theta_star_thn[,1:3])

# - Procedure 1: In terms of the generating process (Just one draw)

## 1 - Get posterior means of the parameters
mu_theta_post_m <- colMeans(mu_theta_thn)
Sigma_theta_post_m <- matrix(colMeans(Sigma_theta_thn),3,3)
phi_thn_post_m <- mean(phi_thn)

## 2 - Generete theta from the DP using the SBP
### - 2.1 Generate psi_k
psi_k_gen <- rep(0, K-1)
psi_k_gen <- rbeta(K-1, 1, phi_thn_post_m)

### - 2.2 Generate pi_k
pi_k <- rep(0, K)
pi_k[1] <- psi_k_gen[1]

for(k in 2:K){
  pi_k[k] <- psi_k_gen[k] * prod(1-psi_k_gen[1:(k-1)])
}

### - 2.3 Generate independent draws of theta
theta_dp_draw <- rep(0,M)
theta_dp_draw <- rmvnorm(1, mean = mu_theta_post_m, sigma = Sigma_theta_post_m, checkSymmetry = TRUE)



# - Procedure 2: In terms of the generating process (Multiple theta draws based on posterior means of the parameters)

## 1 - Get posterior means of the parameters
mu_theta_post_m <- colMeans(mu_theta_thn)
Sigma_theta_post_m <- matrix(colMeans(Sigma_theta_thn),3,3)
phi_thn_post_m <- mean(phi_thn)
alpha_c_post_m <- colMeans(alpha_c_thn)
sigma2_c_post_m <- colMeans(sigma2_c_thn)

n_draws <- 1000

## 2 - Generate theta from the DP using the SBP
### - 2.1 Generate psi_k
psi_k_gen <- rbeta(K, 1, phi_thn_post_m)

### - 2.2 Generate pi_k

pi_k <- rep(0, K) # matrix(0, n_draws, K)
pi_k[1] <- psi_k_gen[1] # pi_k[,1] <- psi_k_gen[,1]

for(k in 2:K){
  pi_k[k] <- psi_k_gen[k] * prod(1-psi_k_gen[1:(k-1)])
}

### - 2.3 Generate independent draws of theta
theta_dp_draw <- matrix(0, M, K)

for(k in 1:K){
theta_dp_draw[,k] <- rmvnorm(1, mean = mu_theta_post_m, sigma = Sigma_theta_post_m, checkSymmetry = TRUE)
}

s_sample <- sample(c(1:K), n_draws, replace=TRUE, prob=pi_k)

# - 3 Draw dependent t=exp(Y)
y_dep_draw <- matrix(NA, n_draws, M)

for(i in 1:n_draws){
  for(c in 1:M){
    y_dep_draw[i,c] <- rnorm(1, mean = theta_dp_draw[c,s_sample[i]], sd=sqrt(sigma2_c_post_m[c]))
  }
}

cor(exp(y_dep_draw))

#================ - Dependent time to event - =========================

pi_bar <- colMeans(pi_mw[seq_bi_thn,])
theta1_bar <- colMeans(theta_star[seq_bi_thn, seq(1,73, 3)])
theta2_bar <- colMeans(theta_star[seq_bi_thn, seq(2,74, 3)])
theta3_bar <- colMeans(theta_star[seq_bi_thn, seq(3,75, 3)])
beta1_bar <- colMeans(beta1[seq_bi_thn,])
beta2_bar <- colMeans(beta2[seq_bi_thn,])
beta3_bar <- colMeans(beta3[seq_bi_thn,])
sigma2_c_bar <- colMeans(sigma2_c[seq_bi_thn,])

n_samples <- 500000
s_sample <- rep(0, n_samples)
y_dep_sample <- matrix(NA, n_samples, M)

for(sample in 1:n_samples){
  s_sample[sample] <- sample(c(1:K), 1, prob=pi_bar)
  y_dep_sample[sample,1] <- rnorm(1, mean=theta1_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[1]))
  y_dep_sample[sample,2] <- rnorm(1, mean=theta2_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[2]))
  y_dep_sample[sample,3] <- rnorm(1, mean=theta3_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[3]))
}

cor(y_dep_sample)
cor(exp(y_dep_sample))


for(sample in 1:n_samples){
  s_sample[sample] <- sample(c(1:K), 1, prob=pi_bar)
  y_dep_sample[sample,1] <- rnorm(1, mean=theta1_bar[s_sample[sample]] + beta1_bar[3], sd=sqrt(sigma2_c_bar[1]))
  y_dep_sample[sample,2] <- rnorm(1, mean=theta2_bar[s_sample[sample]] + beta2_bar[4], sd=sqrt(sigma2_c_bar[2]))
  y_dep_sample[sample,3] <- rnorm(1, mean=theta3_bar[s_sample[sample]] + beta3_bar[3], sd=sqrt(sigma2_c_bar[3]))
}

cor(y_dep_sample)
cor(exp(y_dep_sample))




round(colMeans(pi_mw_thn),4)
round(theta_s1_lp, 4)[c(1:4,15)]



#=============== - Analysis of clusters by Bayes' rule - ===============

beta1_post_m <- colMeans(beta1_thn)
beta2_post_m <- colMeans(beta2_thn)
beta3_post_m <- colMeans(beta3_thn)
theta_star_post_m <- colMeans(theta_star_thn)

pi_post_m <- colMeans(pi_mw_thn)

theta_s1_lp <- theta_star_post_m[seq(1, 73, 3)]
theta_s2_lp <- theta_star_post_m[seq(2, 74, 3)]
theta_s3_lp <- theta_star_post_m[seq(3, 75, 3)]

sigma2_c_post_m <- colMeans(sigma2_c_thn)


post_prob_fc <- function(unit){
  
  # - use log sum exp trick to deal with very small numbers which R inappropriately set to zero
  pf_vector <- rep(0, K)
  th_k_pos <- c(1:M) # position vector for theta
  
  piece1 <- sum(beta1_post_m * X1_cov[unit,])
  piece2 <- sum(beta2_post_m * X2_cov[unit,]) # t_i[unit] * exp(alpha2[i-1,1] + sum(beta2[i-1,] * X2_cov[unit,]))
  piece3 <- sum(beta3_post_m * X3_cov[unit,])
  
  pf_vector <- log(pi_post_m) + d_ci[unit,1] * dnorm(log(t_i[unit]), mean = (piece1 + theta_s1_lp), sd = sqrt(sigma2_c_post_m[1]), log=TRUE) + (1 - d_ci[unit,1]) * pnorm(log(t_i[unit]), mean = (piece1 + theta_s1_lp), sd = sqrt(sigma2_c_post_m[1]), log=TRUE, lower.tail = FALSE) + 
                                d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (piece2 + theta_s2_lp), sd = sqrt(sigma2_c_post_m[2]), log=TRUE) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (piece2 + theta_s2_lp), sd = sqrt(sigma2_c_post_m[2]), log=TRUE, lower.tail = FALSE) + 
                                d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (piece3 + theta_s3_lp), sd = sqrt(sigma2_c_post_m[3]), log=TRUE) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (piece3 + theta_s3_lp), sd = sqrt(sigma2_c_post_m[3]), log=TRUE, lower.tail = FALSE)
  
  max_bayes <- max(pf_vector)
  bayes_alloc <- which(pf_vector==max_bayes)
  
  return(bayes_alloc)
}

Bayes_allocation <- sapply(c(1:sample_size), FUN = post_prob_fc)
Group_analysis <- cbind(X1_cov, d_ci, Bayes_allocation)
colnames(Group_analysis) <- c("AnnualPremium", "DJIA", "AccidDrider", "Gender", "PremFreq", "UWAge", "Surrender", "Death", "Other", "BayesGroup")

classc(100 * nrow(Group_analysis[Group_analysis$BayesGroup==1,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==2,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==3,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==4,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==5,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==6,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==7,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==8,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==9,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==10,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==11,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==12,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==13,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==14,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==15,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==16,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==17,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==18,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==19,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==20,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==21,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==22,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==23,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==24,]) / sample_size,
100 * nrow(Group_analysis[Group_analysis$BayesGroup==25,]) / sample_size)

## - Group selection Analysis
group_select <- 2
Gr_select <- Group_analysis[Group_analysis$BayesGroup==group_select,]

# mean(Gr_select$AnnualPremium)
# mean(Group_analysis$AnnualPremium)

100*nrow(Gr_select) / nrow(Group_analysis)

mean(Gr_select$AnnualPremium) * 526.58 + 560.88

100 * nrow(Gr_select[Gr_select$AccidDrider==1,]) / nrow(Gr_select)
# 100 * nrow(Group_analysis[Group_analysis$AccidDrider==1,]) / nrow(Group_analysis)

100 * nrow(Gr_select[Gr_select$PremFreq==1,]) / nrow(Gr_select)
# 100 * nrow(Group_analysis[Group_analysis$PremFreq==1,]) / nrow(Group_analysis)

100 * nrow(Gr_select[Gr_select$UWAge==1,]) / nrow(Gr_select)
# 100 * nrow(Group_analysis[Group_analysis$UWAge==1,]) / nrow(Group_analysis)

100 * nrow(Gr_select[Gr_select$Surrender==1,]) / nrow(Gr_select)
100 * nrow(Gr_select[Gr_select$Death==1,]) / nrow(Gr_select)
100 * nrow(Gr_select[Gr_select$Other==1,]) / nrow(Gr_select)
# 100 * nrow(Group_analysis[Group_analysis$Surrender==1,]) / nrow(Group_analysis)

c(theta_s1_lp[group_select], theta_s2_lp[group_select], theta_s3_lp[group_select])

cumulative_g <- 0
for(g in 1:K){
  Gr_select <- Group_analysis[Group_analysis$BayesGroup==g,]
  cumulative_g <- cumulative_g + 100*nrow(Gr_select) / nrow(Group_analysis)
  print(paste(g, 100*nrow(Gr_select) / nrow(Group_analysis), cumulative_g))
}



#=============- Cox - Snell residuals - ===============

CS_DPM <- rep(0, sample_size)

log_S <- function(unit){
  Xbeta1 <- beta1_thn[,1] * X1_cov[unit,1] + beta1_thn[,2] * X1_cov[unit,2] + beta1_thn[,3] * X1_cov[unit,3] + beta1_thn[,4] * X1_cov[unit,4] + beta1_thn[,5] * X1_cov[unit,5] + beta1_thn[,6] * X1_cov[unit,6]
  Xbeta2 <- beta2_thn[,1] * X2_cov[unit,1] + beta2_thn[,2] * X2_cov[unit,2] + beta2_thn[,3] * X2_cov[unit,3] + beta2_thn[,4] * X2_cov[unit,4]
  Xbeta3 <- beta3_thn[,1] * X3_cov[unit,1] + beta3_thn[,2] * X3_cov[unit,2] + beta3_thn[,3] * X3_cov[unit,3]
  
  S_t_i <- mean( pi_mw_thn[,1] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,1]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,2]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,2]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,3]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,3]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,2] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,4]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,5]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,5]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,6]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,6]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,3] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,7]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,8]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,8]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,9]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,9]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,4] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,10]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,11]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,11]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,12]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,12]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,5] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,13]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,14]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,14]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,15]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,15]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,6] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,16]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,17]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,17]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,18]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,18]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,7] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,19]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,20]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,20]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,21]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,21]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,8] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,22]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,23]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,23]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,24]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,24]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,9] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,25]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,26]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,26]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,27]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,27]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,10] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,28]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,29]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,29]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,30]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,30]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,11] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,31]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,32]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,32]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,33]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,33]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,12] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,34]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,35]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,35]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,36]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,36]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,13] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,37]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,38]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,38]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,39]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,39]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,14] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,40]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,41]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,41]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,42]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,42]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,15] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,43]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,44]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,44]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,45]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,45]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,16] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,46]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,47]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,47]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,48]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,48]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,17] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,49]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,50]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,50]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,51]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,51]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,18] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,52]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,53]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,53]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,54]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,54]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,19] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,55]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,56]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,56]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,57]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,57]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,20] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,58]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,59]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,59]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,60]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,60]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,21] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,61]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,62]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,62]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,63]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,63]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,22] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,64]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,65]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,65]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,66]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,66]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,23] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,67]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,68]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,68]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,69]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,69]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +   
                 pi_mw_thn[,24] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,70]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,71]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,71]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,72]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,72]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) +
                 pi_mw_thn[,25] * pnorm(log(t_i[unit]), mean = (Xbeta1 + theta_star_thn[,73]), sd=sqrt(sigma2_c_thn[,1]), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,74]), sd=sqrt(sigma2_c_thn[,2])) + (1 - d_ci[unit,2]) * pnorm(log(t_i[unit]), mean = (Xbeta2 + theta_star_thn[,74]), sd=sqrt(sigma2_c_thn[,2]), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,75]), sd=sqrt(sigma2_c_thn[,3])) + (1 - d_ci[unit,3]) * pnorm(log(t_i[unit]), mean = (Xbeta3 + theta_star_thn[,75]), sd=sqrt(sigma2_c_thn[,3]), lower.tail = FALSE)) )
return(-log(S_t_i))
}

i1 <- seq(1, 73, 3)
i2 <- seq(2, 74, 3)
i3 <- seq(3, 75, 3)


log_S <- function(unit){
  Xbeta1 <- beta1_thn[,1] * X1_cov[unit,1] + beta1_thn[,2] * X1_cov[unit,2] + beta1_thn[,3] * X1_cov[unit,3] + beta1_thn[,4] * X1_cov[unit,4] + beta1_thn[,5] * X1_cov[unit,5] + beta1_thn[,6] * X1_cov[unit,6]
  Xbeta2 <- beta2_thn[,1] * X2_cov[unit,1] + beta2_thn[,2] * X2_cov[unit,2] + beta2_thn[,3] * X2_cov[unit,3] + beta2_thn[,4] * X2_cov[unit,4]
  Xbeta3 <- beta3_thn[,1] * X3_cov[unit,1] + beta3_thn[,2] * X3_cov[unit,2] + beta3_thn[,3] * X3_cov[unit,3]
   
  S_t_i <- mean( pi_mw_thn[,1]  * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[1]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[1]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[1]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[1]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[1]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +   
                 pi_mw_thn[,2]  * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[2]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[2]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[2]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[2]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[2]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +   
                 pi_mw_thn[,3]  * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[3]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[3]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[3]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[3]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[3]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +   
                 pi_mw_thn[,4]  * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[4]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[4]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[4]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[4]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[4]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +   
                 pi_mw_thn[,5]  * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[5]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[5]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[5]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[5]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[5]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +    
                 pi_mw_thn[,6]  * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[6]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[6]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[6]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[6]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[6]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +   
                 pi_mw_thn[,7]  * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[7]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[7]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[7]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[7]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[7]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +   
                 pi_mw_thn[,8]  * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[8]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[8]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[8]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[8]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[8]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +     
                 pi_mw_thn[,9]  * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[9]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[9]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[9]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[9]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[9]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +   
                 pi_mw_thn[,10] * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[10]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[10]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[10]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[10]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[10]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +     
                 pi_mw_thn[,11] * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[11]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[11]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[11]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[11]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[11]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +   
                 pi_mw_thn[,12] * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[12]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[12]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[12]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[12]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[12]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +    
                 pi_mw_thn[,13] * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[13]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[13]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[13]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[13]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[13]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +   
                 pi_mw_thn[,14] * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[14]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[14]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[14]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[14]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[14]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +   
                 pi_mw_thn[,15] * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[15]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[15]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[15]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[15]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[15]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +     
                 pi_mw_thn[,16] * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[16]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[16]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[16]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[16]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[16]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +   
                 pi_mw_thn[,17] * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[17]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[17]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[17]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[17]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[17]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +   
                 pi_mw_thn[,18] * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[18]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[18]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[18]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[18]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[18]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +   
                 pi_mw_thn[,19] * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[19]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[19]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[19]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[19]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[19]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +     
                 pi_mw_thn[,20] * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[20]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[20]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[20]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[20]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[20]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +     
                 pi_mw_thn[,21] * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[21]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[21]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[21]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[21]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[21]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +     
                 pi_mw_thn[,22] * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[22]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[22]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[22]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[22]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[22]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +     
                 pi_mw_thn[,23] * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[23]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[23]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[23]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[23]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[23]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +      
                 pi_mw_thn[,24] * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[24]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[24]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[24]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[24]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[24]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE)) +   
                 pi_mw_thn[,25] * pnorm(((log(t_i[unit]) - Xbeta1 - theta_star_thn[,i1[25]]) / sqrt(sigma2_c_thn[,1])), lower.tail = FALSE) * (d_ci[unit,2] * dnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[25]]) / sqrt(sigma2_c_thn[,2]))) + (1 - d_ci[unit,2]) * pnorm(((log(t_i[unit]) - Xbeta2 - theta_star_thn[,i2[25]]) / sqrt(sigma2_c_thn[,2])), lower.tail = FALSE)) * (d_ci[unit,3] * dnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[25]]) / sqrt(sigma2_c_thn[,3]))) + (1 - d_ci[unit,3]) * pnorm(((log(t_i[unit]) - Xbeta3 - theta_star_thn[,i3[25]]) / sqrt(sigma2_c_thn[,3])), lower.tail = FALSE))  )
  return(-log(S_t_i))
}



CS_DPM <- sapply(c(1:sample_size), log_S)

r_CS_DPM <- CS_DPM + 1 - d_ci[,1]

dataset_CS <- cbind(r_CS_DPM, d_ci[,1])
dataset_CS <- as.data.frame(dataset_CS)
KM_DPM_cs <- survfit(Surv(r_CS_DPM, d_ci[,1]) ~ 1,  type="kaplan-meier", conf.type="log")

r_CS_DPM_ord <- sort(r_CS_DPM)
plot(r_CS_DPM_ord[1:21976], -log(KM_DPM_cs$surv), xlim=c(0,2), ylim=c(0,2))
abline(a=0, b=1)


## - Cox Proportional hazard model



# - Posterior summaries of the parameters
