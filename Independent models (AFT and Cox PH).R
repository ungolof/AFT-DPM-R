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



#================ - Simple AFT with independence - ===========================

neg_logL_surr <- function(vdParameters){
  beta <- vdParameters[1:R1]
  l_sigma_surr <- vdParameters[R1+1]
  
  logL <- sum( d_ci[,1] * dnorm(log(t_i), mean = X1_cov %*% beta, sd = exp(l_sigma_surr), log=TRUE) + (1 - d_ci[,1]) * pnorm(log(t_i), mean = X1_cov %*% beta, sd = exp(l_sigma_surr), log=TRUE, lower.tail = FALSE))
  return(-logL)
  
}

neg_logL_d <- function(vdParameters){
  beta <- vdParameters[1:R2]
  l_sigma <- vdParameters[R2+1]
  
  logL <- sum( d_ci[,2] * dnorm(log(t_i), mean = X2_cov %*% beta, sd = exp(l_sigma), log=TRUE) + (1 - d_ci[,2]) * pnorm(log(t_i), mean = X2_cov %*% beta, sd = exp(l_sigma), log=TRUE, lower.tail = FALSE))
  return(-logL)
}

neg_logL_o <- function(vdParameters){
  beta <- vdParameters[1:R3]
  l_sigma <- vdParameters[R3+1]
  
  logL <- sum( d_ci[,3] * dnorm(log(t_i), mean = X3_cov %*% beta, sd = exp(l_sigma), log=TRUE) + (1 - d_ci[,3]) * pnorm(log(t_i), mean = X3_cov %*% beta, sd = exp(l_sigma), log=TRUE, lower.tail = FALSE))
  return(-logL)
}

par_init <- c(1.54962008, -0.04107819, -0.64499319,  0.05340514, 0.03104960, 0.07439069, -0.10359116, -1.35261111)   # c(colMeans(mu_theta_thn)[1], colMeans(beta1_thn), log(sqrt((colMeans(sigma2_c_thn)[1]))))
AFT_surr <- nlm(neg_logL_surr, p=par_init, typsize=par_init, hessian=T, iterlim=1000)

par_init <- c(1.54962008, -0.04107819, -0.64499319,  0.05340514, 0.03104960, -1.35261111)   # c(colMeans(mu_theta_thn)[1], colMeans(beta1_thn), log(sqrt((colMeans(sigma2_c_thn)[1]))))
AFT_d <- nlm(neg_logL_d, p=par_init, typsize=par_init, hessian=T, iterlim=1000)

par_init <- c(1.54962008, -0.04107819, -0.64499319,  0.05340514, -1.35261111)   # c(colMeans(mu_theta_thn)[1], colMeans(beta1_thn), log(sqrt((colMeans(sigma2_c_thn)[1]))))
AFT_o <- nlm(neg_logL_o, p=par_init, typsize=par_init, hessian=T, iterlim=1000)

# - AIC computation

AIC_AFT <- - 2 * (- AFT_surr$minimum) + 2 * (R1 + 1)

# - Using the lognormal distribution
neg_logL_surr_ln <- function(vdParameters){
  beta <- vdParameters[1:R1]
  l_sigma_surr <- vdParameters[R1+1]
  
  logL <- sum( d_ci[,1] * dlnorm(t_i, meanlog = X1_cov %*% beta, sdlog = exp(l_sigma_surr), log=TRUE) + (1 - d_ci[,1]) * plnorm(t_i, meanlog = X1_cov %*% beta, sdlog = exp(l_sigma_surr), log=TRUE, lower.tail = FALSE))
  return(-logL)
  
}

AIC_AFT_ln <- - 2 * (- neg_logL_surr_ln(AFT_surr$estimate)) + 2 * (R1 + 1)

# - Prediction

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


integrand_AFT <- function(t){
  dnorm(t, mean=sum(AFT_surr$estimate[1:R1] * X1_cov[i,]), sd=exp(AFT_surr$estimate[R1+1])) * pnorm(t, mean=sum(AFT_d$estimate[1:R2] * X2_cov[i,]), sd=exp(AFT_d$estimate[R2+1]), lower.tail = FALSE) * pnorm(t, mean=sum(AFT_o$estimate[1:R3] * X3_cov[i,]), sd=exp(AFT_o$estimate[R3+1]), lower.tail = FALSE)
}

F_t_table_AFT <- matrix(NA, sample_size_test, n_quarters) # - takes F(d_q; X)

rates_tb_qt_AFT <- function(j){
  integrate(integrand_AFT, lower=-Inf, upper=log(j))$value
}

for(i in 1:sample_size_test){ # - Calculate the F in formula (9)
  F_t_table_AFT[i,] <- sapply(1:n_quarters, rates_tb_qt_AFT)
}

# \hat{r}_{t,i} - individually predicted rates for computation of eq. (10) in Milhaud-Dutang
r_pred_table_AFT <- matrix(NA, sample_size_test, n_quarters)
r_pred_table_AFT[,1] <- F_t_table_AFT[,1]
r_pred_table_AFT[,2:n_quarters] <- (F_t_table_AFT[,2:n_quarters] - F_t_table_AFT[,1:(n_quarters-1)]) / (1 - F_t_table_AFT[,1:(n_quarters-1)])

r_t_hat_AFT <- colSums(r_pred_table_AFT * risk_set_ind) / R_tq

lines(c(1:n_quarters), r_t_hat_AFT, type='l', lty=4, col='red')
R_RMSEq_AFT <- rep(0, n_quarters)



#===================== - Cox Proportional hazard model - =================================

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


## 1 - Fit a Cox proportional hazard model
Surr_Cox <- coxph(Surv(uslapseagent$duration, uslapseagent$surrender) ~ uslapseagent$annual.premium + uslapseagent$termination.cause + ifelse(uslapseagent$acc.death.rider=="Rider", 1, 0) + ifelse(uslapseagent$gender=="Female", 1, 0) + ifelse(uslapseagent$premium.frequency %in% c("Annual", 'Other'), 1, 0) + ifelse(uslapseagent$underwriting.age %in% c('Young',"Middle"), 1, 0), data=uslapseagent)
Death_Cox <- coxph(Surv(uslapseagent$duration, uslapseagent$death) ~ ifelse(uslapseagent$gender=="Female", 1, 0) + ifelse(uslapseagent$underwriting.age=="Old", 1, 0) + ifelse(uslapseagent$living.place=="Other", 1, 0) + ifelse(uslapseagent$risk.state=="Smoker", 1, 0), data=uslapseagent)
Oth_Cox <- coxph(Surv(uslapseagent$duration, uslapseagent$other) ~ uslapseagent$annual.premium + ifelse(uslapseagent$premium.frequency %in% c("Annual", 'Other'), 1, 0) + ifelse(uslapseagent$acc.death.rider=="Rider", 1, 0), data=uslapseagent)

## 2 - Get a cumulative hazard function
bh_est_sur <- basehaz(Surr_Cox)
bh_est_death <- basehaz(Death_Cox)
bh_est_oth <- basehaz(Oth_Cox)

### 2.1 - Fit a linear model between cumulative hazard and time
model_s <- lm(bh_est_sur$hazard ~ bh_est_sur$time, data = bh_est_sur)
model_d <- lm(bh_est_death$hazard ~ bh_est_death$time, data = bh_est_death)
model_o <- lm(bh_est_oth$hazard ~ bh_est_oth$time, data = bh_est_oth)

### 2.2 - Calculate AIC
log_lik_CoxPH <- rep(0, sample_size)
for(i in 1:sample_size){
  log_lik_CoxPH[i] <- -uslapseagent$duration[i] * model_s$coefficients[2] * exp(sum(Surr_Cox$coefficients * X1_cov[i,])) + uslapseagent$surrender[i] * (log(model_s$coefficients[2]) + sum(Surr_Cox$coefficients * X1_cov[i,]))
}


AIC_is_CoxPH <- - 2 * sum(log_lik_CoxPH) + 2*7

# - Prediction of rates
R1 <- 6
R2 <- 4
R3 <- 3

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


F_t_cox_table <- F_t_table

for(i in 1:sample_size_test){
  for(j in 1:n_quarters){
    F_t_cox_table[i,j] <- (1 - exp(-j * (model_s$coefficients[2] * exp(sum(Surr_Cox$coefficients * X1_cov[i,1:R1])) + model_d$coefficients[2] * exp(sum(Death_Cox$coefficients * X2_cov[i,1:R2])) + model_o$coefficients[2] * exp(sum(Oth_Cox$coefficients * X3_cov[i,1:R3]))))) * (model_s$coefficients[2] * exp(sum(Surr_Cox$coefficients * X1_cov[i,1:R1]))) / (model_s$coefficients[2] * exp(sum(Surr_Cox$coefficients * X1_cov[i,1:R1])) + model_d$coefficients[2] * exp(sum(Death_Cox$coefficients * X2_cov[i,1:R2])) + model_o$coefficients[2] * exp(sum(Oth_Cox$coefficients * X3_cov[i,1:R3])))
  }
}


r_pred_table_cox <- matrix(NA, sample_size_test, n_quarters)
r_pred_table_cox[,1] <- F_t_cox_table[,1]
r_pred_table_cox[,2:n_quarters] <- (F_t_cox_table[,2:n_quarters] - F_t_cox_table[,1:(n_quarters-1)]) / (1 - F_t_cox_table[,1:(n_quarters-1)])

r_t_hat_cox <- colSums(r_pred_table_cox * risk_set_ind) / R_tq

