#=========================== - Data preparation for the test set - ==============================

sample_size_test <- nrow(uslapseagent_test)

#===================== - Final dataset for inference - =====================
R1 <- 6
R2 <- 4
R3 <- 3

# - Covariates storage

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

# - Parameters for analysis - Take the posterior mean
beta1_fin <- colMeans(beta1_thn)
beta2_fin <- colMeans(beta2_thn)
beta3_fin <- colMeans(beta3_thn)
sigma2_1 <- mean(sigma2_c_thn[,1])
sigma2_2 <- mean(sigma2_c_thn[,2])
sigma2_3 <- mean(sigma2_c_thn[,3])
theta_star_fin <- colMeans(theta_star_thn)
pi_mw_m <- colMeans(pi_mw_thn)


# - Support function

integrand <- function(t){
    pi_mw_m[1]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[1] ), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[2] ), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[3] ), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[2]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[4] ), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[5] ), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[6] ), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[3]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[7] ), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[8] ), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[9] ), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[4]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[10]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[11]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[12]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[5]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[13]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[14]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[15]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[6]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[16]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[17]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[18]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[7]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[19]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[20]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[21]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[8]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[22]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[23]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[24]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[9]  * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[25]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[26]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[27]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[10] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[28]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[29]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[30]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[11] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[31]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[32]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[33]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[12] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[34]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[35]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[36]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[13] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[37]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[38]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[39]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[14] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[40]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[41]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[42]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[15] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[43]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[44]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[45]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[16] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[46]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[47]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[48]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[17] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[49]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[50]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[51]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[18] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[52]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[53]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[54]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[19] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[55]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[56]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[57]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[20] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[58]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[59]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[60]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[21] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[61]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[62]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[63]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[22] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[64]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[65]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[66]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[23] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[67]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[68]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[69]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[24] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[70]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[71]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[72]), sd=sqrt(sigma2_3), lower.tail = FALSE) +
    pi_mw_m[25] * dnorm(t, mean=(sum(beta1_fin * X1_cov[i,]) + theta_star_fin[73]), sd=sqrt(sigma2_1)) * pnorm(t, mean=(sum(beta2_fin * X2_cov[i,]) + theta_star_fin[74]), sd=sqrt(sigma2_2), lower.tail = FALSE) * pnorm(t, mean=(sum(beta3_fin * X3_cov[i,]) + theta_star_fin[75]), sd=sqrt(sigma2_3), lower.tail = FALSE)

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

F_t_table <- matrix(NA, sample_size_test, n_quarters) # - takes F(d_q; X)

rates_tb_qt <- function(j){
  integrate(integrand, lower=-Inf, upper=log(j))$value
}

for(i in 1:sample_size_test){ # - Calculate the F in formula (9)
  F_t_table[i,] <- sapply(1:n_quarters, rates_tb_qt)
}

# \hat{r}_{t,i} - individually predicted rates for computation of eq. (10) in Milhaud-Dutang
r_pred_table <- matrix(NA, sample_size_test, n_quarters)
r_pred_table[,1] <- F_t_table[,1]
r_pred_table[,2:n_quarters] <- (F_t_table[,2:n_quarters] - F_t_table[,1:(n_quarters-1)]) / (1 - F_t_table[,1:(n_quarters-1)])

# - individual being at risk in quarter t. First column corresponds to individuals at risk at time 0, 2nd col corresponds to individual
# at the end of the first quarter, and so on
risk_set_ind <- matrix(NA, sample_size_test, n_quarters)
risk_set_ind[,1] <- 1
for(i in 1:sample_size_test){
  for(j in 1:(n_quarters-1)){
    risk_set_ind[i,j+1] <- ifelse(t_i[i]>j,1,0)
  }
}

# - \hat{r}_{t}
r_t_hat <- colSums(r_pred_table * risk_set_ind) / R_tq
## - reformulated \hat{r}_{t}
### r_t_hat <- colSums(r_pred_table) / R_tq # - We shouldn't multiplu r_pred_table by the risk set

plot(c(1:n_quarters), r_t_hat, type='l', ylim=c(0,0.06), lwd=2, xlab="Quarter")
# - plot(c(1:48), r_t_hat[1:48], type='l', ylim=c(0,0.06), lwd=2)

# - Calculation of the indicator variable if $T$ occurs in the quarter
t_i_between <- matrix(0, sample_size_test, n_quarters)

# - Empirical
# t_i between quarters

t_i_between[,1] <- ifelse((t_i<=1) & (d_ci[,1]==1),1,0)
for(i in 2:n_quarters){
  for(j in 1:sample_size_test){
    t_i_between[j,i] <- ifelse(((t_i[j]>(i-1)) & (t_i[j]<= i) & (d_ci[j,1]==1)),1,0)
  }
}

## -  Calculation of the empirical \hat{r}_t
r_t_emp <- colSums(t_i_between) / R_tq

lines(c(1:n_quarters), r_t_emp, lty=2, lwd=2)

## - Lower and upper CI with the sampling distribution of the binomial
sd_bin <- sqrt(colSums(r_pred_table * risk_set_ind) * (R_tq - colSums(r_pred_table * risk_set_ind)) / (R_tq^3))

lines(c(1:60), r_t_hat[1:60] + 2 * sd_bin[1:60], lty=1)
lines(c(1:60), r_t_hat[1:60] - 2 * sd_bin[1:60], lty=1)



# - Prediction using Cox (assume piecewise constant hazard function between two dates)

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

#plot(bh_est_sur$time, bh_est_sur$hazard, type='p', pch=16, cex=0.1, ylim=c(0,0.7))
#abline(a=model_s$coefficients[1], b=model_s$coefficients[2]) # - Check whether the linear fit is reasonable

#plot(bh_est_death$time, bh_est_death$hazard, type='p', pch=16, cex=0.1, ylim=c(0,0.1))
#abline(a=model_d$coefficients[1], b=model_d$coefficients[2]) # - Check whether the linear fit is reasonable

#plot(bh_est_oth$time, bh_est_oth$hazard, type='p', pch=16, cex=0.1, ylim=c(0,0.2))
#abline(a=model_o$coefficients[1], b=model_o$coefficients[2]) # - Check whether the linear fit is reasonable

F_t_cox_table <- F_t_table

for(i in 1:sample_size_test){
  for(j in 1:n_quarters){
    F_t_cox_table[i,j] <- (1 - exp(-j * model_s$coefficients[2] * exp(sum(Surr_Cox$coefficients * X1_cov[i,1:R1])))) * exp(-j * (model_d$coefficients[2] * exp(sum(Death_Cox$coefficients * X2_cov[i,1:R2])) + model_o$coefficients[2] * exp(sum(Oth_Cox$coefficients * X3_cov[i,1:R3]))))  #(model_s$coefficients[2] * exp(sum(Surr_Cox$coefficients * X1_cov[i,2:R1]))) * (1 - exp(-j * (model_s$coefficients[2] * exp(sum(Surr_Cox$coefficients * X1_cov[i,2:R1])) + model_d$coefficients[2] * exp(sum(Death_Cox$coefficients * X2_cov[i,2:R2])) + model_o$coefficients[2] * exp(sum(Oth_Cox$coefficients * X3_cov[i,2:R3]))))) #### 1 - exp(- j * model$coefficients[2] * exp(sum(Surr_Cox$coefficients * X1_cov[i,])))
  }
}

r_pred_table_cox <- matrix(NA, sample_size_test, n_quarters)
r_pred_table_cox[,1] <- F_t_cox_table[,1]
r_pred_table_cox[,2:n_quarters] <- (F_t_cox_table[,2:n_quarters] - F_t_cox_table[,1:(n_quarters-1)]) / (1 - F_t_cox_table[,1:(n_quarters-1)])

r_t_hat_cox <- colSums(r_pred_table_cox * risk_set_ind) / R_tq

lines(c(1:n_quarters), r_t_hat_cox, type='l', lty=3)

# - Prediction using Subdistribution (assume piecewise constant hazard function between two dates)
uslapseagent$durationsbd <- ifelse(uslapseagent$death + uslapseagent$other==1, 70, uslapseagent$duration)
uslapseagent_test$durationsbd <- ifelse(uslapseagent_test$death + uslapseagent_test$other==1, 70, uslapseagent_test$duration)

t_i_sbd <- uslapseagent_test$durationsbd
risk_set_ind_sbd <- matrix(NA, sample_size_test, n_quarters)
risk_set_ind_sbd[,1] <- 1
for(i in 1:sample_size_test){
  for(j in 2:n_quarters){
    risk_set_ind_sbd[i,j] <- ifelse(t_i_sbd[i]>j,1,0)
  }
}

r_tq_sbd <- rep(0, n_quarters)

r_tq_sbd[1] <- sample_size_test
for(q in 2:n_quarters){
  r_tq_sbd[q] <- sum(ifelse(t_i_sbd>q,1,0))
}

Surr_sbd <- coxph(Surv(uslapseagent$durationsbd, uslapseagent$surrender) ~ uslapseagent$annual.premium + uslapseagent$termination.cause + ifelse(uslapseagent$acc.death.rider=="Rider", 1, 0) + ifelse(uslapseagent$gender=="Female", 1, 0) + ifelse(uslapseagent$premium.frequency %in% c("Annual", 'Other'), 1, 0) + ifelse(uslapseagent$underwriting.age %in% c('Young',"Middle"), 1, 0), data=uslapseagent)

Surv_pred_sbd <- basehaz(Surr_sbd, centered=TRUE)

model <- lm(Surv_pred_sbd$hazard ~ Surv_pred_sbd$time + sqrt(Surv_pred_sbd$time), data = Surv_pred_sbd)

#plot(Surv_pred_sbd$time, Surv_pred_sbd$hazard, type='p', pch=16, cex=0.1)
#lines(Surv_pred_sbd$time, model$coefficients[1] + model$coefficients[2] * Surv_pred_sbd$time + model$coefficients[3] * sqrt(Surv_pred_sbd$time))

integrand_sbd <- function(t){
  exp(sum(Surr_sbd$coefficients * X1_cov[i,1:R1])) * exp(-exp(sum(Surr_sbd$coefficients * X1_cov[i,1:R1])) * (model$coefficients[2] * t + model$coefficients[3] * sqrt(t))) * (model$coefficients[2] + 0.5 * model$coefficients[3] / sqrt(t))
}

rates_tb_qt_sbd <- function(j){
  integrate(integrand_sbd, lower=0, upper=j)$value
}

F_t_sbd_table <- F_t_table

for(i in 1:sample_size_test){
  F_t_sbd_table[i,] <- sapply(1:n_quarters, rates_tb_qt_sbd)
}

r_pred_table_sbd <- matrix(NA, sample_size_test, n_quarters)
r_pred_table_sbd[,1] <- F_t_sbd_table[,1]
r_pred_table_sbd[,2:n_quarters] <- (F_t_sbd_table[,2:n_quarters] - F_t_sbd_table[,1:(n_quarters-1)]) / (1 - F_t_sbd_table[,1:(n_quarters-1)])

r_t_hat_sbd <- colSums(r_pred_table_sbd * risk_set_ind_sbd) / r_tq_sbd

lines(c(1:n_quarters), r_t_hat_sbd, type='l', lty=4)

par(mfrow=c(1,1), mar = c(4.5,4,1,2.5))
plot(c(1:60), r_t_hat[1:60], type='l', ylim=c(0,0.06), lwd=2, main="Out-of-sample prediction", ylab="Surrending rates", xlab="Quarter")
lines(c(1:60), r_t_emp[1:60], lty=2, lwd=2)
lines(c(1:60), r_t_hat_cox[1:60], lty=3, lwd=1)
lines(c(1:60), r_t_hat_sbd[1:60], lty=4, lwd=1)
lines(c(1:60), r_t_hat[1:60] + 2 * sd_bin[1:60], lty=1)
lines(c(1:60), r_t_hat[1:60] - 2 * sd_bin[1:60], lty=1)

legend("topright", legend=c("DPM", "Empirical", "Cox PH", 'Subdistribution'), lty=c(1,2,3,4), lwd=c(2,2,1,1), cex=0.75, bty='n')

par(mfrow=c(1,1), mar = c(4.5,4,1,2.5))

par(mfrow=c(1,2), mar = c(4.5,4,1,2.5))
plot(c(1:60), r_t_hat[1:60], type='l', ylim=c(0,0.06), lwd=2, ylab="Surrending rates r_t_hat", xlab="Quarter")
lines(c(1:60), r_t_emp[1:60], lty=2, lwd=2)
lines(c(1:60), r_t_hat_cox[1:60], lty=3, lwd=1)
lines(c(1:60), r_t_hat_sbd[1:60], lty=4, lwd=1)
lines(c(1:60), r_t_hat[1:60] + 2 * sd_bin[1:60], lty=1)
lines(c(1:60), r_t_hat[1:60] - 2 * sd_bin[1:60], lty=1)

legend("topright", legend=c("DPM", "DPM 95% CI", "Empirical", "Cox PH", 'Subdistribution'), lty=c(1,1,2,3,4), lwd=c(2,1,2,1,1), cex=0.75, bty='n')

plot(c(1:60), R_RMSEq_DPM[1:60], type='l', lwd=2, ylim=c(0,0.03), ylab="Rolling RMSE", xlab="Quarter")
lines(c(1:60), R_RMSEq_Cox[1:60], type='l', lty=3)
lines(c(1:60), R_RMSEq_Sbd[1:60], type='l', lty=4)
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=0.75, bty='n')



# - Calculation of the Rolling RMSE
R_RMSEq_DPM <- R_RMSEq_Cox <- R_RMSEq_Sbd <- rep(0, n_quarters)

for(q in 1:n_quarters){
  R_RMSEq_DPM[q] <- sqrt(sum((r_t_hat[1:q]-r_t_emp[1:q])^2)/q)
  R_RMSEq_Cox[q] <- sqrt(sum((r_t_hat_cox[1:q]-r_t_emp[1:q])^2)/q)
  R_RMSEq_Sbd[q] <- sqrt(sum((r_t_hat_sbd[1:q]-r_t_emp[1:q])^2)/q)
}

#par(mfrow=c(1,1), mar = c(4.5,4,1,2.5))
plot(c(1:n_quarters), R_RMSEq_DPM, type='l', lwd=2, ylim=c(0,0.03), main="Rolling RMSE", ylab="Rolling RMSE", xlab="Quarter")
lines(c(1:n_quarters), R_RMSEq_Cox, type='l', lty=3)
lines(c(1:n_quarters), R_RMSEq_Sbd, type='l', lty=4)
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=0.75, bty='n')

# - Prediction by gender
## - Males

R_tq_males <- rep(0, n_quarters)
male_index <- which(uslapseagent_test$gender=="Male")
male_sample_size <- length(male_index)
R_tq_males[1] <- nrow(uslapseagent_test[uslapseagent_test$gender=="Male",])
for(q in 2:n_quarters){
  R_tq_males[q] <- sum(ifelse(t_i[male_index]>q,1,0))
}


# \hat{r}_{t,i} - individually predicted rates for computation of eq. (10) in Milhaud-Dutang
r_pred_table_males <- matrix(NA, nrow(uslapseagent_test[uslapseagent_test$gender=="Male",]), n_quarters)
r_pred_table_males[,1] <- F_t_table[male_index,1]
r_pred_table_males[,2:n_quarters] <- (F_t_table[male_index,2:n_quarters] - F_t_table[male_index,1:(n_quarters-1)]) / (1 - F_t_table[male_index,1:(n_quarters-1)])

# - individual being at risk in quarter t. First column corresponds to individuals at risk at time 0, 2nd col corresponds to individual
# at the end of the first quarter, and so on
risk_set_ind_males <- matrix(NA, nrow(uslapseagent_test[uslapseagent_test$gender=="Male",]), n_quarters)
risk_set_ind_males[,1] <- 1

counter <- 1 
for(i in male_index){
  for(j in 1:(n_quarters-1)){
    risk_set_ind_males[counter,j+1] <- ifelse(t_i[i]>j,1,0)
  }
  counter <- counter + 1
}

# - \hat{r}_{t}
r_t_hat_males <- colSums(r_pred_table_males * risk_set_ind_males) / R_tq_males
## - reformulated \hat{r}_{t}
### r_t_hat <- colSums(r_pred_table) / R_tq # - We shouldn't multiplu r_pred_table by the risk set

plot(c(1:n_quarters), r_t_hat_males, type='l', ylim=c(0,0.06), lwd=2)

# - Calculation of the indicator variable if $T$ occurs in the quarter
t_i_between_males <- matrix(0, nrow(uslapseagent_test[uslapseagent_test$gender=="Male",]), n_quarters)

# - Empirical
# t_i between quarters

t_i_between_males[,1] <- ifelse((t_i[male_index]<=1) & (d_ci[male_index,1]==1),1,0)

for(i in 2:n_quarters){
  counter <- 1
  for(j in male_index){
    t_i_between_males[counter,i] <- ifelse(((t_i[j]>(i-1)) & (t_i[j]<= i) & (d_ci[j,1]==1)),1,0)
    counter <- counter + 1
  }
}

## -  Calculation of the empirical \hat{r}_t
r_t_emp_males <- colSums(t_i_between_males) / R_tq_males

## - Comparison with Cox PH
r_pred_table_cox_males <- matrix(NA, male_sample_size, n_quarters)
r_pred_table_cox_males[,1] <- F_t_cox_table[male_index,1]
r_pred_table_cox_males[,2:n_quarters] <- (F_t_cox_table[male_index,2:n_quarters] - F_t_cox_table[male_index,1:(n_quarters-1)]) / (1 - F_t_cox_table[male_index,1:(n_quarters-1)])

r_t_hat_cox_males <- colSums(r_pred_table_cox_males * risk_set_ind_males) / R_tq_males

## - Comparison with Subdistribution
r_pred_table_sbd_males <- matrix(NA, male_sample_size, n_quarters)
r_pred_table_sbd_males[,1] <- F_t_sbd_table[male_index,1]
r_pred_table_sbd_males[,2:n_quarters] <- (F_t_sbd_table[male_index,2:n_quarters] - F_t_sbd_table[male_index,1:(n_quarters-1)]) / (1 - F_t_sbd_table[male_index,1:(n_quarters-1)])

r_t_hat_sbd_males <- colSums(r_pred_table_sbd_males * risk_set_ind_males) / R_tq_males

## - females

R_tq_females <- rep(0, n_quarters)
female_index <- which(uslapseagent_test$gender=="Female")
female_sample_size <- length(female_index)
R_tq_females[1] <- nrow(uslapseagent_test[uslapseagent_test$gender=="Female",])
for(q in 2:n_quarters){
  R_tq_females[q] <- sum(ifelse(t_i[female_index]>q,1,0))
}


# \hat{r}_{t,i} - individually predicted rates for computation of eq. (10) in Milhaud-Dutang
r_pred_table_females <- matrix(NA, nrow(uslapseagent_test[uslapseagent_test$gender=="Female",]), n_quarters)
r_pred_table_females[,1] <- F_t_table[female_index,1]
r_pred_table_females[,2:n_quarters] <- (F_t_table[female_index,2:n_quarters] - F_t_table[female_index,1:(n_quarters-1)]) / (1 - F_t_table[female_index,1:(n_quarters-1)])

# - individual being at risk in quarter t. First column corresponds to individuals at risk at time 0, 2nd col corresponds to individual
# at the end of the first quarter, and so on
risk_set_ind_females <- matrix(NA, nrow(uslapseagent_test[uslapseagent_test$gender=="Female",]), n_quarters)
risk_set_ind_females[,1] <- 1

counter <- 1 
for(i in female_index){
  for(j in 1:(n_quarters-1)){
    risk_set_ind_females[counter,j+1] <- ifelse(t_i[i]>j,1,0)
  }
  counter <- counter + 1
}

# - \hat{r}_{t}
r_t_hat_females <- colSums(r_pred_table_females * risk_set_ind_females) / R_tq_females
## - reformulated \hat{r}_{t}
### r_t_hat <- colSums(r_pred_table) / R_tq # - We shouldn't multiplu r_pred_table by the risk set

# - Calculation of the indicator variable if $T$ occurs in the quarter
t_i_between_females <- matrix(0, nrow(uslapseagent_test[uslapseagent_test$gender=="Female",]), n_quarters)

# - Empirical
# t_i between quarters

t_i_between_females[,1] <- ifelse((t_i[female_index]<=1) & (d_ci[female_index,1]==1),1,0)

for(i in 2:n_quarters){
  counter <- 1
  for(j in female_index){
    t_i_between_females[counter,i] <- ifelse(((t_i[j]>(i-1)) & (t_i[j]<= i) & (d_ci[j,1]==1)),1,0)
    counter <- counter + 1
  }
}

## -  Calculation of the empirical \hat{r}_t
r_t_emp_females <- colSums(t_i_between_females) / R_tq_females

## - Comparison with Cox PH
r_pred_table_cox_females <- matrix(NA, female_sample_size, n_quarters)
r_pred_table_cox_females[,1] <- F_t_cox_table[female_index,1]
r_pred_table_cox_females[,2:n_quarters] <- (F_t_cox_table[female_index,2:n_quarters] - F_t_cox_table[female_index,1:(n_quarters-1)]) / (1 - F_t_cox_table[female_index,1:(n_quarters-1)])

r_t_hat_cox_females <- colSums(r_pred_table_cox_females * risk_set_ind_females) / R_tq_females


## - Comparison with Subdistribution
r_pred_table_sbd_females <- matrix(NA, female_sample_size, n_quarters)
r_pred_table_sbd_females[,1] <- F_t_sbd_table[female_index,1]
r_pred_table_sbd_females[,2:n_quarters] <- (F_t_sbd_table[female_index,2:n_quarters] - F_t_sbd_table[female_index,1:(n_quarters-1)]) / (1 - F_t_sbd_table[female_index,1:(n_quarters-1)])

r_t_hat_sbd_females <- colSums(r_pred_table_sbd_females * risk_set_ind_females) / R_tq_females



## - Evaluation by payment frequency

## - InfraAnnual

R_tq_infra <- rep(0, n_quarters)
infra_index <- which(uslapseagent_test$premium.frequency=="InfraAnnual")
infra_sample_size <- length(infra_index)
R_tq_infra[1] <- nrow(uslapseagent_test[uslapseagent_test$premium.frequency=="InfraAnnual",])
for(q in 2:n_quarters){
  R_tq_infra[q] <- sum(ifelse(t_i[infra_index]>q-1,1,0))
}


# \hat{r}_{t,i} - individually predicted rates for computation of eq. (10) in Milhaud-Dutang
r_pred_table_infra <- matrix(NA, nrow(uslapseagent_test[uslapseagent_test$premium.frequency=="InfraAnnual",]), n_quarters)
r_pred_table_infra[,1] <- F_t_table[infra_index,1]
r_pred_table_infra[,2:n_quarters] <- (F_t_table[infra_index,2:n_quarters] - F_t_table[infra_index,1:(n_quarters-1)]) / (1 - F_t_table[infra_index,1:(n_quarters-1)])

# - individual being at risk in quarter t. First column corresponds to individuals at risk at time 0, 2nd col corresponds to individual
# at the end of the first quarter, and so on
risk_set_ind_infra <- matrix(NA, nrow(uslapseagent_test[uslapseagent_test$premium.frequency=="InfraAnnual",]), n_quarters)
risk_set_ind_infra[,1] <- 1

counter <- 1 
for(i in infra_index){
  for(j in 1:(n_quarters-1)){
    risk_set_ind_infra[counter,j+1] <- ifelse(t_i[i]>j,1,0)
  }
  counter <- counter + 1
}

# - \hat{r}_{t}
r_t_hat_infra <- colSums(r_pred_table_infra * risk_set_ind_infra) / R_tq_infra
## - reformulated \hat{r}_{t}
### r_t_hat <- colSums(r_pred_table) / R_tq # - We shouldn't multiplu r_pred_table by the risk set

# - Calculation of the indicator variable if $T$ occurs in the quarter
t_i_between_infra <- matrix(0, nrow(uslapseagent_test[uslapseagent_test$premium.frequency=="InfraAnnual",]), n_quarters)

# - Empirical
# t_i between quarters

t_i_between_infra[,1] <- ifelse((t_i[infra_index]<=1) & (d_ci[infra_index,1]==1),1,0)

for(i in 2:n_quarters){
  counter <- 1
  for(j in infra_index){
    t_i_between_infra[counter,i] <- ifelse(((t_i[j]>(i-1)) & (t_i[j]<= i) & (d_ci[j,1]==1)),1,0)
    counter <- counter + 1
  }
}

## -  Calculation of the empirical \hat{r}_t
r_t_emp_infra <- colSums(t_i_between_infra) / R_tq_infra

## - Comparison with Cox PH
r_pred_table_cox_infra <- matrix(NA, infra_sample_size, n_quarters)
r_pred_table_cox_infra[,1] <- F_t_cox_table[infra_index,1]
r_pred_table_cox_infra[,2:n_quarters] <- (F_t_cox_table[infra_index,2:n_quarters] - F_t_cox_table[infra_index,1:(n_quarters-1)]) / (1 - F_t_cox_table[infra_index,1:(n_quarters-1)])

r_t_hat_cox_infra <- colSums(r_pred_table_cox_infra * risk_set_ind_infra) / R_tq_infra

## - Comparison with Subdistribution
r_pred_table_sbd_infra <- matrix(NA, infra_sample_size, n_quarters)
r_pred_table_sbd_infra[,1] <- F_t_sbd_table[infra_index,1]
r_pred_table_sbd_infra[,2:n_quarters] <- (F_t_sbd_table[infra_index,2:n_quarters] - F_t_sbd_table[infra_index,1:(n_quarters-1)]) / (1 - F_t_sbd_table[infra_index,1:(n_quarters-1)])

r_t_hat_sbd_infra <- colSums(r_pred_table_sbd_infra * risk_set_ind_infra) / R_tq_infra


## - Others

R_tq_other <- rep(0, n_quarters)
other_index <- which(uslapseagent_test$premium.frequency!="InfraAnnual")
other_sample_size <- length(other_index)
R_tq_other[1] <- nrow(uslapseagent_test[uslapseagent_test$premium.frequency!="InfraAnnual",])
for(q in 2:n_quarters){
  R_tq_other[q] <- sum(ifelse(t_i[other_index]>q-1,1,0))
}


# \hat{r}_{t,i} - individually predicted rates for computation of eq. (10) in Milhaud-Dutang
r_pred_table_other <- matrix(NA, nrow(uslapseagent_test[uslapseagent_test$premium.frequency!="InfraAnnual",]), n_quarters)
r_pred_table_other[,1] <- F_t_table[other_index,1]
r_pred_table_other[,2:n_quarters] <- (F_t_table[other_index,2:n_quarters] - F_t_table[other_index,1:(n_quarters-1)]) / (1 - F_t_table[other_index,1:(n_quarters-1)])

# - individual being at risk in quarter t. First column corresponds to individuals at risk at time 0, 2nd col corresponds to individual
# at the end of the first quarter, and so on
risk_set_ind_other <- matrix(NA, nrow(uslapseagent_test[uslapseagent_test$premium.frequency!="InfraAnnual",]), n_quarters)
risk_set_ind_other[,1] <- 1

counter <- 1 
for(i in other_index){
  for(j in 1:(n_quarters-1)){
    risk_set_ind_other[counter,j+1] <- ifelse(t_i[i]>j,1,0)
  }
  counter <- counter + 1
}

# - \hat{r}_{t}
r_t_hat_other <- colSums(r_pred_table_other * risk_set_ind_other) / R_tq_other
## - reformulated \hat{r}_{t}
### r_t_hat <- colSums(r_pred_table) / R_tq # - We shouldn't multiplu r_pred_table by the risk set

plot(c(1:n_quarters), r_t_hat_other, type='l', ylim=c(0,0.06), lwd=2)

# - Calculation of the indicator variable if $T$ occurs in the quarter
t_i_between_other <- matrix(0, nrow(uslapseagent_test[uslapseagent_test$premium.frequency!="InfraAnnual",]), n_quarters)

# - Empirical
# t_i between quarters

t_i_between_other[,1] <- ifelse((t_i[other_index]<=1) & (d_ci[other_index,1]==1),1,0)

for(i in 2:n_quarters){
  counter <- 1
  for(j in other_index){
    t_i_between_other[counter,i] <- ifelse(((t_i[j]>(i-1)) & (t_i[j]<= i) & (d_ci[j,1]==1)),1,0)
    counter <- counter + 1
  }
}

## -  Calculation of the empirical \hat{r}_t
r_t_emp_other <- colSums(t_i_between_other) / R_tq_other

## - Comparison with Cox PH
r_pred_table_cox_other <- matrix(NA, other_sample_size, n_quarters)
r_pred_table_cox_other[,1] <- F_t_cox_table[other_index,1]
r_pred_table_cox_other[,2:n_quarters] <- (F_t_cox_table[other_index,2:n_quarters] - F_t_cox_table[other_index,1:(n_quarters-1)]) / (1 - F_t_cox_table[other_index,1:(n_quarters-1)])

r_t_hat_cox_other <- colSums(r_pred_table_cox_other * risk_set_ind_other) / R_tq_other

## - Comparison with Subdistribution
r_pred_table_sbd_other <- matrix(NA, other_sample_size, n_quarters)
r_pred_table_sbd_other[,1] <- F_t_sbd_table[other_index,1]
r_pred_table_sbd_other[,2:n_quarters] <- (F_t_sbd_table[other_index,2:n_quarters] - F_t_sbd_table[other_index,1:(n_quarters-1)]) / (1 - F_t_sbd_table[other_index,1:(n_quarters-1)])

r_t_hat_sbd_other <- colSums(r_pred_table_sbd_other * risk_set_ind_other) / R_tq_other


## - Evaluation by presence of accidental death rider

## - Yes

R_tq_adry <- rep(0, n_quarters)
adry_index <- which(uslapseagent_test$acc.death.rider=="Rider")
adry_sample_size <- length(adry_index)
R_tq_adry[1] <- nrow(uslapseagent_test[uslapseagent_test$acc.death.rider=="Rider",])
for(q in 2:n_quarters){
  R_tq_adry[q] <- sum(ifelse(t_i[adry_index]>q-1,1,0))
}


# \hat{r}_{t,i} - individually predicted rates for computation of eq. (10) in Milhaud-Dutang
r_pred_table_adry <- matrix(NA, nrow(uslapseagent_test[uslapseagent_test$acc.death.rider=="Rider",]), n_quarters)
r_pred_table_adry[,1] <- F_t_table[adry_index,1]
r_pred_table_adry[,2:n_quarters] <- (F_t_table[adry_index,2:n_quarters] - F_t_table[adry_index,1:(n_quarters-1)]) / (1 - F_t_table[adry_index,1:(n_quarters-1)])

# - individual being at risk in quarter t. First column corresponds to individuals at risk at time 0, 2nd col corresponds to individual
# at the end of the first quarter, and so on
risk_set_ind_adry <- matrix(NA, nrow(uslapseagent_test[uslapseagent_test$acc.death.rider=="Rider",]), n_quarters)
risk_set_ind_adry[,1] <- 1

counter <- 1 
for(i in adry_index){
  for(j in 1:(n_quarters-1)){
    risk_set_ind_adry[counter,j+1] <- ifelse(t_i[i]>j,1,0)
  }
  counter <- counter + 1
}

# - \hat{r}_{t}
r_t_hat_adry <- colSums(r_pred_table_adry * risk_set_ind_adry) / R_tq_adry
## - reformulated \hat{r}_{t}
### r_t_hat <- colSums(r_pred_table) / R_tq # - We shouldn't multiplu r_pred_table by the risk set

# - Calculation of the indicator variable if $T$ occurs in the quarter
t_i_between_adry <- matrix(0, nrow(uslapseagent_test[uslapseagent_test$acc.death.rider=="Rider",]), n_quarters)

# - Empirical
# t_i between quarters

t_i_between_adry[,1] <- ifelse((t_i[adry_index]<=1) & (d_ci[adry_index,1]==1),1,0)

for(i in 2:n_quarters){
  counter <- 1
  for(j in adry_index){
    t_i_between_adry[counter,i] <- ifelse(((t_i[j]>(i-1)) & (t_i[j]<= i) & (d_ci[j,1]==1)),1,0)
    counter <- counter + 1
  }
}

## -  Calculation of the empirical \hat{r}_t
r_t_emp_adry <- colSums(t_i_between_adry) / R_tq_adry

## - Comparison with Cox PH
r_pred_table_cox_adry <- matrix(NA, adry_sample_size, n_quarters)
r_pred_table_cox_adry[,1] <- F_t_cox_table[adry_index,1]
r_pred_table_cox_adry[,2:n_quarters] <- (F_t_cox_table[adry_index,2:n_quarters] - F_t_cox_table[adry_index,1:(n_quarters-1)]) / (1 - F_t_cox_table[adry_index,1:(n_quarters-1)])

r_t_hat_cox_adry <- colSums(r_pred_table_cox_adry * risk_set_ind_adry) / R_tq_adry

## - Comparison with Subdistribution
r_pred_table_sbd_adry <- matrix(NA, adry_sample_size, n_quarters)
r_pred_table_sbd_adry[,1] <- F_t_sbd_table[adry_index,1]
r_pred_table_sbd_adry[,2:n_quarters] <- (F_t_sbd_table[adry_index,2:n_quarters] - F_t_sbd_table[adry_index,1:(n_quarters-1)]) / (1 - F_t_sbd_table[adry_index,1:(n_quarters-1)])

r_t_hat_sbd_adry <- colSums(r_pred_table_sbd_adry * risk_set_ind_adry) / R_tq_adry


## - No Rider

R_tq_adrn <- rep(0, n_quarters)
adrn_index <- which(uslapseagent_test$acc.death.rider!="Rider")
adrn_sample_size <- length(adrn_index)
R_tq_adrn[1] <- nrow(uslapseagent_test[uslapseagent_test$acc.death.rider!="Rider",])
for(q in 2:n_quarters){
  R_tq_adrn[q] <- sum(ifelse(t_i[adrn_index]>q-1,1,0))
}


# \hat{r}_{t,i} - individually predicted rates for computation of eq. (10) in Milhaud-Dutang
r_pred_table_adrn <- matrix(NA, nrow(uslapseagent_test[uslapseagent_test$acc.death.rider!="Rider",]), n_quarters)
r_pred_table_adrn[,1] <- F_t_table[adrn_index,1]
r_pred_table_adrn[,2:n_quarters] <- (F_t_table[adrn_index,2:n_quarters] - F_t_table[adrn_index,1:(n_quarters-1)]) / (1 - F_t_table[adrn_index,1:(n_quarters-1)])

# - individual being at risk in quarter t. First column corresponds to individuals at risk at time 0, 2nd col corresponds to individual
# at the end of the first quarter, and so on
risk_set_ind_adrn <- matrix(NA, nrow(uslapseagent_test[uslapseagent_test$acc.death.rider!="Rider",]), n_quarters)
risk_set_ind_adrn[,1] <- 1

counter <- 1 
for(i in adrn_index){
  for(j in 1:(n_quarters-1)){
    risk_set_ind_adrn[counter,j+1] <- ifelse(t_i[i]>j,1,0)
  }
  counter <- counter + 1
}

# - \hat{r}_{t}
r_t_hat_adrn <- colSums(r_pred_table_adrn * risk_set_ind_adrn) / R_tq_adrn
## - reformulated \hat{r}_{t}
### r_t_hat <- colSums(r_pred_table) / R_tq # - We shouldn't multiplu r_pred_table by the risk set

plot(c(1:n_quarters), r_t_hat_adrn, type='l', ylim=c(0,0.06), lwd=2)

# - Calculation of the indicator variable if $T$ occurs in the quarter
t_i_between_adrn <- matrix(0, nrow(uslapseagent_test[uslapseagent_test$acc.death.rider!="Rider",]), n_quarters)

# - Empirical
# t_i between quarters

t_i_between_adrn[,1] <- ifelse((t_i[adrn_index]<=1) & (d_ci[adrn_index,1]==1),1,0)

for(i in 2:n_quarters){
  counter <- 1
  for(j in adrn_index){
    t_i_between_adrn[counter,i] <- ifelse(((t_i[j]>(i-1)) & (t_i[j]<= i) & (d_ci[j,1]==1)),1,0)
    counter <- counter + 1
  }
}

## -  Calculation of the empirical \hat{r}_t
r_t_emp_adrn <- colSums(t_i_between_adrn) / R_tq_adrn

## - Comparison with Cox PH
r_pred_table_cox_adrn <- matrix(NA, adrn_sample_size, n_quarters)
r_pred_table_cox_adrn[,1] <- F_t_cox_table[adrn_index,1]
r_pred_table_cox_adrn[,2:n_quarters] <- (F_t_cox_table[adrn_index,2:n_quarters] - F_t_cox_table[adrn_index,1:(n_quarters-1)]) / (1 - F_t_cox_table[adrn_index,1:(n_quarters-1)])

r_t_hat_cox_adrn <- colSums(r_pred_table_cox_adrn * risk_set_ind_adrn) / R_tq_adrn

## - Comparison with Subdistribution
r_pred_table_sbd_adrn <- matrix(NA, adrn_sample_size, n_quarters)
r_pred_table_sbd_adrn[,1] <- F_t_sbd_table[adrn_index,1]
r_pred_table_sbd_adrn[,2:n_quarters] <- (F_t_sbd_table[adrn_index,2:n_quarters] - F_t_sbd_table[adrn_index,1:(n_quarters-1)]) / (1 - F_t_sbd_table[adrn_index,1:(n_quarters-1)])

r_t_hat_sbd_adrn <- colSums(r_pred_table_sbd_adrn * risk_set_ind_adrn) / R_tq_adrn


# - Plot for final version of the paper
par(mfrow=c(2,3), mar = c(4.5,4,1,2.5))
plot(c(1:60), r_t_hat_adry[1:60], type='l', ylim=c(0,0.06), lwd=2, main="With Acc. death rider", ylab="Surrending rates", xlab="Quarter")
lines(c(1:60), r_t_emp_adry[1:60], lty=2, lwd=2)
lines(c(1:60), r_t_hat_cox_adry[1:60], lty=3, lwd=1)
lines(c(1:60), r_t_hat_sbd_adry[1:60], lty=4, lwd=1)
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=0.75, bty='n')

plot(c(1:60), r_t_hat_adrn[1:60], type='l', ylim=c(0,0.06), lwd=2, main="No Acc. death rider", ylab="Surrending rates", xlab="Quarter")
lines(c(1:60), r_t_emp_adrn[1:60], lty=2, lwd=2)
lines(c(1:60), r_t_hat_cox_adrn[1:60], lty=3, lwd=1)
lines(c(1:60), r_t_hat_sbd_adrn[1:60], lty=4, lwd=1)
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=0.75, bty='n')

plot(c(1:60), R_RMSEq_DPM_adry[1:60], type='l', lwd=2, ylim=c(0,0.03), main="Rolling RMSE by Acc. Death Rid.", ylab="Rolling RMSE", xlab="Quarter")
lines(c(1:60), R_RMSEq_Cox_adry[1:60], type='l', lty=3)
lines(c(1:60), R_RMSEq_Sbd_adry[1:60], type='l', lty=4)
lines(c(1:60), R_RMSEq_DPM_adrn[1:60], type='l', lwd=2, ylim=c(0,0.03), col='grey70')
lines(c(1:60), R_RMSEq_Cox_adrn[1:60], type='l', lty=3, col='grey70')
lines(c(1:60), R_RMSEq_Sbd_adrn[1:60], type='l', lty=4, col='grey70')
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=1, bty='n')
legend("right", legend=c("With", "Without"), pch=16, col=c('black', 'grey70'), cex=1, bty='n')


plot(c(1:60), r_t_hat_infra[1:60], type='l', ylim=c(0,0.06), lwd=2, col='black', main="InfraAnnual payment", ylab="Surrending rates", xlab="Quarter")
lines(c(1:60), r_t_emp_infra[1:60], lty=2, lwd=2, col='black')
lines(c(1:60), r_t_hat_cox_infra[1:60], lty=3, lwd=1, col='black')
lines(c(1:60), r_t_hat_sbd_infra[1:60], lty=4, lwd=1, col='black')
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=0.75, bty='n')

plot(c(1:60), r_t_hat_other[1:60], type='l', ylim=c(0,0.06), lwd=2, main="Annual and Other payment", ylab="Surrending rates", xlab="Quarter")
lines(c(1:60), r_t_emp_other[1:60], lty=2, lwd=2)
lines(c(1:60), r_t_hat_cox_other[1:60], lty=3, lwd=1)
lines(c(1:60), r_t_hat_sbd_other[1:60], lty=4, lwd=1)
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=0.75, bty='n')

plot(c(1:60), R_RMSEq_DPM_infra[1:60], type='l', lwd=2, ylim=c(0,0.03), main="Rolling RMSE by Payment Frequency", ylab="Rolling RMSE", xlab="Quarter", col='black')
lines(c(1:60), R_RMSEq_Cox_infra[1:60], type='l', lty=3, col='black')
lines(c(1:60), R_RMSEq_Sbd_infra[1:60], type='l', lty=4, col='black')
lines(c(1:60), R_RMSEq_DPM_other[1:60], type='l', lwd=2, ylim=c(0,0.03), col='grey70')
lines(c(1:60), R_RMSEq_Cox_other[1:60], type='l', lty=3, col='grey70')
lines(c(1:60), R_RMSEq_Sbd_other[1:60], type='l', lty=4, col='grey70')
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=1, bty='n')
legend("right", legend=c("InfraAnnual", "Annual+Others"), pch=16, col=c('black', 'grey70'), cex=1, bty='n')





par(mfrow=c(2,2), mar = c(4.5,4,1,2.5))
plot(c(1:n_quarters), r_t_hat_males, type='l', ylim=c(0,0.06), lwd=2, col='blue', main="Males", ylab="Surrending rates", xlab="Quarter")
lines(c(1:n_quarters), r_t_emp_males, lty=2, lwd=2, col='blue')
lines(c(1:n_quarters), r_t_hat_cox_males, lty=3, lwd=1, col='blue')
lines(c(1:n_quarters), r_t_hat_sbd_males, lty=4, lwd=1, col='blue')
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=0.75, bty='n')

plot(c(1:n_quarters), r_t_hat_females, type='l', ylim=c(0,0.06), lwd=2, col='red', main="Females", ylab="Surrending rates", xlab="Quarter")
lines(c(1:n_quarters), r_t_emp_females, lty=2, lwd=2, col='red')
lines(c(1:n_quarters), r_t_hat_cox_females, lty=3, lwd=1, col='red')
lines(c(1:n_quarters), r_t_hat_sbd_females, lty=4, lwd=1, col='red')
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=0.75, bty='n')

plot(c(1:n_quarters), r_t_hat_adry, type='l', ylim=c(0,0.06), lwd=2, col='blue', main="Acc. Death Rid. Yes", ylab="Surrending rates", xlab="Quarter")
lines(c(1:n_quarters), r_t_emp_adry, lty=2, lwd=2, col='blue')
lines(c(1:n_quarters), r_t_hat_cox_adry, lty=3, lwd=1, col='blue')
lines(c(1:n_quarters), r_t_hat_sbd_adry, lty=4, lwd=1, col='blue')
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=0.75, bty='n')

plot(c(1:n_quarters), r_t_hat_adrn, type='l', ylim=c(0,0.06), lwd=2, col='red', main="Acc. Death Rid. No", ylab="Surrending rates", xlab="Quarter")
lines(c(1:n_quarters), r_t_emp_adrn, lty=2, lwd=2, col='red')
lines(c(1:n_quarters), r_t_hat_cox_adrn, lty=3, lwd=1, col='red')
lines(c(1:n_quarters), r_t_hat_sbd_adrn, lty=4, lwd=1, col='red')
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=0.75, bty='n')


plot(c(1:n_quarters), r_t_hat_infra, type='l', ylim=c(0,0.06), lwd=2, col='black', main="InfraAnnual payment", ylab="Surrending rates", xlab="Quarter")
lines(c(1:n_quarters), r_t_emp_infra, lty=2, lwd=2, col='black')
lines(c(1:n_quarters), r_t_hat_cox_infra, lty=3, lwd=1, col='black')
lines(c(1:n_quarters), r_t_hat_sbd_infra, lty=4, lwd=1, col='black')
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=0.75, bty='n')

plot(c(1:n_quarters), r_t_hat_other, type='l', ylim=c(0,0.06), lwd=2, col='grey70', main="Annual and Other payment", ylab="Surrending rates", xlab="Quarter")
lines(c(1:n_quarters), r_t_emp_other, lty=2, lwd=2, col='grey70')
lines(c(1:n_quarters), r_t_hat_cox_other, lty=3, lwd=1, col='grey70')
lines(c(1:n_quarters), r_t_hat_sbd_other, lty=4, lwd=1, col='grey70')
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=0.75, bty='n')

# - Calculation of the Rolling RMSE
R_RMSEq_DPM_males <- R_RMSEq_DPM_females <- R_RMSEq_DPM_infra <- R_RMSEq_DPM_other <-  R_RMSEq_Cox_males <- R_RMSEq_Cox_females <- R_RMSEq_Cox_infra <- R_RMSEq_Cox_other <-  R_RMSEq_Sbd_males <- R_RMSEq_Sbd_females <- R_RMSEq_Sbd_infra <- R_RMSEq_Sbd_other <- rep(0, n_quarters)
R_RMSEq_DPM_adry <- R_RMSEq_DPM_adrn <- R_RMSEq_DPM_infra <- R_RMSEq_DPM_other <-  R_RMSEq_Cox_adry <- R_RMSEq_Cox_adrn <- R_RMSEq_Cox_infra <- R_RMSEq_Cox_other <-  R_RMSEq_Sbd_adry <- R_RMSEq_Sbd_adrn <- R_RMSEq_Sbd_infra <- R_RMSEq_Sbd_other <- rep(0, n_quarters)

R_RMSEq_DPM <- R_RMSEq_Sbd <- R_RMSEq_Cox <- rep(0, n_quarters)


for(q in 1:n_quarters){
  R_RMSEq_Cox[q] <- sqrt(sum((r_t_hat_cox[1:q]-r_t_emp[1:q])^2)/q)
  R_RMSEq_Sbd[q] <- sqrt(sum((r_t_hat_sbd[1:q]-r_t_emp[1:q])^2)/q)
  R_RMSEq_DPM[q] <- sqrt(sum((r_t_hat[1:q]-r_t_emp[1:q])^2)/q)
#  R_RMSEq_Cox_males[q] <- sqrt(sum((r_t_hat_cox_males[1:q]-r_t_emp_males[1:q])^2)/q)
#  R_RMSEq_Sbd_males[q] <- sqrt(sum((r_t_hat_sbd_males[1:q]-r_t_emp_males[1:q])^2)/q)
#  R_RMSEq_DPM_females[q] <- sqrt(sum((r_t_hat_females[1:q]-r_t_emp_females[1:q])^2)/q)
#  R_RMSEq_Cox_females[q] <- sqrt(sum((r_t_hat_cox_females[1:q]-r_t_emp_females[1:q])^2)/q)
#  R_RMSEq_Sbd_females[q] <- sqrt(sum((r_t_hat_sbd_females[1:q]-r_t_emp_females[1:q])^2)/q)
  R_RMSEq_DPM_adry[q] <- sqrt(sum((r_t_hat_adry[1:q]-r_t_emp_adry[1:q])^2)/q)
  R_RMSEq_Cox_adry[q] <- sqrt(sum((r_t_hat_cox_adry[1:q]-r_t_emp_adry[1:q])^2)/q)
  R_RMSEq_Sbd_adry[q] <- sqrt(sum((r_t_hat_sbd_adry[1:q]-r_t_emp_adry[1:q])^2)/q)
  R_RMSEq_DPM_adrn[q] <- sqrt(sum((r_t_hat_adrn[1:q]-r_t_emp_adrn[1:q])^2)/q)
  R_RMSEq_Cox_adrn[q] <- sqrt(sum((r_t_hat_cox_adrn[1:q]-r_t_emp_adrn[1:q])^2)/q)
  R_RMSEq_Sbd_adrn[q] <- sqrt(sum((r_t_hat_sbd_adrn[1:q]-r_t_emp_adrn[1:q])^2)/q)
  R_RMSEq_DPM_infra[q] <- sqrt(sum((r_t_hat_infra[1:q]-r_t_emp_infra[1:q])^2)/q)
  R_RMSEq_Cox_infra[q] <- sqrt(sum((r_t_hat_cox_infra[1:q]-r_t_emp_infra[1:q])^2)/q)
  R_RMSEq_Sbd_infra[q] <- sqrt(sum((r_t_hat_sbd_infra[1:q]-r_t_emp_infra[1:q])^2)/q)
  R_RMSEq_DPM_other[q] <- sqrt(sum((r_t_hat_other[1:q]-r_t_emp_other[1:q])^2)/q)
  R_RMSEq_Cox_other[q] <- sqrt(sum((r_t_hat_cox_other[1:q]-r_t_emp_other[1:q])^2)/q)
  R_RMSEq_Sbd_other[q] <- sqrt(sum((r_t_hat_sbd_other[1:q]-r_t_emp_other[1:q])^2)/q)
}

#par(mfrow=c(1,1), mar = c(4.5,4,1,2.5))
plot(c(1:n_quarters), R_RMSEq_DPM, type='l', lwd=2, ylim=c(0,0.03), main="Rolling RMSE", ylab="Rolling RMSE", xlab="Quarter")
lines(c(1:n_quarters), R_RMSEq_Cox, type='l', lty=3)
lines(c(1:n_quarters), R_RMSEq_Sbd, type='l', lty=4)
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=0.75, bty='n')

par(mfrow=c(1,2), mar = c(4.5,4,1,2.5))
plot(c(1:n_quarters), R_RMSEq_DPM_males, type='l', lwd=2, ylim=c(0,0.03), main="Rolling RMSE by gender", ylab="Rolling RMSE", xlab="Quarter", col='blue')
lines(c(1:n_quarters), R_RMSEq_Cox_males, type='l', lty=3, col='blue')
lines(c(1:n_quarters), R_RMSEq_Sbd_males, type='l', lty=4, col='blue')
lines(c(1:n_quarters), R_RMSEq_DPM_females, type='l', lwd=2, ylim=c(0,0.03), col='red')
lines(c(1:n_quarters), R_RMSEq_Cox_females, type='l', lty=3, col='red')
lines(c(1:n_quarters), R_RMSEq_Sbd_females, type='l', lty=4, col='red')
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=1, bty='n')
legend("right", legend=c("Males", "Females"), pch=16, col=c('blue', 'red'), cex=1, bty='n')

plot(c(1:n_quarters), R_RMSEq_DPM_adry, type='l', lwd=2, ylim=c(0,0.03), main="Rolling RMSE by Acc. Death Rid.", ylab="Rolling RMSE", xlab="Quarter", col='blue')
lines(c(1:n_quarters), R_RMSEq_Cox_adry, type='l', lty=3, col='blue')
lines(c(1:n_quarters), R_RMSEq_Sbd_adry, type='l', lty=4, col='blue')
lines(c(1:n_quarters), R_RMSEq_DPM_adrn, type='l', lwd=2, ylim=c(0,0.03), col='red')
lines(c(1:n_quarters), R_RMSEq_Cox_adrn, type='l', lty=3, col='red')
lines(c(1:n_quarters), R_RMSEq_Sbd_adrn, type='l', lty=4, col='red')
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=1, bty='n')
legend("right", legend=c("Yes", "No"), pch=16, col=c('blue', 'red'), cex=1, bty='n')

plot(c(1:n_quarters), R_RMSEq_DPM_infra, type='l', lwd=2, ylim=c(0,0.03), main="Rolling RMSE by Payment Frequency", ylab="Rolling RMSE", xlab="Quarter", col='black')
lines(c(1:n_quarters), R_RMSEq_Cox_infra, type='l', lty=3, col='black')
lines(c(1:n_quarters), R_RMSEq_Sbd_infra, type='l', lty=4, col='black')
lines(c(1:n_quarters), R_RMSEq_DPM_other, type='l', lwd=2, ylim=c(0,0.03), col='grey70')
lines(c(1:n_quarters), R_RMSEq_Cox_other, type='l', lty=3, col='grey70')
lines(c(1:n_quarters), R_RMSEq_Sbd_other, type='l', lty=4, col='grey70')
legend("topright", legend=c("DPM", "Cox PH", 'Subdistribution'), lty=c(1,3,4), lwd=c(2,1,1), cex=1, bty='n')
legend("right", legend=c("InfraAnnual", "Annual+Others"), pch=16, col=c('black', 'grey70'), cex=1, bty='n')


#=================== - Predicted rates by premium amount (standardized) at quarter 10 - ================

ap_points <- length(unique(uslapseagent_test$annual.premium))

r_pred_12_table <- cbind(r_pred_table[,12], uslapseagent_test$annual.premium)
risk_set_ind_12 <- risk_set_ind[,12]
R_t12 <- R_tq[12]



R_t12_ap <- rep(0, ap_points)
ap_position <- matrix(0, sample_size_test, ap_points)

for(i in 1:sample_size_test){
  for(ap in 1:ap_points){
    ap_position[i,ap] <- ifelse(uslapseagent_test$annual.premium[i]==unique(uslapseagent_test$annual.premium)[ap],1,0)
  }
}

adrn_index <- which(uslapseagent_test$acc.death.rider!="Rider")
adrn_sample_size <- length(adrn_index)
R_tq_adrn[1] <- nrow(uslapseagent_test[uslapseagent_test$acc.death.rider!="Rider",])
for(q in 2:n_quarters){
  R_tq_adrn[q] <- sum(ifelse(t_i[adrn_index]>q,1,0))
}


# \hat{r}_{t,i} - individually predicted rates for computation of eq. (10) in Milhaud-Dutang
r_pred_table_adrn <- matrix(NA, nrow(uslapseagent_test[uslapseagent_test$acc.death.rider!="Rider",]), n_quarters)
r_pred_table_adrn[,1] <- F_t_table[adrn_index,1]
r_pred_table_adrn[,2:n_quarters] <- (F_t_table[adrn_index,2:n_quarters] - F_t_table[adrn_index,1:(n_quarters-1)]) / (1 - F_t_table[adrn_index,1:(n_quarters-1)])

# - individual being at risk in quarter t. First column corresponds to individuals at risk at time 0, 2nd col corresponds to individual
# at the end of the first quarter, and so on
risk_set_ind_adrn <- matrix(NA, nrow(uslapseagent_test[uslapseagent_test$acc.death.rider!="Rider",]), n_quarters)
risk_set_ind_adrn[,1] <- 1

counter <- 1 
for(i in adrn_index){
  for(j in 1:(n_quarters-1)){
    risk_set_ind_adrn[counter,j+1] <- ifelse(t_i[i]>j,1,0)
  }
  counter <- counter + 1
}

# - \hat{r}_{t}
r_t_hat_adrn <- colSums(r_pred_table_adrn * risk_set_ind_adrn) / R_tq_adrn
## - reformulated \hat{r}_{t}
### r_t_hat <- colSums(r_pred_table) / R_tq # - We shouldn't multiplu r_pred_table by the risk set

plot(c(1:n_quarters), r_t_hat_adrn, type='l', ylim=c(0,0.06), lwd=2)

# - Calculation of the indicator variable if $T$ occurs in the quarter
t_i_between_adrn <- matrix(0, nrow(uslapseagent_test[uslapseagent_test$acc.death.rider!="Rider",]), n_quarters)

# - Empirical
# t_i between quarters

t_i_between_adrn[,1] <- ifelse((t_i[adrn_index]<=1) & (d_ci[adrn_index,1]==1),1,0)

for(i in 2:n_quarters){
  counter <- 1
  for(j in adrn_index){
    t_i_between_adrn[counter,i] <- ifelse(((t_i[j]>(i-1)) & (t_i[j]<= i) & (d_ci[j,1]==1)),1,0)
    counter <- counter + 1
  }
}

## -  Calculation of the empirical \hat{r}_t
r_t_emp_adrn <- colSums(t_i_between_adrn) / R_tq_adrn

## - Comparison with Cox PH
r_pred_table_cox_adrn <- matrix(NA, adrn_sample_size, n_quarters)
r_pred_table_cox_adrn[,1] <- F_t_cox_table[adrn_index,1]
r_pred_table_cox_adrn[,2:n_quarters] <- (F_t_cox_table[adrn_index,2:n_quarters] - F_t_cox_table[adrn_index,1:(n_quarters-1)]) / (1 - F_t_cox_table[adrn_index,1:(n_quarters-1)])

r_t_hat_cox_adrn <- colSums(r_pred_table_cox_adrn * risk_set_ind_adrn) / R_tq_adrn









#================== - Posterior predictive density - ==============================

n_samples <- 1000000
row_sample <- sample(seq_bi_thn, n_samples, replace=TRUE)
pi_sample <- pi_mw[row_sample,]
theta1_sample1 <- theta_star[row_sample, seq(1, 73, 3)]
beta1_sample <- beta1[row_sample,]
sigma2_1_sample <- sigma2_c[row_sample, 1]

s_sample <- rep(0, n_samples)
theta1_sample2 <- rep(0, n_samples) #matrix(NA, nrow=n_samples, col=K)
y_sample <- matrix(NA, nrow=n_samples, ncol=8)

for(sample in 1:n_samples){
  s_sample[sample] <- sample(c(1:K), 1, prob=pi_sample[sample,])
  theta1_sample2[sample] <- theta1_sample1[sample, s_sample[sample]]
  
  # - UW Age: Young-Middle
  ## - Acc. Death Rid: NO, PayFreq: I
  y_sample[sample,1] <- rnorm(1, mean = (beta1_sample[sample,6] + theta1_sample2[sample]), sd = sqrt(sigma2_1_sample[sample]))
  ## - Acc. Death Rid: NO, PayFreq: Annual+Other
  y_sample[sample,2] <- rnorm(1, mean = (beta1_sample[sample,6] + beta1_sample[sample,5] + theta1_sample2[sample]), sd = sqrt(sigma2_1_sample[sample]))
  ## - Acc. Death Rid: Y, PayFreq: I
  y_sample[sample,3] <- rnorm(1, mean = (beta1_sample[sample,6] + beta1_sample[sample,3] + theta1_sample2[sample]), sd = sqrt(sigma2_1_sample[sample]))
  ## - Acc. Death Rid: Y, PayFreq: Annual+Other
  y_sample[sample,4] <- rnorm(1, mean = (beta1_sample[sample,6] + beta1_sample[sample,5] + beta1_sample[sample,3] + theta1_sample2[sample]), sd = sqrt(sigma2_1_sample[sample]))

    # - UW Age: Old
  ## - Acc. Death Rid: NO, PayFreq: I
  y_sample[sample,5] <- rnorm(1, mean = (theta1_sample2[sample]), sd = sqrt(sigma2_1_sample[sample]))
  ## - Acc. Death Rid: NO, PayFreq: Annual+Other
  y_sample[sample,6] <- rnorm(1, mean = (beta1_sample[sample,5] + theta1_sample2[sample]), sd = sqrt(sigma2_1_sample[sample]))
  ## - Acc. Death Rid: Y, PayFreq: I
  y_sample[sample,7] <- rnorm(1, mean = (beta1_sample[sample,3] + theta1_sample2[sample]), sd = sqrt(sigma2_1_sample[sample]))
  ## - Acc. Death Rid: Y, PayFreq: Annual+Other
  y_sample[sample,8] <- rnorm(1, mean = (beta1_sample[sample,5] + beta1_sample[sample,3] + theta1_sample2[sample]), sd = sqrt(sigma2_1_sample[sample]))
  
}

t_sample <- exp(y_sample)

par(mfrow=c(2,2), mar = c(4.5,4,1,2.5))

## - Plots for slides
par(mfrow=c(1,2), mar = c(4.5,4,1,2.5))
plot(density(t_sample[,1]), main="UW Age: 0-54", xlab="T (quarters)", xlim=c(0,250), lwd=2)
lines(density(t_sample[,2]), lty=2, lwd=2, col="grey70")
lines(density(t_sample[,3]), lty=3, lwd=2, col="red")
lines(density(t_sample[,4]), lty=4, lwd=2, col="blue")

plot(density(t_sample[,5]), main="UW Age: 55 and over", xlab="T (quarters)", xlim=c(0,250), lwd=2)
lines(density(t_sample[,6]), lty=2, lwd=2, col="grey70")
lines(density(t_sample[,7]), lty=3, lwd=2, col="red")
lines(density(t_sample[,8]), lty=4, lwd=2, col="blue")
legend("topright", legend=c("ADR: N; PF: Infrannual", "ADR: N; PF: Ann.+Ot.", "ADR: Y; PF: Infrannual", "ADR: Y; PF: Ann.+Ot."), lty=c(1,2,3,4), lwd=2, col=c("black", "grey70", "red", "blue"), cex=0.75, bty='n')


par(mfrow=c(1,2), mar = c(4.5,4,1,2.5))
plot(AnnPrem_seq, exp(E_y_sample[1] + mean(beta1[seq_bi_thn,1]) * AnnPrem_seq), type='l', lty=1, ylim=c(30,45), xlab="Ann. Prem. (std)", ylab="E(T) (quarters)", lwd=2)
lines(AnnPrem_seq, exp(E_y_sample[2] + mean(beta1[seq_bi_thn,1]) * AnnPrem_seq), lty=2, lwd=2, col='grey70')
lines(AnnPrem_seq, exp(E_y_sample[3] + mean(beta1[seq_bi_thn,1]) * AnnPrem_seq), lty=3, lwd=2, col='red')
lines(AnnPrem_seq, exp(E_y_sample[4] + mean(beta1[seq_bi_thn,1]) * AnnPrem_seq), lty=4, lwd=2, col='blue')

plot(AnnPrem_seq, exp(E_y_sample[5] + mean(beta1[seq_bi_thn,1]) * AnnPrem_seq), type='l', lty=1, ylim=c(32,50), xlab="Ann. Prem. (std)", ylab="E(T) (quarters)", lwd=2)
lines(AnnPrem_seq, exp(E_y_sample[6] + mean(beta1[seq_bi_thn,1]) * AnnPrem_seq), lty=2, lwd=2, col='grey70')
lines(AnnPrem_seq, exp(E_y_sample[7] + mean(beta1[seq_bi_thn,1]) * AnnPrem_seq), lty=3, lwd=2, col='red')
lines(AnnPrem_seq, exp(E_y_sample[8] + mean(beta1[seq_bi_thn,1]) * AnnPrem_seq), lty=4, lwd=2, col='blue')
legend("topright", legend=c("ADR: N; PF: Infrannual", "ADR: N; PF: Ann.+Ot.", "ADR: Y; PF: Infrannual", "ADR: Y; PF: Ann.+Ot."), lty=c(1,2,3,4), lwd=2, col=c("black", "grey70", "red", "blue"), cex=0.75, bty='n')

## - Plots for paper

plot(density(t_sample[,1]), main="UW Age: 0-54", xlab="T (quarters)", xlim=c(0,250), lwd=1.5)
lines(density(t_sample[,2]), lty=2, lwd=1.5)
lines(density(t_sample[,3]), lty=3, lwd=1.5)
lines(density(t_sample[,4]), lty=4, lwd=1.5)

plot(density(t_sample[,5]), main="UW Age: 55 and over", xlab="T (quarters)", xlim=c(0,250), lwd=1.5)
lines(density(t_sample[,6]), lty=2, lwd=1.5)
lines(density(t_sample[,7]), lty=3, lwd=1.5)
lines(density(t_sample[,8]), lty=4, lwd=1.5)
legend("topright", legend=c("ADR: N; PF: Infrannual", "ADR: N; PF: Ann.+Ot.", "ADR: Y; PF: Infrannual", "ADR: Y; PF: Ann.+Ot."), lty=c(1,2,3,4), lwd=c(1.5,1.5,1.5,1.5), cex=0.75, bty='n')


plot(density(y_sample[,1]), main="", xlab="Time to surrending")
lines(density(y_sample[,2]), lty=2)
lines(density(y_sample[,3]), lty=3)
lines(density(y_sample[,4]), lty=4)

plot(density(y_sample[,5]), main="", xlab="Time to surrending")
lines(density(y_sample[,6]), lty=2)
lines(density(y_sample[,7]), lty=3)
lines(density(y_sample[,8]), lty=4)

#========================== - Posterior predictive mean - =============================

E_y_sample <- rep(NA, 8)

  # - UW Age: Young-Middle
  ## - Acc. Death Rid: NO, PayFreq: I
  E_y_sample[1] <- mean(beta1[seq_bi_thn,6]) + sum(colMeans(pi_mw[seq_bi_thn,]) * colMeans(theta_star[seq_bi_thn,seq(1,73,3)]))
  ## - Acc. Death Rid: NO, PayFreq: Annual+Other
  E_y_sample[2] <- mean(beta1[seq_bi_thn,6]) + mean(beta1[seq_bi_thn,5]) + sum(colMeans(pi_mw[seq_bi_thn,]) * colMeans(theta_star[seq_bi_thn,seq(1,73,3)]))
  ## - Acc. Death Rid: Y, PayFreq: I
  E_y_sample[3] <- mean(beta1[seq_bi_thn,6]) + mean(beta1[seq_bi_thn,3]) + sum(colMeans(pi_mw[seq_bi_thn,]) * colMeans(theta_star[seq_bi_thn,seq(1,73,3)]))
  ## - Acc. Death Rid: Y, PayFreq: Annual+Other
  E_y_sample[4] <- mean(beta1[seq_bi_thn,6]) + mean(beta1[seq_bi_thn,5]) + mean(beta1[seq_bi_thn,3]) + sum(colMeans(pi_mw[seq_bi_thn,]) * colMeans(theta_star[seq_bi_thn,seq(1,73,3)]))

  # - UW Age: Old
  ## - Acc. Death Rid: NO, PayFreq: I
  E_y_sample[5] <- sum(colMeans(pi_mw[seq_bi_thn,]) * colMeans(theta_star[seq_bi_thn,seq(1,73,3)]))
  ## - Acc. Death Rid: NO, PayFreq: Annual+Other
  E_y_sample[6] <- mean(beta1[seq_bi_thn,5]) + sum(colMeans(pi_mw[seq_bi_thn,]) * colMeans(theta_star[seq_bi_thn,seq(1,73,3)]))
  ## - Acc. Death Rid: Y, PayFreq: I
  E_y_sample[7] <- mean(beta1[seq_bi_thn,3]) + sum(colMeans(pi_mw[seq_bi_thn,]) * colMeans(theta_star[seq_bi_thn,seq(1,73,3)]))
  ## - Acc. Death Rid: Y, PayFreq: Annual+Other
  E_y_sample[8] <- mean(beta1[seq_bi_thn,5]) + mean(beta1[seq_bi_thn,3]) + sum(colMeans(pi_mw[seq_bi_thn,]) * colMeans(theta_star[seq_bi_thn,seq(1,73,3)]))
  

AnnPrem_seq <- seq(-3,3, 0.1)
par(mfrow=c(1,2), mar = c(4.5,4,1,2.5))
plot(AnnPrem_seq, exp(E_y_sample[1] + mean(beta1[seq_bi_thn,1]) * AnnPrem_seq), type='l', lty=1, ylim=c(30,45), xlab="Ann. Prem. (std)", ylab="E(T) (quarters)")
lines(AnnPrem_seq, exp(E_y_sample[2] + mean(beta1[seq_bi_thn,1]) * AnnPrem_seq), lty=2)
lines(AnnPrem_seq, exp(E_y_sample[3] + mean(beta1[seq_bi_thn,1]) * AnnPrem_seq), lty=3)
lines(AnnPrem_seq, exp(E_y_sample[4] + mean(beta1[seq_bi_thn,1]) * AnnPrem_seq), lty=4)

plot(AnnPrem_seq, exp(E_y_sample[5] + mean(beta1[seq_bi_thn,1]) * AnnPrem_seq), type='l', lty=1, ylim=c(32,50), xlab="Ann. Prem. (std)", ylab="E(T) (quarters)")
lines(AnnPrem_seq, exp(E_y_sample[6] + mean(beta1[seq_bi_thn,1]) * AnnPrem_seq), lty=2)
lines(AnnPrem_seq, exp(E_y_sample[7] + mean(beta1[seq_bi_thn,1]) * AnnPrem_seq), lty=3)
lines(AnnPrem_seq, exp(E_y_sample[8] + mean(beta1[seq_bi_thn,1]) * AnnPrem_seq), lty=4)


