library(ggplot2)
library(tidyr)
library(dplyr)

#==========
#QUESTION 1
#==========
#observation times in days
t.obs <- c(0,7,14,30*(1:6),30*seq(12,30,by=6))

df <- read.table("Weight.txt", header = TRUE)
colnames(df) <- c("Boy", "0d", "7d", "14d", "1m", "2m", "3m", "4m", "5m", "6m", "12m", "18m", "24m", "30m")


mean_boys  <- colMeans(df[df[,1] == 1, 2:14])
mean_girls <- colMeans(df[df[,1] == 0, 2:14])

plot(t.obs, mean_boys, type = "n", 
     ylim = c(min(df[, 2:14]), max(df[, 2:14])),
     xlab = "Days since birth", 
     ylab = "Weight (kg)",
     main = "Growth Evolution: Boys vs Girls")

grid()

#Boys in Blue, Girls in orange
#individual observations
for(i in 1:nrow(df)) {
  lines(t.obs, df[i, 2:14], 
        col = ifelse(df[i,1] == 1, rgb(0,0,1,0.05), rgb(1,1,0,0.05))) 
}

lines(t.obs, mean_boys,  col = "blue",  lwd = 3, type = "o", pch = 16)
lines(t.obs, mean_girls, col = "orange", lwd = 3, type = "o", pch = 16)

legend("bottomright", 
       legend = c("Boys (Mean)", "Girls (Mean)"), 
       col = c("blue", "orange"), 
       lwd = 3, pch = 16, bty = "n")


#=========
#QUESTION 2
#==========


#QUESTION (d)

estimate_betas <- function(means, times) {
  b0 <- means[1] #Initial weight
  b1 <- (means[13] - means[11]) / (times[13] - times[11]) #Three last observations
  intercept_longterm <- means[13] - (b1 * times[13]) #beta0 + Beta2
  b2 <- intercept_longterm - b0
  v_initial <- (means[2] - means[1]) / (times[2] - times[1]) #First 7 days velocity
  b3 <- v_initial - b1 #See Taylor expansion of a)
  return(c(Beta0 = b0, Beta1 = b1, Beta2 = b2, Beta3 = b3))
}

cat("Rough Estimates for Boys:\n")
print(estimate_betas(mean_boys, t.obs))

cat("\nRough Estimates for Girls:\n")
print(estimate_betas(mean_girls, t.obs))


#QUESTION (e)
is_boy <- df[, 1] == 1
weights <- as.matrix(df[, 2:14])

#Apply range function at every age
range_boys  <- apply(weights[is_boy, ], 2, range)
range_girls <- apply(weights[!is_boy, ], 2, range)

#Estimate sigma by max-min at every age
sigmas_boys  <- (range_boys[2, ] - range_boys[1, ]) / 6
sigmas_girls <- (range_girls[2, ] - range_girls[1, ]) / 6

sigmas <- data.frame(
  Time = colnames(df)[2:14],
  Days = t.obs,
  Sigmas_Boys = sigmas_boys,
  Sigmas_Girls = sigmas_girls
)

print(sigmas)


#=========
#QUESTION 3
#==========


#QUESTION (b)

# For Boys
boys_data <- df[df$Boy == 1, 2:14]
Y_boys <- as.matrix(boys_data)

mu_function <- function(times, gamma) {
  gamma0 <- gamma[1]
  gamma1 <- gamma[2]
  gamma2 <- gamma[3]
  gamma3 <- gamma[4]
  
  beta0 <- exp(gamma0)
  beta1 <- exp(-gamma1)
  beta2 <- exp(gamma2)
  beta3 <- exp(-gamma3)
  mu <- beta0 +
    beta1 * times +
    beta2 * (1 - exp(-(beta3 / beta2) * times))
  return(mu)
}

# Log-posterior
logposterior <- function(theta, Y, times) {
  gamma <- theta[1:4]
  tau   <- theta[5]
  if (tau <= 0) {
    return(-Inf)
  }
  mu <- mu_function(times, gamma)
  RSS <- sum(
    (Y - matrix(mu,
                nrow = nrow(Y),
                ncol = length(mu),
                byrow = TRUE))^2
  )
  n <- nrow(Y)
  m <- ncol(Y)
  logpost <- (n * m / 2 - 1) * log(tau) -
    (tau / 2) * RSS
  return(logpost)
}


#QUESTION (c)

# Starting values from Question 2(d), for boys
beta0_init <- 2.950
beta1_init <- 0.0058
beta2_init <- 5.032
beta3_init <- 0.0370

# Rough estimate from Question 2(e)
sigma_init <- 0.25
start <- c(
  log(beta0_init),
  -log(beta1_init),
  log(beta2_init),
  -log(beta3_init),
  1 / sigma_init^2
)

neglogpost <- function(theta) {
  return(-logposterior(theta, Y_boys, t.obs))
}

# Laplace approximation
fit <- nlm(
  f = neglogpost,
  p = start,
  hessian = TRUE
)
theta_hat <- fit$estimate
H <- fit$hessian
Sigma_theta <- solve(H)

cat("Posterior mode:\n")
print(theta_hat)

cat("\nCovariance matrix:\n")
print(Sigma_theta)


#=========
#QUESTION 4
#==========

#QUESTION (b)(i)

library(MASS)
set.seed(123)
M <- 10000

# Initial values of gamma and tau
gamma_current <- theta_hat[1:4]
tau_current   <- theta_hat[5]

# Covariance matrix for gamma from the Laplace approximation
Sigma_gamma <- Sigma_theta[1:4, 1:4]

c_tune <- 1.5
gamma_chain <- matrix(NA, nrow = M, ncol = 4)
tau_chain   <- numeric(M)
colnames(gamma_chain) <- c("gamma0", "gamma1", "gamma2", "gamma3")
accept <- 0

# Compute S(gamma)
compute_RSS <- function(gamma, Y, times) {
  mu <- mu_function(times, gamma)
  RSS <- sum(
    (Y - matrix(mu,
                nrow = nrow(Y),
                ncol = length(mu),
                byrow = TRUE))^2
  )
  
  return(RSS)
}

for (s in 1:M) {
  # Metropolis step for gamma
  gamma_proposal <- as.numeric(
    MASS::mvrnorm(
      n = 1,
      mu = gamma_current,
      Sigma = c_tune^2 * Sigma_gamma
    )
  )
  theta_current  <- c(gamma_current, tau_current)
  theta_proposal <- c(gamma_proposal, tau_current)
  log_alpha <- logposterior(theta_proposal, Y_boys, t.obs) -
    logposterior(theta_current, Y_boys, t.obs)
  if (log(runif(1)) < log_alpha) {
    gamma_current <- gamma_proposal
    accept <- accept + 1
  }
  
  # Gibbs step for tau
  RSS_current <- compute_RSS(gamma_current, Y_boys, t.obs)
  n <- nrow(Y_boys)
  m <- ncol(Y_boys)
  shape_tau <- n * m / 2
  rate_tau  <- RSS_current / 2
  tau_current <- rgamma(1, shape = shape_tau, rate = rate_tau)

  gamma_chain[s, ] <- gamma_current
  tau_chain[s] <- tau_current
}

acceptance_rate <- accept / M
acceptance_rate   # = 0.3674 pour un ctune de 1, 0.2548 pour c_tune de 1.3, 0.2126 pour c_tune de 1.5, 0.1903 pour c_tune de 1.6


#QUESTION (b)(ii)
par(mfrow = c(2,2))

#par(mfrow = c(3, 2))

plot(gamma_chain[, 1], type = "l",
     main = "Traceplot gamma0", ylab = "gamma0")

plot(gamma_chain[, 2], type = "l",
     main = "Traceplot gamma1", ylab = "gamma1")

plot(gamma_chain[, 3], type = "l",
     main = "Traceplot gamma2", ylab = "gamma2")

plot(gamma_chain[, 4], type = "l",
     main = "Traceplot gamma3", ylab = "gamma3")

plot(tau_chain, type = "l",
     main = "Traceplot tau", ylab = "tau")


par(mfrow = c(2, 2))

acf(gamma_chain[, 1], main = "ACF gamma0")
acf(gamma_chain[, 2], main = "ACF gamma1")
acf(gamma_chain[, 3], main = "ACF gamma2")
acf(gamma_chain[, 4], main = "ACF gamma3")
acf(tau_chain, main = "ACF tau")


#QUESTION (b)(iii)

library(coda)

mcmc_chain <- mcmc(
  cbind(
    gamma_chain,
    tau = tau_chain
  )
)

effectiveSize(mcmc_chain)



#QUESTION (b)(iv)

beta0_chain <- exp(gamma_chain[, 1])
beta1_chain <- exp(-gamma_chain[, 2])
beta2_chain <- exp(gamma_chain[, 3])
beta3_chain <- exp(-gamma_chain[, 4])

sigma_chain <- 1 / sqrt(tau_chain)

posterior_summary <- function(x) {
  c(
    Mean = mean(x),
    Median = median(x),
    Q2.5 = quantile(x, 0.025),
    Q97.5 = quantile(x, 0.975)
  )
}

results <- rbind(
  beta0 = posterior_summary(beta0_chain),
  beta1 = posterior_summary(beta1_chain),
  beta2 = posterior_summary(beta2_chain),
  beta3 = posterior_summary(beta3_chain),
  sigma = posterior_summary(sigma_chain)
)

results

# QUESTION 5

# Setup data matrices and sample sizes specifically for boys
N_boys <- nrow(Y_boys)
M_times <- length(t.obs)

jags_data_q5 <- list(
  Y = Y_boys, 
  t = t.obs, 
  N = N_boys, 
  M = M_times
)

jags_model_q5_string <- "
model {
  for (i in 1:N) {
    for (j in 1:M) {
      mu[i,j] <- beta0 + beta1 * t[j] + beta2 * (1 - exp(-(beta3 / beta2) * t[j]))
      Y[i,j] ~ dnorm(mu[i,j], tau)
    }
  }
  
  beta0 ~ dnorm(0, 1.0E-6) T(0, )
  beta1 ~ dnorm(0, 1.0E-6) T(0, )
  beta2 ~ dnorm(0, 1.0E-6) T(0, )
  beta3 ~ dnorm(0, 1.0E-6) T(0, )
  tau   ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
}
"

# Set initial parameter values across three distinct chains using Q3 optimization outputs
inits_q5 <- list(
  list(beta0 = 2.96, beta1 = 0.00550, beta2 = 5.29, beta3 = 0.0374, tau = 9.72),
  list(beta0 = 2.90, beta1 = 0.00500, beta2 = 5.10, beta3 = 0.0360, tau = 9.00),
  list(beta0 = 3.02, beta1 = 0.00600, beta2 = 5.50, beta3 = 0.0390, tau = 10.50)
)

set.seed(123)
model_q5 <- jags.model(textConnection(jags_model_q5_string), data = jags_data_q5, inits = inits_q5, n.chains = 3, n.adapt = 1000)
update(model_q5, 2000)

params_q5 <- c("beta0", "beta1", "beta2", "beta3", "sigma")
samples_q5 <- coda.samples(model_q5, variable.names = params_q5, n.iter = 10000)

cat("\n--- Q5(a): Summary Results ---\n")
print(summary(samples_q5))
print(gelman.diag(samples_q5))
print(effectiveSize(samples_q5))

png("q5a_traceplots.png", width = 900, height = 600)
par(mfrow = c(3, 2), mar = c(2, 2, 2, 1))
for (p in params_q5) { 
  traceplot(samples_q5[, p], main = paste("Traceplot", p)) 
}
dev.off()


# Q5(b) - Compare analytical M-W-G posterior with JAGS output
q4_mean  <- c(beta0=2.958,   beta1=0.00553, beta2=5.287,  beta3=0.0374,  sigma=0.321)
q4_lower <- c(beta0=2.925,   beta1=0.00530, beta2=5.131,  beta3=0.0365,  sigma=0.310)
q4_upper <- c(beta0=2.992,   beta1=0.00575, beta2=5.455,  beta3=0.0384,  sigma=0.332)

q5_stats  <- summary(samples_q5)$statistics
q5_quant  <- summary(samples_q5)$quantiles

cat("\n--- Q5(b): Comparing Custom M-W-G vs JAGS Sampler Calculations ---\n")
cat(sprintf("%-8s  %9s  %20s    %9s  %20s\n", "Param", "Q4 Mean", "Q4 95% CI", "Q5 Mean", "Q5 95% CI"))
for (p in params_q5) {
  cat(sprintf("%-8s  %9.5f  [%8.5f, %8.5f]    %9.5f  [%8.5f, %8.5f]\n",
              p, q4_mean[p], q4_lower[p], q4_upper[p], q5_stats[p, "Mean"], q5_quant[p, "2.5%"], q5_quant[p, "97.5%"]))
}


# Q5(c) - Out of sample weight simulation at t = 90
chains_q5 <- as.matrix(samples_q5)
n_sims    <- nrow(chains_q5)
t_pred    <- 90
y_pred_q5 <- numeric(n_sims)

for (s in 1:n_sims) {
  mu_90        <- chains_q5[s, "beta0"] + chains_q5[s, "beta1"] * t_pred + chains_q5[s, "beta2"] * (1 - exp(-(chains_q5[s, "beta3"] / chains_q5[s, "beta2"]) * t_pred))
  y_pred_q5[s] <- rnorm(1, mean = mu_90, sd = chains_q5[s, "sigma"])
}

cat("\n--- Q5(c): Fixed Model Predictions at Day 90 ---\n")
print(quantile(y_pred_q5, probs = c(0.025, 0.5, 0.975)))

png("q5c_histogram.png", width = 800, height = 500)
hist(y_pred_q5, breaks = 60, col = "steelblue", border = "white", prob = TRUE, main = "Q5c: Prediction Profile (Boys Model Only)", xlab = "Weight (kg)")
lines(density(y_pred_q5), col = "black", lwd = 2)
abline(v = quantile(y_pred_q5, c(0.025, 0.975)), col = "red", lwd = 2, lty = 2)
dev.off()




# QUESTION 6


# Fit alternative full-population multi-level model formulation
Y_all <- as.matrix(df[, 2:14])
boy   <- df$Boy
N_all <- nrow(Y_all)

jags_data_q6 <- list(
  Y = Y_all, 
  t = t.obs, 
  boy = boy, 
  N = N_all, 
  M = M_times
)

jags_model_q6_string <- "
model {
  for (i in 1:N) {
    beta0i[i] ~ dnorm(beta0, tau0)
    
    b0_eff[i] <- beta0i[i] + delta0 * boy[i]
    b1_eff[i] <- beta1     + delta1 * boy[i]
    b2_eff[i] <- beta2     + delta2 * boy[i]
    b3_eff[i] <- beta3     + delta3 * boy[i]

    for (j in 1:M) {
      mu[i,j] <- b0_eff[i] + b1_eff[i] * t[j] + b2_eff[i] * (1 - exp(-(b3_eff[i] / b2_eff[i]) * t[j]))
      Y[i,j] ~ dnorm(mu[i,j], tau)
    }
  }
  
  beta0 ~ dnorm(0, 1.0E-6) T(0, )
  beta1 ~ dnorm(0, 1.0E-6) T(0, )
  beta2 ~ dnorm(0, 1.0E-6) T(0, )
  beta3 ~ dnorm(0, 1.0E-6) T(0, )
  
  delta0 ~ dnorm(0, 1.0E-6)
  delta1 ~ dnorm(0, 1.0E-6)
  delta2 ~ dnorm(0, 1.0E-6) T(-beta2, )
  delta3 ~ dnorm(0, 1.0E-6) T(-beta3, )
  
  tau    ~ dgamma(0.001, 0.001)
  tau0   ~ dgamma(0.001, 0.001)
  sigma  <- 1 / sqrt(tau)
  sigma0 <- 1 / sqrt(tau0)
}
"

inits_q6 <- list(
  list(beta0=2.84, beta1=0.00590, beta2=4.31, beta3=0.03330, delta0=0.11, delta1=0.0, delta2=0.97, delta3=0.004, tau=10.0, tau0=17.0),
  list(beta0=2.80, beta1=0.00550, beta2=4.10, beta3=0.03200, delta0=0.09, delta1=0.0, delta2=0.85, delta3=0.003, tau=9.5,  tau0=15.0),
  list(beta0=2.90, beta1=0.00630, beta2=4.50, beta3=0.03450, delta0=0.13, delta1=0.0, delta2=1.10, delta3=0.005, tau=10.5, tau0=19.0)
)

set.seed(123)
model_q6 <- jags.model(textConnection(jags_model_q6_string), data = jags_data_q6, inits = inits_q6, n.chains = 3, n.adapt = 2000)
update(model_q6, 5000)

params_q6 <- c("beta0", "beta1", "beta2", "beta3", "delta0", "delta1", "delta2", "delta3", "sigma", "sigma0")
samples_q6 <- coda.samples(model_q6, variable.names = params_q6, n.iter = 20000)

cat("\n--- Q6: Mixed-Effects Performance Evaluation ---\n")
print(summary(samples_q6))
print(gelman.diag(samples_q6))
print(effectiveSize(samples_q6))

png("q6_traceplots.png", width = 1200, height = 1200)
par(mfrow = c(5, 2))
for (p in params_q6) { 
  traceplot(samples_q6[, p], main = paste("Traceplot", p)) 
}
dev.off()


# Q6(d) - Population profile projection vs subject level intercept uncertainty
chains_q6 <- as.matrix(samples_q6)
n_sims_q6 <- nrow(chains_q6)
y_pred_q6 <- numeric(n_sims_q6)

for (s in 1:n_sims_q6) {
  b0   <- chains_q6[s, "beta0"]
  b1   <- chains_q6[s, "beta1"] + chains_q6[s, "delta1"]
  b2   <- chains_q6[s, "beta2"] + chains_q6[s, "delta2"]
  b3   <- chains_q6[s, "beta3"] + chains_q6[s, "delta3"]
  d0   <- chains_q6[s, "delta0"]
  sig  <- chains_q6[s, "sigma"]
  sig0 <- chains_q6[s, "sigma0"]
  
  # Account for individual variations by including random intercept noise
  beta0i        <- rnorm(1, mean = b0 + d0, sd = sig0)
  mu_90         <- beta0i + b1 * t_pred + b2 * (1 - exp(-(b3 / b2) * t_pred))
  y_pred_q6[s]  <- rnorm(1, mean = mu_90, sd = sig)
}

cat("\n--- Q6(d): Day 90 Prediction with Individual Variation Included ---\n")
print(quantile(y_pred_q6, probs = c(0.025, 0.5, 0.975)))
cat("Interval Width Q6d:", diff(quantile(y_pred_q6, c(0.025, 0.975))), "\n")
cat("Interval Width Q5c:", diff(quantile(y_pred_q5, c(0.025, 0.975))), "\n")

png("q6d_histogram.png", width = 800, height = 500)
hist(y_pred_q6, breaks = 60, col = "steelblue", border = "white", prob = TRUE, main = "Q6d: Predictive Density (Random Boy Model)", xlab = "Weight (kg)")
lines(density(y_pred_q6), col = "black", lwd = 2)
abline(v = quantile(y_pred_q6, c(0.025, 0.975)), col = "red", lwd = 2, lty = 2)
dev.off()

# Final comparative layout generation
png("q5c_q6d_comparison.png", width = 900, height = 500)
dens_q5 <- density(y_pred_q5)
dens_q6 <- density(y_pred_q6)
plot(dens_q6, col = "darkorange", lwd = 2, main = "Day 90 Posterior Density Shift: Q5c vs Q6d", xlab = "Weight (kg)", xlim = c(4.5, 7.5), ylim = c(0, max(dens_q5$y, dens_q6$y) * 1.1))
lines(dens_q5, col = "steelblue", lwd = 2)
abline(v = quantile(y_pred_q5, c(0.025, 0.975)), col = "steelblue", lty = 2, lwd = 1.5)
abline(v = quantile(y_pred_q6, c(0.025, 0.975)), col = "darkorange", lty = 2, lwd = 1.5)
legend("topright", legend = c("Q5c Fixed-Effects Profile", "Q6d Mixed-Effects Profile"), col = c("steelblue", "darkorange"), lwd = 2, bty = "n")
dev.off()



