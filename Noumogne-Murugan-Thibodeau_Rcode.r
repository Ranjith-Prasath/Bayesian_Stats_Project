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

# Boys only
boys_data <- df[df$Boy == 1, 2:14]
Y_boys <- as.matrix(boys_data)

# Mean growth function
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

# Transformation from beta to gamma
start <- c(
  log(beta0_init),
  -log(beta1_init),
  log(beta2_init),
  -log(beta3_init),
  1 / sigma_init^2
)

# Negative log-posterior
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

# Initial values from Question 3
gamma_current <- theta_hat[1:4]
tau_current   <- theta_hat[5]

# Covariance matrix for gamma from the Laplace approximation
Sigma_gamma <- Sigma_theta[1:4, 1:4]

# Tuning parameter
# Since the proposal is multivariate, we target an acceptance rate around 20%
c_tune <- 1

# Storage
gamma_chain <- matrix(NA, nrow = M, ncol = 4)
tau_chain   <- numeric(M)

colnames(gamma_chain) <- c("gamma0", "gamma1", "gamma2", "gamma3")

accept <- 0

# Function computing S(gamma)
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
  
  # -------------------------
  # Metropolis step for gamma
  # -------------------------
  
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
  
  # -------------------------
  # Gibbs step for tau
  # -------------------------
  
  RSS_current <- compute_RSS(gamma_current, Y_boys, t.obs)
  
  n <- nrow(Y_boys)
  m <- ncol(Y_boys)
  
  shape_tau <- n * m / 2
  rate_tau  <- RSS_current / 2
  
  tau_current <- rgamma(1, shape = shape_tau, rate = rate_tau)
  
  # Store values
  gamma_chain[s, ] <- gamma_current
  tau_chain[s] <- tau_current
}

acceptance_rate <- accept / M
acceptance_rate


#QUESTION (b)(ii)

par(mfrow = c(3, 2))

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


par(mfrow = c(3, 2))

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




