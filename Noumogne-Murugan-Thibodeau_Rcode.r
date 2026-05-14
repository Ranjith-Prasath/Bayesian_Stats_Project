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
