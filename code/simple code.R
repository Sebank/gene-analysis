library(lme4)
library(insight)
library(performance)
library(ggplot2)

set.seed(123)

beta = matrix(c(2,0),ncol=1)

tau = 0.2

# this alpha is extracted from the DESeq2 data for the simulations to be more comparable to real data
alpha = 0.25
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.01168  0.10547  0.25233  0.57844  0.63406 56.06058

n = 200

# conditional model
# data.frame("E" = rbind(c(1, 0), c(1, 1)) %*% beta, "Var" = rbind(c(1, 0), c(1, 1)) %*% beta + alpha * (rbind(c(1, 0), c(1, 1)) %*% beta)^2)

# marginal model
marginal = data.frame("E" = exp(rbind(c(1, 0), c(1, 1)) %*% beta) * exp(tau^2/2),
           "Var" = exp(rbind(c(1, 0), c(1, 1)) %*% beta) * exp(tau^2/2) +
             alpha * exp(rbind(c(1, 0), c(1, 1)) %*% beta)^2 * exp(2*tau^2) +
             exp(rbind(c(1, 0), c(1, 1)) %*% beta)^2 * exp(tau^2) * (exp(tau^2) - 1))

# two observations from the same tissue from the same patient (want from different tissue)
# x=cbind(rep(1,n),c(rep(0,n/2),rep(1,n/2)))

# Intercept for all and alternating effect, as that captures each pair having one observations with and without effect
x = cbind(rep(1, n), rep(0:1, n/2))

# beware that the function utilizes x, beta and tau, which are global variables
create_dataset = function(group = 100, n = 200){
  # input: 
  #   - group (int)
  #     number of grouped pairs for the simulated dataset
  #   - n (int)
  #     number of total samples (needs to be larger than 2*group)
  # output:
  #   - y (list, float)
  #     simulated negative binomial data
  #
  # Function that creates a dataset based on a negative binomial GLMM
  
  pers = rep(1:group, each = 2)
  gamma = rep(rnorm(group, mean = 0, sd = tau), each = 2)
  
  # give all single samples their own random effect and identifier
  if(2 * group < n){
    pers = c(pers, (group + 1):(n - group))
    gamma = c(gamma, rnorm(n - 2*group, mean = 0, sd = tau))
  }
  
  eta = x %*% beta + gamma
  
  y = rnbinom(n = n, size = 1/alpha, mu = exp(eta))
  
  return(
    list("y" = y, "gamma" = gamma, "eta" = eta, "pers" = pers)
    )
}

make_conditional_plot = function(conditional){
  # input:
  #   - conditional (double)
  #     expected value for conditional model for a random person of the GLMM
  # output: 
  #   - obj (ggplot2)
  #     Bar plot of the conditional distribution mass function for a random person of the GLMM
  #
  # Function that creates a bar plot of the conditional distribiution mass function for a random person of the GLMM 
  # and plots it with an observation conditional on the effect being zero
  
  obj = ggplot() +
    geom_col(aes(x = 0:n_plot, y = dnbinom(x = 0:n_plot, size = 1/alpha, mu = exp(beta[1]))), fill = "#CC0000") +
    geom_col(aes(x = 0:n_plot, y = dnbinom(x = 0:n_plot, size = 1/alpha, mu = exp(conditional))), fill = "#000099", alpha = 0.6) +
    labs(x = "", y = "")
  return(obj)
}

# collecting relevant data to initialize the random variables before launching a function that depends on them (not good coding practice)
n=200
tau = 0.2
x = cbind(rep(1, n), rep(0:1, n/2))
beta = matrix(c(2,0),ncol=1)
result = create_dataset(group = 100, n = 200)

y = result$y
eta = result$eta
gamma = result$gamma
pers = result$pers

# Examine model itself and compare to marginal

# For the time being both 1 and 2 are the same, so we only consider 1
# Due to the variance being so large it is hard to get a close estimate to the theoretical value (particularly with high alpha or mu)
mean(y)/marginal$E[1]
var(y)/marginal$Var[1]

# define endpoint of probability mass function
n_plot = 35

q = sample(eta, 9)
#
a1 = make_conditional_plot(q[1])
a2 = make_conditional_plot(q[2])
a3 = make_conditional_plot(q[3])
a4 = make_conditional_plot(q[4])
a5 = make_conditional_plot(q[5])
a6 = make_conditional_plot(q[6])
a7 = make_conditional_plot(q[7])
a8 = make_conditional_plot(q[8])
a9 = make_conditional_plot(q[9])
library(gridExtra)
grid.arrange(a1, a2, a3, a4, a5, a6, a7, a8, a9, nrow = 3)

# ggplot() + 
#   geom_col(aes(x = 0:n_plot, y = dnbinom(x = 0:n_plot, size = 1/alpha, mu = exp(beta[1])), fill = "GLM equivalent")) +
#   geom_col(aes(x = 0:n_plot, y = dnbinom(x = 0:n_plot, size = 1/alpha, mu = exp(q[i])), fill = "Conditional"), alpha = 0.6) +
#   labs(x = "x", y = "y", title = "Probability mass function")
# # 3x3 plot of conditional

# plot realization
pos = max(table(y))/length(y)*1.1
ggplot() + 
#  geom_col(aes(x = 0:n_plot, y = dnbinom(x = 0:n_plot, size = 1/alpha, mu = exp(beta[1])), fill = "GLM equivalent")) +
  geom_bar(aes(x = y, y = ..count../sum(..count..), fill = "realization"), alpha = 0.7) + 
  geom_segment(aes(x = mean(y) - sqrt(var(y)), xend = mean(y) + sqrt(var(y)), y = pos, yend = pos)) + 
  geom_point(aes(x = mean(y), y = pos)) + 
  annotate("label", x = mean(y), y = pos + 0.01, label = "mean +/- sd")
  labs(x = "x", y = "y", title = "observed draws from negbin GLMM")

fit = glm(y ~ x-1, family=MASS::negative.binomial(link = "log", theta = 1/alpha))

# summary(fit, dispersion=1)
# 
# vcov(fit, dispersion = 1)
# 
# sqrt(diag(vcov(summary(fit, dispersion = 1))))

fit2 = glmer(formula = y ~ x-1 + (1|pers), family = MASS::negative.binomial(link = "log", theta = 1/alpha))

# summary(fit2)
# 
# vcov(fit2)
# 
# sqrt(diag(vcov(summary(fit2))))
# 
# icc(fit2)

# potential fail check
summary(fit2)$optinfo$conv$lme4$code

m = 10000
results = data.frame("intercept" = rep(NA, m), "effect" = rep(NA, m), 
                     "intercept.se" = rep(NA, m), "effect.se" = rep(NA, m))
results_GLMM = data.frame("intercept" = rep(NA, m), "effect" = rep(NA, m), 
                          "intercept.se" = rep(NA, m), "effect.se" = rep(NA, m), 
                          "ICC" = rep(NA, m), "ICC_conditional" = rep(NA, m))
n=200
tau = 1
x = cbind(rep(1, n), rep(0:1, n/2))
beta = matrix(c(2,0),ncol=1)
pers = create_dataset(group = 100, n = 200)$pers

# simulate different datasets and fit GLM and GLMM models to them, then store relevant values for comparison
for(i in 1:m){
  tryCatch(
    expr = {
      if(i %% (m/100) == 0) print(paste("Progression:", 100*i/m, "%"))
      temp = create_dataset(group = 100, n = 200)$y
      temp_fit = glm(temp ~ x-1,family=MASS::negative.binomial(link = "log", theta = 1/alpha))
      temp_GLMM = glmer(formula = temp ~ x-1 + (1|pers), family = MASS::negative.binomial(link = "log", theta = 1/alpha))
      
      results$intercept[i] = coef(summary(temp_fit, dispersion = 1))[1, 1]
      results$intercept.se[i] = coef(summary(temp_fit, dispersion = 1))[1, 2]
      results$effect[i] = coef(summary(temp_fit, dispersion = 1))[2, 1]
      results$effect.se[i] = coef(summary(temp_fit, dispersion = 1))[2, 2]
      
      results_GLMM$intercept[i] = coef(summary(temp_GLMM))[1, 1]
      results_GLMM$intercept.se[i] = coef(summary(temp_GLMM))[1, 2]
      results_GLMM$effect[i] = coef(summary(temp_GLMM))[2, 1]
      results_GLMM$effect.se[i] = coef(summary(temp_GLMM))[2, 2]
      
      temp_icc = icc(temp_GLMM)
      results_GLMM$ICC[i] = temp_icc$ICC_adjusted
      results_GLMM$ICC_conditional[i] = temp_icc$ICC_conditional
    }, 
    warning = function(cond){
    }
  )
}

# plot counts of P values for GLMM and GLM simulations
ggplot() + 
  geom_histogram(aes(x = 2 * pnorm(abs(results$effect/results$effect.se), lower.tail = FALSE), 
                     y = ..count../sum(..count..), fill = "GLM"), breaks = 0:25/25) + 
  geom_histogram(aes(x = 2 * pnorm(abs(results_GLMM$effect/results_GLMM$effect.se), lower.tail = FALSE), 
                     y = ..count../sum(..count..), fill = "GLMM"), alpha = 0.5, breaks = 0:25/25) +
  labs(x = "P value for models", y = "proportion", title = "counts of P values for GLMM and GLM simulations")

# plot distribution of ICC when simulating a correlated GLMM
ggplot() + 
  geom_histogram(aes(x = results_GLMM$ICC, y = ..count../sum(..count..), fill = "adjusted ICC"), breaks = 0:100/100) + 
  geom_histogram(aes(x = results_GLMM$ICC_conditional, y = ..count../sum(..count..), fill = "conditional ICC"), alpha = 0.5, breaks = 0:100/100) + 
  labs(x = "estimated ICC", title = "distribution of ICC when simulating a correlated GLMM")

# plot comparison of p values for ICC
ggplot() + 
  geom_point(aes(x = 2 * pnorm(abs(results$effect/results$effect.se), lower.tail = FALSE), 
                 y = 2 * pnorm(abs(results_GLMM$effect/results_GLMM$effect.se), lower.tail = FALSE)), alpha = 0.5) +
  labs(x = "P value GLM", y = "P value GLMM", title = "comparison of p values for ICC")

# plot difference in P value between GLMM against GLM as reference
ggplot() +
  geom_point(aes(x = results_GLMM$ICC,
                 y = 2 * pnorm(abs(results_GLMM$effect/results_GLMM$effect.se), lower.tail = FALSE) -
                 2 * pnorm(abs(results$effect/results$effect.se), lower.tail = FALSE) ), alpha = 0.5) +
  labs(x = "ICC adjusted", y = "difference in P value between GLMM against GLM as reference")

remove.na = !is.na(results_GLMM$intercept)
y = list("GLMM" = rep(0, 101), "GLM" = rep(0, 101))
for(i in 1:101){
  y$GLMM[i] = sum(2 * pnorm(abs(results_GLMM$effect/results_GLMM$effect.se)[remove.na], 
                         lower.tail = FALSE) <= (i-1)/100)/length(results_GLMM$intercept[remove.na])
  y$GLM[i] =  sum(2 * pnorm(abs(results$effect/results$effect.se)[remove.na], 
                         lower.tail = FALSE) <= (i-1)/100)/length(results$intercept[remove.na])
}
ggplot() +
  geom_line(aes(x = c(0, 1), y = c(0, 1), col = "truth")) + 
  geom_line(aes(x = 0:100/100, 
                y = y$GLMM, col = "GLMM")) +
  geom_line(aes(x = 0:100/100, 
                y = y$GLM, col = "GLM")) +
  labs(x = "p values", y = "proportion", title = "cumulative p values")
