library(lme4)
library(insight)
library(performance)
library(ggplot2)
library(MASS)
library(ggpubr)

set.seed(123)

# define global variables once at top
beta = matrix(c(3,0),ncol=1)
tau = 0.4

n = 200

# Intercept for all and alternating effect, as that captures each pair having one observations with and without effect
x = cbind(rep(1, n), rep(0:1, n/2))

# this alpha is extracted from the DESeq2 data for the simulations to be more comparable to real data
alpha = 0.25
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.01168  0.10547  0.25233  0.57844  0.63406 56.06058

# define endpoint of probability mass function
n_plot = 100

group = 20

# conditional model
# data.frame("E" = rbind(c(1, 0), c(1, 1)) %*% beta, "Var" = rbind(c(1, 0), c(1, 1)) %*% beta + alpha * (rbind(c(1, 0), c(1, 1)) %*% beta)^2)

# marginal model
marginal = data.frame("E" = exp(rbind(c(1, 0), c(1, 1)) %*% beta) * exp(tau^2/2),
           "Var" = exp(rbind(c(1, 0), c(1, 1)) %*% beta) * exp(tau^2/2) +
             alpha * exp(rbind(c(1, 0), c(1, 1)) %*% beta)^2 * exp(2*tau^2) +
             exp(rbind(c(1, 0), c(1, 1)) %*% beta)^2 * exp(tau^2) * (exp(tau^2) - 1))

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
  
  if(group != 0){
    pers = rep(1:group, each = 2)
    gamma = rep(rnorm(group, mean = 0, sd = tau), each = 2)
  }else{
    pers = c()
    gamma = c()
  }
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

correlation = function(mu, alpha, tau){
  # input:
  #   - mu (double)
  #     expected value negative binomial
  #   - alpha (double)
  #     dispersion parameter in the negative binomial
  #   - tau (double)
  #     standard deviation of random effects, or estimate of
  # output: 
  #   - corr (double)
  #
  # Function that estimates the correlation between two observatios assumed to follow a negative binomial GLMM
  
  # covariance between two observations from the negative binomial
  corr = prod(mu)*exp(tau^2)*(exp(tau^2) - 1)
  # divide it by the root of the variance for one observation
  corr = corr/sqrt(mu[1] + exp(tau^2/2) +                # poisson part
                     alpha*mu[1]^2*exp(2*tau^2) +        # negative binomial part
                     mu[1]*exp(tau^2)*(exp(tau^2) - 1))  # mixed effect part
  # and the same for the other observation
  corr = corr/sqrt(mu[2] + exp(tau^2/2) +                # poisson part
                     alpha*mu[2]^2*exp(2*tau^2) +        # negative binomial part
                     mu[2]*exp(tau^2)*(exp(tau^2) - 1))  # mixed effect part
  return(corr)
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
    geom_col(aes(x = 0:n_plot, y = dnbinom(x = 0:n_plot, size = 1/alpha, mu = exp(beta[1])), fill = "gamma = 0")) +
    geom_col(aes(x = 0:n_plot, y = dnbinom(x = 0:n_plot, size = 1/alpha, mu = exp(conditional)), fill = "gamma != 0"), alpha = 0.6) +
    annotate("label", x = n_plot - n_plot/4, y = 0.02, label = paste("gamma =", round(conditional - beta[1], 3))) +
    labs(x = "", y = "", fill ="") + xlim(-0.5, n_plot) + ylim(0, 0.07)
  return(obj)
}

# # plot a single realization
# pos = max(table(y))/length(y)*1.1
# ggplot() + 
#   #  geom_col(aes(x = 0:n_plot, y = dnbinom(x = 0:n_plot, size = 1/alpha, mu = exp(beta[1])), fill = "GLM equivalent")) +
#   geom_bar(aes(x = y, y = ..count../sum(..count..), fill = "realization"), alpha = 0.7) + 
#   geom_segment(aes(x = mean(y) - sqrt(var(y)), xend = mean(y) + sqrt(var(y)), y = pos, yend = pos)) + 
#   geom_point(aes(x = mean(y), y = pos)) + 
#   annotate("label", x = mean(y), y = pos + 0.01, label = "mean +/- sd") + 
#   labs(x = "x", y = "y", title = "observed draws from negbin GLMM")

# Plot marginal distribution
marginal_distribution = c()
for(i in 1:100){
  marginal_distribution = c(marginal_distribution, create_dataset(group = 0)$y)
}
y_pos = max(c(max(table(marginal_distribution))/length(marginal_distribution), dnbinom(x = 0:n_plot, size = 1/alpha, mu = exp(beta[1]))))*1.1
ggplot() + 
  geom_histogram(aes(x = marginal_distribution, y = ..count../sum(..count..), fill = "GLMM", col = "GLMM"), breaks = (-1/2):(n_plot + 1/2)) +
  geom_col(aes(x = 0:n_plot, y = dnbinom(x = 0:n_plot, size = 1/alpha, mu = exp(beta[1])), fill = "GLM", col = "GLM"), alpha = 0.4) +
  geom_segment(aes(x = mean(marginal_distribution) - sd(marginal_distribution), 
                   xend = mean(marginal_distribution) + sd(marginal_distribution), 
                   y = y_pos, yend = y_pos, col = "GLMM")) + 
  geom_point(aes(x = mean(marginal_distribution), y = y_pos, col = "GLMM")) +
  geom_segment(aes(x = exp(beta[1]) - sqrt(exp(beta[1]) + alpha*exp(2*beta[1])), 
                                           xend = exp(beta[1]) + sqrt(exp(beta[1]) + alpha*exp(2*beta[1])),
                   y = y_pos-0.001, yend = y_pos-0.001, col = "GLM")) +
  geom_point(aes(x = exp(beta[1]), y = y_pos - 0.001, col = "GLM")) +
  annotate("label", x = mean(marginal_distribution), y = y_pos + 0.005, label = "mean +/- sd") +
  xlim(-0.5, n_plot) + guides(col = "none") + 
  labs(x = "x", y = "Proportion", fill = "")
cat(paste(" \t\t E\t\t\t Var\n Theoretical:\t", marginal$E[1], "\t", sqrt(marginal$Var[1]), "\n Estimated:\t", mean(marginal_distribution), "\t\t", sd(marginal_distribution)))


# collecting relevant data to initialize the random variables before launching a function that depends on them (not good coding practice)
result = create_dataset(group = group, n = n)

data_frame = data.frame("y" = result$y, "x" = x, "pers" = result$pers)
y = result$y
eta = result$eta
gamma = result$gamma
pers = result$pers

# Examine model itself and compare to marginal

# For the time being both 1 and 2 are the same, so we only consider 1
# Due to the variance being so large it is hard to get a close estimate to the theoretical value (particularly with high alpha or mu)
mean(y)/marginal$E[1]
var(y)/marginal$Var[1]

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

ggarrange(a1, a2, a3, a4, a5, a6, a7, a8, a9, ncol = 3, nrow = 3, common.legend = TRUE)

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
  annotate("label", x = mean(y), y = pos + 0.01, label = "mean +/- sd") + 
  labs(x = "x", y = "y", title = "observed draws from negbin GLMM")

fit = glm.nb(y ~ x-1, link = "log", data = data_frame)

# summary(fit, dispersion=1)
# 
# vcov(fit, dispersion = 1)
# 
# sqrt(diag(vcov(summary(fit, dispersion = 1))))

fit2 = glmer.nb(formula = y ~ x-1 + (1|pers), data = data_frame)

# summary(fit2)
# 
# vcov(fit2)
# 
# sqrt(diag(vcov(summary(fit2))))
# 
# icc(fit2)

# potential fail check
summary(fit2)$optinfo$conv$lme4$code

# check correlation for specific set of settings
correlation(
  mu = exp(coef(summary(fit2))[1, 1]) * c(1, exp(coef(summary(fit2))[2, 1])), 
  alpha = 1/getME(fit2, "glmer.nb.theta"),
  tau = fit2@pp$theta
)

m = 10000
results = data.frame("intercept" = rep(NA, m), "effect" = rep(NA, m), 
                     "intercept.se" = rep(NA, m), "effect.se" = rep(NA, m), "alpha" = rep(NA, m))
results_GLMM = data.frame("intercept" = rep(NA, m), "effect" = rep(NA, m), 
                          "intercept.se" = rep(NA, m), "effect.se" = rep(NA, m), 
                          "ICC" = rep(NA, m), "ICC_conditional" = rep(NA, m), "alpha" = rep(NA, m), "corre" = rep(NA, m), "tau" = rep(NA, m))

time = Sys.time()

# simulate different datasets and fit GLM and GLMM models to them, then store relevant values for comparison
for(i in 1:m){
  tryCatch(
    expr = {
      if(i %% (m/100) == 0){
        print(paste("Progression:", 100*i/m, "%."), sep = "")
        print(Sys.time() - time)
        }
      temp = create_dataset(group = group, n = n)
      temp_fit = glm.nb(y ~ x-1, link = "log", data = temp)
      temp_GLMM = glmer.nb(formula = y ~ x-1 + (1|pers), data = temp)
      
      coef_summmary_temp_fit = coef(summary(temp_fit, dispersion = 1))
      results$intercept[i] = coef_summmary_temp_fit[1, 1]
      results$intercept.se[i] = coef_summmary_temp_fit[1, 2]
      results$effect[i] = coef_summmary_temp_fit[2, 1]
      results$effect.se[i] = coef_summmary_temp_fit[2, 2]
      results$alpha[i] = 1/temp_fit$theta
      
      coef_summmary_temp_GLMM = coef(summary(temp_GLMM))
      results_GLMM$intercept[i] = coef_summmary_temp_GLMM[1, 1]
      results_GLMM$intercept.se[i] = coef_summmary_temp_GLMM[1, 2]
      results_GLMM$effect[i] = coef_summmary_temp_GLMM[2, 1]
      results_GLMM$effect.se[i] = coef_summmary_temp_GLMM[2, 2]
      results_GLMM$alpha[i] = 1/getME(temp_GLMM, "glmer.nb.theta")
      results_GLMM$tau[i] = temp_GLMM@pp$theta
      
      results_GLMM$corre[i] = correlation(
        mu = exp(results_GLMM$intercept[i]) * c(1, exp(results$effect[i])), 
        alpha = results_GLMM$alpha[i],
        tau = results_GLMM$tau[i]
      )
      # temp_icc = icc(temp_GLMM)
      # results_GLMM$ICC[i] = temp_icc$ICC_adjusted
      # results_GLMM$ICC_conditional[i] = temp_icc$ICC_conditional
    }, 
    warning = function(cond){
    }, 
    error = function(cond){
    }
  )
}

# equivalent to pchisq((results$effect/results$effect.se)^2, df = 1, lower.tail = FALSE)
P = list("GLM" = 2 * pnorm(abs(results$effect/results$effect.se), lower.tail = FALSE), 
         "GLM.adj" = c(), 
         "GLMM" = 2 * pnorm(abs(results_GLMM$effect/results_GLMM$effect.se), lower.tail = FALSE), 
         "GLMM.adj" = c())

P$GLM.adj = sort(P$GLM[!is.na(P$GLM)])
P$GLMM.adj = sort(P$GLMM[!is.na(P$GLMM)])

for(i in (length(P$GLM.adj) - 1):1){
  P$GLM.adj[i] = min(P$GLM.adj[i + 1], P$GLM.adj[i]*length(P$GLM.adj)/i)
  P$GLMM.adj[i] = min(P$GLMM.adj[i + 1], P$GLMM.adj[i]*length(P$GLMM.adj)/i)
}

# plot counts of P values for GLMM and GLM simulations
ggplot() + 
  geom_histogram(aes(x = P$GLM, 
                     y = ..count../sum(..count..), fill = "GLM"), breaks = 0:25/25) + 
  geom_histogram(aes(x = P$GLMM, 
                     y = ..count../sum(..count..), fill = "GLMM"), alpha = 0.5, breaks = 0:25/25) +
  labs(x = "P value for models", y = "proportion", title = "counts of P values for GLMM and GLM simulations", fill = "Type of model") +
  xlim(0, 1) + ylim(0, 0.12)
ggsave(paste("P count tau = ", tau, ", group = ", group, ", number of simuations = ", m, ".pdf"), width = 15, height = 10, units = "cm")

# plot distribution of ICC when simulating a correlated GLMM
ggplot() + 
  geom_histogram(aes(x = results_GLMM$corre, y = ..count../sum(..count..)), breaks = 0:100/100) + 
  labs(x = "estimated correlation", title = "distribution of correlation when simulating a correlated GLMM", y = "proportion", fill = "") + 
  xlim(0, 1) + ylim(0, 0.07)
ggsave(paste("correlation tau = ", tau, ", group = ", group, ", number of simuations = ", m, ".pdf"), width = 15, height = 10, units = "cm")

# plot comparison of p values for ICC
ggplot() + 
  geom_point(aes(x = P$GLM, 
                 y = P$GLMM), alpha = 0.1) +
  labs(x = "P value GLM", y = "P value GLMM", title = "comparison of p values for correlated data") +
  xlim(0, 0.07) + ylim(0, 0.07)
ggsave(paste("P compare restricted tau = ", tau, ", group = ", group, ", number of simuations = ", m, ".pdf"), width = 15, height = 15, units = "cm")

# plot comparison of p values for ICC
ggplot() + 
  geom_point(aes(x = P$GLM, 
                 y = P$GLMM), alpha = 0.1) +
  labs(x = "P value GLM", y = "P value GLMM", title = "comparison of p values for correlated data") +
  xlim(0, 1) + ylim(0, 1)
ggsave(paste("P compare tau = ", tau, ", group = ", group, ", number of simuations = ", m, ".pdf"), width = 15, height = 15, units = "cm")

# plot difference in P value between GLMM against GLM as reference
ggplot() +
  geom_point(aes(x = results_GLMM$corre,
                 y = P$GLMM -
                   P$GLM ), alpha = 0.06) +
  labs(x = "estimated correlation", y = "difference in P value of GLMM subtracted by GLM") +
  xlim(0, 1) + ylim(-1, 1)
ggsave(paste("P difference tau = ", tau, ", group = ", group, ", number of simuations = ", m, ".pdf"), width = 15, height = 10, units = "cm")

# Does this do anything, adjused P values
remove.na = !is.na(results_GLMM$intercept)
x_cumulative = list("GLM" = c(0, sort(P$GLM[remove.na])), "GLMM" = c(0, sort(P$GLMM[remove.na])))
y = list("GLM" = (0:(length(x_cumulative$GLM) - 1))/(length(x_cumulative$GLM) - 1), "GLMM" = (0:(length(x_cumulative$GLMM) - 1))/(length(x_cumulative$GLM) - 1))
# y = list("GLMM" = rep(0, 101), "GLM" = rep(0, 101))
#  for(i in 1:legnth(x_cumulative)){
#   y$GLMM[i] = sum(P$GLMM[remove.na] <= (i-1)/100)/length(results_GLMM$intercept[remove.na])
#   y$GLM[i] =  sum(P$GLM[remove.na] <= (i-1)/100)/length(results$intercept[remove.na])
# }
ggplot() +
  geom_line(aes(x = c(0, 1), y = c(0, 1), col = "exact")) + 
  geom_line(aes(x = x_cumulative$GLMM, 
                y = y$GLMM, col = "GLMM")) +
  geom_line(aes(x = x_cumulative$GLM, 
                y = y$GLM, col = "GLM")) +
  labs(x = "p values", y = "cumulative proportion", title = "cumulative p values", col = "")
ggsave(paste("P cumulative tau = ", tau, ", group = ", group, ", number of simuations = ", m, ".pdf"), width = 15, height = 15, units = "cm")

# distribution of alpha
ggplot() + geom_histogram(aes(x = results$alpha, fill = "GLM", y = ..count../sum(..count..)), alpha = 0.8) + 
  geom_histogram(aes(x = results_GLMM$alpha, fill = "GLMM", y = ..count../sum(..count..)), alpha = 0.4) + 
  geom_vline(xintercept = mean(results$alpha[!is.na(results$alpha)])) + 
  geom_vline(xintercept = mean(results_GLMM$alpha[!is.na(results$alpha)])) + 
  labs(x = "alpha", y = "proportion", fill = "", title = "Distribution of alpha for estimated GLM and GLMM models")
ggsave(paste("alpha distribution tau = ", tau, ", group = ", group, ", number of simuations = ", m, ".pdf"), width = 15, height = 15, units = "cm")

x_seb = list("GLM" = y$GLM, "GLMM" = y$GLMM)
y_seb = list("GLM" = c(0), "GLMM" = c(0))

for(i in 2:(length(x_seb$GLM))){
  y_seb$GLM = c(y_seb$GLM, y_seb$GLM[i - 1] + y$GLM[i]*(x_cumulative$GLM[i] - x_cumulative$GLM[i - 1]))
  y_seb$GLMM = c(y_seb$GLMM, y_seb$GLMM[i - 1] + y$GLMM[i]*(x_cumulative$GLMM[i] - x_cumulative$GLMM[i - 1]))
}

ggplot() + 
  geom_line(aes(x = x_seb$GLMM, y = y_seb$GLMM, col = "GLMM")) +
  geom_line(aes(x = x_seb$GLM, y = y_seb$GLM, col = "GLM"), alpha = 0.8) +
  lims(x = c(0, 0.05), y = c(0, 0.05^2*3/4))

ggplot() + 
  geom_line(aes(x = x_seb$GLMM, y = y_seb$GLMM, col = "GLMM")) +
  geom_line(aes(x = x_seb$GLM, y = y_seb$GLM, col = "GLM"), alpha = 0.8)

num = list("GLM" = tail(which(x_seb$GLM < 0.05), 1), "GLMM" = tail(which(x_seb$GLMM < 0.05), 1))
y_num = list("GLM" = y_seb$GLM[num$GLM] + y$GLM[num$GLM + 1]*(0.05 - x_cumulative$GLM[num$GLM]), 
             "GLMM" = y_seb$GLM[num$GLMM] + y$GLM[num$GLMM + 1]*(0.05 - x_cumulative$GLMM[num$GLMM]))

compare.power = list("GLM.seb" = tail(y_seb$GLM, 1), "GLMM.seb" = tail(y_seb$GLMM, 1), 
     "GLM.seb.0.05" = y_num$GLM,"GLMM.seb.0.05" = y_num$GLMM, 
     "GLM.ecdf" = 1 - x_cumulative$GLM[head(which(y$GLM > 0.05), 1)], 
     "GLMM.ecdf" = 1 - x_cumulative$GLMM[head(which(y$GLMM > 0.05), 1)])
compare.power = read.csv("D:\\Essential folders\\Documents\\RStudio\\master\\gene-analysis\\code\\compare.power.csv", header = TRUE, row.names = 1)
compare.power = rbind(compare.power, list("GLM.seb" = tail(y_seb$GLM, 1), "GLMM.seb" = tail(y_seb$GLMM, 1), 
                          "GLM.seb.0.05" = y_num$GLM,"GLMM.seb.0.05" = y_num$GLMM, 
                          "GLM.ecdf" = 1 - x_cumulative$GLM[head(which(y$GLM > 0.05), 1)], 
                          "GLMM.ecdf" = 1 - x_cumulative$GLMM[head(which(y$GLM > 0.05), 1)]))
rownames(compare.power) = c(head(rownames(compare.power), -1), paste("tau", tau, "m*", 2*group/n))
write.csv(compare.power, "D:\\Essential folders\\Documents\\RStudio\\master\\gene-analysis\\code\\compare.power.csv")
