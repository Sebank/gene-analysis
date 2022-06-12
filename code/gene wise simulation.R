library(lme4)
library(insight)
library(performance)
library(ggplot2)
library(MASS)
library(ggpubr)

set.seed(123)

# define global variables once at top
beta = matrix(c(5.3, 2),ncol=1)
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

group = 100

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
  # Function that estimates the correlation between two observations assumed to follow a negative binomial GLMM
  
  # covariance between two observations from the negative binomial
  corr = prod(mu)*exp(tau^2)*(exp(tau^2) - 1)
  # divide it by the root of the variance for one observation
  corr = corr/sqrt(mu[1]*exp(tau^2/2) +                # poisson part
                     alpha*mu[1]^2*exp(2*tau^2) +        # negative binomial part
                     mu[1]^2*exp(tau^2)*(exp(tau^2) - 1))  # mixed effect part
  # and the same for the other observation
  corr = corr/sqrt(mu[2]*exp(tau^2/2) +                # poisson part
                     alpha*mu[2]^2*exp(2*tau^2) +        # negative binomial part
                     mu[2]^2*exp(tau^2)*(exp(tau^2) - 1))  # mixed effect part
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
    labs(x = "", y = "", fill ="") + xlim(-0.5, n_plot) + ylim(0, 0.06) + theme_bw()
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
  labs(x = "x", y = "Proportion", fill = "") + 
  theme_bw()
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
# #
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
# pos = max(table(y))/length(y)*1.1
# ggplot() + 
# #  geom_col(aes(x = 0:n_plot, y = dnbinom(x = 0:n_plot, size = 1/alpha, mu = exp(beta[1])), fill = "GLM equivalent")) +
#   geom_bar(aes(x = y, y = ..count../sum(..count..), fill = "realization"), alpha = 0.7) + 
#   geom_segment(aes(x = mean(y) - sqrt(var(y)), xend = mean(y) + sqrt(var(y)), y = pos, yend = pos)) + 
#   geom_point(aes(x = mean(y), y = pos)) + 
#   annotate("label", x = mean(y), y = pos + 0.01, label = "mean +/- sd") + 
#   labs(x = "x", y = "y", title = "observed draws from negbin GLMM")

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

m = 1000
results = data.frame("intercept" = rep(NA, m), "effect" = rep(NA, m), 
                     "intercept.se" = rep(NA, m), "effect.se" = rep(NA, m), "alpha" = rep(NA, m))
results_GLMM = data.frame("intercept" = rep(NA, m), "effect" = rep(NA, m), 
                          "intercept.se" = rep(NA, m), "effect.se" = rep(NA, m), 
                          "alpha" = rep(NA, m), "corre" = rep(NA, m), "tau" = rep(NA, m))

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

# equivalent to 2 * pnorm(abs(results_GLMM$effect/results_GLMM$effect.se), lower.tail = FALSE)
P = list("GLM" = pchisq((results$effect/results$effect.se)^2, df = 1, lower.tail = FALSE), 
         "GLM.adj" = c(), 
         "GLMM" = 2 * pnorm(abs(results_GLMM$effect/results_GLMM$effect.se), lower.tail = FALSE), 
         "GLMM.adj" = c()
    )

write.csv(list("GLM" = results, "GLMM" = results_GLMM), 
          paste("D:\\Essential folders\\Documents\\RStudio\\master\\gene-analysis\\code\\csv\\results tau", 
                tau, "group", group, "beta", beta[2], ".csv"))
# read.csv("D:\\Essential folders\\Documents\\RStudio\\master\\gene-analysis\\code\\csv\\results tau  0.15 group 100 beta 0 .csv", header = TRUE, row.names = 1)

m = 1000
results = list("tau" = rep(c(0.2, 0.3, 0.4), 3),
            "group" = rep(c(30, 60, 100), each = 3), 
            "GLM" = list("intercept" = array(NA, c(m, 9)), "effect" = array(NA, c(m, 9)), 
                         "intercept.se" = array(NA, c(m, 9)), "effect.se" = array(NA, c(m, 9)), 
                         "alpha" = array(NA, c(m, 9)), "P" = array(NA, c(m, 9))
              ),
            "GLMM" = list("intercept" = array(NA, c(m, 9)), "effect" = array(NA, c(m, 9)), 
                          "intercept.se" = array(NA, c(m, 9)), "effect.se" = array(NA, c(m, 9)), 
                          "alpha" = array(NA, c(m, 9)), "corre" = array(NA, c(m, 9)), 
                          "tau" = array(NA, c(m, 9)), "P" = array(NA, c(m, 9))
              )
            )

# make data easier accessible
# read data from 9 different simulations, all need to have been run before this can be run
for(i in 1:length(results$group)){
  resultsorary = read.csv(paste("D:\\Essential folders\\Documents\\RStudio\\master\\gene-analysis\\code\\csv\\results tau", 
                                   results$tau[i] , "group", results$group[i], ".csv"), header = TRUE, row.names = 1)
  results$GLM$intercept[, i] = resultsorary$GLM.intercept
  results$GLM$effect[, i] = resultsorary$GLM.effect
  results$GLM$intercept.se[, i] = resultsorary$GLM.intercept.se
  results$GLM$effect.se[, i] = resultsorary$GLM.effect.se
  results$GLM$alpha[, i] = resultsorary$GLM.alpha
  
  results$GLMM$intercept[, i] = resultsorary$GLMM.intercept
  results$GLMM$effect[, i] = resultsorary$GLMM.effect
  results$GLMM$intercept.se[, i] = resultsorary$GLMM.intercept.se
  results$GLMM$effect.se[, i] = resultsorary$GLMM.effect.se
  results$GLMM$alpha[, i] = resultsorary$GLMM.alpha
  results$GLMM$corre[, i] = resultsorary$GLMM.corre
  results$GLMM$tau[, i] = resultsorary$GLMM.tau
}

results$GLM$P = pchisq((results$GLM$effect/results$GLM$effect.se)^2, df = 1, lower.tail = FALSE)
results$GLMM$P = pchisq((results$GLMM$effect/results$GLMM$effect.se)^2, df = 1, lower.tail = FALSE)

make_cumulative_plot = function(results, i){
  # input:
  #   - results (list)
  #     list of lists with specific organization
  #   - i (int)
  #     integer informing of what plot that should be plotted and informs of which observation that should be plotted
  # output: 
  #   - obj (ggplot2)
  #     Bar plot of the conditional distribution mass function for a random person of the GLMM
  #
  # Function that creates a bar plot of the conditional distribution mass function for a random person of the GLMM 
  # and plots it with an observation conditional on the effect being zero
  
  GLM.P = results$GLM$P[, i]
  GLMM.P = results$GLMM$P[, i]
  
  GLM.P = GLM.P[!is.na(GLM.P)]
  GLMM.P = GLMM.P[!is.na(GLMM.P)]
  GLM.P = sort(GLM.P)
  GLMM.P = sort(GLMM.P)
  
  
  #mess with theme to remove whitespace between figures
  obj = ggplot() +
    geom_line(aes(x = c(0, 1), y = c(0, 1), col = "exact")) +
    geom_line(aes(x = GLM.P, 
                  y = (1:length(GLM.P))/length(GLM.P), col = "GLM")) + 
    geom_line(aes(x = GLMM.P, 
                  y = (1:length(GLMM.P))/length(GLMM.P), col = "GLMM")) +
    #labs(x = "p values", y = "cumulative proportion", title = "cumulative p values", col = "") +
    labs(x = "", y = "", col = "Model") + 
    annotate("label", x = 0.25, y = 0.75, 
             label = paste("cor", round(correlation(mu = rep(exp(beta[1]), 2), 
                                              alpha = alpha, tau = results$tau[i]), 2), 
                           "m*", results$group[i]/100)) + 
    # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    #       plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
    #       panel.background = element_rect(fill = "white"),
    #       panel.grid.major = element_line(colour = "grey"),
    #       panel.grid.minor = element_line(colour = "grey", size = 0.25)) +
    # theme_bw() + 
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1))
  # make border plots include x and y notations
  # should use cases instead, this is slower
  if(i == 1 || i == 4){
    obj = obj + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey"),
          panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else if(i == 8 || i == 9){
    obj = obj + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(colour = "grey"),
                panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else if(i == 7){
    obj = obj + theme(plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(colour = "grey"),
                panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else{
    obj = obj +  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                          plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                          panel.background = element_rect(fill = "white"),
                          panel.grid.major = element_line(colour = "grey"),
                          panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }
  return(obj)
}

b1 = make_cumulative_plot(results, 1)
b2 = make_cumulative_plot(results, 2)
b3 = make_cumulative_plot(results, 3)
b4 = make_cumulative_plot(results, 4)
b5 = make_cumulative_plot(results, 5)
b6 = make_cumulative_plot(results, 6)
b7 = make_cumulative_plot(results, 7)
b8 = make_cumulative_plot(results, 8)
b9 = make_cumulative_plot(results, 9)

ggarrange(b1, b2, b3, b4, b5, b6, b7, b8, b9, ncol = 3, nrow = 3, common.legend = TRUE, align = "hv")

# # exploring facet wrap (would recommend usage of data frame, and has large borders)
# ggplot() +
#   scale_x_continuous(limits = c(0, 1), breaks = c(0, 1)) + 
#   scale_y_continuous(limits = c(0, 1), breaks = c(0, 1)) + 
#   facet_wrap(~ rep(results$tau, each = m) * rep(results$group, each = m), nrow = 3) + theme_bw()


make_effect_plot = function(results, i){
  # input:
  #   - results (list)
  #     list of lists with specific organization
  #   - i (int)
  #     integer informing of what plot that should be plotted and informs of which observation that should be plotted
  # output: 
  #   - obj (ggplot2)
  #     Bar plot of the conditional distribution mass function for a random person of the GLMM
  #
  # Function that creates a bar plot of the conditional distribution mass function for a random person of the GLMM 
  # and plots it with an observation conditional on the effect being zero
  
  
  
  GLM.effect = results$GLM$effect[, i]
  GLMM.effect = results$GLMM$effect[, i]
  
  GLM.effect = GLM.effect[!is.na(GLM.effect)]
  GLMM.effect = GLMM.effect[!is.na(GLMM.effect)]
  
  
  #mess with theme to remove whitespace between figures
  obj = ggplot() +
    geom_histogram(aes(x = GLM.effect, y = ..count../sum(..count..), fill = "GLM"), binwidth = 0.01) + 
    geom_histogram(aes(x = GLMM.effect, y = ..count../sum(..count..), fill = "GLMM"), alpha = 0.5, binwidth = 0.01) + 
    geom_vline(aes(xintercept = mean(GLM.effect), col = "GLM")) + geom_vline(aes(xintercept = mean(GLMM.effect), col = "GLMM")) +
    #labs(x = "p values", y = "cumulative proportion", title = "cumulative p values", col = "") +
    labs(x = "", y = "", fill = "Model") + guides(col = "none") +
    annotate("label", x = -0.2, y = 0.06,
             label = paste("cor", round(correlation(mu = rep(exp(beta[1]), 2),
                                                    alpha = alpha, tau = results$tau[i]), 2),
                           "m*", results$group[i]/100)) +
    # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    #       plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
    #       panel.background = element_rect(fill = "white"),
    #       panel.grid.major = element_line(colour = "grey"),
    #       panel.grid.minor = element_line(colour = "grey", size = 0.25)) +
    # theme_bw() +
    scale_x_continuous(limits = c(-0.3, 0.3), breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3)) +
    scale_y_continuous(limits = c(0, 0.07), breaks = c(0, 0.07))
  # make border plots include x and y notations
  # should use cases instead, this is slower
  if(i == 1 || i == 4){
    obj = obj + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                      plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else if(i == 8 || i == 9){
    obj = obj + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                      plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else if(i == 7){
    obj = obj + theme(plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else{
    obj = obj +  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                       axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                       panel.background = element_rect(fill = "white"),
                       panel.grid.major = element_line(colour = "grey"),
                       panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }
  return(obj)
}

b1 = make_effect_plot(results, 1)
b2 = make_effect_plot(results, 2)
b3 = make_effect_plot(results, 3)
b4 = make_effect_plot(results, 4)
b5 = make_effect_plot(results, 5)
b6 = make_effect_plot(results, 6)
b7 = make_effect_plot(results, 7)
b8 = make_effect_plot(results, 8)
b9 = make_effect_plot(results, 9)

ggarrange(b1, b2, b3, b4, b5, b6, b7, b8, b9, ncol = 3, nrow = 3, common.legend = TRUE, align = "hv")

make_dispersion_plot = function(results, i){
  # input:
  #   - results (list)
  #     list of lists with specific organization
  #   - i (int)
  #     integer informing of what plot that should be plotted and informs of which observation that should be plotted
  # output: 
  #   - obj (ggplot2)
  #     Bar plot of the conditional distribution mass function for a random person of the GLMM
  #
  # Function that creates a bar plot of the conditional distribution mass function for a random person of the GLMM 
  # and plots it with an observation conditional on the effect being zero
  
  
  
  GLM.dipsersion = results$GLM$alpha[, i]
  GLMM.dipsersion = results$GLMM$alpha[, i]
  
  GLM.dipsersion = GLM.dipsersion[!is.na(GLM.dipsersion)]
  GLMM.dipsersion = GLMM.dipsersion[!is.na(GLMM.dipsersion)]
  
  
  #mess with theme to remove whitespace between figures
  obj = ggplot() +
    geom_histogram(aes(x = GLM.dipsersion, y = ..count../sum(..count..), fill = "GLM"), alpha = 0.7, binwidth = 0.01) + 
    geom_histogram(aes(x = GLMM.dipsersion, y = ..count../sum(..count..), fill = "GLMM"), alpha = 0.5, binwidth = 0.01) + 
    geom_vline(aes(xintercept = mean(GLM.dipsersion), col = "GLM")) + geom_vline(aes(xintercept = mean(GLMM.dipsersion), col = "GLMM")) +
    #labs(x = "p values", y = "cumulative proportion", title = "cumulative p values", col = "") +
    labs(x = "", y = "", fill = "Model") + guides(col = "none") +
    annotate("label", x = 0.2, y = 0.15,
             label = paste("cor", round(correlation(mu = rep(exp(beta[1]), 2),
                                                    alpha = alpha, tau = results$tau[i]), 2),
                           "m*", results$group[i]/100)) +
    # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    #       plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
    #       panel.background = element_rect(fill = "white"),
    #       panel.grid.major = element_line(colour = "grey"),
    #       panel.grid.minor = element_line(colour = "grey", size = 0.25)) +
    # theme_bw() +
    scale_x_continuous(limits = c(0.1, 0.6), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) +
    scale_y_continuous(limits = c(0, 0.18), breaks = c(0, 0.07, 0.14))
  # make border plots include x and y notations
  # should use cases instead, this is slower
  if(i == 1 || i == 4){
    obj = obj + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                      plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else if(i == 8 || i == 9){
    obj = obj + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                      plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else if(i == 7){
    obj = obj + theme(plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else{
    obj = obj +  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                       axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                       panel.background = element_rect(fill = "white"),
                       panel.grid.major = element_line(colour = "grey"),
                       panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }
  return(obj)
}

b1 = make_dispersion_plot(results, 1)
b2 = make_dispersion_plot(results, 2)
b3 = make_dispersion_plot(results, 3)
b4 = make_dispersion_plot(results, 4)
b5 = make_dispersion_plot(results, 5)
b6 = make_dispersion_plot(results, 6)
b7 = make_dispersion_plot(results, 7)
b8 = make_dispersion_plot(results, 8)
b9 = make_dispersion_plot(results, 9)

ggarrange(b1, b2, b3, b4, b5, b6, b7, b8, b9, ncol = 3, nrow = 3, common.legend = TRUE, align = "hv")

make_intercept_plot = function(results, i){
  # input:
  #   - results (list)
  #     list of lists with specific organization
  #   - i (int)
  #     integer informing of what plot that should be plotted and informs of which observation that should be plotted
  # output: 
  #   - obj (ggplot2)
  #     Bar plot of the conditional distribution mass function for a random person of the GLMM
  #
  # Function that creates a bar plot of the conditional distribution mass function for a random person of the GLMM 
  # and plots it with an observation conditional on the intercept being zero
  
  
  
  GLM.intercept = results$GLM$intercept[, i]
  GLMM.intercept = results$GLMM$intercept[, i]
  
  GLM.intercept = GLM.intercept[!is.na(GLM.intercept)]
  GLMM.intercept = GLMM.intercept[!is.na(GLMM.intercept)]
  
  
  #mess with theme to remove whitespace between figures
  obj = ggplot() +
    geom_histogram(aes(x = GLM.intercept, y = ..count../sum(..count..), fill = "GLM"), binwidth = 0.01, alpha = 0.8) + 
    geom_histogram(aes(x = GLMM.intercept, y = ..count../sum(..count..), fill = "GLMM"), alpha = 0.5, binwidth = 0.01) + 
    geom_vline(aes(xintercept = mean(GLM.intercept), col = "GLM")) + geom_vline(aes(xintercept = mean(GLMM.intercept), col = "GLMM")) +
    #labs(x = "p values", y = "cumulative proportion", title = "cumulative p values", col = "") +
    labs(x = "", y = "", fill = "Model") + guides(col = "none") +
    annotate("label", x = 5.15, y = 0.075,
             label = paste("cor", round(correlation(mu = rep(exp(beta[1]), 2),
                                                    alpha = alpha, tau = results$tau[i]), 2),
                           "m*", results$group[i]/100)) +
    # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    #       plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
    #       panel.background = element_rect(fill = "white"),
    #       panel.grid.major = element_line(colour = "grey"),
    #       panel.grid.minor = element_line(colour = "grey", size = 0.25)) +
    # theme_bw() +
    scale_x_continuous(limits = c(5.05, 5.6), breaks = c(5.1, 5.2, 5.3, 5.4, 5.5, 5.6)) +
    scale_y_continuous(limits = c(0, 0.09), breaks = c(0, 0.03, 0.06, 0.09))
  # make border plots include x and y notations
  # should use cases instead, this is slower
  if(i == 1 || i == 4){
    obj = obj + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                      plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else if(i == 8 || i == 9){
    obj = obj + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                      plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else if(i == 7){
    obj = obj + theme(plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else{
    obj = obj +  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                       axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                       panel.background = element_rect(fill = "white"),
                       panel.grid.major = element_line(colour = "grey"),
                       panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }
  return(obj)
}

b1 = make_intercept_plot(results, 1)
b2 = make_intercept_plot(results, 2)
b3 = make_intercept_plot(results, 3)
b4 = make_intercept_plot(results, 4)
b5 = make_intercept_plot(results, 5)
b6 = make_intercept_plot(results, 6)
b7 = make_intercept_plot(results, 7)
b8 = make_intercept_plot(results, 8)
b9 = make_intercept_plot(results, 9)

ggarrange(b1, b2, b3, b4, b5, b6, b7, b8, b9, ncol = 3, nrow = 3, common.legend = TRUE, align = "hv")

make_intercept.se_plot = function(results, i){
  # input:
  #   - results (list)
  #     list of lists with specific organization
  #   - i (int)
  #     integer informing of what plot that should be plotted and informs of which observation that should be plotted
  # output: 
  #   - obj (ggplot2)
  #     Bar plot of the conditional distribution mass function for a random person of the GLMM
  #
  # Function that creates a bar plot of the conditional distribution mass function for a random person of the GLMM 
  # and plots it with an observation conditional on the intercept.se being zero
  
  
  
  GLM.intercept.se = results$GLM$intercept.se[, i]
  GLMM.intercept.se = results$GLMM$intercept.se[, i]
  
  GLM.intercept.se = GLM.intercept.se[!is.na(GLM.intercept.se)]
  GLMM.intercept.se = GLMM.intercept.se[!is.na(GLMM.intercept.se)]
  
  
  #mess with theme to remove whitespace between figures
  obj = ggplot() +
    geom_histogram(aes(x = GLM.intercept.se, y = ..count../sum(..count..), fill = "GLM"), binwidth = 0.001, alpha = 0.8) + 
    geom_histogram(aes(x = GLMM.intercept.se, y = ..count../sum(..count..), fill = "GLMM"), alpha = 0.5, binwidth = 0.001) + 
    geom_vline(aes(xintercept = mean(GLM.intercept.se), col = "GLM")) + geom_vline(aes(xintercept = mean(GLMM.intercept.se), col = "GLMM")) +
    #labs(x = "p values", y = "cumulative proportion", title = "cumulative p values", col = "") +
    labs(x = "", y = "", fill = "Model") + guides(col = "none") + 
    annotate("label", x = 0.045, y = 0.15,
             label = paste("cor", round(correlation(mu = rep(exp(beta[1]), 2),
                                                    alpha = alpha, tau = results$tau[i]), 2),
                           "m*", results$group[i]/100)) +
    # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    #       plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
    #       panel.background = element_rect(fill = "white"),
    #       panel.grid.major = element_line(colour = "grey"),
    #       panel.grid.minor = element_line(colour = "grey", size = 0.25)) +
    # theme_bw() +
    scale_x_continuous(limits = c(0.04, 0.08), breaks = c(0.04, 0.05, 0.06, 0.07, 0.08)) +
    scale_y_continuous(limits = c(0, 0.18), breaks = c(0, 0.03, 0.06, 0.09, 0.12, 0.15))
  # make border plots include x and y notations
  # should use cases instead, this is slower
  if(i == 1 || i == 4){
    obj = obj + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                      plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else if(i == 8 || i == 9){
    obj = obj + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                      plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else if(i == 7){
    obj = obj + theme(plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else{
    obj = obj +  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                       axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                       panel.background = element_rect(fill = "white"),
                       panel.grid.major = element_line(colour = "grey"),
                       panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }
  return(obj)
}

b1 = make_intercept.se_plot(results, 1)
b2 = make_intercept.se_plot(results, 2)
b3 = make_intercept.se_plot(results, 3)
b4 = make_intercept.se_plot(results, 4)
b5 = make_intercept.se_plot(results, 5)
b6 = make_intercept.se_plot(results, 6)
b7 = make_intercept.se_plot(results, 7)
b8 = make_intercept.se_plot(results, 8)
b9 = make_intercept.se_plot(results, 9)

ggarrange(b1, b2, b3, b4, b5, b6, b7, b8, b9, ncol = 3, nrow = 3, common.legend = TRUE, align = "hv")

make_effect.se_plot = function(results, i){
  # input:
  #   - results (list)
  #     list of lists with specific organization
  #   - i (int)
  #     integer informing of what plot that should be plotted and informs of which observation that should be plotted
  # output: 
  #   - obj (ggplot2)
  #     Bar plot of the conditional distribution mass function for a random person of the GLMM
  #
  # Function that creates a bar plot of the conditional distribution mass function for a random person of the GLMM 
  # and plots it with an observation conditional on the effect.se being zero
  
  
  
  GLM.effect.se = results$GLM$effect.se[, i]
  GLMM.effect.se = results$GLMM$effect.se[, i]
  
  GLM.effect.se = GLM.effect.se[!is.na(GLM.effect.se)]
  GLMM.effect.se = GLMM.effect.se[!is.na(GLMM.effect.se)]
  
  
  #mess with theme to remove whitespace between figures
  obj = ggplot() +
    geom_histogram(aes(x = GLM.effect.se, y = ..count../sum(..count..), fill = "GLM"), binwidth = 0.001, alpha = 0.8) + 
    geom_histogram(aes(x = GLMM.effect.se, y = ..count../sum(..count..), fill = "GLMM"), alpha = 0.5, binwidth = 0.001) + 
    geom_vline(aes(xintercept = mean(GLM.effect.se), col = "GLM")) + geom_vline(aes(xintercept = mean(GLMM.effect.se), col = "GLMM")) +
    #labs(x = "p values", y = "cumulative proportion", title = "cumulative p values", col = "") +
    labs(x = "", y = "", fill = "Model") + guides(col = "none") +
    annotate("label", x = 0.06, y = 0.1,
             label = paste("cor", round(correlation(mu = rep(exp(beta[1]), 2),
                                                    alpha = alpha, tau = results$tau[i]), 2),
                           "m*", results$group[i]/100)) +
    # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    #       plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
    #       panel.background = element_rect(fill = "white"),
    #       panel.grid.major = element_line(colour = "grey"),
    #       panel.grid.minor = element_line(colour = "grey", size = 0.25)) +
    # theme_bw() +
    scale_x_continuous(limits = c(0.05, 0.11), breaks = c(0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11)) +
    scale_y_continuous(limits = c(0, 0.125), breaks = c(0, 0.03, 0.06, 0.09, 0.12))
  # make border plots include x and y notations
  # should use cases instead, this is slower
  if(i == 1 || i == 4){
    obj = obj + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                      plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else if(i == 8 || i == 9){
    obj = obj + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                      plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else if(i == 7){
    obj = obj + theme(plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else{
    obj = obj +  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                       axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                       panel.background = element_rect(fill = "white"),
                       panel.grid.major = element_line(colour = "grey"),
                       panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }
  return(obj)
}

b1 = make_effect.se_plot(results, 1)
b2 = make_effect.se_plot(results, 2)
b3 = make_effect.se_plot(results, 3)
b4 = make_effect.se_plot(results, 4)
b5 = make_effect.se_plot(results, 5)
b6 = make_effect.se_plot(results, 6)
b7 = make_effect.se_plot(results, 7)
b8 = make_effect.se_plot(results, 8)
b9 = make_effect.se_plot(results, 9)

ggarrange(b1, b2, b3, b4, b5, b6, b7, b8, b9, ncol = 3, nrow = 3, common.legend = TRUE, align = "hv")

make_tau_plot = function(results, i){
  # input:
  #   - results (list)
  #     list of lists with specific organization
  #   - i (int)
  #     integer informing of what plot that should be plotted and informs of which observation that should be plotted
  # output: 
  #   - obj (ggplot2)
  #     Bar plot of the conditional distribution mass function for a random person of the GLMM
  #
  # Function that creates a bar plot of the conditional distribution mass function for a random person of the GLMM 
  # and plots it with an observation conditional on the tau being zero
  
  
  
  GLMM.tau = results$GLMM$tau[, i]
  
  GLMM.tau = GLMM.tau[!is.na(GLMM.tau)]
  
  
  #mess with theme to remove whitespace between figures
  obj = ggplot() +
    geom_histogram(aes(x = GLMM.tau, y = ..count../sum(..count..)), alpha = 0.5, breaks = 0:30/50) + 
    geom_vline(aes(xintercept = mean(GLMM.tau), col = "GLMM")) +
    #labs(x = "p values", y = "cumulative proportion", title = "cumulative p values", col = "") +
    labs(x = "", y = "") + guides(col = "none") +
    annotate("label", x = 0.1, y = 0.125,
             label = paste("cor", round(correlation(mu = rep(exp(beta[1]), 2),
                                                    alpha = alpha, tau = results$tau[i]), 2),
                           "m*", results$group[i]/100)) +
    # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    #       plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
    #       panel.background = element_rect(fill = "white"),
    #       panel.grid.major = element_line(colour = "grey"),
    #       panel.grid.minor = element_line(colour = "grey", size = 0.25)) +
    # theme_bw() +
    scale_x_continuous(limits = c(0, 0.65), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) +
    scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15))
  # make border plots include x and y notations
  # should use cases instead, this is slower
  if(i == 1 || i == 4){
    obj = obj + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                      plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else if(i == 8 || i == 9){
    obj = obj + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                      plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else if(i == 7){
    obj = obj + theme(plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }else{
    obj = obj +  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                       axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       plot.margin = unit(c(-0, -0.5, -0.5, -0.5), "cm"),
                       panel.background = element_rect(fill = "white"),
                       panel.grid.major = element_line(colour = "grey"),
                       panel.grid.minor = element_line(colour = "grey", size = 0.25))
  }
  return(obj)
}

b1 = make_tau_plot(results, 1)
b2 = make_tau_plot(results, 2)
b3 = make_tau_plot(results, 3)
b4 = make_tau_plot(results, 4)
b5 = make_tau_plot(results, 5)
b6 = make_tau_plot(results, 6)
b7 = make_tau_plot(results, 7)
b8 = make_tau_plot(results, 8)
b9 = make_tau_plot(results, 9)

ggarrange(b1, b2, b3, b4, b5, b6, b7, b8, b9, ncol = 3, nrow = 3, common.legend = TRUE, align = "hv")
