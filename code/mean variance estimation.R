library(lme4)
library(insight)
library(performance)
library(ggplot2)
library(MASS)
library(ggpubr)
library(DESeq2)

set.seed(1)

G = 1000
# mean(log(rowMeans(count)[rowMeans(count) != 0]))
# 1.959193
# sd(log(rowMeans(count)[rowMeans(count) != 0]))
# 3.420544
# ggplot() + geom_histogram(aes(x = log(rowMeans(count)), y = ..count../sum(..count..)))
#
# with only keep
# mean(log(rowMeans(count[keep, ])))
# 4.491542
# sd(log(rowMeans(count)[keep]))
# 1.912126
beta = cbind(rnorm(G, mean = 4.5, sd = 2), rep(0, G))

# alpha = alpha_0/mu + alpha_1
mean_var = c(10, 0.15)

alpha = rnorm(n = G, mean = mean_var[1]/exp(beta[, 1] + beta[, 2]/2) + mean_var[2], sd = 0.1)
# check of plot, compare to plotDispEsts(dds) # full for completeness
ggplot() + geom_line(aes(y = mean_var[1]/(exp(1:8)) + mean_var[2], x = exp(1:8))) + scale_x_log10() + scale_y_log10()
tau = 0.4
group = 100

# define functions for later use
create_dataset = function(group = 100, n = 200, beta = c(1, 0), alpha = 0.25){
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

alpha_est = list("GLM" = rep(NA, G), "GLMM" = rep(NA, G))

m = 200
x = cbind(rep(1, m), rep(c(0, 1), m/2))

counts = matrix(0, nrow = G, ncol = m)
iter = rep(0, G)
P = list("GLM" = rep(NA, G), "GLMM" = rep(NA, G), "DESeq" = rep(NA, G))

for(g in 1:G){
  tryCatch(expr = {
  temp = create_dataset(group = group, n = m, beta = beta[g, ], alpha = alpha[g])
  counts[g, ] = temp$y
  temp_fit = glm.nb(y ~ x-1, link = "log", data = temp)
  temp_GLMM = glmer.nb(formula = y ~ x-1 + (1|pers), data = temp)
  alpha_est$GLM[g] = 1/temp_fit$theta
  alpha_est$GLMM[g] = 1/getME(temp_GLMM, "glmer.nb.theta")
  
  iter[g] = temp_GLMM@nevals
  P$GLM[g] = coef(summary(temp_fit))[2, 4]
  P$GLMM[g] = coef(summary(temp_GLMM))[2, 4]
  }, 
  warning = function(cond){
  }, 
  error = function(cond){
  }
  )
}

counts = counts[!is.na(alpha_est$GLM), ]
iter = iter[!is.na(alpha_est$GLM)]
alpha_est$GLM = alpha_est$GLM[!is.na(alpha_est$GLM)]
alpha_est$GLMM = alpha_est$GLMM[!is.na(alpha_est$GLMM)]

info = data.frame("id" = temp$pers, "effect" = factor(x[, 2]))

dds = DESeqDataSetFromMatrix(countData = counts, colData = info, 
                             design = ~ effect)
dds = estimateSizeFactors(dds)

# three level estimate dispersion

dds = estimateDispersionsGeneEst(dds)
# just alter the relevant estimate?
# dds@rowRanges@elementMetadata@listData$dispGeneEst = alpha_est$GLMM
# dds@rowRanges@elementMetadata@listData$dispGeneIter = iter

dds = estimateDispersionsFit(dds)
dds = estimateDispersionsMAP(dds)
plotDispEsts(dds[1:10, ])

dds = nbinomWaldTest(dds)

# plot under H_0 of p values with mean variance relationship implemented from DESeq2, data originating from GLMM data
ggplot() + geom_histogram(aes(x =results(dds)[, 5], fill = "DESeq2"), breaks = 0:20/20, alpha = 0.5) + 
  geom_histogram(aes(x = P$GLM, fill = "GLM"), breaks = 0:20/20, alpha = 0.5) + 
  geom_histogram(aes(x = P$GLMM, fill = "GLMM"), breaks = 0:20/20, alpha = 0.5) + 
  theme_bw()

#' relevant changes:
#'   
#'   dds@rowRanges@elementMetadata
#' @listData 
#' - base Mean, keep, mean(counts[i, ] / dds@colData@listData$sizeFactor)
#' - base Var, keep, var(counts[i, ] / dds@colData@listData$sizeFactor)
#' - allZero, keep
#' - dispGeneEst, easy change
#' - dispGeneIter, easy change
#' @elementMetadata - descriptions, keep
#' 
#' dds@assays@data@listData$mu
#' ?, note that this is not mean(counts[i, ] / dds@colData@listData$sizeFactor) * dds@colData@listData$sizeFactor[j]
#' No idea, but is observation and gene specific
#' 
#' Does DESeq2 utilitze base var for estimates of dispersion?
#'   As it assumes GLM this might not translate to GLMM. 
