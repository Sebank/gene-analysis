library(lme4)
library(insight)
library(performance)

set.seed(123)

beta=matrix(c(10,0),ncol=1)

tau= 1

alpha=1.1

n=200

x=cbind(rep(1,n),c(rep(0,n/2),rep(1,n/2)))

pers=rep(1:100,each=2)

pers

gamma=rep(rnorm(n/2,0,tau),each=2)

gamma

eta=x%*%beta+gamma

y=rnbinom(n=n,size=1/alpha,mu=exp(eta))

y

fit=glm(y~x-1,family=MASS::negative.binomial(link = "log", theta = 1/alpha))

summary(fit,dispersion=1)

vcov(fit, dispersion = 1)

sqrt(diag(vcov(summary(fit, dispersion = 1))))

fit2 = glmer(formula = y~x-1 + (1|pers), family = MASS::negative.binomial(link = "log", theta = 1/alpha))

summary(fit2)

vcov(fit2)

sqrt(diag(vcov(summary(fit2))))

icc(fit2)

# potential fail check
summary(fit2)$optinfo$conv$lme4$code

m = 10000
results = data.frame("intercept" = rep(NA, m), "effect" = rep(NA, m), "intercept.se" = rep(NA, m), "effect.se" = rep(NA, m))
results_GLMM = data.frame("intercept" = rep(NA, m), "effect" = rep(NA, m), "intercept.se" = rep(NA, m), "effect.se" = rep(NA, m), "ICC" = rep(NA, m), "ICC_conditional" = rep(NA, m))

for(i in 1:m){
  tryCatch(
    expr = {
      temp = rnbinom(n=n,size=1/alpha,mu=exp(eta))
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
      print(cond)
    }
  )
}

ggplot() + 
  geom_histogram(aes(x = 2 * pnorm(abs(results$effect/results$effect.se), lower.tail = FALSE), y = ..count../sum(..count..), fill = "GLM"), alpha = 1) + 
  geom_histogram(aes(x = 2 * pnorm(abs(results_GLMM$effect/results_GLMM$effect.se), lower.tail = FALSE), y = ..count../sum(..count..), fill = "GLMM"), alpha = 0.5) +
  labs(x = "P value for ICC", title = "counts of P values for GLMM and GLM simulations")
ggplot() + 
  geom_histogram(aes(x = results_GLMM$ICC, y = ..count../sum(..count..), fill = "adjusted ICC"), binwidth = 0.01) + 
  geom_histogram(aes(x = results_GLMM$ICC_conditional, y = ..count../sum(..count..), fill = "conditional ICC"), alpha = 0.5, binwidth = 0.01) + 
  labs(x = "estimated ICC", title = "distribution of ICC when simulating a correlated GLMM")

ggplot() + 
  geom_point(aes(x = 2 * pnorm(abs(results$effect/results$effect.se), lower.tail = FALSE), 
                 y = 2 * pnorm(abs(results_GLMM$effect/results_GLMM$effect.se), lower.tail = FALSE)), alpha = 0.05) +
  labs(x = "P value GLM", y = "P value GLMM", title = "comparison of p values for ICC")

ggplot() +
  geom_point(aes(x = results_GLMM$ICC,
                 y = 2 * pnorm(abs(results_GLMM$effect/results_GLMM$effect.se), lower.tail = FALSE) -
                 2 * pnorm(abs(results$effect/results$effect.se), lower.tail = FALSE) ), alpha = 0.05) +
  labs(x = "ICC adjusted", y = "difference in P value between GLM against GLMM as reference")
