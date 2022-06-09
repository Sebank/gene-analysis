library(stringr)
library(ggvenn)
library(edgeR)
library(limma)
library(DESeq2)
library(lme4)
library(hexbin) 
library(apeglm)
library(readr)
library(ggplot2)

#wald.test
library(aod)

#file location
fileloc = "Z:\\Atle van Beelen Granlund/"
git = "D:\\Essential folders\\Documents\\RStudio\\master\\gene-analysis\\code"


###########################################################################
#pre-processing
###########################################################################

# info about the samples to connect to disease, tissue, inflammation, pasid

pas_info_tot = read.delim(paste(fileloc,'sample_info.txt',sep=""), 
                          sep="\t", header = TRUE, row.names = 1)
count_tot = read.delim(paste(fileloc, "quasi.gene.quant", sep = ""), 
                       sep = "\t", check.names = FALSE)

#remove data that is considered bad
removal = c(which(pas_info_tot$Sample_ID == "3115i"), 
            which(pas_info_tot$Sample_ID == "3073i_new"))

#remove outliers
removal = c(removal, 
            which(pas_info_tot$Sample_ID == "3009F"), #the one with 3 observations
            which(pas_info_tot$Sample_ID == "TT27F_t")) #only singular observation

#this has been run twice for some reason
removal = c(removal, which(pas_info_tot$Sample_ID == "3002i"))

# only sample with less than a million counts (306740)
removal = c(removal, which(pas_info_tot$Sample_ID == "TT38_t"))

#Corrections from Atle
pas_info_tot$Sample_Biosource[pas_info_tot$Sample_ID == "3114i_t"] = "Colon tissue"
pas_info_tot$Sample_Biosource[pas_info_tot$Sample_ID == "3151i"] = "Ileum tissue"

pas_info_tot = pas_info_tot[-removal, ]
count_tot = count_tot[, pas_info_tot$Sample_ID]

#extract relevant columns
pas_info <- data.frame("sampleid" = pas_info_tot$Sample_ID, 
                       "diseaseinflammation" = pas_info_tot$Sample_Group, 
                       "tissue" = pas_info_tot$Sample_Biosource)
#remove tissue from tissue
pas_info$tissue = substr(pas_info$tissue, 1, 5)

#using grep to remove final segments of string ending with irrelevant information
pas_info$sampleid = sub("[iFS]", "", pas_info$sampleid)
pas_info$sampleid = sub("_t", "", pas_info$sampleid)

#Assuming we follow tissue instead of inflammation
#remove i from disease information, as it should be part of tissue as well
pas_info$diseaseinflammation = sub("i", "", pas_info$diseaseinflammation)
pas_info$diseaseinflammation = sub("F", "H", pas_info$diseaseinflammation)

#more a double check than anything else
#unique(pas_info$diseaseinflammation)

#make the variable disease inflammation two variables
#replaced F with H for "healthy"
AU = grep("..[AU]", pas_info$diseaseinflammation)
pas_info$disease = rep("H", length(pas_info$diseaseinflammation))
pas_info$disease[AU] = substr(pas_info$diseaseinflammation[AU], 1, 
                              nchar(pas_info$diseaseinflammation[AU]) - 1)


A = grep("A", pas_info$diseaseinflammation)
H = grep("H", pas_info$diseaseinflammation)
pas_info$inflammation = rep("U", length(pas_info$diseaseinflammation))
pas_info$inflammation[A] = "A"
pas_info$inflammation[H] = "H"

# finds tables for comparing how much data we have for each group
table(pas_info$disease)
table(pas_info$inflammation)
table(pas_info$tissue,pas_info$inflammation,pas_info$disease)


#Only examine CD and H, not UC
count = count_tot[, !pas_info$disease == "UC"]
pas_info = pas_info[!pas_info$disease == "UC", ]

# Check distribution of tissues and inflammation, and correlated data
table(pas_info$tissue,pas_info$inflammation)
table(table(pas_info$sampleid)) 

# change the names of count to the modified names from pas_info
# observe that the count was sorted based on the names of pas_info before manipulation
colnames(count) = pas_info$sampleid

#change F to H in diseaseinflammation and make variables factors
pas_info$diseaseinflammation[pas_info$diseaseinflammation == "F"] = "H"
pas_info$diseaseinflammation = factor(pas_info$diseaseinflammation)
pas_info$inflammation = factor(pas_info$inflammation, c("H", "U", "A"))
pas_info$tissue = factor(pas_info$tissue)

###########################################################################
#limma-voom (model)
###########################################################################

#convert to relevant setup for limma-voom

dge = DGEList(count)

#make as factors and order factors
pas_info$inflammation = pas_info$inflammation
pas_info$tissue = pas_info$tissue

#"remove counts that have zero or very low counts"
design = model.matrix(~ inflammation*tissue, data = pas_info)
keep = filterByExpr(round(count), design)
dge = dge[keep, , keep.lib.sizes=FALSE]

dge = calcNormFactors(dge)

#voom
v = voom(dge, design = design, plot=TRUE)

# plot of mean trend in s^{1/2}
plotMDS(v, labels = substr(pas_info$tissue, 1, 1), 
        col = as.numeric(pas_info$diseaseinflammation))
legend("top", legend = levels(pas_info$diseaseinflammation), col = 1:3, pch = 15)

# compute consensus correlation for pipeline that assumes correlated data
voom_cor = duplicateCorrelation(v,design,block=pas_info$sampleid)

# make model with and without considering correlation
vfit = lmFit(v, design, weights = v$weights)
vfit_cor = lmFit(v, design, block = pas_info$sampleid, correlation = voom_cor$consensus, weights = v$weights)

# create contrast matrix for consideration, mentioned in master
C = rbind("IU"=c(0,0,0,0,1,0), #(CU-CH)-(IU-IH)
          "IA"=c(0,0,0,0,0,1), #(CA-CH)-(IA-IH)
          "IA minus IU"=c(0,0,0,0,-1,1)) #(CA-CU)-(IA-IU)

# fit specific contrasts and consider the empirical Bayes model for limma-voom
vfit2 = contrasts.fit(vfit, t(C))
vfit2 = eBayes(vfit2, trend = FALSE)

vfit2_cor = contrasts.fit(vfit_cor, t(C))
vfit2_cor = eBayes(vfit2_cor, trend = FALSE)

#for specific coef use coef = 1, 2 or 3, as contrast has three coefs
table_vfit2 = topTable(vfit2, number = 10, coef = 3)
table_vfit2_cor = topTable(vfit2_cor, number = 10, coef = 3)

# plot residual variance after transformation with respect to weights
plotSA(vfit2)
plotSA(vfit2_cor)

nr = 3
# create indicator variable to identify differentially expressed genes for the different contrasts
indicator = rep(0, length(rownames(dge$count)))
indicator = indicator + (rownames(dge$count) %in% rownames(topTable(vfit2, number = 10000, coef = nr, 
                                                                    p.value = 0.05))[topTable(vfit2, coef = nr, number = 1000,
                                                                                              p.value = 0.05)[, 5] < 0.05])
# plot MA plot with highlighted significant genes
limma::plotMA(vfit2, coef = nr, status = indicator, hl.cex = c(0.6, 0.1))

# Similarly for the correlated data 
indicator_cor = rep(0, length(rownames(dge$count)))
indicator_cor = indicator_cor + (rownames(dge$count) %in% rownames(topTable(vfit2_cor, number = 10000, coef = nr, 
                                                                        p.value = 0.05))[topTable(vfit2_cor, coef = nr, number = 1000,
                                                                                                  p.value = 0.05)[, 5] < 0.05])
limma::plotMA(vfit2_cor, coef = nr, status = indicator_cor, hl.cex = c(0.6, 0.1))

# combine histogram of p values into one plot for reduction of total number of plots
ggplot() + geom_histogram(aes(vfit2$p.value[, nr], y = ..count../sum(..count..), fill = "no"), breaks = 0:20/20, alpha = 0.6, col = "white") +
  geom_histogram(aes(vfit2_cor$p.value[, nr], y = ..count../sum(..count..), fill = "yes"), breaks = 0:20/20, alpha = 0.3, col = "white") +
  labs(x = "p-values", y = "proportion", fill = "correlation") + lims(y = c(0, 0.195)) + 
  theme_bw()

# plot distribution of atanh transform for correlation
G = length(voom_cor$atanh.correlations)
voom_cor$atanh.correlations = sort(voom_cor$atanh.correlations)
ggplot() + 
  geom_histogram(aes(x = voom_cor$atanh.correlations, y = ..count../sum(..count..)), binwidth = 0.02) + 
  geom_vline(aes(xintercept = c(voom_cor$atanh.correlations[round(0.15*G)], 
                                voom_cor$atanh.correlations[round(0.85*G)]))) +
  labs(x = "atanh correlation", y = "proportion") + theme_bw()

###########################################################################
#deseq (model)
###########################################################################

# create deseq2 compatible dds
dds = DESeqDataSetFromMatrix(countData = round(count), colData = pas_info, 
                             design = ~ inflammation*tissue)

#pre processing, this is that same as is done for limma-voom
keep_DESeq = filterByExpr(counts(dds), design)
dds = dds[keep_DESeq, ]


# equal to three steps beneath
# dds <- DESeq(dds)

# normalization step
dds = estimateSizeFactors(dds)
# mean variance estimation step
dds = estimateDispersions(dds)
# effect estimation and hypotesis testing step
dds = nbinomWaldTest(dds) # , betaPrior = TRUE


contrast.UI <- results(dds, contrast = C[1, ], alpha = 0.05)
# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="Intercept", type="apeglm")

#for consideration
removal.UI = contrast.UI[is.na(contrast.UI$padj), ]
contrast.UI = contrast.UI[!is.na(contrast.UI$padj), ]
contrast.UI = contrast.UI[order(contrast.UI$padj), ]

contrast.AI <- results(dds, contrast = C[2, ], alpha = 0.05)

#for consideration
removal.AI = contrast.AI[is.na(contrast.AI$padj), ]
contrast.AI = contrast.AI[!is.na(contrast.AI$padj), ]
contrast.AI = contrast.AI[order(contrast.AI$padj), ]


contrast.complex <- results(dds, contrast = C[3, ], alpha = 0.05)

#for consideration
removal.complex = contrast.complex[is.na(contrast.complex$padj), ]
contrast.complex = contrast.complex[!is.na(contrast.complex$padj), ]
contrast.complex = contrast.complex[order(contrast.complex$padj), ]

# plot mean variance relationship
plotDispEsts(dds[sample(1:dim(dds)[1], 1000), ]) # smaller for visibility
plotDispEsts(dds) # full for completeness

# MA plot
DESeq2::plotMA(contrast.UI, ylim = c(-4.5, 4.5))
#DESeq2::plotPCA(rlogTransformation(dds), intgroup = c("tissue", "diseaseinflammation"))

#plots of potential interest
ggplot() + geom_point(aes(estimateSizeFactorsForMatrix(count), colSums(count)))
#observation 88 is the smallest

###########################################################################
# GLMM
###########################################################################

# computes the correlation between two observations based on offset + estimated coefficients (mu),
# standard deviation of mixed effect (tau) and dispersion estimate (alpha) for negative binomial distribution
correlation = function(mu, alpha, tau){
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

# find correlated observations and group them together
library(MASS)
name = names(table(pas_info$sampleid))[ table(pas_info$sampleid) > 1]
group = matrix(0, length(pas_info$sampleid), length(name))
colnames(group) = name

for(i in 1:length(name)){
  group[which(pas_info$sampleid == name[i]), i] = 1
}

# define filtered counts and extract GLMM dispersion
# dispersion = 1/theta 
genes = dds@assays@data@listData$counts

# create checkContrast (long so import instead)
n = length(rownames(genes))

checkContrast = data.frame("UI.glmm" = rep(NA, n), "UI.se" = rep(NA, n), "UI.P" = rep(NA, n),
                           "AI.glmm" = rep(NA, n), "AI.se" = rep(NA, n), "AI.P" = rep(NA, n),
                           "complex.glmm" = rep(NA, n), "complex.se" = rep(NA, n), "complex.P" = rep(NA, n),
                           "mixed.sd" = rep(NA, n), "alpha" = rep(NA, n))
rownames(checkContrast) = rownames(genes)
cova = vector(mode = "double", length = dim(group)[2])
corre = matrix(nrow = n, ncol = dim(group)[2])

# If continue old run
# time = Sys.time()
# for(i in 1:dim(genes)[1]){
#   tryCatch(expr = {
#     if(i %% 275 == 1){
#       print(paste("Progression:", 100*i/27500, "%."), sep = "")
#       print(Sys.time() - time)
#     }
#     temp = glmer.nb(
#       formula = genes[i, ] ~ tissue * inflammation + (1|sampleid),
#       family = MASS::negative.binomial(link = "log"),
#       offset = log(dds@colData@listData$sizeFactor), data = pas_info
#     )
#     checkContrast$UI.glmm[i] = fixef(temp)[5]*log2(exp(1))
#     checkContrast$AI.glmm[i] = fixef(temp)[6]*log2(exp(1))
#     checkContrast$complex.glmm[i] = ((C[3, ] %*% fixef(temp)) * log2(exp(1)))[1, 1]
#     
#     checkContrast$UI.se[i] = sqrt(t(C[1, ]) %*% vcov(temp) %*% C[1, ])[1, 1] * log2(exp(1))
#     checkContrast$AI.se[i] = sqrt(t(C[2, ]) %*% vcov(temp) %*% C[2, ])[1, 1] * log2(exp(1))
#     checkContrast$complex.se[i] = sqrt(t(C[3, ]) %*% vcov(temp) %*% C[3, ])[1, 1] * log2(exp(1))
#     
#     
#     checkContrast$mixed.sd[i] = temp@pp$theta
#     checkContrast$alpha[i] = 1/getME(temp, "glmer.nb.theta")
#     
#     P = coef(summary(temp))[, 4]
#     
#     checkContrast$UI.P[i] = P[5]
#     checkContrast$AI.P[i] = P[6]
#     checkContrast$complex.P[i] = pchisq( (checkContrast$complex.glmm[i]/checkContrast$complex.se[i])^2 , df = 1, lower.tail = FALSE)
#     
#     for(l in 1:dim(group)[2]){
#       j = which(group[, l] == 1)[1]
#       k = which(group[, l] == 1)[2]
#       
#       tau = temp@pp$theta
#       mu_cov = fitted(temp)[c(j, k)]
#       corre[i, l] = correlation(mu = mu_cov, alpha = dispersion[i], tau = tau)
#     }
#   }, warning = function(cond){
#     print(cond)
#   }, error = function(cond){
#     print(cond)
#   },
#   silent = TRUE
#   )
# }

# write.csv(checkContrast, "D:\\Essential folders\\Documents\\RStudio\\master\\gene-analysis\\code\\checkContrast.csv")
# write.csv(corre, "D:\\Essential folders\\Documents\\RStudio\\master\\gene-analysis\\code\\corre.csv")

# These require the location of the git directory
corre = read.csv(paste(git, "\\corre.csv", sep = ""), header = TRUE, row.names = 1)
checkContrast = read.csv(paste(git, "\\checkContrast.csv", sep = ""), header = TRUE, row.names = 1)

# adjusted p-values
checkContrast = checkContrast[!is.na(checkContrast$UI.P), ]
complex.P.adj = order(checkContrast$complex.P)
UI.P.adj = order(checkContrast$UI.P[])
AI.P.adj = order(checkContrast$AI.P)
checkContrast$complex.P.adj = checkContrast$complex.P
checkContrast$UI.P.adj = checkContrast$UI.P
checkContrast$AI.P.adj = checkContrast$AI.P

relevant_genes = dim(counts(dds))[1]

i = 10944

for(i in (dim(checkContrast)[1] - 1):1){
  checkContrast$complex.P.adj[complex.P.adj[i]] = min(checkContrast$complex.P.adj[complex.P.adj[i + 1]], 
                                                      checkContrast$complex.P.adj[complex.P.adj[i]] * relevant_genes /i)
  checkContrast$UI.P.adj[UI.P.adj[i]] =  min(checkContrast$UI.P.adj[UI.P.adj[i + 1]], 
                                             checkContrast$UI.P.adj[UI.P.adj[i]] * relevant_genes /i)
  checkContrast$AI.P.adj[AI.P.adj[i]] =  min(checkContrast$AI.P.adj[AI.P.adj[i + 1]], 
                                             checkContrast$AI.P.adj[AI.P.adj[i]] * relevant_genes /i)
}

# contrast names are UI, AI and complex (which has the alternative interpretation of the comparison of the comparison)
# custom MA plot
found_significant = checkContrast$AI.P.adj < 0.05
ggplot() + 
  geom_point(aes(x = rowMeans(t(t(count[rownames(checkContrast[!found_significant, ]), ])/dds@colData@listData$sizeFactor)), 
                 y = checkContrast$AI.glmm[!found_significant], col = "no"), size = 0.4) +
  geom_point(aes(x = rowMeans(t(t(count[rownames(checkContrast[found_significant, ]), ])/dds@colData@listData$sizeFactor)), 
                 y = checkContrast$AI.glmm[found_significant], col = "yes"), size = 0.6) + 
  scale_x_log10() +
  labs(x = "mean of normalized counts", y = "log fold change", col = "Found significant") +
  theme_bw() + scale_y_continuous(limits = c(-4.5, 4.5), breaks = c(-4, -2, 0, 2, 4))

# mean variance
consider = which(!checkContrast$alpha > 10)
ggplot() + 
  geom_point(aes(x = rowMeans(t(t(count[rownames(checkContrast[consider, ]), ])/dds@colData@listData$sizeFactor)),
                 y = checkContrast$alpha[consider]), size = 0.5, alpha = 0.3) + 
  scale_x_log10() + scale_y_log10() + 
  labs(x = "mean of normalized counts", y = "dispersion parameter estimate") + 
  theme_bw()

tanh(mean(atanh(corre[, 1]), na.rm = TRUE, trim = 0.15))
# correlation plot for GLMM
ggplot() + geom_histogram(aes(x = atanh(corre[, 1]), y = ..count../sum(..count..)), breaks = 0:35/20) +
  geom_vline(aes(xintercept = c(quantile(corre[, 1], 0.15, na.rm = TRUE), quantile(corre[, 1], 0.85, na.rm = TRUE)))) + 
  labs(x = "arctanh of correlation", y = "proportion") + theme_bw()

# Plot of distribution of p values
ggplot() + 
  geom_histogram(aes(x = contrast.AI$pvalue, y = ..count../sum(..count..), fill = "no"), breaks = 0:20/20, alpha = 0.6, col = "white") + 
  geom_histogram(aes(x = checkContrast$AI.P, y = ..count../sum(..count..), fill = "yes"), breaks = 0:20/20, alpha = 0.3, col = "white") + 
  theme_bw() + 
  labs(y = "proportion", x = "p-values", fill = "correlation")

###########################################################################
#analysis of results
###########################################################################

#general venn diagram

# compare relevant genes for the different models
number = 10000
venn.1 = list(
  "voom_cor" = rownames(topTable(vfit2_cor, number = number, coef = 1, p.value = 0.05)),
  "voom" = rownames(topTable(vfit2, number = number, coef = 1, p.value = 0.05)),
  # "DESeq2" = rownames(contrast.UI[1:number, ]))
  "DESeq2" = rownames(contrast.UI[contrast.UI$padj < 0.05, ]), 
  "GLMM" = rownames(checkContrast)[checkContrast$UI.P.adj < 0.05])
ggvenn(venn.1, show_percentage = FALSE)

venn.2 = list(
  "voom_cor" = rownames(topTable(vfit2_cor, number = number, coef = 2, p.value = 0.05)),  
  "voom" = rownames(topTable(vfit2, number = number, coef = 2, p.value = 0.05)),
  # "DESeq2" = rownames(contrast.AI[1:number, ]))
  "DESeq2" = rownames(contrast.AI[contrast.AI$padj < 0.05, ]), 
  "GLMM" = rownames(checkContrast)[checkContrast$AI.P.adj < 0.05])
ggvenn(venn.2, show_percentage = FALSE)

venn.3 = list(
  "voom_cor" = rownames(topTable(vfit2_cor, number = number, coef = 3, p.value = 0.05)),  
  "voom" = rownames(topTable(vfit2, number = number, coef = 3, p.value = 0.05)),
  # "DESeq2" = rownames(contrast.complex[1:number, ]))
  "DESeq2" = rownames(contrast.complex[contrast.complex$padj < 0.05, ]), 
  "GLMM" = rownames(checkContrast)[checkContrast$complex.P.adj < 0.05])
ggvenn(venn.3, show_percentage = FALSE)

###########################################################################
# Residual plots
###########################################################################
genes = dds@assays@data@listData$counts
dispersion = dds@rowRanges@elementMetadata@listData$dispersion

i = which(rownames(genes) == "ENSG00000155380")
GLMM = glmer.nb(
  formula = genes[i, ] ~ -1 + design + (1|sampleid),
  family = MASS::negative.binomial(link = "log"),
  offset = log(dds@colData@listData$sizeFactor), data = pas_info
)
sf = estimateSizeFactorsForMatrix(dds@assays@data$counts)
ncounts = t(t(dds@assays@data$counts[i, ])/sf)

fitted.common.scale = t(t(assays(dds)[["mu"]])/sizeFactors(dds))

# implemented quantile residual, observe size = 1/alpha (alpha is our dispersion)
residual = list("voom" = residuals(vfit[i, ], v[i, ])*sqrt(vfit2_cor$s2.post[i]), 
                      "voom.corr" = residuals(vfit_cor[i, ], v[i, ])*sqrt(vfit2_cor$s2.post[i]), 
                      "DESeq" = qnorm(
                          runif(dim(count)[2], 
                                min = pnbinom(round(counts(dds[i ], normalized=TRUE), 0), 
                                          size = 1/dispersion[i], mu = fitted.common.scale[i, ]), 
                                max = pnbinom(round(counts(dds[i ], normalized=TRUE) + 1, 0), 
                                              size = 1/dispersion[i], mu = fitted.common.scale[i, ])
                          )
                        ),
                      "GLMM" = qnorm(
                          runif(dim(count)[2], 
                                min = pnbinom(round(counts(dds[i ], normalized=TRUE), 0), 
                                         size = getME(GLMM, "glmer.nb.theta"), 
                                         mu = exp(design %*% fixef(GLMM))), 
                                max = pnbinom(round(counts(dds[i ], normalized=TRUE) + 1, 0), 
                                              size = getME(GLMM, "glmer.nb.theta"), 
                                              mu = exp(design %*% fixef(GLMM)))
                          )
                        )
)

g1 = ggplot() + 
  stat_qq(aes(sample = sort(residual$voom))) + stat_qq_line(aes(sample = sort(residual$voom))) + 
  labs(x = "theoretical", y = "observed") + 
  theme_bw()

g2 = ggplot() + 
  stat_qq(aes(sample = sort(residual$voom.cor))) + stat_qq_line(aes(sample = sort(residual$voom.cor))) + 
  labs(x = "theoretical", y = "observed") + 
  theme_bw()

g3 = ggplot() + 
  stat_qq(aes(sample = sort(residual$DESeq))) + stat_qq_line(aes(sample = sort(residual$DESeq))) + 
  labs(x = "theoretical", y = "observed") + 
  theme_bw()

g4 = ggplot() + 
  stat_qq(aes(sample = sort(residual$GLMM))) + stat_qq_line(aes(sample = sort(residual$GLMM))) + 
  labs(x = "theoretical", y = "observed") + 
  theme_bw()

ggarrange(g1, g2, g3, g4, nrow = 2, ncol = 2, labels = c("voom", "voom_cor", "DESeq2", "GLMM"), label.x = 0.1)

###########################################################################
# More DESeq2 plots
###########################################################################

# Variance stabilizing transform
vsd = vst(dds)
DESeq2::plotPCA(vsd, intgroup = c("tissue", "inflammation"))

res = results(dds, contrast = C[3, ], alpha = 0.05)#LFCshrink, type= "apeglm")
summary(res)
res[order(res$padj), ] %>% head()
# results are considered by four plots, histogram of p values, MA plot, ordination plot and a heatmap

ggplot(as(res, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)

DESeq2::plotMA(res, ylim = c(-5, 5))

# plot rlog transform as PCA to get an idea of variance
# Slow
# pas_rlog = rlogTransformation(dds)
# DESeq2::plotPCA(pas_rlog, intgroup = c("tissue", "inflammation"))

library("pheatmap")
select = order(rowMeans(assay(vsd)), decreasing = TRUE)[1:200]
pheatmap( assay(vsd)[select, ],
          scale = "row",
          annotation_col = as.data.frame(
            colData(vsd)[, c("tissue", "inflammation")]
          ), show_rownames =  FALSE, show_colnames = FALSE, 
          clustering_distance_cols = "correlation", 
          clustering_distance_rows = "correlation"
)
select = sample(1:dim(count)[1], 10)
pheatmap( assay(vsd)[select, ],
          scale = "row",
          annotation_col = as.data.frame(
            colData(vsd)[, c("tissue", "inflammation")]
          ), show_rownames =  FALSE, show_colnames = FALSE
)

# sanity checks
# plot of library sizes for our observations
cs <- colSums(assay(dds, "counts"))
hist(cs/1e6, col="grey", border="white",
     main="", xlab="column sums (per million)")


###########################################################################
# Results tables and plots
###########################################################################

#analyze data
ggplot(as(contrast.complex, "data.frame"), aes(x = pvalue)) + 
  geom_histogram(binwidth = 0.01, fill = "Blue", boundary = 0)

DESeq2::plotMA(contrast.complex, ylim = c(-5, 5))
# DESeq2::plotPCA(pas_rlog, intgroup = c("tissue", "diseaseinflammation"))
#plots of potential interest
ggplot() + geom_point(aes(estimateSizeFactorsForMatrix(count), colSums(count)))
#observation 88 is the smallest
sf = estimateSizeFactorsForMatrix(dds@assays@data$counts)
ncounts = t(t(dds@assays@data$counts)/sf)
ggplot() + geom_hex(aes(x = log(rowMeans(ncounts)), y = log(rowVars(ncounts)))) + 
  coord_fixed() + theme(legend.position = "none") + 
  geom_abline(slope = 1:2, color = c("forestgreen", "red"))
#green is variance for Poisson, red is for gamma


C.coef.DESeq2 = rbind("CH"=c(1,0,0,0,0,0), 
                      "CU"=c(1,0,1,0,0,0),
                      "CA"=c(1,0,0,1,0,0),
                      "IH"=c(1,1,0,0,0,0), 
                      "IU"=c(1,1,1,0,1,0),
                      "IA"=c(1,1,0,1,0,1))
C.coef.voom = rbind("CH"=c(1,0,0,0,0,0), 
                    "CU"=c(1,1,0,0,0,0),
                    "CA"=c(1,0,1,0,0,0),
                    "IH"=c(1,0,0,1,0,0), 
                    "IU"=c(1,1,0,1,1,0),
                    "IA"=c(1,0,1,1,0,1))
C.coef.names = factor(rownames(C.coef.DESeq2), levels = unique(rownames(C.coef.DESeq2)))

#plot of all observations
pas_info$comment = paste(substr(pas_info$tissue, 1, 1), pas_info$inflammation, sep = "")

# "ENSG00000155380"
name = "ENSG00000280390" #rownames(contrast.complex[1, ])
y_DESeq2 = coef(dds[name])
ggplot()  + 
  geom_point(aes(x = factor(pas_info$comment, levels = c("CH", "CU", "CA", "IH", "IU", "IA")), y = unlist(count[name, ])), alpha = 0.3) + 
  geom_point(aes(x = C.coef.names, y = 2^(C.coef.DESeq2 %*% t(y_DESeq2))), col = "red")
#plot of all observations


limma_baseline = eBayes(vfit, trend = FALSE)
y_limma = coef(limma_baseline[name, ])
ggplot() + geom_point(aes(x = C.coef.names, y = C.coef.voom %*% t(y_limma)))

# "ENSG00000109099"
name = rownames(topTable(vfit2, number = 1, coef = 3))
y_DESeq2 = coef(dds[name])
ggplot() + geom_point(aes(x = C.coef.names, y = C.coef.DESeq2 %*% t(y_DESeq2)))

limma_baseline = eBayes(vfit, trend = FALSE)
y_limma = coef(limma_baseline[name, ])
ggplot() + geom_point(aes(x = C.coef.names, y = C.coef.voom %*% t(y_limma)))

number = 1e5
voom_top = topTable(vfit2, number = number, coef = 3, p.value = 0.05)
voom_cor_top = topTable(vfit2_cor, number = number, coef = 3, p.value = 0.05)
DESeq_top = contrast.complex[contrast.complex$padj < 0.05, ]
GLMM_top = checkContrast[checkContrast$complex.P.adj < 0.05, ]
# Final result table
# soring gives same order for output from biomaRt as in dataFrame
name = sort(unique(c(rownames(voom_cor_top), 
                 rownames(voom_top), 
                 rownames(DESeq_top), 
                 rownames(GLMM_top))))


# extract gene names based on ensembl gene id
library('biomaRt')

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) #Hent det humane datasettet fra biomaRt

genes <- name # df$genes #Her er ID’ene du har, og som du ønsker å oversette

G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)

# add values not found
not_found = which(!(name %in% G_list$ensembl_gene_id))
G_list = rbind(G_list, data.frame("ensembl_gene_id" = name[not_found], "hgnc_symbol" = rep(NA, length(not_found))))
rownames(G_list) = G_list$ensembl_gene_id


# General base mean is extracted from DESeq2
# Order is decied by sum of adjusted p values, where model did not converge (NA) ajd.p-value is set to 1
# Alternative is to consider order (where not significant would be total genes in list)
compare = data.frame("names" = G_list[name, ]$hgnc_symbol, 
                     "base.mean" = contrast.complex[name, ]$baseMean, 
# voom without cor
           "voom.LFC" = round(topTable(vfit2, number = number, coef = 3)[name, ]$logF, 3), 
           "voom.adj.p" = topTable(vfit2, number = number, coef = 3)[name, ]$adj.P.Val, 
           "voom.place" = cbind(voom_top, 
                                "number" = 1:dim(voom_top)[1])[name, ]$number,
# voom with cor
           "voom.cor.LFC" = round(topTable(vfit2_cor, number = number, coef = 3)[name, ]$logFC, 3), 
           "voom.cor.adj.p" = topTable(vfit2_cor, number = number, coef = 3)[name, ]$adj.P.Val,  
           "voom.cor.place" = cbind(voom_cor_top, 
                                        "number" = 1:dim(voom_cor_top)[1])[name, ]$number,
# DESeq2
           "DESeq2.LFC" = round(contrast.complex[name, ]$log2FoldChange, 3),  
           "DESeq2.adj.p" = contrast.complex[name, ]$padj,
           "DESeq2.place" = cbind(DESeq_top, 
                                  "number" = 1:dim(DESeq_top)[1])[name, ]$number, 
# GLMM
           "GLMM.LFC" = round(checkContrast[name, ]$complex.glmm, 3), 
           "GLMM.adj.p" = checkContrast[name, ]$complex.P.adj, 
           "GLMM.place" = cbind(GLMM_top, 
                                "number" = 1:dim(GLMM_top)[1])[name, ]$number
           )
rownames(compare) = name
p.val = c("voom.adj.p", "voom.cor.adj.p", "DESeq2.adj.p", "GLMM.adj.p")
compare[p.val][is.na(compare[p.val])] = 1
compare = compare[order(rowSums(compare[p.val])), ]
compare[p.val][compare[p.val] == 1] = NA

write.csv(compare[c("names", "base.mean", "voom.LFC", "voom.cor.LFC", "DESeq2.LFC", "GLMM.LFC")], 
          "D:\\Essential folders\\Documents\\RStudio\\master\\gene-analysis\\code\\compare.csv")