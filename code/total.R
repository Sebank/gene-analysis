library(stringr)
library(ggvenn)
library(edgeR)
library(DESeq2)
library(lme4)
library(hexbin) 
library(apeglm)
library(readr)

#wald.test
library(aod)

#file location
fileloc = "Z:/Atle van Beelen Granlund/"
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

#finds indexes for New RNA samples, these might be bad
#New_RNA = grep("New RNA", pas_info_tot$Comment)
#names_for_referance = pas_info_tot$Sample_ID[New_RNA]

pas_info_tot = pas_info_tot[-removal, ]
count_tot = count_tot[, pas_info_tot$Sample_ID]

#extract relevant columns
pas_info <- data.frame("sampleid" = pas_info_tot$Sample_ID, 
                       "diseaseinflammation" = pas_info_tot$Sample_Group, 
                       "tissue" = pas_info_tot$Sample_Biosource)
#remove tissue from tissue
pas_info$tissue = substr(pas_info$tissue, 1, 5)

#using grep to remove final segments of string ending with irrelevant information
firstpart = str_split_fixed(pas_info$sampleid, "_", 2)[,1]
iFS = grep("*[iFS]", firstpart)
firstpart[iFS] = substr(firstpart[iFS], 1, nchar(firstpart[iFS]) - 1)
pas_info$sampleid = firstpart

#Assuming we follow tissue instead of inflammation
#remove i from disease information, as it should be part of tissue as well
i = grepl("*i", pas_info$diseaseinflammation)
pas_info$diseaseinflammation[i] = substr(pas_info$diseaseinflammation[i], 1, 
                                         nchar(pas_info$diseaseinflammation[i]) - 1)

#more a double check than anything else
#unique(pas_info$diseaseinflammation)

#make the variable disease inflammation two variables
#replaced F with H for "healthy"
AU = grep("..[AU]", pas_info$diseaseinflammation)
pas_info$disease = rep("H", length(pas_info$diseaseinflammation))
pas_info$disease[AU] = substr(pas_info$diseaseinflammation[AU], 1, 
                              nchar(pas_info$diseaseinflammation[AU]) - 1)


A = grep("A", pas_info$diseaseinflammation)
H = grep("F", pas_info$diseaseinflammation)
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

dim(count)
dim(pas_info)
table(pas_info$tissue,pas_info$inflammation)
table(table(pas_info$sampleid)) # 13 med 2 og 1 med 3 prøver

#change the names of count to the modified names from pas_info
#observe that the count was sorted based on the names of pas_info before manipulation
colnames(count) = pas_info$sampleid

#change F to H in diseaseinflammation and make variables factors
pas_info$diseaseinflammation[pas_info$diseaseinflammation == "F"] = "H"
pas_info$diseaseinflammation = factor(pas_info$diseaseinflammation)
pas_info$inflammation = factor(pas_info$inflammation)
pas_info$tissue = factor(pas_info$tissue)

###########################################################################
#limma/voom (model)
###########################################################################

#convert to relevant setup for limma

dge = DGEList(count)

#make as factors and order factors
pas_info$inflammation = factor(pas_info$inflammation, c("H", "U", "A"))
pas_info$tissue = factor(pas_info$tissue, c("Colon", "Ileum"))
table(pas_info$inflammation,pas_info$tissue)
#"remove counts that have zero or very low counts"
design = model.matrix(~ inflammation*tissue, data = pas_info)
keep = filterByExpr(round(count), design)
dge = dge[keep, , keep.lib.sizes=FALSE]
dim(dge$counts)
# 27534, before 58003 
#use TMM normalization method
dge = calcNormFactors(dge)



#limma vs voom

#limma
#"prior.count damp down the variances of logarithms of low counts"
logCPM = cpm(dge, log=TRUE)
dim(logCPM)

#plotMDS
plotMDS(logCPM, labels = substr(pas_info$tissue, 1, 1), 
        col = as.numeric(pas_info$diseaseinflammation))

C=rbind("(CU-CH)-(IU-IH)"=c(0,0,0,0,1,0), 
        "(CA-CH)-(IA-IH)"=c(0,0,0,0,0,1),
        "(CA-CU)-(IA-IU)"=c(0,0,0,0,1,-1))

#implement correlation
corfit = duplicateCorrelation(logCPM, design, block = pas_info$sampleid)

fit = lmFit(logCPM, design, block = pas_info$sampleid, correlation = corfit$consensus)

fit2=contrasts.fit(fit,t(C))

fit3 = eBayes(fit2, trend=TRUE)

table_fit3 = topTable(fit3, number = 10)

limma::plotMD(fit3, coef = 3)
limma::plotSA(fit3, col = "green")

#voom
v = voom(dge, design, plot=TRUE) 
plotMDS(v, labels = substr(pas_info$tissue, 1, 1), 
        col = as.numeric(pas_info$diseaseinflammation))
legend("top", legend = levels(pas_info$diseaseinflammation), col = 1:3, pch = 15)

voom_cor = duplicateCorrelation(v,design,block=pas_info$sampleid)

vfit = lmFit(v, design, block = pas_info$sampleid, correlation = voom_cor$consensus, weights = v$weights)
# inne i v er både en ny respons og vekter inn i GLS

C = rbind("IU"=c(0,0,0,0,1,0), #(CU-CH)-(IU-IH)
          "IA"=c(0,0,0,0,0,1), #(CA-CH)-(IA-IH)
          "IA minus IU"=c(0,0,0,0,-1,1)) #(CA-CU)-(IA-IU)
vfit2 = contrasts.fit(vfit, t(C))
vfit2 = eBayes(vfit2, trend = FALSE)
#for specific coef use coef = 1, 2 or 3, as contrast has three coefs
table_vfit2 = topTable(vfit2, number = 10)

ggplot() + geom_histogram(aes(vfit2$p.value[, 3]), binwidth = 0.01, fill = "blue")

###########################################################################
#deseq (model)
###########################################################################

dds = DESeqDataSetFromMatrix(countData = round(count), colData = pas_info, 
                             design = ~ inflammation*tissue)
dds$tissue
#pre processing
length(dds[, 1])
keep_DESeq = filterByExpr(counts(dds), design)
dds = dds[keep_DESeq, ]
length(dds[, 1])
#58003 vs 43358
#vs 27708

# equal to three part beneath
# dds <- DESeq(dds)

dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)
dds = nbinomWaldTest(dds)

attributes(dds@modelMatrix)$dimnames[[2]]
resultsNames(dds) # lists the coefficients

contrast.UI <- results(dds, contrast = C[1, ])
# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="Intercept", type="apeglm")

#for consideration
removal.UI = contrast.UI[is.na(contrast.UI$padj), ]
contrast.UI = contrast.UI[!is.na(contrast.UI$padj), ]
contrast.UI = contrast.UI[order(contrast.UI$padj), ]

contrast.AI <- results(dds, contrast = C[2, ])

#for consideration
removal.AI = contrast.AI[is.na(contrast.AI$padj), ]
contrast.AI = contrast.AI[!is.na(contrast.AI$padj), ]
contrast.AI = contrast.AI[order(contrast.AI$padj), ]


contrast.complex <- results(dds, contrast = C[3, ])

#for consideration
removal.complex = contrast.complex[is.na(contrast.complex$padj), ]
contrast.complex = contrast.complex[!is.na(contrast.complex$padj), ]
contrast.complex = contrast.complex[order(contrast.complex$padj), ]

#analyze data
ggplot(as(contrast.complex, "data.frame"), aes(x = pvalue)) + 
  geom_histogram(binwidth = 0.01, fill = "Blue", boundary = 0)

DESeq2::plotMA(contrast.complex, ylim = c(-10, 10))
#DESeq2::plotPCA(rlogTransformation(dds), intgroup = c("tissue", "diseaseinflammation"))
#plots of potential interest
ggplot() + geom_point(aes(estimateSizeFactorsForMatrix(count), colSums(count)))
#observation 88 is the smallest
sf = estimateSizeFactorsForMatrix(dds@assays@data$counts)
ncounts = t(t(dds@assays@data$counts)/sf)
ggplot() + geom_hex(aes(x = log(rowMeans(ncounts)), y = log(rowVars(ncounts)))) + 
  coord_fixed() + theme(legend.position = "none") + 
  geom_abline(slope = 1:2, color = c("forestgreen", "red"))
#green is variance for Poisson, red is for gamma

###########################################################################
#analysis of results
###########################################################################

#general venn diagram

#check coef with DESeq2, only considering with correlation
number = 1000000
venn.1 = list(#"limma" = rownames(topTable(fit3, number = number, coef = 1, p.value = 0.05)), 
             "voom" = rownames(topTable(vfit2, number = number, coef = 1, p.value = 0.05)),
            # "DESeq2" = rownames(contrast.UI[1:number, ]))
            "DESeq2" = rownames(contrast.UI[contrast.UI$padj < 0.05, ]))
ggvenn(venn.1, show_percentage = FALSE)

venn.2 = list("limma" = rownames(topTable(fit3, number = number, coef = 2, p.value = 0.05)), 
              # "voom" = rownames(topTable(vfit2, number = number, coef = 2, p.value = 0.05)),
              # "DESeq2" = rownames(contrast.AI[1:number, ]))
              "DESeq2" = rownames(contrast.AI[contrast.AI$padj < 0.05, ]))
ggvenn(venn.2, show_percentage = FALSE)

venn.3 = list(# "limma" = rownames(topTable(fit3, number = number, coef = 3, p.value = 0.05)), 
               "voom" = rownames(topTable(vfit2, number = number, coef = 3, p.value = 0.05)),
              # "DESeq2" = rownames(contrast.complex[1:number, ]))
            "DESeq2" = rownames(contrast.complex[contrast.complex$padj < 0.05, ]))
ggvenn(venn.3, show_percentage = FALSE)


###########################################################################
#Differences between limma and DESeq
###########################################################################

#pre processing
ggvenn(list("limma" = names(keep), "DESeq" = names(keep_DESeq)), show_percentage = FALSE)

###########################################################################
#personal comments, code to quality check changes
###########################################################################

# plot av estimerte og selv estimert mean variance relationship
ggplot() + 
  geom_smooth(aes(x = fit3_alt$Amean, y = sqrt(fit3_alt$sigma)), method = "gam", formula = y ~ s(x, bs = "cs")) + 
  geom_point(aes(x = fit3_alt$Amean, y = (fit3_alt$s2.prior)^(1/4)))

#f is a measurement of how much of the data is relevant for each function value, just random guess so far
curve = lowess(x = fit3_alt$Amean, y = sqrt(fit3_alt$sigma), f = 1/4)

ggplot() + 
  geom_point(aes(curve$x, curve$y), col = "red") + 
  geom_point(aes(x = fit3$Amean, y = (fit3$s2.prior)^(1/4) ), col = "blue")



###########################################################################
#Messing with DESeq2 intro
###########################################################################
# create a new DESeq2 object and manipulate it for further examination
dataSet = DESeqDataSetFromMatrix(
  countData = round(count), 
  colData = pas_info,
  design = ~ tissue * inflammation
)

#pre processing
length(dataSet[, 1])
keep_DESeq = filterByExpr(counts(dataSet), design)
dataSet = dataSet[keep_DESeq, ]
length(dataSet[, 1])
#58003 vs 43358
#vs 27708


#full, slow
dataSet = estimateSizeFactors(dataSet)
dataSet = estimateDispersions(dataSet)
dataSet = nbinomWaldTest(dataSet) # betaPrior = TRUE

attributes(dataSet@modelMatrix)$dimnames[[2]]
resultsNames(dataSet) # lists the coefficients

res = results(dataSet, contrast = C[1, ])#LFCshrink, type= "apeglm")
summary(res)
res[order(res$padj), ] %>% head()
# results are considered by four plots, histogram of p values, MA plot, ordination plot and a heatmap

ggplot(as(res, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)

DESeq2::plotMA(res, ylim = c(-10, 10))

# plot rlog transform as PCA to get an idea of variance
# # # Slow
# pas_rlog = rlogTransformation(dataSet)
# DESeq2::plotPCA(pas_rlog, intgroup = c("tissue", "inflammation"))
# 
# library("pheatmap")
# select = order(rowMeans(assay(pas_rlog)), decreasing = TRUE)[1:30]
# pheatmap( assay(pas_rlog)[select, ], 
#           scale = "row",
#           annotation_col = as.data.frame(
#             colData(pas_rlog)[, c("tissue", "inflammation")]
#             )
#           )

# sanity checks
# plot of library sizes for our observations
cs <- colSums(assay(dataSet, "counts"))
hist(cs/1e6, col="grey", border="white",
     main="", xlab="column sums (per million)")
# why are there counts with less than 5? There is just 1, TT38, colsum 306740

###########################################################################
#deseq (model)
###########################################################################

contrast.UI <- results(dds, contrast = C[1, ])
# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="Intercept", type="apeglm")

#for consideration
removal.UI = contrast.UI[is.na(contrast.UI$padj), ]
contrast.UI = contrast.UI[!is.na(contrast.UI$padj), ]
contrast.UI = contrast.UI[order(contrast.UI$padj), ]

contrast.AI <- results(dds, contrast = C[2, ])

#for consideration
removal.AI = contrast.AI[is.na(contrast.AI$padj), ]
contrast.AI = contrast.AI[!is.na(contrast.AI$padj), ]
contrast.AI = contrast.AI[order(contrast.AI$padj), ]


contrast.complex <- results(dds, contrast = C[3, ])

#for consideration
removal.complex = contrast.complex[is.na(contrast.complex$padj), ]
contrast.complex = contrast.complex[!is.na(contrast.complex$padj), ]
contrast.complex = contrast.complex[order(contrast.complex$padj), ]

#analyze data
ggplot(as(contrast.complex, "data.frame"), aes(x = pvalue)) + 
  geom_histogram(binwidth = 0.01, fill = "Blue", boundary = 0)

DESeq2::plotMA(contrast.complex)
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


# code for Atle meeting
number = 10000
venn.limma = list("IU" = rownames(topTable(vfit2, number = number, coef = 1, p.value = 0.05)), 
                  "IA" = rownames(topTable(vfit2, number = number, coef = 2, p.value = 0.05)),
                  "IA - IU" = rownames(topTable(vfit2, number = number, coef = 3, p.value = 0.05))
                  )
ggvenn(venn.limma, show_percentage = FALSE)

venn.DESeq2 = list("IU" = rownames(contrast.UI[contrast.UI$padj < 0.05, ]),
                   "IA" = rownames(contrast.AI[contrast.AI$padj < 0.05, ]),
                   "IA - IU" = rownames(contrast.complex[contrast.complex$padj < 0.05, ])
                   )
ggvenn(venn.DESeq2, show_percentage = FALSE)


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
y_DESeq2 = coef(dataSet[name])
ggplot()  + 
  geom_point(aes(x = factor(pas_info$comment, levels = c("CH", "CU", "CA", "IH", "IU", "IA")), y = unlist(count[name, ])), alpha = 0.3) + 
  geom_point(aes(x = C.coef.names, y = 2^(C.coef.DESeq2 %*% t(y_DESeq2))), col = "red")
#plot of all observations


limma_baseline = eBayes(vfit, trend = FALSE)
y_limma = coef(limma_baseline[name, ])
ggplot() + geom_point(aes(x = C.coef.names, y = C.coef.voom %*% t(y_limma)))

# "ENSG00000109099"
name = rownames(topTable(vfit2, number = 1, coef = 3))
y_DESeq2 = coef(dataSet[name])
ggplot() + geom_point(aes(x = C.coef.names, y = C.coef.DESeq2 %*% t(y_DESeq2)))

limma_baseline = eBayes(vfit, trend = FALSE)
y_limma = coef(limma_baseline[name, ])
ggplot() + geom_point(aes(x = C.coef.names, y = C.coef.voom %*% t(y_limma)))

# biomart extract names of genes
# biomart, find gen
# library('biomaRt')
# 
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) #Hent det humane datasettet fra biomaRt
# 
# genes <- rownames(temp_limma) # df$genes #Her er ID’ene du har, og som du ønsker å oversette
# 
# G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
# 
# merge(df,G_list,by.x="gene",by.y="ensembl_gene_id") #Slå opp i biomaRt objekt for å hente ut gensymbol som tilsvarer ID.

#messing with GLMM
library(MASS)
name = names(table(pas_info$sampleid))[ table(pas_info$sampleid) > 1]
group = matrix(0, length(pas_info$sampleid), length(name))
colnames(group) = name

for(i in 1:length(name)){
  group[which(pas_info$sampleid == name[i]), i] = 1
}

# dispersion = 1/theta 
dispersion = dataSet@rowRanges@elementMetadata@listData$dispersion
genes = dataSet@assays@data@listData$counts

# 1189 is first error
i = 5
GLMM = glmer(
  formula = genes[i, ] ~ - 1 + design + (1|pas_info$sampleid), 
  family = MASS::negative.binomial(link = "log", theta=1/dispersion[i]), 
  offset = log(dataSet@colData@listData$sizeFactor)
  )
fisher = -GLMM@optinfo$derivs$Hessian

GLM = glm(
  formula = genes[i, ] ~ -1 + design, 
  family = MASS::negative.binomial(link = "log", theta=1/dispersion[i]), 
  offset = log(dataSet@colData@listData$sizeFactor)
)

contrast.UI[rownames(genes)[i], ]
contrast.AI[rownames(genes)[i], ]

# computes the correlation between two observations based on offset + estimated coefficients (mu),
# standard deviation of mixed effect (tau) and dispersion estimate (alpha) for negative binomial distribution
correlation = function(mu, alpha, tau){
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

# using well established estimates instead of manual

# create checkContrast (long so import instead)
n = length(rownames(genes))

checkContrast = data.frame("UI" = contrast.UI[rownames(genes), ][, 2], "UI.glmm" = rep(NA, n), "UI.GLM" = rep(NA, n), "UI.P" = rep(NA, n), "UI.GLM.P" = rep(NA, n),
                           "AI" = contrast.AI[rownames(genes), ][, 2], "AI.glmm" = rep(NA, n), "AI.GLM" = rep(NA, n), "AI.P" = rep(NA, n), "AI.GLM.P" = rep(NA, n),
                           "mixed.sd" = rep(NA, n))
rownames(checkContrast) = rownames(genes)
cova = vector(mode = "double", length = dim(group)[2])
corre = matrix(nrow = n, ncol = dim(group)[2])

for(i in 1:dim(genes)[1]){
  if(!is.na(checkContrast$UI[i]) & !is.na(checkContrast$AI[i])){
    tryCatch(expr = {
      temp = glmer(
        formula = genes[i, ] ~ -1 + design + (1|pas_info$sampleid),
        family = MASS::negative.binomial(link = "log", theta=1/dispersion[i]),
        offset = log(dataSet@colData@listData$sizeFactor)
      )
      checkContrast$UI.glmm[i] = temp@beta[5]*log2(exp(1))
      checkContrast$AI.glmm[i] = temp@beta[6]*log2(exp(1))
      checkContrast$mixed.sd[i] = temp@optinfo$val[1]

      P = coef(summary(temp))[, 4]

      checkContrast$UI.P[i] = P[5]
      checkContrast$AI.P[i] = P[6]

      for(l in 1:dim(group)[2]){
        j = which(group[, l] == 1)[1]
        k = which(group[, l] == 1)[2]

        tau = temp@pp$theta
        mu_cov = fitted(temp)[c(j, k)]
        corre[i, l] = correlation(mu = mu_cov, alpha = dispersion[i], tau = tau)
      }
    }, warning = function(cond){
      print(cond)
    }, error = function(cond){
      print(cond)
    },
    silent = TRUE
    )


    tryCatch(expr = {temp = glm(
      formula = genes[i, ] ~ -1 + design,
      family = MASS::negative.binomial(link = "log", theta=1/dispersion[i]),
      offset = log(dataSet@colData@listData$sizeFactor)
    )

    P = coef(summary(temp, dispersion = 1))[, 4]

    checkContrast$UI.GLM.P[i] = P[5]
    checkContrast$AI.GLM.P[i] = P[6]

    checkContrast$UI.GLM[i] = temp$coefficients[5]*log2(exp(1))
    checkContrast$AI.GLM[i] = temp$coefficients[6]*log2(exp(1))
    }, warning = function(cond){
      print(cond)
    })
  }
}

# write.csv(checkContrast, "D:\\Essential folders\\Documents\\RStudio\\master\\gene-analysis\\code\\checkContrast.csv")
# write.csv(corre, "D:\\Essential folders\\Documents\\RStudio\\master\\gene-analysis\\code\\corre.csv")

# These require the location of the git directory
corre = read.csv(paste(git, "\\corre.csv", sep = ""), header = TRUE, row.names = 1)
checkContrast = read.csv(paste(git, "\\checkContrast.csv", sep = ""), header = TRUE, row.names = 1)
# checkContrast = checkContrast[order(checkContrast$mixed.sd, decreasing = TRUE), ]

#make comparable to contrast from DESeq2
contrast.UI.glmm = data.frame("UI" = checkContrast$UI.glmm[checkContrast$UI.P != 0], "P" = checkContrast$UI.P[checkContrast$UI.P != 0])
rownames(contrast.UI.glmm) = rownames(checkContrast[checkContrast$UI.P != 0, ])
contrast.UI.glmm = contrast.UI.glmm[order(contrast.UI.glmm$P), ]

# correlation of mixed effects 
summary(GLMM)

# plot mean correlation for a gene against standard deviation of the random effect for a GLMM
ggplot() + geom_point(aes(x = checkContrast$mixed.sd, y = rowMeans(corre)), alpha = 0.1)

# Initial idea, requires more explanation than it is worth
# ggplot() + geom_point(aes(x = checkContrast$mixed.sd[!(checkContrast$mixed.sd < 0.3 & 
#                                                          abs(checkContrast$UI.glmm - checkContrast$UI.GLM) > 1)], 
#                           y = (checkContrast$UI.glmm - checkContrast$UI.GLM)[!(checkContrast$mixed.sd < 0.3 & 
#                                                                                  abs(checkContrast$UI.glmm - checkContrast$UI.GLM) > 1)]
#                           ), alpha = 0.1)
# checks of convergence
# 9, 12, 22, 23 is safe
# 11635 was observed to also not converge for GLM, the only one where sd was estimated to a larger value
# 21681, sd = 0
# 23430, GLMM does converge, but no GLM
# 10643, neither converges
# which((checkContrast$UI.glmm - checkContrast$UI.GLM) > 1)
# i = which((checkContrast$UI.glmm - checkContrast$UI.GLM) > 5)[1]
i = which(rownames(genes) == "ENSG00000164035")
GLMM = glmer(
  formula = genes[i, ] ~ 1 + design[, -1] + (1|pas_info$sampleid),
  family = MASS::negative.binomial(link = "log", theta=1/dispersion[i]),
  offset = log(dataSet@colData@listData$sizeFactor)
)

GLM = glm(
  formula = genes[i, ] ~ -1 + design,
  family = MASS::negative.binomial(link = "log", theta=1/dispersion[i]),
  offset = log(dataSet@colData@listData$sizeFactor)
)

summary(GLM, dispersion = 1)
summary(GLMM)

# remove samples that did not converge and perform checks
include = checkContrast$AI.GLM != 0 & !is.na(corre[, 1])
# final estimate for difference between coef when GLMM converged (there might be some GLM which did not converge, but also might not, should be fixed)
ggplot() + geom_point(aes(x = checkContrast$mixed.sd[include], 
                          y = rowMeans(corre[include, ])), alpha = 0.1) + 
  labs(x = "standard deviation of mixed effects", y = "mean of correlations for a gene", title = "correlation against standard deviation")
ggplot() + geom_point(aes(x = rowMeans(corre[include, ]), 
                          y = (checkContrast$UI.glmm - checkContrast$UI.GLM)[include]), alpha = 0.1) +
  labs(x = "Average correlation", y = "GLMM effect - GLM effect", title = "UI plot for converged genes")
ggplot() + geom_point(aes(x = checkContrast$UI.P[include], 
                          y = checkContrast$UI.GLM.P[include]), alpha = 0.1) +
  labs(x = "P from GLMM", y = "P from GLM", title = "comparison of estimated P values")

# alt checks
ggplot() + geom_point(aes(x = rowMeans(corre[include, ]), 
                          y = abs(checkContrast$UI.glmm - checkContrast$UI.GLM)[include]), alpha = 0.1) +
  labs(x = "mixed effect average correlation", y = "absolute value of GLMM effect - GLM effect", title = "UI plot for converged genes")
# ENSG00000152977 fails for some reason
ggplot() + geom_point(aes(x = rowMeans(corre[include, ]), 
                          y = (checkContrast$UI.glmm/checkContrast$UI.GLM)[include]), alpha = 0.1) +
  labs(x = "mixed effect average correlation", y = "GLMM effect divided by GLM effect", title = "UI plot for converged genes")
# extract names to get contrast from DESeq2, then consider what is NA for both
names = rownames(checkContrast[include, ])
ggplot() + geom_point(aes(x = checkContrast[names, ]$UI.P, 
                          y = contrast.UI[names, ]$pvalue), alpha = 0.1) + 
  labs(x = "p value for GLMM", y = "p value for DESeq2", title = "compare P values between models and DESeq2")
ggplot() + geom_point(aes(x = checkContrast[names, ]$UI.GLM.P, 
                          y = contrast.UI[names, ]$pvalue), alpha = 0.1)+ 
  labs(x = "p value for GLM", y = "p value for DESeq2", title = "compare P values between models and DESeq2")
ggplot() + geom_histogram(aes(rowMeans(corre[include, ]), y = ..count../sum(..count..)), binwidth = 0.01)
ggplot() + geom_point(aes(x = rowMeans(corre[include, ]), 
                          y = (checkContrast$UI.P[include] - checkContrast$UI.GLM.P[include])), alpha = 0.3) +
  labs(x = "mean of correlation", y = "difference in P value between GLMM and GLM", title = "difference of P value against correlation for UI")

# find percent for samples with estimated correlation greeater than x
temp = rowMeans(corre[include, ])
length(temp[temp >= 0.1])/length(temp)
length(temp[temp >= 0.05])/length(temp)

# simulation and comparison between GLMM and GLM
g = 5
GLMM = glmer(
  formula = genes[g, ] ~ - 1 + design + (1|pas_info$sampleid), 
  family = MASS::negative.binomial(link = "log", theta=1/dispersion[g]), 
  offset = log(dataSet@colData@listData$sizeFactor)
)
fisher = -GLMM@optinfo$derivs$Hessian

GLM = glm(
  formula = genes[g, ] ~ -1 + design, 
  family = MASS::negative.binomial(link = "log", theta=1/dispersion[g]), 
  offset = log(dataSet@colData@listData$sizeFactor)
)

n = 100
beta = GLMM@beta
tau = GLMM@pp$theta


tissue = c()
inflammation = c()
pers = c()

for(i in 1:n){
  pers = c(pers, toString(i), toString(i))
  u = runif(1)
  if(u <= 1/7){
    tissue = c(tissue, "Colon")
    inflammation = c(inflammation, "H")
    tissue = c(tissue, "Ileum")
    inflammation = c(inflammation, "H")
  }
  else if(u <= 2/7){
    tissue = c(tissue, "Colon")
    inflammation = c(inflammation, "U")
    tissue = c(tissue, "Colon")
    inflammation = c(inflammation, "A")
  }
  else if(u <= 3/7){
    tissue = c(tissue, "Colon")
    inflammation = c(inflammation, "U")
    tissue = c(tissue, "Ileum")
    inflammation = c(inflammation, "U")
  }
  else if(u <= 4/7){
    tissue = c(tissue, "Colon")
    inflammation = c(inflammation, "U")
    tissue = c(tissue, "Ileum")
    inflammation = c(inflammation, "A")
  }
  else if(u <= 5/7){
    tissue = c(tissue, "Colon")
    inflammation = c(inflammation, "A")
    tissue = c(tissue, "Ileum")
    inflammation = c(inflammation, "U")
  }
  else if(u <= 6/7){
    tissue = c(tissue, "Colon")
    inflammation = c(inflammation, "A")
    tissue = c(tissue, "Ileum")
    inflammation = c(inflammation, "A")
  }
  else{
    tissue = c(tissue, "Ileum")
    inflammation = c(inflammation, "U")
    tissue = c(tissue, "Ileum")
    inflammation = c(inflammation, "A")
  }
}

tissue = factor(tissue, levels = c("Colon", "Ileum"))
inflammation = factor(inflammation, c("H", "U", "A"))

m = 10
est = data.frame("GLMM" = rep(NA, m), "GLM" = rep(NA, m), "se" = rep(NA, m), "GLM.P" = rep(NA, m), "GLMM.P" = rep(NA, m))

for(i in 1:m){
  u = rnorm(100, sd = tau)
  Y = rnbinom(n = 2*n, size = 1/dispersion[g], mu = exp(model.matrix(~tissue * inflammation) %*% beta + model.matrix(~ -1 + pers) %*% u))
  
  GLMM_compare = glmer(formula = Y ~ - 1 + model.matrix(~tissue * inflammation) + (1|pers), 
                       family = MASS::negative.binomial(link = "log", theta = 1/dispersion[g]))
  GLM_compare = glm(formula = Y ~ - 1 + model.matrix(~tissue * inflammation), 
                    family = MASS::negative.binomial(link = "log", theta = 1/dispersion[g]))
  
  est$GLM[i] = coef(GLM_compare)[5]
  est$GLMM[i] = fixef(GLMM_compare)[5]
  est$se[i] = GLMM_compare@pp$theta
  # est$GLM.P[i] = coef(summary(GLM_compare, dispersion = 1))[6, 4]
  # est$GLMM.P[i] = coef(summary(GLMM_compare))[6, 4]
  est$GLM.P[i] = wald.test(vcov(summary(GLM_compare, dispersion = 1)), coef(GLM_compare), Terms = 6, H0 = beta[6])$result$chi2[3]
  est$GLMM.P[i] = wald.test(vcov(summary(GLMM_compare)), GLMM_compare@beta, Terms = 6, H0 = beta[6])$result$chi2[3]
}
ggplot(aes(y = ..count../sum(..count..)), data = est) +
  geom_histogram(aes(x = GLM.P, fill = "GLM"), alpha = 0.7, binwidth = 0.01) +
  geom_histogram(aes(x = GLMM.P, fill = "GLMM"), alpha = 0.3, binwidth = 0.01)

# read values
# est = read.csv(paste(git, "\\est.csv", sep = ""), header = TRUE, row.names = 1)
wald.test(vcov(summary(GLM, dispersion = 1)), coef(GLM), Terms = 6, H0 = 1)$result$chi2[3]

# messing with se for negative binomial
# Mike love uses e and then converts restults to 2
i = 1
GLM = glm(
  formula = genes[i, ] ~ -1 + design,
  family = MASS::negative.binomial(link = "log", theta=1/dispersion[i]),
  offset = log(dataSet@colData@listData$sizeFactor)
)
vcov(summary(GLM, dispersion = 1))
X = design
mu_hat = fitted(GLM)
W = mu_hat/(1 + dispersion[i] * mu_hat)
sigma = solve(t(X) %*% diag(W) %*% X)
sigma
vcov(summary(GLM, dispersion = 1))/sigma

GLMM = glmer(
  formula = genes[i, ] ~ - 1 + design + (1|pas_info$sampleid), 
  family = MASS::negative.binomial(link = "log", theta=1/dispersion[i]), 
  offset = log(dataSet@colData@listData$sizeFactor)
)
vcov(GLMM)
tau = GLMM@pp$theta
mu_hat = fitted(GLMM)
# initial thought
# W = mu_hat^2/(mu_hat + dispersion[i] * mu_hat^2 * exp(tau^2) + mu_hat^2 * (exp(tau^2) - 1))
# in reality W = mu_hat^2/(mu_hat + dispersion[i] * mu_hat^2), same as for GLM, and consider a seperate one for random effects
W = weights(GLMM, type = "working")


sigma = solve(t(X) %*% diag(W) %*% X)
sigma
vcov(GLMM)/sigma

res.ape = lfcShrink(dds, coef = 5, type = "apeglm")
library(insight)
library(performance)

ICC = c()
ICC_conditional = c()
performance::icc(GLMM)
for(i in which(rowMeans(corre) > 1e-6)){
  GLMM = glmer(
    formula = genes[i, ] ~ -1 + design + (1|pas_info$sampleid),
    family = MASS::negative.binomial(link = "log", theta=1/dispersion[i]),
    offset = log(dataSet@colData@listData$sizeFactor)
  )
  ICC = c(ICC, performance::icc(GLMM)$ICC_adjusted)
  ICC_conditional = c(ICC_conditional, performance::icc(GLMM)$ICC_conditional)
}

ggplot() + geom_point(aes(x = ICC, y = rowMeans(corre)[which(rowMeans(corre) > 1e-6)]), alpha = 0.1) +
  labs(x = "Adjusted ICC", y = "mean of correlation (self computed)")
