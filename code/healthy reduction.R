library(stringr)
library(ggvenn)
library(edgeR)
library(DESeq2)
library(lme4)
library(hexbin) 
library(apeglm)

#file location
fileloc = "../Atle van Beelen Granlund/"


###########################################################################
#pre-processing
###########################################################################

# info about the samples to connect to disease, tissue, inflammation, pasid

pas_info_tot <- read.delim(paste(fileloc,'sample_info.txt',sep=""), 
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

# Remove H altogether for comparison
count = count[, pas_info$disease != "H"]
pas_info = pas_info[pas_info$disease != "H", ]


#change F to H in diseaseinflammation and make variables factors
pas_info$diseaseinflammation = factor(pas_info$diseaseinflammation)
pas_info$inflammation = factor(pas_info$inflammation)
pas_info$tissue = factor(pas_info$tissue)



###########################################################################
#limma/voom (model)
###########################################################################

#convert to relevant setup for limma

dge = DGEList(count)

#make as factors and order factors
pas_info$inflammation = factor(pas_info$inflammation, c("U", "A"))
pas_info$tissue = factor(pas_info$tissue, c("Colon", "Ileum"))
table(pas_info$inflammation,pas_info$tissue)
#"remove counts that have zero or very low counts"
design = model.matrix(~ inflammation*tissue, data = pas_info)
keep = filterByExpr(dge, design)
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

C=rbind("(CA-CU)-(IA-IU)"=c(0,0,0,1))

#implement correlation
corfit = duplicateCorrelation(logCPM, design, block = pas_info$sampleid)

fit = lmFit(logCPM, design, block = pas_info$sampleid, correlation = corfit$consensus)

fit_alt = lmFit(logCPM, design)

fit2=contrasts.fit(fit,t(C))
fit2_alt = contrasts.fit(fit_alt, t(C))

fit3 = eBayes(fit2, trend=TRUE)
fit3_alt = eBayes(fit2_alt, trend = TRUE)

table_fit3 = topTable(fit3, number = 10)
table_fit3_alt = topTable(fit3_alt, number = 10)

limma::plotMD(fit3, coef = 1)
limma::plotSA(fit3, col = "green")

#voom
v = voom(dge, design, plot=TRUE) 
plotMDS(v, labels = substr(pas_info$tissue, 1, 1), 
        col = as.numeric(pas_info$diseaseinflammation))
legend("top", legend = levels(pas_info$diseaseinflammation), col = 1:3, pch = 15)

voom_cor = duplicateCorrelation(v,design,block=pas_info$sampleid)

vfit = lmFit(v, design, block = pas_info$sampleid, correlation = voom_cor$consensus, weights = v$weights)
vfit_alt = lmFit(v, design, weights = v$weights)
# inne i v er både en ny respons og vekter inn i GLS

vfit2 = contrasts.fit(vfit, t(C))
vfit2_alt = contrasts.fit(vfit_alt, t(C))
vfit2 = eBayes(vfit2, trend = FALSE)
vfit2_alt = eBayes(vfit2_alt, trend = FALSE)
#for specific coef use coef = 1, 2 or 3, as contrast has three coefs
table_vfit2 = topTable(vfit2, number = 10)
table_vfit2_alt = topTable(vfit2_alt, number = 10)

ggplot() + geom_histogram(aes(vfit2$p.value[, 1]), binwidth = 0.01, fill = "blue")

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

# beneath
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

#analyze data
ggplot(as(contrast.UI, "data.frame"), aes(x = pvalue)) + 
  geom_histogram(binwidth = 0.01, fill = "Blue", boundary = 0)

DESeq2::plotMA(contrast.UI)
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

#double check that adj.P.val is used
tail(topTable(vfit2_alt, number = 10000, p.value = 0.05))


#check coef with DESeq2, only considering with correlation
number = 1000000
venn.1 = list(#"limma" = rownames(topTable(fit3, number = number, coef = 1, p.value = 0.05)), 
               "voom" = rownames(topTable(vfit2, number = number, coef = 1, p.value = 0.05)),
              # "DESeq2" = rownames(contrast.UI[1:number, ]))
              "DESeq2" = rownames(contrast.UI[contrast.UI$padj < 0.05, ]))
ggvenn(venn.1, show_percentage = FALSE)

###########################################################################
#Differences between limma and DESeq
###########################################################################

#pre processing
ggvenn(list("limma" = names(keep), "DESeq" = names(keep_DESeq)), show_percentage = FALSE)

###########################################################################
#personal comments, code to quality check changes
###########################################################################

#double check if changes made above
#unique(pas_info$diseaseinflammation[pas_info$tissue=="Ileum tissue"])
#pas_info$sampleid[pas_info$tissue=="Ileum tissue" & pas_info$diseaseinflammation == "CDA"]
#pas_info$sampleid[pas_info$tissue=="Ileum tissue" & pas_info$diseaseinflammation == "CDU"]

#unique(pas_info$diseaseinflammation[pas_info$tissue=="Colon tissue"])
#pas_info$sampleid[pas_info$tissue=="Colon tissue" & pas_info$diseaseinflammation == "CDAi"]
#pas_info$sampleid[pas_info$tissue=="Colon tissue" & pas_info$diseaseinflammation == "CDUi"]

#finds indexes for New RNA samples, these might be bad
#New_RNA = grep("New RNA", pas_info_tot$Comment)
#names_for_referance = pas_info_tot$Sample_ID[New_RNA]

#more a double check than anything else
#unique(pas_info$diseaseinflammation)

# #finds tables for comparing how much data we have for each group
# table(pas_info$disease)
# table(pas_info$inflammation)
# table(pas_info$tissue,pas_info$inflammation)
# 
# table(table(pas_info$sampleid))


# cbind(pas_info$sampleid, colnames(count))
# 
# #finne konsepter fra data
# 
# summary(dge)
# dge
# dim(dge)
# # 27672 110
# dim(dge$counts)
# dim(dge$samples)
# dge$samples
# 
# 
# fit0=lm(logCPM[1,]~design-1)
# fit0_2 = lmer(logCPM[1,]~design-1 + (1 | pas_info$sampleid))
# dim(fit0$coefficients)
# dim(C)
# fit0.C=C%*%matrix(fit0$coefficients,ncol=1)
# fit0.C
# fit0.Ccov=C%*%vcov(fit0)%*%t(C)
# fit0.Ccov
# fcontrast3=fit0.C[3,1]^2/fit0.Ccov[3,3]
# fcontrast2=fit0.C[2,1]^2/fit0.Ccov[2,2]
# pf(fcontrast2,1,fit0$df.residual,lower.tail=FALSE)
# pf(fcontrast3,1,fit0$df.residual,lower.tail=FALSE)
# 
# dim(fit$coefficients)
# fit$coefficients[1,]
# summary(fit0)
# 
# str(fit2$contrasts)
# names(fit2)
# names(fit)
# dim(fit2$coefficients)
# fit$coefficients[1,]
# 
# # leter etter pverdier og testobs
# summary(fit0)
# 
# fit3$t[1,2]*sqrt(fit3$s2.post[1])*fit3$stdev.unscaled[1,2]
# fit0$coefficients[6]
# fit3$s2.prior[1]
# fit3$s2.post[1]
# 
# plot(fit3$Amean,fit3$s2.prior)
# plot(fit3$Amean,sqrt(fit3$s2.post),pch=".")
# lines(lowess(fit3$Amean,sqrt(fit3$s2.post)),col=2)
# 
# names(fit3)
# dim(fit3$p.value)
# fit3$p.value[1,]
# fit3$stdev.unscaled[1,]
# fit3$sigma[1]
# fit3$Amean[1]
# 
# str(v)
# hj=getEAWP(v)
# 
# 

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
#vs 27708 (27581 after low count removal)


#full, slow
dataSet = estimateSizeFactors(dataSet)
dataSet = estimateDispersions(dataSet)
dataSet = nbinomWaldTest(dataSet)

attributes(dataSet@modelMatrix)$dimnames[[2]]
resultsNames(dataSet) # lists the coefficients

res = results(dataSet, contrast = C[1, ])#LFCshrink, type= "apeglm")
summary(res)
res[order(res$padj), ] %>% head()
# results are considered by four plots, histogram of p values, MA plot, ordination plot and a heatmap

ggplot(as(res, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)

DESeq2::plotMA(dataSet)

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
cs <- colSums(assay(dataSet, "counts"))
hist(cs/1e6, col="grey", border="white",
     main="", xlab="column sums (per million)")
# why are there counts with less than 5? There is just 1, TT38, colsum 306740

# copy paste
cts <- assay(dataSet, "counts")[,c(1,4)] # choose random samples from different groups and select genes that have more than 5 as count (to give meaningful results)
idx <- rowSums(cts >= 5) == 2
cts <- cts[idx,]
p <- sweep(cts, 2, colSums(cts), "/")

mean.log.p <- rowMeans(log10(p))
hist(mean.log.p,  col="grey", border="white",
     main="", xlab="log10 proportion (geometric mean)")
#observe that most genes fall in the range -7 to -4 (log 10 proportion)

lfc <- log2(p[,2]) - log2(p[,1])
hist(lfc[between(mean.log.p, -6, -4)], 
     breaks=seq(-10,10,by=.1),
     xlim=c(-5,5),
     col="grey50", border="white",
     main="", xlab="log2 fold change of proportion")

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

# "ENSG00000155380"
y_DESeq2 = coef(dataSet[rownames(contrast.complex[1, ])])
ggplot() + geom_point(aes(x = C.coef.names, y = C.coef.DESeq2 %*% t(y_DESeq2)))

limma_baseline = eBayes(vfit, trend = TRUE)
y_limma = coef(limma_baseline[rownames(contrast.complex[1, ]), ])
ggplot() + geom_point(aes(x = C.coef.names, y = C.coef.voom %*% t(y_limma)))

# "ENSG00000109099"
y_DESeq2 = coef(dataSet[rownames(topTable(vfit2, number = 1, coef = 3))])
ggplot() + geom_point(aes(x = C.coef.names, y = C.coef.DESeq2 %*% t(y_DESeq2)))

limma_baseline = eBayes(vfit, trend = TRUE)
y_limma = coef(limma_baseline[rownames(topTable(vfit2, number = 1, coef = 3)), ])
ggplot() + geom_point(aes(x = C.coef.names, y = C.coef.voom %*% t(y_limma)))