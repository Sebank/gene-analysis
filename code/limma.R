library(stringr)
library(ggvenn)
library(edgeR)
library(DESeq2)
library(lme4)

#file location
fileloc = "Z:/Atle van Beelen Granlund/"


###########################################################################
#pre-processing
###########################################################################

# info about the samples to connect to disease, tissue, inflammation, pasid

pas_info_tot <- read.delim(paste(fileloc,'sample_info.txt',sep=""),sep="\t", header = TRUE, row.names = 1)
count_tot = read.delim(paste(fileloc, "quasi.gene.quant", sep = ""), sep = "\t", check.names = FALSE)
#remove data that is considered bad
removal = c(which(pas_info_tot$Sample_ID == "3115i"), which(pas_info_tot$Sample_ID == "3073i_new"))
#this has been run twice for some reason
removal = c(removal, which(pas_info_tot$Sample_ID == "3002i"))

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
pas_info$diseaseinflammation[i] = substr(pas_info$diseaseinflammation[i], 1, nchar(pas_info$diseaseinflammation[i]) - 1)

#more a double check than anything else
#unique(pas_info$diseaseinflammation)

#make the variable disease inflammation two variables
#replaced F with H for "healthy"
AU = grep("..[AU]", pas_info$diseaseinflammation)
pas_info$disease = rep("H", length(pas_info$diseaseinflammation))
pas_info$disease[AU] = substr(pas_info$diseaseinflammation[AU], 1, nchar(pas_info$diseaseinflammation[AU]) - 1)


A = grep("A", pas_info$diseaseinflammation)
H = grep("F", pas_info$diseaseinflammation)
pas_info$inflammation = rep("U", length(pas_info$diseaseinflammation))
pas_info$inflammation[A] = "A"
pas_info$inflammation[H] = "H"

#remove diseaseinflamation (never run this twice!!!)
# kanskje lurt å beholde den kolonnen?
di=pas_info[,2]
pas_info = pas_info[-2]

# #finds tables for comparing how much data we have for each group
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
keep = filterByExpr(dge, design)
dge = dge[keep, , keep.lib.sizes=FALSE]
dim(dge$counts)
# 27672, before 58003 
#use TMM normalization method
dge = calcNormFactors(dge)



#limma vs voom

#limma
#"prior.count damp down the variances of logarithms of low counts"
logCPM = cpm(dge, log=TRUE, prior.count=3)
dim(logCPM)

C=rbind("(CU-CH)-(IU-IH)"=c(0,0,0,0,1,0),"(CA-CH)-(IA-IH)"=c(0,0,0,0,0,1),"(CA-CU)-(IA-IU)"=c(0,0,0,0,1,-1))

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

#voom
v = voom(dge, design, plot=TRUE) 
plotMDS(v)

res=duplicateCorrelation(v,design,block=pas_info$sampleid)

vfit = lmFit(v, design, block = pas_info$sampleid, correlation = res$consensus)
vfit_alt = lmFit(v, design)
# inne i v er både en ny respons og vekter inn i GLS

vfit = eBayes(vfit)
vfit_alt = eBayes(vfit_alt)
topTable(vfit, coef=ncol(design))
#alternative, gives slightly different results
#vfit = treat(vfit, lfc=log2(1.2))
#topTreat(vfit, coef=ncol(design))

C = rbind("IU"=c(0,0,0,0,1,0), #(CU-CH)-(IU-IH)
          "IA"=c(0,0,0,0,0,1), #(CA-CH)-(IA-IH)
          "IA minus IU"=c(0,0,0,0,1,-1)) #(CA-CU)-(IA-IU)
vfit2 = contrasts.fit(vfit, t(C))
vfit2_alt = contrasts.fit(vfit_alt, t(C))
vfit2 = eBayes(vfit2, trend = TRUE)
vfit2_alt = eBayes(vfit2_alt, trend = TRUE)
#for specific coef use coef = 1, 2 or 3, as contrast has three coefs
table_vfit2 = topTable(vfit2, number = 10)
table_vfit2_alt = topTable(vfit2_alt, number = 10)

###########################################################################
#deseq (model)
###########################################################################

dds = DESeqDataSetFromMatrix(countData = round(count), colData = pas_info, design = ~ inflammation*tissue)
dds$tissue
#pre processing
length(dds[, 1])
keep = filterByExpr(counts(dds), design)
dds = dds[keep, ]
length(dds[, 1])
#58003 vs 43358
#vs 19668

dds <- DESeq(dds)
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
ggplot(as(contrast.UI, "data.frame"), aes(x = pvalue)) + 
  geom_histogram(binwidth = 0.01, fill = "Blue", boundary = 0)

DESeq2::plotMA(contrast.UI, ylim = c(-2, 2))

###########################################################################
#analysis of results
###########################################################################

#general venn diagram

venn = list("limma" = rownames(table_fit3_alt),
            "limma with correlation" = rownames(table_fit3), 
            "voom with correlation" = rownames(table_vfit2),
            "voom" = rownames(table_vfit2_alt))
ggvenn(venn, show_percentage = FALSE)

venn.p = list("limma" = rownames(topTable(fit3_alt, number = 10000, p.value = 0.05)),
            "limma with correlation" = rownames(topTable(fit3, number = 10000, p.value = 0.05)), 
            "voom with correlation" = rownames(topTable(vfit2, number = 10000, p.value = 0.05)),
            "voom" = rownames(topTable(vfit2_alt, number = 10000, p.value = 0.05)))
ggvenn(venn.p, show_percentage = FALSE)

#double check that adj.P.val is used
tail(topTable(vfit2_alt, number = 10000, p.value = 0.05))


#check coef with DESeq2, only considering with correlation
number = 10000
venn.1 = list("limma" = rownames(topTable(fit3, number = number, coef = 1, p.value = 0.05)), 
            "voom" = rownames(topTable(vfit2, number = number, coef = 1, p.value = 0.05)),
            "DESeq2" = rownames(contrast.UI[1:number, ]))
            # "DESeq2" = rownames(contrast.UI[contrast.UI$padj < 0.05, ]))
ggvenn(venn.1, show_percentage = FALSE)

venn.2 = list("limma" = rownames(topTable(fit3, number = number, coef = 2, p.value = 0.05)), 
              "voom" = rownames(topTable(vfit2, number = number, coef = 2, p.value = 0.05)),
              "DESeq2" = rownames(contrast.AI[1:number, ]))
              # "DESeq2" = rownames(contrast.AI[contrast.AI$padj < 0.05, ]))
ggvenn(venn.2, show_percentage = FALSE)

venn.3 = list("limma" = rownames(topTable(fit3, number = number, coef = 2, p.value = 0.05)), 
              "voom" = rownames(topTable(vfit2, number = number, coef = 2, p.value = 0.05)),
              "DESeq2" = rownames(contrast.complex[1:number, ]))
            # "DESeq2" = rownames(contrast.complex[contrast.complex$padj < 0.05, ]))
ggvenn(venn.3, show_percentage = FALSE)

#personal comments, code to quality check changes

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
# table(pas_info$tissue,pas_info$inflammation,pas_info$disease)
# 
# #not perfect way to do this
# n_occur = data.frame(table(id = pas_info$sampleid))
# #this gives us levels
# #ids with more than one sample
# sum(n_occur$Freq > 1)
# sum(n_occur$Freq == 2)
# sum(n_occur$Freq == 3)
# #ids with one sample
# sum(n_occur$Freq == 1)


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
# #lmFit(logCPM, design, block = pas_info$sampleid, correlation = corfit$consensus)
# # lmFit
# function (object, design = NULL, ndups = NULL, spacing = NULL, 
#           block = NULL, correlation, weights = NULL, method = "ls", 
#           ...) 
# {
#   y <- getEAWP(object)
#   design <- as.matrix(design)
#   
#   ne <- nonEstimable(design)
#   ndups <- 1
#   spacing <- 1
#   
#   method <- match.arg(method, c("ls", "robust"))
#   
#   fit <- gls.series(y$exprs, design = design, ndups = ndups, 
#                       spacing = spacing, block = block, correlation = correlation, 
#                       weights = weights, ...)
#   
#     n <- rowSums(is.na(fit$coefficients))
#     n <- sum(n > 0 & n < NCOL(fit$coefficients))
# 
#   fit$genes <- y$probes
#   fit$Amean <- y$Amean
#   fit$method <- method
#   fit$design <- design
#   new("MArrayLM", fit)
# }
# 
# 
# #gls.series(getEAWP(logCPM)$exprs, design = design, ndups = 1, 
# #spacing = 1, block = pas_info$sampleid, correlation = corfit$consensus, 
# #weights = NULL, ...)
# #gls.series
# function (M, design = NULL, ndups = 2, spacing = 1, block = NULL, 
#           correlation = NULL, weights = NULL, ...) 
# {
#   M <- as.matrix(M)
#   ngenes <- nrow(M)
#   narrays <- ncol(M)
#   design <- as.matrix(design)
#   nbeta <- ncol(design)
#   coef.names <- colnames(design)
# 
#   ndups <- spacing <- 1
#   block <- as.vector(block)
# 
#   ub <- unique(block)
#   nblocks <- length(ub)
#   Z <- matrix(block, narrays, nblocks) == matrix(ub, narrays, 
#                                                  nblocks, byrow = TRUE)
#   cormatrix <- Z %*% (correlation * t(Z))
#     
#   diag(cormatrix) <- 1
#   stdev.unscaled <- matrix(NA, ngenes, nbeta, dimnames = list(rownames(M), 
#                                                               coef.names))
#   NoProbeWts <- all(is.finite(M)) && (is.null(weights) || !is.null(attr(weights, 
#                                                                         "arrayweights")))
#   V <- cormatrix
# 
#   cholV <- chol(V)
#   y <- backsolve(cholV, t(M), transpose = TRUE)
#   dimnames(y) <- rev(dimnames(M))
#   X <- backsolve(cholV, design, transpose = TRUE)
#   dimnames(X) <- dimnames(design)
#   fit <- lm.fit(X, y)
#   if (fit$df.residual > 0) {
#     if (is.matrix(fit$effects)) 
#       fit$sigma <- sqrt(colMeans(fit$effects[-(1:fit$rank), 
#                                              , drop = FALSE]^2))
#     else fit$sigma <- sqrt(mean(fit$effects[-(1:fit$rank)]^2))
#   }
#   else fit$sigma <- rep_len(NA_real_, ngenes)
#   fit$fitted.values <- fit$residuals <- fit$effects <- NULL
#   fit$coefficients <- t(fit$coefficients)
#   fit$cov.coefficients <- chol2inv(fit$qr$qr, size = fit$qr$rank)
#   est <- fit$qr$pivot[1:fit$qr$rank]
#   dimnames(fit$cov.coefficients) <- list(coef.names[est], 
#                                          coef.names[est])
#   stdev.unscaled[, est] <- matrix(sqrt(diag(fit$cov.coefficients)), 
#                                   ngenes, fit$qr$rank, byrow = TRUE)
#   fit$stdev.unscaled <- stdev.unscaled
#   fit$df.residual <- rep_len(fit$df.residual, length.out = ngenes)
#   dimnames(fit$stdev.unscaled) <- dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
#   fit$pivot <- fit$qr$pivot
#   fit$ndups <- ndups
#   fit$spacing <- spacing
#   fit$block <- block
#   fit$correlation <- correlation
#   return(fit)
# }
# 
# 
# 
# #corfit = duplicateCorrelation(logCPM, design, block = pas_info$sampleid)
# #duplicate correlation
# function (object, design = NULL, ndups = 2L, spacing = 1L, block = NULL, 
#           trim = 0.15, weights = NULL) 
# {
#   y <- getEAWP(object)
#   M <- y$exprs
#   ngenes <- nrow(M)
#   narrays <- ncol(M)
#   design <- as.matrix(design)
#   nbeta <- ncol(design)
#   QR <- qr(design)
#   
#   MaxBlockSize <- max(table(block))
#   design.block <- model.matrix(~factor(block))
#   design.block <- design.block[, -1, drop = FALSE]
#   QtBlock <- qr.qty(QR, design.block)
#   
#   ndups <- 1L
#   nspacing <- 1L
#   Array <- block
#   
#   rho <- rep_len(NA_real_, ngenes)
#   nafun <- function(e) NA
#   for (i in 1:ngenes) {
#     y <- drop(M[i, ])
#     o <- is.finite(y)
#     A <- factor(Array[o])
#     nobs <- sum(o)
#     nblocks <- length(levels(A))
#     if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
#         (nobs - 1L)) {
#       y <- y[o]
#       X <- design[o, , drop = FALSE]
#       Z <- model.matrix(~0 + A)
#       if (!is.null(weights)) {
#         w <- drop(weights[i, ])[o]
#         s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
#                                                                X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
#                       error = nafun)
#       }
#       else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
#                                                                   X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
#                          error = nafun)
#       if (!is.na(s[1])) 
#         rho[i] <- s[2]/sum(s)
#     }
#   }
#   rhomax <- 0.99
#   
#   rhomin <- 1/(1 - MaxBlockSize) + 0.01
#   
#   m <- min(rho, 0, na.rm = TRUE)
#   if (m < rhomin) 
#     rho[rho < rhomin] <- rhomin
#   m <- max(rho, 0, na.rm = TRUE)
#   if (m > rhomax) 
#     rho[rho > rhomax] <- rhomax
#   
#   arho <- atanh(rho)
#   mrho <- tanh(mean(arho, trim = trim, na.rm = TRUE))
#   list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
# }