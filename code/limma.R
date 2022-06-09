library(stringr)
library(ggvenn)
library(edgeR)
library(DESeq2)
library(lme4)
library(hexbin)

#file location
fileloc = "Z:/Atle van Beelen Granlund/"


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
#this has been run twice for some reason
removal = c(removal, which(pas_info_tot$Sample_ID == "3002i"))

#Corrections from Atle
pas_info_tot$Sample_Biosource[pas_info_tot$Sample_ID == "3114i_t"] = "Colon tissue"
pas_info_tot$Sample_Biosource[pas_info_tot$Sample_ID == "3151i"] = "Ileum tissue"

#finds indexes for New RNA samples, these might be bad
#New_RNA = grep("New RNA", pas_info_tot$Comment)
#names_for_reference = pas_info_tot$Sample_ID[New_RNA]

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

#make as factors and order factors
pas_info$inflammation = factor(pas_info$inflammation, c("H", "U", "A"))
pas_info$tissue = factor(pas_info$tissue, c("Colon", "Ileum"))
table(pas_info$inflammation,pas_info$tissue)

###########################################################################
#limma/voom (model)
###########################################################################

#convert to relevant setup for limma

dge = DGEList(count)


#"remove counts that have zero or very low counts"
design = model.matrix(~ inflammation*tissue, data = pas_info)
keep = filterByExpr(dge, design)
dge = dge[keep, , keep.lib.sizes=FALSE]
dim(dge$counts)
# 27534, before 58003 
#use TMM normalization method
dge = calcNormFactors(dge)

#check largest library size against smallest
max(colSums(dge$counts))/min(colSums(dge$counts))

#limma vs voom

#limma
#"prior.count damp down the variances of logarithms of low counts"
logCPM = cpm(dge, log=TRUE)
dim(logCPM)

#plotMDS
plotMDS(logCPM, labels = substr(pas_info$tissue, 1, 1), 
        col = as.numeric(pas_info$diseaseinflammation))

C=rbind("(IU-IH)-(CU-CH)"=c(0,0,0,0,1,0), 
        "(IA-IH)-(CA-CH)"=c(0,0,0,0,0,1),
        "(IA-IU)-(CA-CU)"=c(0,0,0,0,-1,1))

#implement correlation
corfit = duplicateCorrelation(logCPM, design, block = pas_info$sampleid)

fit = lmFit(logCPM, design, block = pas_info$sampleid, correlation = corfit$consensus)


fit2 = contrasts.fit(fit,t(C))

fit3 = eBayes(fit2, trend=TRUE)

#consensus correlation
corfit$consensus.correlation
#plots and tables
limma::plotMD(fit3, coef = 1)
limma::plotMD(fit3, coef = 2)
limma::plotMD(fit3, coef = 3)
limma::plotSA(fit3, col = "green")
ggplot() + geom_histogram(aes(fit3$p.value[, 1]), binwidth = 0.01, fill = "blue")
ggplot() + geom_histogram(aes(fit3$p.value[, 2]), binwidth = 0.01, fill = "blue")
ggplot() + geom_histogram(aes(fit3$p.value[, 3]), binwidth = 0.01, fill = "blue")
topTable(fit3, coef = 1, number = 20)
topTable(fit3, coef = 2, number = 20)
topTable(fit3, coef = 3, number = 20)


#convert to relevant setup for voom

dge = DGEList(count)

#"remove counts that have zero or very low counts"
design = model.matrix(~ inflammation*tissue, data = pas_info)
keep = filterByExpr(dge, design)
dge = dge[keep, , keep.lib.sizes=FALSE]
dim(dge$counts)
# 27534, before 58003 
#use TMM normalization method
dge = calcNormFactors(dge)

#voom
v = voom(dge, design, plot=TRUE) 
plotMDS(v, labels = substr(pas_info$tissue, 1, 1), 
        col = as.numeric(pas_info$diseaseinflammation))
legend("top", legend = levels(pas_info$diseaseinflammation), col = 1:3, pch = 15)

res = duplicateCorrelation(v, design, block=pas_info$sampleid)

vfit = lmFit(v, design, block = pas_info$sampleid, correlation = res$consensus, weights = v$weights)
# inside of v is both a new response and weights, which are put into GLS

C=rbind("(IU-IH)-(CU-CH)"=c(0,0,0,0,1,0), 
        "(IA-IH)-(CA-CH)"=c(0,0,0,0,0,1),
        "(IA-IU)-(CA-CU)"=c(0,0,0,0,-1,1))
vfit2 = contrasts.fit(vfit, t(C))
vfit2 = eBayes(vfit2, trend = FALSE)

# consensus correlation
res$consensus.correlation
#plots and tables
limma::plotSA(vfit2)
limma::plotMD(vfit2, coef = 1)
limma::plotMD(vfit2, coef = 2)
limma::plotMD(vfit2, coef = 3)
ggplot() + geom_histogram(aes(vfit2$p.value[, 1]), binwidth = 0.01, fill = "blue")
ggplot() + geom_histogram(aes(vfit2$p.value[, 2]), binwidth = 0.01, fill = "blue")
ggplot() + geom_histogram(aes(vfit2$p.value[, 3]), binwidth = 0.01, fill = "blue")
topTable(vfit2, coef = 1, number = 20)
topTable(vfit2, coef = 2, number = 20)
topTable(vfit2, coef = 3, number = 20)


###########################################################################
#analysis of results
###########################################################################


#only considering with correlation
number = 1000000
venn.1 = list("limma" = rownames(topTable(fit3, number = number, coef = 1, p.value = 0.05)), 
              "voom" = rownames(topTable(vfit2, number = number, coef = 1, p.value = 0.05)))

ggvenn(venn.1, show_percentage = FALSE)

venn.2 = list("limma" = rownames(topTable(fit3, number = number, coef = 2, p.value = 0.05)), 
              "voom" = rownames(topTable(vfit2, number = number, coef = 2, p.value = 0.05)))

ggvenn(venn.2, show_percentage = FALSE)

venn.3 = list("limma" = rownames(topTable(fit3, number = number, coef = 3, p.value = 0.05)), 
              "voom" = rownames(topTable(vfit2, number = number, coef = 3, p.value = 0.05)))

ggvenn(venn.3, show_percentage = FALSE)

G = length(voom_cor$atanh.correlations)
voom_cor$atanh.correlations = sort(voom_cor$atanh.correlations)
ggplot() + 
        geom_histogram(aes(x = voom_cor$atanh.correlations, y = ..count../sum(..count..)), binwidth = 0.02) + 
                geom_vline(aes(xintercept = c(voom_cor$atanh.correlations[round(0.15*G)], 
                                      voom_cor$atanh.correlations[round(0.85*G)]))) +
        labs(x = "atanh correlation", y = "proportion")
