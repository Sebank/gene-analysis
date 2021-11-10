library(stringr)

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
pas_info <- data.frame("sampleid" = pas_info_tot$Sample_ID, "diseaseinflammation" = pas_info_tot$Sample_Group, "tissue" = pas_info_tot$Sample_Biosource)

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



###########################################################################
#limma/voom (model)
###########################################################################

#convert to relevant setup for limma
library(edgeR)
dge = DGEList(count)

#make as factors and order factors
pas_info$inflammation = factor(pas_info$inflammation, c("H", "U", "A"))
pas_info$tissue = factor(pas_info$tissue, c("Colon tissue", "Ileum tissue"))
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
#analysis of results
###########################################################################

#venn diagram
library(ggvenn)
venn = list("limma with correlation" = rownames(table_fit3), 
            "limma" = rownames(table_fit3_alt),
            "voom with correlation" = rownames(table_vfit2),
            "voom" = rownames(table_vfit2_alt))
ggvenn(venn, show_percentage = FALSE)

venn.p = list("limma with correlation" = rownames(topTable(fit3, number = 10000, p.value = 0.05)), 
            "limma" = rownames(topTable(fit3_alt, number = 10000, p.value = 0.05)),
            "voom with correlation" = rownames(topTable(vfit2, number = 10000, p.value = 0.05)),
            "voom" = rownames(rownames(topTable(vfit2_alt, number = 10000, p.value = 0.05))))
ggvenn(venn.p, show_percentage = FALSE)

#double check that adj.P.val is used
tail(topTable(vfit2_alt, number = 10000, p.value = 0.05))

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


cbind(pas_info$sampleid,colnames(count_tot))

#finne konsepter fra data

summary(dge)
dge
dim(dge)
# 27672 110
dim(dge$counts)
dim(dge$samples)
dge$samples


fit0=lm(logCPM[1,]~design-1)
dim(fit0$coefficients)
dim(C)
fit0.C=C%*%matrix(fit0$coefficients,ncol=1)
fit0.C
fit0.Ccov=C%*%vcov(fit0)%*%t(C)
fit0.Ccov
fcontrast3=fit0.C[3,1]^2/fit0.Ccov[3,3]
fcontrast2=fit0.C[2,1]^2/fit0.Ccov[2,2]
pf(fcontrast2,1,fit0$df.residual,lower.tail=FALSE)
pf(fcontrast3,1,fit0$df.residual,lower.tail=FALSE)

dim(fit$coefficients)
fit$coefficients[1,]
summary(fit0)

str(fit2$contrasts)
names(fit2)
names(fit)
dim(fit2$coefficients)
fit$coefficients[1,]

# leter etter pverdier og testobs
summary(fit0)

fit3$t[1,2]*sqrt(fit3$s2.post[1])*fit3$stdev.unscaled[1,2]
fit0$coefficients[6]
fit3$s2.prior[1]
fit3$s2.post[1]

plot(fit3$Amean,fit3$s2.prior)
plot(fit3$Amean,sqrt(fit3$s2.post),pch=".")
lines(lowess(fit3$Amean,sqrt(fit3$s2.post)),col=2)

names(fit3)
dim(fit3$p.value)
fit3$p.value[1,]
fit3$stdev.unscaled[1,]
fit3$sigma[1]
fit3$Amean[1]

str(v)
hj=getEAWP(v)