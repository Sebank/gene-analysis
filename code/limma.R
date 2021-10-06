library(stringr)

#file location
fileloc = "Z:/Atle van Beelen Granlund/"

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

#double check if changes made above
#unique(pas_info$diseaseinflammation[pas_info$tissue=="Ileum tissue"])
#pas_info$sampleid[pas_info$tissue=="Ileum tissue" & pas_info$diseaseinflammation == "CDA"]
#pas_info$sampleid[pas_info$tissue=="Ileum tissue" & pas_info$diseaseinflammation == "CDU"]

#unique(pas_info$diseaseinflammation[pas_info$tissue=="Colon tissue"])
#pas_info$sampleid[pas_info$tissue=="Colon tissue" & pas_info$diseaseinflammation == "CDAi"]
#pas_info$sampleid[pas_info$tissue=="Colon tissue" & pas_info$diseaseinflammation == "CDUi"]

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
pas_info = pas_info[-2]

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

#Only examine CD and H, not UC
count = count_tot[, !pas_info$disease == "UC"]
pas_info = pas_info[!pas_info$disease == "UC", ]

#convert to relevant setup for limma
library(edgeR)
dge = DGEList(count)

#make as factors and order factors
pas_info$inflammation = factor(pas_info$inflammation, c("H", "U", "A"))
pas_info$tissue = factor(pas_info$tissue, c("Colon tissue", "Ileum tissue"))

#"remove counts that have zero or very low counts"
design = model.matrix(~ inflammation*tissue, data = pas_info)
keep = filterByExpr(dge, design)
dge = dge[keep, , keep.lib.sizes=FALSE]

#use TMM normalization method
dge = calcNormFactors(dge)

#limma vs voom
#limma
#"prior.count damp down the variances of logarithms of low counts"
logCPM = cpm(dge, log=TRUE, prior.count=3)

fit = lmFit(logCPM, design)
fit = eBayes(fit, trend=TRUE)
topTable(fit, coef=ncol(design))

#voom
v = voom(dge, design, plot=TRUE) 
plotMDS(v)
#alternatively don't need normalization for voom 
#v = voom(cont, design, plot=TRUE)

vfit = lmFit(v, design)
vfit = eBayes(fit)
topTable(vfit, coef=ncol(design))
#alternative, gives slightly different results
#vfit = treat(vfit, lfc=log2(1.2))
#topTreat(vfit, coef=ncol(design))
