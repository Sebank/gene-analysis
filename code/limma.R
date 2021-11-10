library(stringr)

#file location
fileloc = "Z:/Atle van Beelen Granlund/"

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

#remove bad samples
pas_info_tot = pas_info_tot[-removal, ]
count_tot = count_tot[, pas_info_tot$Sample_ID]

#extract relevant columns
pas_info <- data.frame("sampleid" = pas_info_tot$Sample_ID, 
                       "diseaseinflammation" = pas_info_tot$Sample_Group, 
                       "tissue" = pas_info_tot$Sample_Biosource)

#using grep to remove final segments of string ending with irrelevant information
firstpart = str_split_fixed(pas_info$sampleid, "_", 2)[,1]
iFS = grep("*[iFS]", firstpart)
firstpart[iFS] = substr(firstpart[iFS], 1, nchar(firstpart[iFS]) - 1)
pas_info$sampleid = firstpart

#remove i from disease information, as it should be part of tissue as well
i = grepl("*i", pas_info$diseaseinflammation)
pas_info$diseaseinflammation[i] = 
  substr(pas_info$diseaseinflammation[i], 1, nchar(pas_info$diseaseinflammation[i]) - 1)

#make the variable disease inflammation two variables
#replaced F with H for "healthy"
AU = grep("..[AU]", pas_info$diseaseinflammation)
pas_info$disease = rep("H", length(pas_info$diseaseinflammation))
pas_info$disease[AU] = 
  substr(pas_info$diseaseinflammation[AU], 1, nchar(pas_info$diseaseinflammation[AU]) - 1)

#set inflamation variable
A = grep("A", pas_info$diseaseinflammation)
H = grep("F", pas_info$diseaseinflammation)
pas_info$inflammation = rep("U", length(pas_info$diseaseinflammation))
pas_info$inflammation[A] = "A"
pas_info$inflammation[H] = "H"
pas_info$inflammation = factor(pas_info$inflammation, c("H", "U", "A"))

#remove diseaseinflamation (never run this twice!!!)
pas_info = pas_info[-2]

#Only examine CD and H, not UC
count = count_tot[, !pas_info$disease == "UC"]
pas_info = pas_info[!pas_info$disease == "UC", ]

#examines how many patients in each group
table(pas_info$tissue, pas_info$disease)

#convert to relevant setup for limma
library(edgeR)
dge = DGEList(count)

#make as factors and order factors
pas_info$inflammation = factor(pas_info$inflammation, c("H", "U", "A"))
pas_info$tissue = factor(pas_info$tissue, c("Colon tissue", "Ileum tissue"))

#"remove counts that have zero or very low counts"
design = model.matrix(~ inflammation*tissue, data = pas_info)#[, -1]
keep = filterByExpr(dge, design)
dge = dge[keep, , keep.lib.sizes=FALSE]

#use TMM normalization method
dge = calcNormFactors(dge)

#limma vs voom
#limma
#"prior.count damp down the variances of logarithms of low counts"
logCPM = cpm(dge, log=TRUE, prior.count=3)

fit = lmFit(logCPM, design)
#trend - should an intensity-trend be allowed for prior value
fit = eBayes(fit, trend=TRUE)
#topTable(fit)

#using the contrast matrix
C = rbind("IU"=c(0,0,0,0,1,0), #(CU-CH)-(IU-IH)
          "IA"=c(0,0,0,0,0,1), #(CA-CH)-(IA-IH)
          "IA minus IU"=c(0,0,0,0,1,-1)) #(CA-CU)-(IA-IU)
fit2 = contrasts.fit(fit,t(C))
fit2 = eBayes(fit2, trend = TRUE)
#for specific coef use coef = 1, 2 or 3, as contrast has three coefs
table_fit2 = topTable(fit2)

#voom
v = voom(dge, design, plot=TRUE) 
plotMDS(v)
#alternatively don't need normalization for voom 
#v = voom(count, design, plot=TRUE, normalize = "quantile")

vfit = lmFit(v, design)
vfit = eBayes(vfit)
#topTable(vfit)

#alternative, gives slightly different results
#vfit = treat(vfit, lfc=log2(1.2))
#topTreat(vfit)

vfit2 = contrasts.fit(vfit, t(C))
vfit2 = eBayes(vfit2, trend = TRUE)
#for specific coef use coef = 1, 2 or 3, as contrast has three coefs
table_vfit2 = topTable(vfit2)

#venn diagram
library(ggvenn)
venn = list("limma" = rownames(table_fit2), 
            "voom" = rownames(table_vfit2))
ggvenn(venn, show_percentage = FALSE)


#split count from filter by expr
selected_count = count[rownames(dge), ]

R = vector("numeric", dim(selected_count)[2])

log_selected_count = count[rownames(dge), ]

for(i in 1:dim(selected_count)[2]){
  R[i] = sum(selected_count[, i])
  log_selected_count[, i] = log2((selected_count[, i] + 0.5)/(R[i] + 1)*10^6)
}

beta = matrix(0, dim(selected_count)[1], dim(design[, -1])[2])

for(i in 1:dim(selected_count)[1]){
  beta[i, ] = lm(unlist(log_selected_count[i, ]) ~ design[, -1])$coef[-1]
}
rownames(beta) = rownames(selected_count)
colnames(beta) = colnames(design)[-1]

beta_C = beta %*% t(C[, -1])

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