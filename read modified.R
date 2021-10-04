library(stringr)

#file location
fileloc = "Z:/Atle van Beelen Granlund/"

# info about the samples to connect to disease, tissue, inflammation, pasid

pas_info_tot <- read.delim(paste(fileloc,'sample_info.txt',sep=""),sep="\t", header = TRUE, row.names = 1)
count_tot = read.delim(paste(fileloc, "quasi.gene.quant", sep = ""), sep = "\t", check.names = FALSE)
#removes New RNA samples
New_RNA = grep("New RNA", pas_info_tot$Comment)
#add random comment from Atle
New_RNA = c(New_RNA, which(pas_info_tot$Sample_ID == "3115i"))
#names_for_referance = pas_info_tot$Sample_ID[New_RNA]
pas_info_tot = pas_info_tot[-New_RNA, ]
count_tot = count_tot[, pas_info_tot$Sample_ID]

#extract relevant columns
pas_info <- data.frame("sampleid" = pas_info_tot$Sample_ID, "diseaseinflammation" = pas_info_tot$Sample_Group, "tissue" = pas_info_tot$Sample_Biosource)

#using grep to remove final segments of string ending with irrelevant information
firstpart = str_split_fixed(pas_info$sampleid, "_", 2)[,1]
iFS = grep("*[iFS]", firstpart)
firstpart[iFS] = substr(firstpart[iFS], 1, nchar(firstpart[iFS]) - 1)
pas_info$sampleid = firstpart

#comments of how to find something, more than relevant code
unique(pas_info$diseaseinflammation[pas_info$tissue=="Ileum tissue"])
pas_info$sampleid[pas_info$tissue=="Ileum tissue" & pas_info$diseaseinflammation == "CDA"]
pas_info$sampleid[pas_info$tissue=="Ileum tissue" & pas_info$diseaseinflammation == "CDU"]

unique(pas_info$diseaseinflammation[pas_info$tissue=="Colon tissue"])
pas_info$sampleid[pas_info$tissue=="Colon tissue" & pas_info$diseaseinflammation == "CDAi"]
pas_info$sampleid[pas_info$tissue=="Colon tissue" & pas_info$diseaseinflammation == "CDUi"]

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

#finds tables for comparing how much data we have for each group
table(pas_info$disease)
table(pas_info$inflammation)
table(pas_info$tissue,pas_info$inflammation,pas_info$disease)

#not perfect way to do this
n_occur = data.frame(table(id = pas_info$sampleid))
#this gives us levels
#ids with more than one sample
sum(n_occur$Freq > 1)
sum(n_occur$Freq == 2)
sum(n_occur$Freq == 3)
#ids with one sample
sum(n_occur$Freq == 1)

