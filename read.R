library(stringr)

#fileloc="/Volumes/nice/p180/Atle van Beelen Granlund/"
fileloc = "Z:/Atle van Beelen Granlund/"

file <- read.csv(paste(fileloc,"colonCounts.csv",sep=""))
ccounts <- file
ccounts <- file[,-1]
rownames(ccounts) <- file[,1]
colnames(ccounts)

file <- read.csv(paste(fileloc,"ileumCounts.csv",sep=""))
icounts <- file
icounts <- file[,-1]
rownames(icounts) <- file[,1]
colnames(icounts)

# info about the samples to connect to disease, tissue, inflammation, pasid

#count.fields(paste(fileloc,'sample_info.txt',sep=""))
pas_info_tot <- read.table(paste(fileloc,'sample_info.txt',sep=""),sep="\t",skip=1)
dim(pas_info_tot)
colnames(pas_info_tot)=scan(paste(fileloc,'sample_info.txt',sep=""),what="s",n=24)
pas_info <- data.frame("sampleid" = pas_info_tot$Sample_ID, "diseaseinflammation" = pas_info_tot$Sample_Group, "tissue" = pas_info_tot$Sample_Biosource)
#extract relevant columns

# neste: fikse på diseaseinflammation og på pasid

# først pasid

firstpart=str_split_fixed(pas_info$sampleid, "_", 2)[,1]
# for alle unntatt disse som starter på TT skal vi fjerne siste char - evt kan vi beholde numeriske verdier?

# denne virker for alle unntatt de TTene som får ett siffer for lite
firstpartU=substr(str_split_fixed(pas_info$sampleid, "_", 2)[,1], 1, nchar(str_split_fixed(pas_info$sampleid, "_", 2)[,1])-1)


# innføre T/F for dette
nn=nchar(str_split_fixed(pas_info$sampleid, "_", 2)[,1])
# sjekker siste char
substring(firstpart[1],nn[1])

isnum=rep(0,length(nn))
for (i in 1:length(nn))
{
  if(!is.na(as.numeric(substring(firstpart[i],nn[i])))) isnum[i]=1
}

pas_info$pasid=firstpart
pas_info$pasid[isnum==0]=firstpartU[isnum==0]
pas_info$pasid

# fikse disease and inflammation
# først - ta bort i for de som er ileum - det bør henge sammen med tissue 
pas_info$diseaseinflammation[pas_info$tissue=="Ileum tissue"]
pas_info$diseaseinflammation[pas_info$tissue=="Colon tissue"]
# stemmer ikke helt, en med i på Colon og en uten i på Ileum

nn=nchar(pas_info$diseaseinflammation)
nn

hj=pas_info$diseaseinflammation
for (i in 1:length(nn))
{
  if (substring(pas_info$diseaseinflammation[i],nn[i])=="i") hj[i]=substring(pas_info$diseaseinflammation[i],1,nn[i]-1)
}
table(hj)

disease=rep(NA, length(nn))
inflammation=rep(NA,length(nn))
for (i in 1:length(nn))
{
  if (hj[i]=="CDA") 
  {
    disease[i]="CD"; inflammation[i]="A"
  }
  if (hj[i]=="CDU") 
  {
    disease[i]="CD"; inflammation[i]="U"
  }
  if (hj[i]=="UCA") 
  {
    disease[i]="UC"; inflammation[i]="A"
  }
  if (hj[i]=="UCU") 
  {
    disease[i]="UC"; inflammation[i]="U"
  }
  if (hj[i]=="F") 
  {
    disease[i]="H"; inflammation[i]="H"
  }
}
  
pas_info$disease=disease
pas_info$inflammation=inflammation

table(pas_info$disease)
table(pas_info$inflammation)
table(pas_info$tissue,pas_info$inflammation,pas_info$disease)

