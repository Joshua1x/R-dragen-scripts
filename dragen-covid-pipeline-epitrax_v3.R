#!/usr/bin/env Rscript
# Author: OLP/01/2022
# Usage: Getting collection dates and other metadata from EpiTrax for samples with LabIds and Customer Ids
# Parameters: Only one: Run Name i.e. Rscript / . . . ./dragen-covid-pipeline-epitrax-v2.R UT-A01290-220113
# Note: include - in the filtering
library("dplyr", warn.conflicts = FALSE)
library("grid")
library("rbin")#
library("tidyverse")
library("lubridate", warn.conflicts = FALSE)
library("tibble")
library("plyr")
#install.packages("sos")
#library("sos")

# Step one: date from argument
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
args = commandArgs(trailingOnly=TRUE)
args[1]<-c('UT-A01290-220301'); # only to works on script preparation
s<- paste('/Volumes/NGS/Analysis/covidseq/',args[1],sep="/"); head(s)
ss<-gsub("\\.*.*/","", s);date<-ymd(substr(ss, 11, 16)); head(ss); head(date)

# Step 2: Reading covidseq_summary
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
loc<-paste(paste('/Volumes/NGS/Analysis/covidseq',args[1],sep="/"),"sumcovidseq",sep ="/"); head(loc)
covidseq_summary<-read.csv(paste(loc,"covidseq_summary.csv", sep="/"));
covidseq_summary<-covidseq_summary[,c(1,9)]; colnames(covidseq_summary)[2] <- 'Test result';#  sample and pangolin_lineage


# Step 3: Selecting ng only PatientSamples from the SampleSheet
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# filtering Empy weels. 
# The samplesheet was modified on 2022-02-26 where Empty well is denoted as Empty but if this
# script is used for data before that the Empy wells entries are filtered here
sample<-covidseq_summary$sample
sample<-gsub("\\Emp.*","", sample);
covidseq_summary$sample<-sample

# Now filtering only PatientSamples
loc2 <-paste("/Volumes/NGS/Analysis/covidseq",args[1],sep="/"); head(loc2)
SampleSheet<-read.csv(paste(loc2,"covidseq_output/SampleSheet_Intermediate.csv", sep="/"));
ss1<-tail(SampleSheet,-13);#ss : SampleSheet
ss<-ss1[,c(1,3)]; # sampleNumber and PatientSample
colnames(ss) <- c("Sample_Accession","Sample_Type");
ss<-filter(ss, Sample_Type == "PatientSample"); 
covidseq_summary_filtered<-merge(ss,covidseq_summary, by.x="Sample_Accession", by.y= "sample")
covidseq_summary_filtered<-covidseq_summary_filtered[,c(1,3)]


write.csv(covidseq_summary_filtered,'/Users/olinto-linaresperdomo/Desktop/A020301/R-Results/covidseq_summary_filtered100.csv')

#Step 4: Detecting last files from NAS: NGS_Covid_...*csv & All_ncovid_NoNGS...*csv files
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
tmpshot <- fileSnapshot("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/", pattern=glob2rx("All*.csv"))#
nNGS<-rownames(tmpshot$info[which.max(tmpshot$info$mtime),]);
tmpshot <- fileSnapshot("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/", pattern=glob2rx("NGS*.csv"))#
NGS<-rownames(tmpshot$info[which.max(tmpshot$info$mtime),]);
cat("2.- These files have been detected (IDGenomics_NAS) as the last files created: ",nNGS, " and ",NGS,"\n")


#Step 5:Preparing all data from NAS into a file: LIMS (NGS+nNGS)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
covidseq_summaryIdLab<-covidseq_summary_filtered[c(grep("^[[:digit:]]", covidseq_summary_filtered$Sample_Accession)), ]; # only samples with LabId
NGS<-read.csv(paste("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/",NGS,sep="")); 
NGSS<-NGS[,c(3,4,5,7,9,12)];head(NGSS);
nNGS<-read.csv(paste("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/",nNGS,sep=""));
nNGSS<-nNGS[,c(3,4,5,7,9,11)];head(nNGSS);
write.csv(NGSS,'/Users/olinto-linaresperdomo/Desktop/A020301/R-Results/NGS.csv')
write.csv(nNGS,'/Users/olinto-linaresperdomo/Desktop/A020301/R-Results/nNGS.csv')
LIMS<-bind_rows(NGSS, nNGSS);head(LIMS);
write.csv(LIMS,'/Users/olinto-linaresperdomo/Desktop/A020301/R-Results/LIMS.csv')

#Step 6:matrix1: Collecting demographic data for those samples with labID (left join tables)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
matrix1<-merge(covidseq_summaryIdLab,LIMS, by.x="Sample_Accession", by.y= "sampleNumber", all.x=TRUE);head(matrix1); # left join
matrix1<-matrix1[,c(3,4,5,7,1,2)]; head(matrix1)
colnames(matrix1)[5] <- 'sampleNumber';
matrix1<-distinct(matrix1)
write.csv(matrix1,'/Users/olinto-linaresperdomo/Desktop/A020301/R-Results/matrix1.csv')

# Step 7: matrix2: Demographic data for samples with IHC (Exxxxxxx/Xxxxxx) left join tables)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
LIMS2<-(LIMS %>% filter(grepl('/',submitterId))); head(LIMS2)
IHC1<-gsub("\\.*.*/","", (LIMS2 %>% filter(grepl('/',submitterId)))$submitterId );
LIMS2$IHC1<-IHC1
head(LIMS2)
LIMS2$submitterId<-NULL
write.csv(LIMS2,'/Users/olinto-linaresperdomo/Desktop/A020301/R-Results/LIMS2.csv')
covidseq_summaryIHC1<-covidseq_summary_filtered[-c(grep("^[[:digit:]]", covidseq_summary_filtered$Sample_Accession)), ]; # samples with customer Ids
head(covidseq_summaryIHC1)
matrix2<-merge(covidseq_summaryIHC1,LIMS2, by.x="Sample_Accession", by.y= "IHC1");# inner join
head(matrix2)
matrix2<-distinct(matrix2)
matrix22<-matrix2; 
head(matrix22); # to same the IHC variable
head(covidseq_summaryIHC1);head(LIMS2)
matrix2<-matrix2[,c(3,4,5,7,6,2)]; head(matrix2)
head(matrix2)
lapply(matrix2,class)
matrix2$sampleNumber<-as.character(matrix2$sampleNumber)
matrix2<-distinct(matrix2)
write.csv(matrix2,'/Users/olinto-linaresperdomo/Desktop/A020301/R-Results/matrix2.csv')
write.csv(matrix22,'/Users/olinto-linaresperdomo/Desktop/A020301/R-Results/matrix22.csv')

# Step 8: matrix3: Demographic data for samples with any customer but no with xxxxxxx/xxxxxxxx
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
LIMS3<-dplyr::filter(LIMS, !grepl("/",submitterId)); head(LIMS3);
write.csv(LIMS3,'/Users/olinto-linaresperdomo/Desktop/A020301/R-Results/LIMS3.csv')
head(matrix22)
L1<-matrix22$Sample_Accession;
L2<-covidseq_summaryIHC1$Sample_Accession; head(L2); # 
L3<-L2[-pmatch(L1,L2)]; 
covidseq_summaryOthers<-merge(covidseq_summaryIHC1,L3,by.x="Sample_Accession", by.y= "y",all.y = TRUE); # right join)
matrix3<-merge(covidseq_summaryOthers,LIMS3, by.x="Sample_Accession", by.y= "submitterId",all.x=TRUE); # 
matrix3$sampleNumber<-as.character(matrix3$sampleNumber)


########coalesce
x1<-matrix3$Sample_Accession;
x2<-matrix3$sampleNumber;
cc<-coalesce(x1,x2);
matrix3$sampleNumber<-cc
matrix3<-matrix3[,c(3,4,5,6,7,2)];
write.csv(matrix3,'/Users/olinto-linaresperdomo/Desktop/A020301/R-Results/matrix3.csv')
# step 9:  preparing ngs file
ngs_file<-bind_rows(matrix1,matrix2, matrix3);
write.csv(ngs_file,'/Users/olinto-linaresperdomo/Desktop/A020301/R-Results/ns_first.csv')
ngs_file$'Test result'[is.na(ngs_file$'Test result')]= "No able to be sequenced";
ngs_file2<-cbind(ngs_file, 'Name of facility performing the test'='UPHL', 'ordering facility' ='UHPL',
                 'Name or description of the test ordered'="SARS-CoV-2 Whole Genome Sequencing",
                 'lab test date'=ymd(date),'lab test status' ="preliminary",'ordering facility' ="UPHL")
ngs_file2$'lab test date'<-as.character(ngs_file2$'lab test date')

ngs_file2[is.na(ngs_file2)]= ""
ngs_file<-ngs_file2[,c(1,2,3,7,9,4,5,10,11,6,8)];head(ngs_file)
lapply(ngs_file, class)

write.csv(ngs_file,'/Users/olinto-linaresperdomo/Desktop/A020301/R-Results/ngs_final_R.csv',row.names = FALSE)
ngs_file_name<-paste(paste('ngs',args[1],sep="_"),'csv',sep=".");
print(ngs_file_name)
print(ngs_file)
head(ngs_file)
lapply(ngs_file, class)
# data<-read.csv('/Users/olinto-linaresperdomo/Desktop/A020301/R-Results/these two files are differents/ngs_UT-A01290-220301.csv')

print(data)
lapply(data, class)
lapply(ngs_file, class)

print(getwd())
print(setwd('/Users/olinto-linaresperdomo/Desktop/ngsP'))
setwd('/Users/olinto-linaresperdomo/Desktop/ngsP')

#dat<-read.csv('/Users/olinto-linaresperdomo/Desktop/A220224/ngs_UT-A01290-220224.csv')
#print(dat)
# write results in covidseq_summary folder
# write.csv(ngs_file,paste(loc,ngs_file_name,sep="/"),row.names=FALSE)
cat("The file",ngs_file_name, "should be completed and saved at" ,loc ,"\n")
# sending data to EpiTrax-personmatching
# write.csv(dat,'/Volumes/DDCP/Division Shared Files/DCPIP EDX/UPHL/COVID/WGS Results/PersonMatching/ngs_UT-A01290-220224x.csv',row.names=FALSE)


