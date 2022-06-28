#!/usr/bin/env Rscript
# Author: OLP/01/2022
# Usage: Getting metadata for ngs file
# Parameters: Run Name 
#             i.e. Rscript / . . . ./ngs.R UT-A01290-220528

library("easypackages");libraries("tidyverse","lubridate","grid","rbin","data.table")

args = commandArgs(trailingOnly=T)

date<-ymd(substr(args[1], 11, 16));runPath <-paste('/Volumes/NGS/Analysis/covidseq',args[1],sep="/")
df1 = read.csv(paste(runPath,"covidseq_output/SampleSheet_Intermediate.csv",sep ="/"))
df1<-tail(df1,-13);df1<-df1[,c(1,3)];
colnames(df1) <- c("Sample_Accession","Sample_Type")
M<-filter(df1, Sample_Type == "PatientSample")
M<-M["Sample_Accession"]
cat(" Detected ", length(M$Sample_Accession), " patient's samples");

# Getting Pangolin lineage results
pango<-read.csv(paste(runPath,"covidseq_output/lineage/combined_lineage_report.csv",sep ="/"))[,c("taxon","lineage")]
var <-paste('',args[1],sep='-'); head(var)
pango<-pango %>% mutate_at("taxon", str_replace, var, "")
M<-merge(M,pango, by.x="Sample_Accession", by.y= "taxon", all.x=T)
head(M, n = 200)
# Searching EpiTrax;
cat(" Searching for patient's info on EpiTrax")
EpiTrax_Export<-rownames(fileSnapshot("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/", 
                                      pattern=glob2rx("export_1551*.csv"))$info[which.max(fileSnapshot("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/",
                                                                                                       pattern=glob2rx("export_1551*.csv"))$info$mtime),]);
M_epitrax <- fread(paste('/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/',EpiTrax_Export,sep=""), 
                   select = c("lab_accession_no","person_last_name","person_first_name", "patient_birth_date","lab_collection_date"))   

M_epi<-merge(M,M_epitrax,by.x="Sample_Accession", by.y= "lab_accession_no", all.x=T)
colnames(M_epi) <- c("Sample_Accession","Test result","lastName","firstName","DOB","collectionDate") 
M_epi$DOB<-as.Date(M_epi$DOB);M_epi$collectionDate<-as.Date(M_epi$collectionDate)
M_epi<- M_epi %>% mutate(across(, ~as.character(.)));M_epi<-unique(M_epi)
Mx<-M_epi$Sample_Accession[which(!is.na(M_epi$collectionDate))]
Matrix1<-merge(Mx,M_epi,by.x ="x", by.y ="Sample_Accession", all.x = T)
head(Matrix1, n = 200)

# Searching LIMS
cat("Searching for patient's info on LIMS")
Mx0<-M_epi$Sample_Accession[which(is.na(M_epi$collectionDate))]; 
Mx0<-merge(Mx0,M, by.x = "x", by.y ="Sample_Accession", all.x  = T)

NGS_covid<-rownames(fileSnapshot("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/",
                                 pattern=glob2rx("NGS*.csv"))$info[which.max(fileSnapshot("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/",                                                                                           pattern=glob2rx("NGS*.csv"))$info$mtime),])
NGS_covid <- fread(paste('/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/',NGS_covid,sep=""), 
                   select = c("submitterId","sampleNumber","lastName","firstName","DOB","collectionDate")) 
NGS_covid<- NGS_covid %>% mutate(across(, ~as.character(.)))

head(NGS_covid, n = 200)
NGS_covid1<-select(NGS_covid,c("sampleNumber","lastName","firstName","DOB","collectionDate"))
NGS_covid2<-select(NGS_covid,c("submitterId","lastName","firstName","DOB","collectionDate"))

All_nocovid<-rownames(fileSnapshot("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/",
                                 pattern=glob2rx("All*.csv"))$info[which.max(fileSnapshot("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/", 
                                                                                          pattern=glob2rx("All*.csv"))$info$mtime),])
All_nocovid <- fread(paste('/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/',All_nocovid,sep=""), 
                   select = c("sampleNumber","lastName","firstName","DOB","collectionDate")) 
All_nocovid<- All_nocovid %>% mutate(across(, ~as.character(.)))
head(All_nocovid, n = 200)
colnames(NGS_covid2)[1]<-("sampleNumber")
LIMS <-bind_rows(NGS_covid1, NGS_covid2, All_nocovid);

Matrix2<-merge(Mx0, LIMS, by.x = "x", by.y ="sampleNumber"); 
colnames(Matrix2)[2]<-("Test result")
ngs_file<-bind_rows(Matrix1, Matrix2)
colnames(ngs_file)[1]<-("sampleNumber")  
ngs_file2<-cbind(ngs_file,'Name of facility performing the test'='UPHL', 'ordering facility' ='UHPL',
               'Name or description of the test ordered'="SARS-CoV-2 Whole Genome Sequencing",
               'lab test date'=ymd(date),'lab test status' ="preliminary")

ngs_file3<-select(ngs_file2,"lastName",	"firstName",	"DOB",	"Name of facility performing the test","Name or description of the test ordered",
                  "collectionDate","sampleNumber","lab test date","lab test status",	"Test result",	"ordering facility")

# saving results
head(ngs_file3, n = 200)
arg<-args[1]; ngs_file_name<-paste('ngs',arg,sep="_");
d=paste(runPath,'R-sumcovidseq',sep="/");
dir.create(file.path(d,"" ), recursive = T); 
LL<-paste(substring(Sys.time(), 1,10),substring(Sys.time(), 12,13),substring(Sys.time(), 15,16),substring(Sys.time(), 18,19),sep = "-")
dx<- paste(ngs_file_name,LL,".csv",sep="");
dx<-paste(d,dx, sep ='/');write.csv(ngs_file3,dx,row.names=F)
cat("ngs file is completed and can be found at: ",d)


