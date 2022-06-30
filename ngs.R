#!/usr/bin/env Rscript
# Author: OLP/01/2022
# Usage: Getting metadata for ngs file
# Parameters: Run Name 
#             i.e. Rscript / . . . ./ngs.R UT-A01290-220528
library("easypackages");libraries("tidyverse","lubridate","grid","rbin","data.table","progress")

  args = commandArgs(trailingOnly=T)
  
  date<-ymd(substr(args[1], 11, 16));runPath <-paste('/Volumes/NGS/Analysis/covidseq',args[1],sep="/")
  df1 = read.csv(paste(runPath,"covidseq_output/SampleSheet_Intermediate.csv",sep ="/"))
  df1<-tail(df1,-13);df1<-df1[,c(1,3)];
  colnames(df1) <- c("Sample_Accession","Sample_Type")
  M<-filter(df1, Sample_Type == "PatientSample")
  M<-M["Sample_Accession"]
  cat(" Detected ", length(M$Sample_Accession), " patient's samples")
# Getting Pangolin lineage 
  pango<-read.csv(paste(runPath,"covidseq_output/lineage/combined_lineage_report.csv",sep ="/"))[,c("taxon","lineage")]
  var <-paste('',args[1],sep='-'); head(var)
  pango<-pango %>% mutate_at("taxon", str_replace, var, "")
  M<-merge(M,pango, by.x="Sample_Accession", by.y= "taxon", all.x=T)
  M$lineage[is.na(M$lineage)]= "Not able to be sequenced"
  head(M, n = 20)
# From LIMS
  cat("Searching the LIMS")
  NGS_covid<-rownames(fileSnapshot("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/",
                                   pattern=glob2rx("NGS*.csv"))$info[which.max(fileSnapshot("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/",                                                                                           pattern=glob2rx("NGS*.csv"))$info$mtime),])
  NGS_covid <- fread(paste('/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/',NGS_covid,sep=""), 
                     select = c("submitterId","sampleNumber","lastName","firstName","DOB","collectionDate")) 
  NGS_covid<- NGS_covid %>% mutate(across(, ~as.character(.)))
  
  NGS_covid1<-select(NGS_covid,c("sampleNumber","lastName","firstName","DOB","collectionDate"))
  NGS_covid2<-select(NGS_covid,c("submitterId","lastName","firstName","DOB","collectionDate"))
  head(NGS_covid1, n = 200);  head(NGS_covid2, n = 200)
  All_nocovid<-rownames(fileSnapshot("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/",
                                     pattern=glob2rx("All*.csv"))$info[which.max(fileSnapshot("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/", 
                                                                                              pattern=glob2rx("All*.csv"))$info$mtime),])
  All_nocovid <- fread(paste('/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/',All_nocovid,sep=""), 
                       select = c("sampleNumber","lastName","firstName","DOB","collectionDate")) 
  All_nocovid<- All_nocovid %>% mutate(across(, ~as.character(.)))
  head(All_nocovid, n = 200)
  colnames(NGS_covid2)[1]<-("sampleNumber")
  LIMS <-bind_rows(NGS_covid1, NGS_covid2, All_nocovid)
  Matrix<-merge(M,LIMS, by.x = "Sample_Accession", by.y = "sampleNumber", all.x = T)
  Matrix1<-merge(Matrix$Sample_Accession[which(!is.na(Matrix$collectionDate))],Matrix, by.x = "x", by.y ="Sample_Accession", all.x = T)
  colnames(Matrix1) <- c("sampleNumber","Test result","lastName","firstName","DOB","collectionDate")
# From EpiTrax
  EpiTrax_Export<-rownames(fileSnapshot("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/", 
                                        pattern=glob2rx("export_1551*.csv"))$info[which.max(fileSnapshot("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/",
                                                                                                         pattern=glob2rx("export_1551*.csv"))$info$mtime),])
  cat("Please wait, searching EpiTrax","\n")
  M_epitrax <- fread(paste('/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/',EpiTrax_Export,sep=""), 
                       select = c("lab_accession_no","person_last_name","person_first_name", "patient_birth_date","lab_collection_date"))   
 
  Mx<-merge(Matrix$Sample_Accession[which(is.na(Matrix$collectionDate))],M, by.x ="x", by.y ="Sample_Accession", all.x = T)
  Matrix2<-merge(Mx,M_epitrax, by.x = "x", by.y ="lab_accession_no", all.x = T)
  
  Matrix2$lab_collection_date<-as.Date(Matrix2$lab_collection_date)
  Matrix2<-unique(Matrix2); Matrix2<- Matrix2 %>% mutate(across(, ~as.character(.)))
  colnames(Matrix2) <- c("sampleNumber","Test result","lastName","firstName","DOB","collectionDate") 
  head(Matrix2, n = 200)
  ngs_file<-bind_rows(Matrix1, Matrix2)
  ngs_file<-cbind(ngs_file,'Name of facility performing the test'='UPHL', 'ordering facility' ='UHPL',
                   'Name or description of the test ordered'="SARS-CoV-2 Whole Genome Sequencing",
                   'lab test date'=ymd(date),'lab test status' ="preliminary")
  
  ngs_file<-select(ngs_file,"lastName",	"firstName",	"DOB",	"Name of facility performing the test","Name or description of the test ordered",
                    "collectionDate","sampleNumber","lab test date","lab test status",	"Test result",	"ordering facility")
  ngs_file[is.na(ngs_file)] <- ""
# Saving results
  d=paste(runPath,'R-sumcovidseq',sep="/")
  dir.create(file.path(d,"" ), recursive = T)
  LL<-paste(substring(Sys.time(), 1,10),substring(Sys.time(), 12,13),substring(Sys.time(), 15,16),substring(Sys.time(), 18,19),sep = "-")
  dx<- paste(paste('ngs',args[1],sep="_"),LL,".csv",sep="")
  dx<-paste(d,dx, sep ='/');write.csv(ngs_file,dx,row.names=F)
  cat("ngs file is completed and can be found at: ",d)




