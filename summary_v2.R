#!/usr/bin/env Rscript
# Author: OLP/2022
# Usage: summary for PatientSamples
library("plyr", warn.conflicts = F)
library("dplyr", warn.conflicts = F)
library("rbin")#
library("tidyverse")
library("lubridate", warn.conflicts = F)
library("tibble")
library("stringr") 
require("seqinr"); 
library("xml2")

args = commandArgs(trailingOnly=T)
date<-ymd(substr(args[1], 11, 16));runPath <-paste('/Volumes/NGS/Analysis/covidseq',args[1],sep="/"); head(date)
runPath


# Read Type & Read Length from output/ . . . ./RunParameters
fileName<-"/Volumes/NGS/Output/A01290/220416_A01290_0116_BHLMVHDSX3/RunParameters.xml"
ReadType<-xmlToDataFrame(nodes=getNodeSet(xmlParse(read_xml(fileName)),"//ReadType"));ReadType
ReadLength<-xmlToDataFrame(nodes=getNodeSet(xmlParse(read_xml(fileName)),"//Read1NumberOfCycles"));ReadLength 
#=============================================
# Paths to Analysis
ss1 = paste(runPath,"covidseq_output/SampleSheet_Intermediate.csv",sep ="/"); 
ss2 = paste(runPath,"covidseq_output/lineage/combined_lineage_report.csv",sep ="/"); 
ss3 = paste(runPath,"covidseq_output/Logs_Intermediates/VirusDetection/kmer/Summary.tsv",sep ="/");
#=============================================
df1 = read.csv(ss1); # SampleSheet
df2 = read.csv(ss2); # combined_llineage_report
df3 = read.table(file = ss3, sep = '\t', header = T); # sc2_amplicons 
#=============================================
# Type of Sequencer
Sequencer <-substr(args[1], 4, 4); Sequencer
if (Sequencer =="A"){ TypeOfSequencer<-"NovaSeq";Pt<-"A01290"} else {TypeOfSequencer<-"NextSeq";Pt<-"NB551133"};Pt; 
t=.5;Sys.sleep(t)
#===========================================
# Searching for basic info from samplesheet
instrument_type<-df1[4,2];cat("instrument_type:", instrument_type,"\n");Sys.sleep(t)
wetlab_assay<-df1[5,2]; cat("wetlab_assay:", wetlab_assay,"\n");Sys.sleep(t)
pangolin_version<-df2[1,9];cat("pangolin_version :",pangolin_version,"\n");Sys.sleep(t)
pangoLEARNversion<-df2[1,8];cat("pangoLEARNversion :",pangoLEARNversion,"\n");Sys.sleep(t)
software = "Illumina DRAGEN COVID pipeline (R.U.O.), vs. 1.0.1"; cat("software: ",software,"\n");Sys.sleep(t)
lineage_tools = "Dragen-covid-lineage-tools, vs. 1.1.0.0"; cat ("Lineage tools: ",lineage_tools,"\n");Sys.sleep(t)
#=============================================
# Read Type & Read Length from output/ . . . ./RunParameters
fileName<-"/Volumes/NGS/Output/A01290/220416_A01290_0116_BHLMVHDSX3/RunParameters.xml"
ReadType<-xmlToDataFrame(nodes=getNodeSet(xmlParse(read_xml(fileName)),"//ReadType"));ReadType
ReadLength<-xmlToDataFrame(nodes=getNodeSet(xmlParse(read_xml(fileName)),"//Read1NumberOfCycles"));ReadLength 
#=============================================
# n of PatientSample. 
df1<-tail(df1,-13);
df1<-df1[,c(1,3)]; 
colnames(df1) <- c("Sample_Accession","Sample_Type");
df1<-filter(df1, Sample_Type == "PatientSample"); 
n<-length(df1$Sample_Type); cat("Number of PatientSamples",n,"\n");Sys.sleep(t)
#=========================================
# Random Number Generator for Sample_ID
RNG<-paste(floor(runif(n, min=99, max=1000)),floor(runif(n, min=99, max=1000)),sep =""); 
# Checking duplicates of random numbers
while(length(RNG)!= length(unique(RNG)))
     { RNG<-paste(floor(runif(n, min=99, max=1000)),floor(runif(n, min=99, max=1000)),sep ="");}
IDD<-paste("UT-UPHL",substr(args[1], 11, 16),sep="-");
Sample_ID<-paste(IDD,RNG,sep=""); 
summary1<-cbind(df1,Sample_ID);
#============================================
#left join: SampleSheet & Combined Lineage Report files
var <-paste('',args[1],sep='-'); head(var)
df2 <-df2 %>% mutate_at("taxon", str_replace, var, "")
df2<-df2[,c("taxon","lineage","scorpio_call")]
head(df2); #Taxon, lineage, scorpio_call
df3 = read.table(file = ss3, sep = '\t', header = TRUE); # 
df3<-df3 %>% mutate_at("Sample", str_replace, var, ""); head(df3)
df3<-df3[,c("Sample","nTargetsDetected.SARS.CoV2")]; 
colnames(df3)[2] <- 'sc2_amplicons';#
summary1<-merge(summary1,df2, by.x="Sample_Accession", by.y= "taxon", all.x=TRUE); 
head(summary1)
summary1<-merge(summary1,df3, by.x="Sample_Accession", by.y= "Sample", all.x=TRUE);
head(summary1)
#==========================================
# adding basic information: 
summary1<-cbind(summary1,'wetlab_assay'= wetlab_assay,"instrument_type"=instrument_type,"software"=software,"lineage_tools"=lineage_tools,
                "pangolin_version"=pangolin_version,'read_Type' = ReadType,"readLength" = ReadLength)
                colnames(summary1)[12]<-("Read_Type");
                colnames(summary1)[13]<-("Read_Length");
                head(summary1)
#========================================
#Nucleotides:num of actg & n
num_actgx<-vector();num_nx<-vector(); 
runPath1 = paste(runPath,"covidseq_output/Sample_Analysis/",sep ="/"); 
num_actgx = 0; num_nx = 0;
for (i in 1:n){
     sample=summary1$Sample_Accession[i]; 
     amplicons=summary1$sc2_amplicons[i];
      if(amplicons >=50)
      { 
        Lyapunovv<- Sys.glob(paste(runPath1,sample,"*/*.fasta",sep = ""),dirmark = TRUE); 
        a =  str_count(read.fasta(Lyapunovv, as.string = T, forceDNAtolower = T), pattern = "a|c|t|g");num_actgx[i] = a;
        b =  str_count(read.fasta(Lyapunovv, as.string = T, forceDNAtolower = T))-a; num_nx[i] = b;
      }
      else {# AmpliconThereshold <50, no fasta file
             a = b =0; num_actgx[i] = 0; num_nx[i] = 0; 
           }
     cat(sample," \t",a,"\t",b,"\n")
      }
df_new <- cbind(as.character(num_actgx), as.character(num_nx)); df_new
summary<-cbind(summary1,df_new); 
colnames(summary)[14]<-("num_actg");colnames(summary)[15]<-("num_n"); 
#========================================
# save results
d=paste(runPath,'Rsumcovidseq',sep="/"); 
dir.create(file.path(d,"" ), recursive = TRUE)
write.csv(summary,paste(d,paste("Rcovidseq_summary",Sys.time(),".csv",sep="-"),sep="/"),row.names=FALSE)
cat("Summary file is complete and can be found at ",paste(d,'Rcovidseq_summary.csv',sep="/"))



