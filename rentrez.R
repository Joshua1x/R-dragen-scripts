#/usr/bin/env Rscript
# Author: JLinares/08/2022
#       Usage: Getting fasta files from NCBI/GenBank/nuccore for a specific time window
#       Parameters: time-window. ie: Rscript /. . . ./GenBank_fasta_rentrez.R 220705 220708
library("easypackages");libraries("ape","rbin","rentrez","lubridate","seqinr","tidyverse")
options(timeout = 700, "width"=400)
args = commandArgs(trailingOnly=T);date1<-ymd(substr(args[1], 1, 6));date2<-ymd(substr(args[2], 1, 6))
ff<-!is.na(list(args[1],args[2]));
if(ff[1]!=T |ff[2]!=T) {cat(crayon::bold("Arguments entry error, Try again1","\n"))
} else {

  if (date2 > date1){
      base<-"HUMAN AND UT NOT UPHL AND SARS-CoV-2 AND 25000:30000[SLEN] AND (COLLECTION_DATE="
      # sequences length: 25k to 30K, inclusive.
      a<-seq(ymd(date1), ymd(date2), by="days");a
      cat(crayon::bold("Please wait, searching the NCBI/GenBank/nuccore DB for fasta sequences between 25k - 30k nucleotides long ","\n"))
      Records<-vector()
      for (i in 1:length(a)){
          tt<-a[i];locus= paste(base,tt,sep = "")
          r_search <- entrez_search(db="nuccore", term=locus)
          Records<-c(Records, r_search$ids)
          cat(r_search$ids,"\n")
          }
    cat(crayon::bold("Please wait, downloading/storing", length(Records), "fasta sequences","\n"))
    # NCBI restrictions: For long list of IDs: Use web history or split it in blocks
    bigqueryindex <- split(Records,ceiling(seq_along(Records)/300))
    SARS2_seq<-vector()

    for (i in bigqueryindex){
         part_seq<-read.GenBank(i,as.character = T , chunk.size = 300)
         SARS2_seq<-c(SARS2_seq,part_seq)
         }
    setwd("/Volumes/NGS/Bioinformatics/Lin")
    write.dna(SARS2_seq,paste("GenBank_",Sys.Date(),".fasta",sep =""), format = "fasta")
    cat(crayon::bold("sequences stored at ",getwd(), "as: ", paste("GenBank_",Sys.Date(),".fasta",sep =""),"\n"))

  } else { cat(crayon::bold("Arguments entry error, Try again!","\n"))}
}

