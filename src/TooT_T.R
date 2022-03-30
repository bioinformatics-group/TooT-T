#! /usr/bin/Rscript

###################################################
## Name: TooT_T
## input: fasta file containing the unknown/testing protein sequences
## output: The predicted class of each protein sequence, 1= transporter, 0=non-transporter in csv format
## author: Munira Alballa
## reference:
##################################################

standardized<-function(x,rmean,rsd){((x-rmean)/rsd)}
pop.sd<-function(x){sqrt(sum((x-mean(x))^2)/length(x))}
normalize<- function(matrix){
  data=matrix
  standardizedData<- matrix
  for( i in 1:length(data[,1]) )#until L
  {
    #compute mean across the 20 aa
    rmean= mean(data[i,])
    rsd=pop.sd(data[i,])
    standardizedData[i,]<- standardized(data[i,],rmean,rsd)
  }
  return(standardizedData)
  
  
}


getpredictionsATH<- function(fastafile)
{
  # perform blast against TCDB with different threasholds: TCDB_Exact, TCDB_High, TCDB_Med
  querySeq=fastafile
  tcdbpath=paste0(dbpath,"/TCDB.fasta")
  blastpSeq(seq=fastafile, database.path=tcdbpath,output.path=intermediateFiles)

  results<-read.table(paste0(intermediateFiles,"out.txt"))
  colnames(results)<- c("qseqid","sseqid","pident", "length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qcovs")
  results$qseqid= sub(".+?\\.","", results$qseqid)
  tcdbseq<- read.fasta(tcdbpath,seqtype = "AA",as.string=T)
  tcdbnames<-sub("\\|.*","",sub(".+?\\|","",names(tcdbseq)))
  tcdblengths<- nchar(tcdbseq)
  QuerySeq= read.fasta(querySeq,seqtype = "AA",as.string=T)
  Quernames<-names(QuerySeq)
  Qlengths<- nchar(QuerySeq)

  Qlen<- c()
  Slen<- c()
  for(i in 1:length(results[,1]))
  {
    #query length
    leng<- unname(Qlengths[which(Quernames==results$qseqid[i])])
    if(length(leng)==0)
      print(results$qseqid[i])
    Qlen<-c(Qlen,leng)
    #subject length
    len<- unname(tcdblengths[which(tcdbnames==results$sseqid[i])])
    Slen<-c(Slen,len)
  }

  results2<- cbind(results, Querylength=Qlen, Slength=Slen, diffPer= (abs(Qlen-Slen)/Qlen) )

  HighThreasholdIDS<-results$qseqid[which((results$pident>=40 )& (results$evalue<= 1e-20) &(results$qcovs >= 70) & (results2$diffPer <= .1)) ]
  ExactMatchesID<- results$qseqid[which(results$pident==100.00 & results$evalue<0.000000000000001)]
  MedIDs<-results$qseqid[which(results$evalue<=1e-8 )]
  
  seqs<- readFASTA(test_fasta)
  ATH_predictions<- matrix(data=0,nrow=length(names(seqs)), ncol=3)
  row.names(ATH_predictions)<-names(seqs)
 
  
  ATH_predictions[which( names(seqs) %in% ExactMatchesID),1]<- 1
  ATH_predictions[which( names(seqs) %in% HighThreasholdIDS),2]<- 1
  ATH_predictions[which( names(seqs) %in% MedIDs),3]<- 1
  
  return(ATH_predictions)
    
    
}


getpredictions_psibasedmodels<-function(fastafile)
{
  
  files<- c("MSAAACpsi3.csv","MSADCpsi3.csv","MSAPAACpsi3.csv")
  Texnames<-  c("psiAAC","psiPAAC","psiPseAAC")
  
  seqs<- readFASTA(fastafile)
  names(seqs)<- sub("\\|.*","",sub(".+?\\|","", names(seqs)))
  psi_predictions<- matrix(data=0,nrow=length(names(seqs)), ncol=3)
  row.names(psi_predictions)<-names(seqs)
  
  for (z in 1:length(files))
  {
    t1=read.csv(paste0(compostions, files[z]), header = TRUE,sep=",", stringsAsFactors = FALSE)
    
    testpredictors<- as.matrix(t1[,c(-1,-2)])
    testpredictors<- normalize(testpredictors)
    load( paste0(TooTTdir,"/models/",Texnames[z],"_TvdNT.rda"))
    colnames(testpredictors)<- attr(svm.fit$terms, "term.labels")
    row.names(testpredictors)<-names(seqs)
    p1<- predict(svm.fit,testpredictors,probability=TRUE)
    psi_predictions[rownames(attr(p1,"probabilities")),z]<- as.vector(p1)
  }
  
  return(psi_predictions)
}

get_gbm_ensemble<-function(metadata)
{
  tes<- metadata
  tes<- mapply(tes, FUN=as.vector)
  tes<- mapply(tes, FUN=as.numeric)
  
  tes<-matrix(data=tes, ncol=length(metadata[1,]), nrow=length(metadata[,1]))
  
  load(paste0(paste0(TooTTdir,"/models/","meta_gbm.rda")))
  pred <- predict(xgb.fit.final,tes )
  results<-ifelse(pred>.5,1,0)
  
  return(results)
  
}
args <- commandArgs(trailingOnly=TRUE)

terminate <- FALSE

out <- "."
TooTTdir <- "."
db<-"./db/"
for(i in args){
  arg = strsplit(i, "=")[[1]];
  
  switch(arg[1],
         "-query"={
           query <- arg[2]
         },
         "-out"={
           out <- arg[2]
         },
         "-TooTT"={
           TooTTdir <- normalizePath(arg[2])
         },
         "-db"={
           db <- normalizePath(arg[2])
         },
         "-help"={
           cat("TooTT v1.0 (Oct. 2019)\n")
           cat("\n")
           cat("Usage: TooTT -query=<input> [-TooTT=<TooTTdir>] [-out=<outdir>] [-db=<database path>]\n")
           cat("\n")
           cat("\t<input> is your sequence input file in fasta format\n")
           cat("\t<out> is the output directory where you want the predicted results, formatted as csv\n")
           cat("\t\t<out> defaults to '",out,"'\n")
           cat("\t<TooTTdir> is the directory where the base TooT-T files are located")
           cat("\t\t<TooTTdir> defaults to '",TooTTdir,"'\n")
           cat("\t <database path> is the relative path to the database, it should include TCDB for retrieving ATH prediction in addition to 
               the choice of homology database for psi-compositions, tested using Swiss-Port databses (2018) \n")
           cat("\t\t<database path> defaults to",TooTTdir,"/db/\n")
           cat("\n")
           terminate <- TRUE
           break
         }
  )
}

if(!terminate) {
  
  if(!exists("query")) {
    stop("-query has not been passed")
  }
  #TooTTdir="./TooT-T/"

  test_fasta <- normalizePath(path.expand(query))
  #test_fasta<- "./TooT-T/Transportertest_short.fasta"
  resultspath <- paste0(normalizePath(path.expand(out)),"/")
  #resultspath<-"./TooT-T/TooT-T-output/"
  dbpath=paste0(TooTTdir, "/db/")
  require(seqinr)
  library("Biostrings")
  library("stringr")
  require(protr)
  library(ISLR)
  library(e1071)
  library(caret)
  library(R.utils)
  wd=normalizePath(path.expand(".")) # change the the tool directory
  
  if (isAbsolutePath(db)){
    dbpath <- db
  }else{
    dbpath <- paste0(TooTTdir, db)
  }
  
  compostions=paste0(TooTTdir,"/intermediate_files/Compositions/")
  intermediateFiles=paste0(TooTTdir,"/intermediate_files/")

  #testing data with unknown substrates
  source(paste0(TooTTdir,"/src/psi_compostions.R"))
  
  psi_compostions(test_fasta)
  ATH<-getpredictionsATH(test_fasta)
  ML<-getpredictions_psibasedmodels(test_fasta)
  #meta takes following order:
  #("tcdb_exact","tcdb_high","tcdb_med", "psiAAC", "psiPAAC", "psiPseAAC")
  predictions<- get_gbm_ensemble(cbind.data.frame(ATH,ML))
  
  
  
  # write results
  seqs<- readFASTA(test_fasta)
  names(seqs)<- sub("\\|.*","",sub(".+?\\|","", names(seqs)))
  print(paste0( "TooT-Toutput is found at: ", resultspath, "TooTTout.csv"))
  write.csv(cbind(UniProtID=names(seqs),predictions ),paste0(resultspath,"TooTTout.csv"))
  
  
}
