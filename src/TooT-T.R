#! /usr/bin/Rscript

###################################################
## Name: TooT_T
## input: fasta file containing the unknown/testing protein sequences
## output: The predicted class of each protein sequence, 1= transporter, 0=non-transporter in csv format
## author: Munira Alballa
## reference:
##################################################
  
suppressMessages(suppressWarnings(library(seqinr)))
suppressMessages(suppressWarnings(library("Biostrings")))
suppressMessages(suppressWarnings(library("stringr")))
suppressMessages(suppressWarnings(library(protr)))
suppressMessages(suppressWarnings(library(ISLR)))
suppressMessages(suppressWarnings(library(e1071)))
suppressMessages(suppressWarnings(library(caret)))
suppressMessages(suppressWarnings(library(R.utils)))
suppressMessages(suppressWarnings(library("funr")))

validateTool <- function(toolName, note) {
  tool <- Sys.which(toolName)
  if(tool == '')
    stop(cat("Please install ", toolName, " and make sure it is visiable in the path (", note, ").\n", sep=""))
  return(tool)
}

standardized<-function(x,rmean,rsd){((x-rmean)/rsd)}
pop.sd<-function(x){sqrt(sum((x-mean(x))^2)/length(x))}
normalize<- function(matrix){
  data=matrix
  standardizedData<- matrix
  for( i in 1:length(data[,1]) ) { #until L
    #compute mean across the 20 aa
    rmean= mean(data[i,])
    rsd=pop.sd(data[i,])
    standardizedData[i,]<- standardized(data[i,],rmean,rsd)
  }
  return(standardizedData)  
  
}


getpredictionsATH<- function(fastafile) {
  # perform blast against TCDB with different threasholds: TCDB_Exact, TCDB_High, TCDB_Med
  blastpSeq(seq=fastafile, output.path=intermediateFiles)

  results<-read.table(file.path(intermediateFiles,"out.txt"))
  colnames(results)<- c("qseqid","sseqid","pident", "length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qcovs")
  results$qseqid= sub(".+?\\.","", results$qseqid)
  tcdbseq<- read.fasta(tcdb,seqtype = "AA",as.string=T)
  tcdbnames<-sub("\\|.*","",sub(".+?\\|","",names(tcdbseq)))
  tcdblengths<- nchar(tcdbseq)
  currentSeqs= read.fasta(fastafile, seqtype = "AA",as.string=T)
  qNames<-names(currentSeqs)
  qLengths<- nchar(currentSeqs)

  Qlen<- c()
  Slen<- c()
  for(i in 1:length(results[,1])) {
    #query length
    leng<- unname(qLengths[which(qNames==results$qseqid[i])])
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
  
  currentSeqs<- readFASTA(fastafile)
  ATH_predictions<- matrix(data=0,nrow=length(names(currentSeqs)), ncol=3)
  row.names(ATH_predictions)<-names(currentSeqs)
 
  ATH_predictions[which( names(currentSeqs) %in% ExactMatchesID),1]<- 1
  ATH_predictions[which( names(currentSeqs) %in% HighThreasholdIDS),2]<- 1
  ATH_predictions[which( names(currentSeqs) %in% MedIDs),3]<- 1
  
  return(ATH_predictions)
       
}


getpredictions_psibasedmodels<-function(fastafile) {
  
  files<- c("MSAAACpsi3.csv","MSADCpsi3.csv","MSAPAACpsi3.csv")
  Texnames<-  c("psiAAC","psiPAAC","psiPseAAC")
  
  currentSeqs<- readFASTA(fastafile)
  names(currentSeqs)<- sub("\\|.*","",sub(".+?\\|","", names(currentSeqs)))
  psi_predictions<- matrix(data=0,nrow=length(names(currentSeqs)), ncol=3)
  row.names(psi_predictions)<-names(currentSeqs)
  
  for(z in 1:length(files)) {
    t1=read.csv(file.path(compositions, files[z]), header = TRUE,sep=",", stringsAsFactors = FALSE)
    
    testpredictors<- as.matrix(t1[,c(-1,-2)])
    testpredictors<- normalize(testpredictors)
    currentTvdNT <- file.path(TooTTdir, "models", paste0(Texnames[z], "_TvdNT.rda"))
    load(currentTvdNT)
    colnames(testpredictors)<- attr(svm.fit$terms, "term.labels")
    row.names(testpredictors)<-names(currentSeqs)
    p1<- predict(svm.fit,testpredictors,probability=TRUE)
    psi_predictions[rownames(attr(p1,"probabilities")),z]<- as.vector(p1)
  }
  
  return(psi_predictions)
}

get_gbm_ensemble<-function(metadata) {
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

out <- normalizePath(".")
TooTTdir <- normalizePath(file.path(funr::get_script_path(), ".."))
db <- normalizePath(file.path(TooTTdir, "db"), mustWork = FALSE)
work <- normalizePath(".")
for(i in args){
  arg = strsplit(i, "=")[[1]];
  
  switch(arg[1],
         "-query"={
           toottquery <- normalizePath(arg[2])
         },
         "-out"={
           out <- normalizePath(arg[2])
         },
         "-TooTT"={
           TooTTdir <- normalizePath(normalizePath(arg[2]))
         },
         "-db"={
           db <- normalizePath(normalizePath(arg[2]))
         },
         "-help"={
           cat("TooTT v1.1 (Apr. 2022)\n")
           cat("\n")
           cat("Usage: TooTT -query=<input> [-out=<outdir>] [-db=<database path>] [-work=<work path>] [-TooTT=<TooTTdir>]\n")
           cat("\n")
           cat("\t<input> is your sequence input file in fasta format\n")
           cat("\t<out> is the output directory where you want the predicted results, formatted as csv\n")
           cat("\t\t<out> defaults to '.' ('",out,"')\n", sep="")
           cat("\t <database path> is the relative path to the database, it should include TCDB for retrieving ATH prediction in addition to
               the choice of homology database for psi-compositions, tested using Swiss-Port databses (2018) \n")
           cat("\t\t<database path> defaults to '",db,"'\n", sep="")
           cat("\t<work path> is the path to the working directory for intermediate files. It will be created as needed.\n")
           cat("\t\t<database path> defaults to '.' ('",work,"')\n", sep="")
           cat("\t<TooTTdir> is the directory where the base TooT-T files are located")
           cat("\t\t<TooTTdir> defaults to '",TooTTdir,"'\n", sep="")
           cat("\n")
           terminate <- TRUE
           break
         }
  )
}

if(!terminate) {
  
  if(!exists("toottquery")) {
    stop("-query has not been passed")
  }

  if(!file.exists(toottquery)) {
    stop("The specified query file does not exist: '", toottquery,"'", sep="")
  }

#
# Validate that the db directory and required db files exist
#
if(!file.exists(db)) {
   stop("The specified database directory does not exist: '", db,"'", sep="")
}

dbFiles <- c("SwissOct18.fasta.psi", "SwissOct18.fasta.psd", "SwissOct18.fasta.pog", "SwissOct18.fasta.psq", "SwissOct18.fasta.pin", "SwissOct18.fasta", "SwissOct18.fasta.phr")
missingFiles <- list()
for(file in dbFiles) {
   if(!file.exists(file.path(db, file))) {
      missingFiles <- append(missingFiles, file)
   }
}

if(length(missingFiles) > 0) {
   stop("Unable to find some files in your db directory ('", db,"')\n", paste(missingFiles, collapse=", "))
}

swissprotdb <- file.path(db, "SwissOct18.fasta");

dbFiles <- c("TCDB.fasta.psi", "TCDB.fasta.psd", "TCDB.fasta.pog", "TCDB.fasta.psq", "TCDB.fasta.pin", "TCDB.fasta", "TCDB.fasta.phr")
missingFiles <- list()
for(file in dbFiles) {
   if(!file.exists(file.path(db, file))) {
      missingFiles <- append(missingFiles, file)
   }
}

if(length(missingFiles) > 0) {
   stop("Unable to find some files in your db directory ('", db,"')\n", paste(missingFiles, collapse=", "))
}

tcdb <- file.path(db, "TCDB.fasta");



#
# Validate that tools exist: psiblast and mview, specifically
#
blastp.path <- validateTool('blastp', 'included in the NCBI BLAST+')
psiblast.path <- validateTool('psiblast', 'included in the NCBI BLAST+')

mview.path <- validateTool('mview', 'found at https://desmid.github.io/mview/')

#
# Validate the outpit dir
#
if(!file.exists(out)) {
   stop("The specified output directory does not exist: '", out,"'", sep="")
}

#
# Validate and set up the working directory
#
if(!file.exists(work)) {
   stop("The specified base for your working directory does not exist: '", work,"'", sep="")
}

if(!file.exists(file.path(work, "work"))) {
   dir.create(file.path(work, "work"))
}

intermediateFiles = file.path(work, "work", "TooT-T")
if(!file.exists(intermediateFiles)) {
   dir.create(intermediateFiles)
}

compositions = file.path(intermediateFiles, "Compositions")
if(!file.exists(compositions)) {
   dir.create(compositions)
}

#
# Lastly, validate the TooT-T Dir that they might have passed... I hate this param. We should also validate that all other required pieces are there
# which should catch broken installs
#
if(!file.exists(TooTTdir)) {
   stop("The specified base for your application does not exist does not exist (you may want to consider using the dynamically generated one): '", TooTTdir,"'", sep="")
}

psi_compositionsSource <- file.path(TooTTdir, "src", "psi_compositions.R")
if(!file.exists(psi_compositionsSource)) {
   stop("A required source file to run TooTT does not exist, though it should be located next to this script: '", psi_compositionsSource,"'", sep="")
}

  #testing data with unknown substrates
  source(psi_compositionsSource)
  
  psi_compositions(toottquery)
  ATH<-getpredictionsATH(toottquery)
  ML<-getpredictions_psibasedmodels(toottquery)
  #meta takes following order:
  #("tcdb_exact","tcdb_high","tcdb_med", "psiAAC", "psiPAAC", "psiPseAAC")
  predictions<- get_gbm_ensemble(cbind.data.frame(ATH,ML))
  
  # write results
  currentSeqs<- readFASTA(toottquery)
  names(currentSeqs)<- sub("\\|.*","",sub(".+?\\|","", names(currentSeqs)))
  print(paste0("TooT-Toutput is found at: ", file.path(out,"TooTTout.csv")))
  write.csv(cbind(UniProtID=names(currentSeqs), predictions), file.path(out,"TooTTout.csv"))  
  
}
