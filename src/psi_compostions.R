firstrun=F
blastpSeq<- function(seq, start.pos = 1L, end.pos = nchar(seq), 
                     blastp.path = NULL, makeblastdb.path = NULL, 
                     database.path = NULL, silent = TRUE, 
                     evalue = 10L, output.path=resultspath){
  
  if (Sys.which('makeblastdb') == '' & is.null(makeblastdb.path))
    stop('Please install makeblastdb (included in the NCBI BLAST+) or specify makeblastdb.path.')
  
  if (Sys.which('blastp') == '' & is.null(blastp.path))
    stop('Please install blastp (included in the NCBI BLAST+) or specify blastp.path.')
  
  makeblastdb.path = if (!is.null(makeblastdb.path)) makeblastdb.path else Sys.which('makeblastdb')
  blastp.path = if (!is.null(blastp.path)) blastp.path else Sys.which('blastp')
  
  if (is.null(database.path)) stop('Must specify the database (a FASTA file) path')
  if (is.null(output.path)) stop('Must specify the output path')
  
  N = end.pos - start.pos + 1L
  
  # Prepare data for Blastp
  cmddb = paste0(shQuote(makeblastdb.path), ' -dbtype prot -in ', 
                 shQuote(database.path),' -parse_seqids')        
  if(firstrun)
  {
    print("performing make blastdb")
    # print( cmddb)
    if (silent == TRUE) system(cmddb, ignore.stdout = TRUE) else system(cmddb)
  }
  
  # Basic parameters for Blastp
  tmp = tempfile('Blastp')
      
  # Additional parameters for Blastp
  if (!is.null(evalue)) {
    if (evalue <= 0L) {
      stop('evalue must be > 0')
    }
  }
  # Run Blastp
  cmdblastp = paste(
    paste0(shQuote(blastp.path),
           ' -comp_based_stats 1 -db ', shQuote(database.path),
           ' -query ', shQuote(seq),  ' -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs\"  -num_alignments 1',' -out ', paste0(shQuote(output.path),"out.txt")))
  
  print("******************************")
  print(cmdblastp)
  if (silent == TRUE) system(cmdblastp, ignore.stdout = F) else system(cmdblastp)      
  #get the hit sequences Id
  data = read.table(paste0(output.path,"out.txt"))
 
  
}


psiblastSeq <- function(seq,
                        start.pos = 1L,
                        end.pos = nchar(seq),
                        psiblast.path = NULL,
                        database.path = NULL,
                        silent = TRUE,
                        evalue = 10L,
                        numOfIterations = 3,
                        output.path = path) {
  if (Sys.which('psiblast') == '' & is.null(psiblast.path))
    stop('Please install psiblast (included in the NCBI BLAST+) or specify psiblast.path.')
  
  blastp.path = if (!is.null(psiblast.path))
    psiblast.path
  else
    Sys.which('psiblast')
  
  if (is.null(database.path))
    stop('Must specify the database (a FASTA file) path')
  if (is.null(output.path))
    stop('Must specify the output path')
  
  N = end.pos - start.pos + 1L
  
  
  # Basic parameters for Blastp
  tmp = tempfile('Blastp2')
  queryFasta = paste0(tmp, '.fasta')
  tabularfile = paste0(tmp, '.txt')
  querySeq = Biostrings::AAStringSet(as.character(seq))
  Biostrings::writeXStringSet(querySeq, queryFasta)
  # Additional parameters for Blastp
  if (!is.null(evalue)) {
    if (evalue <= 0L) {
      stop('evalue must be > 0')
    }
  }
  # Run Blastp
  cmdblastp = paste(
    paste0(
      shQuote(blastp.path),
      ' -db ',
      shQuote(database.path),
      ' -query ',
      shQuote(queryFasta),
      ' -num_iterations ',
      numOfIterations ,
      ' -evalue ',
      evalue,
      
      " -out psiblastDsbAOut -out_pssm PSSMDsbA -out_ascii_pssm ascii_mtx_file"
    )
  )
  
  print(cmdblastp)
  if (silent == TRUE)
    system(cmdblastp, ignore.stdout = TRUE)
  else
    system(cmdblastp)
  
  
}
ClaculateCompostions_return<- function(subdirName, filename)
{ 
  
  con= paste0(subdirName,filename)
  print(con)
  se<- read.fasta(con,seqtype = "AA",as.string=T)
  
  seqlist<- unlist(lapply(se,as.character))
  seqlist<-sapply(seqlist,gsub,pattern="[^GPAVLIMCFYWHKRQNEDST]",replacement="") # getting rid of gaps? not the smartest thing to do, need more work
  names(seqlist)<- sub(".+\\:","",names(se))
  names(seqlist)<- sub("\\|.*","",sub(".+?\\|","", names(seqlist)))
  seqlist=seqlist[which(nchar(seqlist)>30)]
  #print(names(seqlist))
  
  MSAAAC<- lapply(seqlist, extractAAC)
  outputAAC <- matrix(unlist(MSAAAC), ncol = 20, byrow = TRUE)
  rownames(outputAAC)<- names(seqlist)
  write.csv(outputAAC, paste0(subdirName,"AAC3",".csv"))
  MSAPseAAC<- lapply(seqlist, extractPAAC)
  
  outputPseAAC <- matrix(unlist(MSAPseAAC), ncol = 50, byrow = TRUE)
  rownames(outputPseAAC)<- names(seqlist)
  write.csv(outputAAC, paste0(subdirName,"PseAAC3",".csv"))
  
  MSADC<- lapply(seqlist, extractDC)
  outputDC <- matrix(unlist(MSADC), ncol = 400, byrow = TRUE)
  rownames(outputDC)<- names(seqlist)
  write.csv( outputDC,paste0(subdirName,"DC3",".csv"))
  return(list(AAC=outputAAC, DC=outputDC, PseAAC=outputPseAAC))
  
}

PreparedataforPSI = function(seq, start.pos = 1L, end.pos = nchar(seq), 
                                blastp.path = NULL, makeblastdb.path = NULL, 
                                database.path = NULL, iter = 5, silent = TRUE, 
                                evalue = 10L, output.path=resultspath) {
  #1- run blastp on datafiles
  seqname=names(seq)
  SeqDirectory=paste0(output.path,names(seq),"/")
  
  dir.create(SeqDirectory, showWarnings = FALSE, recursive = FALSE, mode = "0777")   
  setwd(SeqDirectory)
  
  psiblastSeq(
    seq,
    start.pos ,
    end.pos ,
    psiblast.path = NULL,
    database.path ,
    silent ,
    evalue = 0.001,
    numOfIterations = 3,
    SeqDirectory
  )
  
  
  
  
}

converttofasta<- function(seq, output.path)
{
  dirName<-paste0(output.path,names(seq),"/")
  
  for (iteration in 1:3)
  {
    name = paste0("psi", iteration, ".fasta")
    command = paste0(
      MviewPath,
      " -in blast -cycle ",
      iteration,
      " ",
      dirName,
      "psiblastDsbAOut -out fasta> ",
      dirName,
      name
    )
    print(command)
    system(command)
  }
  
}

ClaculateCompostions<- function(subdirName, filename)
{ 
  
  con= paste0(subdirName,filename)
  print(con)
  se<- read.fasta(con,seqtype = "AA",as.string=T)
  
  seqlist<- unlist(lapply(se,as.character))
  seqlist<-sapply(seqlist,gsub,pattern="[^GPAVLIMCFYWHKRQNEDST]",replacement="") # getting rid of gaps? not the smartest thing to do, need more work
  names(seqlist)<- sub(".+\\:","",names(se))
  names(seqlist)<- sub("\\|.*","",sub(".+?\\|","", names(seqlist)))
  seqlist=seqlist[which(nchar(seqlist)>30)]
  if(length(seqlist)==0)
    return(NA)
  #print(names(seqlist))
  
  MSAAAC<- lapply(seqlist, extractAAC)
  outputAAC <- matrix(unlist(MSAAAC), ncol = 20, byrow = TRUE)
  rownames(outputAAC)<- names(seqlist)
  write.csv(outputAAC, paste0(subdirName,"AAC3",".csv"))
  
  
  MSAPseAAC<- lapply(seqlist, extractPAAC)
  outputPseAAC <- matrix(unlist(MSAPseAAC), ncol = 50, byrow = TRUE)
  rownames(outputPseAAC)<- names(seqlist)
  write.csv(outputPseAAC, paste0(subdirName,"PseAAC3",".csv"))
  
  MSADC<- lapply(seqlist, extractDC)
  outputDC <- matrix(unlist(MSADC), ncol = 400, byrow = TRUE)
  rownames(outputDC)<- names(seqlist)
  write.csv( outputDC,paste0(subdirName,"DC3",".csv"))
  
}



#MviewPath = "~/mview-1.64/bin/./mview"
ClaculateEachCompostions<- function(seq, output.path)
{
  subdirName= paste0(intermediateFiles,names(seq),"/")
  
  
  if(readLines(paste0(subdirName,"psi3.fasta"),n=1)!="mview: no alignments found")
  {
    ClaculateCompostions(subdirName,"psi3.fasta")
  }else if (readLines(paste0(subdirName,"psi2.fasta"),n=1)!="mview: no alignments found"){
    print("no psi3.fasta but Pse2")
    ClaculateCompostions(subdirName,"psi2.fasta")
  }else if (readLines(paste0(subdirName,"psi1.fasta"),n=1)!="mview: no alignments found"){
    print("no psi2.fasta but Pse1")
    ClaculateCompostions(subdirName,"psi1.fasta")
  }
  else {
    print("nomatch")
  }
  
  closeAllConnections()
  
  
}

psi_compostions<- function(fastafile)
{
  seqs<- readFASTA(fastafile)
  names(seqs)<- sub("\\|.*","",sub(".+?\\|","", names(seqs)))
  for(j in c(1:length(seqs)))
  {
    x<- seqs[j]
    PreparedataforPSI(seq= x,database.path=paste0(dbpath,"/SwissOct18.fasta"),output.path=intermediateFiles, evalue=.001)
    converttofasta(seq= x,intermediateFiles) #write function rather than using Mview # todo
    ClaculateEachCompostions(seq= x,intermediateFiles)
  }
  
  
  AAfiles<-  names(seqs)
  ListofSubsetDetailedMSAAACStatistics<- matrix(ncol = 20+1,nrow = length(AAfiles)) #contrains AAC for each seq in the subset
  ListofSubsetDetailedMSADCStatistics <- matrix(ncol=400+1,nrow=length(AAfiles))
  ListofSubsetDetailedMSAPAACStatistics <- matrix(ncol=(50+1),nrow=length(AAfiles))
  
  
  for(i in c(1:length(AAfiles)))
  {
    
    subdirName<- paste0(intermediateFiles,AAfiles[i],"/")
    
    number=3
    if(file.exists(paste0(subdirName,"/AAC",number,".csv")))
    {
      outputAAC <- read.csv( paste0(subdirName,"/AAC",number,".csv"), stringsAsFactors = FALSE)
      Ids<- outputAAC[,1]
      
      outputAAC <-outputAAC[,-1]
      AAC<- apply(outputAAC, 2, mean)
      
      outputDC <- read.csv( paste0(subdirName,"/DC",number,".csv"), stringsAsFactors = FALSE)
      outputDC <-outputDC[,-1]
      DC<- apply(outputDC, 2, mean)
      
      outputPseAAC <- read.csv( paste0(subdirName,"/PseAAC",number,".csv"), stringsAsFactors = FALSE)
      outputPseAAC <-outputPseAAC[,-1]
      PseAAC<- apply(outputPseAAC, 2, mean)
      
      ListofSubsetDetailedMSAAACStatistics[i,]<-  c(NA,AAC)
      ListofSubsetDetailedMSADCStatistics[i,]<-   c(NA,DC)
      ListofSubsetDetailedMSAPAACStatistics[i,]<-  c(NA,PseAAC)    
      
    }else{
    }
    
  }
  
  
  ListofSubsetDetailedMSAAACStatistics<-cbind(UniprotID=AAfiles,ListofSubsetDetailedMSAAACStatistics)
  ListofSubsetDetailedMSADCStatistics<- cbind(UniprotID=AAfiles,ListofSubsetDetailedMSADCStatistics)
  ListofSubsetDetailedMSAPAACStatistics<-cbind(UniprotID=AAfiles,ListofSubsetDetailedMSAPAACStatistics)
  
  write.csv(ListofSubsetDetailedMSAAACStatistics, file = paste0(compostions,"MSAAACpsi3.csv"), row.names = F)
  write.csv(ListofSubsetDetailedMSADCStatistics, file = paste0(compostions,"MSADCpsi3.csv"), row.names = F)
  write.csv(ListofSubsetDetailedMSAPAACStatistics, file =  paste0(compostions,"MSAPAACpsi3.csv"), row.names = F)
}





