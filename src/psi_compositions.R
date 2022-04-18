blastpSeq<- function(seq, start.pos = 1L, end.pos = nchar(seq), 
                     silent = TRUE, 
                     evalue = 10L, output.path=resultspath){
  
  if (is.null(output.path)) stop('Must specify the output path')
  
  N = end.pos - start.pos + 1L
  
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
           ' -comp_based_stats 1 -db ', shQuote(tcdb),
           ' -query ', shQuote(seq),  ' -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs\"  -num_alignments 1',' -out ', 
           file.path(output.path,"out.txt")))
  
  print("******************************")
  print(cmdblastp)
  if (silent == TRUE) system(cmdblastp, ignore.stdout = F) else system(cmdblastp)      
  #get the hit sequences Id
  data = read.table(file.path(output.path,"out.txt"))
  
}


psiblastSeq <- function(seq,
                        start.pos = 1L,
                        end.pos = nchar(seq),
                        silent = TRUE,
                        evalue = 10L,
                        numOfIterations = 3,
                        output.path = path) {

  if (is.null(output.path))
    stop('Must specify the output path')
  
  N = end.pos - start.pos + 1L
  
  
  # Basic parameters for Blastp/psiblast
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
  # Run Psiblast
  cmdpsiblast = paste(
    paste0(
      shQuote(psiblast.path),
      ' -db ',
      shQuote(swissprotdb),
      ' -query ',
      shQuote(queryFasta),
      ' -num_iterations ',
      numOfIterations ,
      ' -evalue ',
      evalue,
      " -out psiblastDsbAOut -out_pssm PSSMDsbA -out_ascii_pssm ascii_mtx_file"
    )
  )
  
  print("******************************")
  print(cmdpsiblast)
  if (silent == TRUE)
    system(cmdpsiblast, ignore.stdout = TRUE)
  else
    system(cmdpsiblast)
  
  
}
ClaculateCompositions_return<- function(subdirName, filename) { 
  
  con= paste0(subdirName,filename)
  print(con)
  se<- read.fasta(con,seqtype = "AA",as.string=T)
  
  seqlist<- unlist(lapply(se,as.character))
  seqlist<-sapply(seqlist,gsub,pattern="[^GPAVLIMCFYWHKRQNEDST]",replacement="") # getting rid of gaps? not the smartest thing to do, need more work
  names(seqlist)<- sub(".+\\:","",names(se))
  names(seqlist)<- sub("\\|.*","",sub(".+?\\|","", names(seqlist)))
  seqlist=seqlist[which(nchar(seqlist)>30)]
  
  MSAAAC<- lapply(seqlist, extractAAC)
  outputAAC <- matrix(unlist(MSAAAC), ncol = 20, byrow = TRUE)
  rownames(outputAAC)<- names(seqlist)
  write.csv(outputAAC, file.path(subdirName, "AAC3.csv"))
  MSAPseAAC<- lapply(seqlist, extractPAAC)
  
  outputPseAAC <- matrix(unlist(MSAPseAAC), ncol = 50, byrow = TRUE)
  rownames(outputPseAAC)<- names(seqlist)
  write.csv(outputAAC, file.path(subdirName,"PseAAC3.csv"))
  
  MSADC<- lapply(seqlist, extractDC)
  outputDC <- matrix(unlist(MSADC), ncol = 400, byrow = TRUE)
  rownames(outputDC)<- names(seqlist)
  write.csv( outputDC, file.path(subdirName,"DC3.csv"))
  return(list(AAC=outputAAC, DC=outputDC, PseAAC=outputPseAAC))
  
}

PreparedataforPSI = function(seq, start.pos = 1L, end.pos = nchar(seq), 
                                iter = 5, silent = TRUE, 
                                evalue = 10L, output.path=resultspath) {
  #1- run blastp on datafiles
  seqname=names(seq)
  SeqDirectory=file.path(output.path,names(seq))
  
  dir.create(SeqDirectory, showWarnings = FALSE, recursive = FALSE, mode = "0777")   
  setwd(SeqDirectory)
  
  psiblastSeq(
    seq,
    start.pos ,
    end.pos ,
    silent ,
    evalue = 0.001,
    numOfIterations = 3,
    SeqDirectory
  )
}

converttofasta<- function(seq, output.path) {
  dirName <- file.path(output.path,names(seq))
  
  for (iteration in 1:3) {
    name = paste0("psi", iteration, ".fasta")
    command = paste0(
      mview.path,
      " -in blast -cycle ",
      iteration,
      " ",
      file.path(dirName, "psiblastDsbAOut"),
      " -out fasta> ",
      file.path(dirName, name)
    )
    print(command)
    system(command)
  }
}

ClaculateCompositions<- function(subdirName, filename) { 
  
  con= file.path(subdirName,filename)
  print(con)
  se<- read.fasta(con,seqtype = "AA",as.string=T)
  
  seqlist<- unlist(lapply(se,as.character))
  seqlist<-sapply(seqlist,gsub,pattern="[^GPAVLIMCFYWHKRQNEDST]",replacement="") # getting rid of gaps? not the smartest thing to do, need more work
  names(seqlist)<- sub(".+\\:","",names(se))
  names(seqlist)<- sub("\\|.*","",sub(".+?\\|","", names(seqlist)))
  seqlist=seqlist[which(nchar(seqlist)>30)]
  if(length(seqlist)==0)
    return(NA)
  
  MSAAAC<- lapply(seqlist, extractAAC)
  outputAAC <- matrix(unlist(MSAAAC), ncol = 20, byrow = TRUE)
  rownames(outputAAC)<- names(seqlist)
  write.csv(outputAAC, file.path(subdirName, "AAC3.csv"))
  
  MSAPseAAC<- lapply(seqlist, extractPAAC)
  outputPseAAC <- matrix(unlist(MSAPseAAC), ncol = 50, byrow = TRUE)
  rownames(outputPseAAC)<- names(seqlist)
  write.csv(outputPseAAC, file.path(subdirName, "PseAAC3.csv"))
  
  MSADC<- lapply(seqlist, extractDC)
  outputDC <- matrix(unlist(MSADC), ncol = 400, byrow = TRUE)
  rownames(outputDC)<- names(seqlist)
  write.csv( outputDC, file.path(subdirName, "DC3.csv"))
  
}

ClaculateEachCompositions<- function(seq, output.path) {
  subdirName= file.path(intermediateFiles,names(seq))
  
  if(readLines(file.path(subdirName,"psi3.fasta"),n=1)!="mview: no alignments found") {
    ClaculateCompositions(subdirName,"psi3.fasta")
  }else if (readLines(file.path(subdirName,"psi2.fasta"),n=1)!="mview: no alignments found"){
    print("no psi3.fasta but Pse2")
    ClaculateCompositions(subdirName,"psi2.fasta")
  }else if (readLines(file.path(subdirName,"psi1.fasta"),n=1)!="mview: no alignments found"){
    print("no psi2.fasta but Pse1")
    ClaculateCompositions(subdirName,"psi1.fasta")
  } else {
    print("nomatch")
  }
  
  closeAllConnections()
  
}

psi_compositions<- function(fastafile) {
  seqs<- readFASTA(fastafile)
  names(seqs)<- sub("\\|.*","",sub(".+?\\|","", names(seqs)))
  for(j in c(1:length(seqs))) {
    x<- seqs[j]
    PreparedataforPSI(seq=x,output.path=intermediateFiles, evalue=.001)
    converttofasta(seq= x,intermediateFiles) #write function rather than using Mview # todo
    ClaculateEachCompositions(seq= x,intermediateFiles)
  }
  
  AAfiles<-  names(seqs)
  ListofSubsetDetailedMSAAACStatistics<- matrix(ncol = 20+1,nrow = length(AAfiles)) #contrains AAC for each seq in the subset
  ListofSubsetDetailedMSADCStatistics <- matrix(ncol=400+1,nrow=length(AAfiles))
  ListofSubsetDetailedMSAPAACStatistics <- matrix(ncol=(50+1),nrow=length(AAfiles))
  
  for(i in c(1:length(AAfiles))) {
    subdirName<- file.path(intermediateFiles,AAfiles[i])
    
    number=3
    if(file.exists(paste0(subdirName,"/AAC",number,".csv"))) {
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
      
    } else { }
  }
  
  ListofSubsetDetailedMSAAACStatistics<-cbind(UniprotID=AAfiles,ListofSubsetDetailedMSAAACStatistics)
  ListofSubsetDetailedMSADCStatistics<- cbind(UniprotID=AAfiles,ListofSubsetDetailedMSADCStatistics)
  ListofSubsetDetailedMSAPAACStatistics<-cbind(UniprotID=AAfiles,ListofSubsetDetailedMSAPAACStatistics)
  
  write.csv(ListofSubsetDetailedMSAAACStatistics, file = file.path(compositions,"MSAAACpsi3.csv"), row.names = F)
  write.csv(ListofSubsetDetailedMSADCStatistics, file = file.path(compositions,"MSADCpsi3.csv"), row.names = F)
  write.csv(ListofSubsetDetailedMSAPAACStatistics, file = file.path(compositions,"MSAPAACpsi3.csv"), row.names = F)
}
