# needleman alignment of two para seqs
#library(Biostrings)
#library(seqinr)

NdlPos<-function(paraseq) {
  
  seqa<-paste0(paraseq[[1]],collapse = '')
  
  seqb<-paste0(paraseq[[2]],collapse = '')
  
  needle<-pairwiseAlignment(pattern=seqb,subject = seqa,type='global',
                            substitutionMatrix = 'BLOSUM62',
                            gapOpening=10,gapExtension=0.5, scoreOnly=F)
  
  alia<-unlist(strsplit(as.character(alignedSubject(needle)),split=''))
  alib<-unlist(strsplit(as.character(alignedPattern(needle)),split=''))
  
  needlealign<-c(alia,alib)
  
  ndl.pos<-c(c(1:length(alia)),c(1:length(alib)))
  
  ORF<-c(rep(names(paraseq)[1],times=length(alia)),rep(names(paraseq)[2],times=length(alib)))
  
  ndlali<-as.data.frame(cbind(ORF,needlealign,ndl.pos))
  
  ndlali$ori.pos<-NA
  
  ndlali$ori.pos[which(ndlali$ORF==names(paraseq)[1] & ndlali$needlealign!='-')]<-c(1:length(paraseq[[1]]))
  ndlali$ori.pos[which(ndlali$ORF==names(paraseq)[2] & ndlali$needlealign!='-')]<-c(1:length(paraseq[[2]]))
  
  
  ndlali$label<-as.character(ndlali$ori.pos)
  
  ndlali$label[is.na(ndlali$label)]<-'-'
  
  return(ndlali)
  
}
