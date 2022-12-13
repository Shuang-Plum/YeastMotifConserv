# return original aas seq without MSA alignment gap -
#co_filep<-'E:/Yeast1011Data/datasets/CO_aligned'
#paralogs<-'YAL007C_YOR016C'

paraGetSeq<-function(paralogs,co_filep){
  
  paralogs<-unlist(strsplit(paralogs,'_'))
  p1<-paralogs[1]
  p2<-paralogs[2]
  
  co1<-read.fasta(file=paste0(co_filep,'/',p1,'_aligned.fasta'),seqtype = 'AA',
                  as.string = F,strip.desc = T)
  
  co2<-read.fasta(file=paste0(co_filep,'/',p2,'_aligned.fasta'),seqtype = 'AA',
                  as.string = F,strip.desc = T)
  
  snames1<-getName(co1)
  snames2<-getName(co2)
  
  p1seq<-getSequence(co1[pmatch('S288C',snames1)],as.string = F)
  p2seq<-getSequence(co2[pmatch('S288C',snames2)],as.string = F)
  
  p1seq<-unlist(p1seq)
  p2seq<-unlist(p2seq)
  
  # remove the - in MSA alignment to get original aas seq
  p1seq<-p1seq[which(p1seq!='-')]
  p2seq<-p2seq[which(p2seq!='-')]
  
  paraseq<-list(p1seq,p2seq)
  
  names(paraseq)<-c(p1,p2)
  
  return(paraseq)
  
  
}
