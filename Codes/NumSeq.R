# number the chr vector of protein seq 
# only number the aas, skip count the '-'

NumSeq<-function(refseq) {
  
  trulen<-sum(refseq!='-')
  
  numvec<-vector(mode='character',length=length(refseq))
  
  numvec[which(refseq=='-')]<-'-'
  
  numvec[which(numvec!='-')]<-as.character(c(1:trulen))
  
  return(numvec)
  
}
