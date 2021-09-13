# get sequence
GetSeq<-function(gene,refstr,co_filep){
  
  co<-read.fasta(file=paste0(co_filep,'/',gene,'_aligned.fasta'),seqtype = 'AA',
                 as.string = F,strip.desc = T,)
  
  snames<-getName(co)
  snames<-strsplit(snames,split=paste0('_',gene,'_'),fixed=T)
  # get the strain name
  snames<-sapply(snames,'[[',1)
  # replace S288C to STC to be consistent
  snames<-gsub('S288C','STC',snames)
  
  refseq<-getSequence(co[which(snames==refstr)],as.string = F)
  
  refseq<-unlist(refseq)
  
  return(refseq)
  
  
}