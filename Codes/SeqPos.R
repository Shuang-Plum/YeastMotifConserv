# plot refstrain seq for the gene with match motif highlighted

# take in searchMotif results
# generate table with xpos and ypos for each matched position
# library(seqinr)
# library(ggplot2)
# co_filep<-'E:/SeqMoid/datasets/CO_aligned'

# gene<-'YOR267C'
# 
# motif<-'KK'
# 
# refstr<-'STC'

SeqPos<-function(refseq,sm.res) {
  # co<-read.fasta(file=paste0(co_filep,'/',gene,'_aligned.fasta'),seqtype = 'AA',
  #                as.string = F,strip.desc = T,)
  # 
  # snames<-getName(co)
  # snames<-strsplit(snames,split=paste0('_',gene,'_'),fixed=T)
  # # get the strain name
  # snames<-sapply(snames,'[[',1)
  # # replace S288C to STC to be consistent
  # snames<-gsub('S288C','STC',snames)
  # 
  # refseq<-getSequence(co[which(snames==refstr)],as.string = F)
  # 
  # refseq<-unlist(refseq)
  
  # split refseq into 30 aas a line
  #refseq.spl<-strsplit(refseq, "(?<=.{30})", perl = TRUE)[[1]]
  
  # define x and y pos for each aas in refseq
  if (length(sm.res)>0) {

    naas<-50
    
    xstep=3
    ystep=(-5)
    
    seqlen<-length(refseq)
    nline<-ceiling(seqlen/naas)
    
    xmax<-xstep*(naas+1)
    ymax<-(-nline*ystep)
    
    if (nline<2) {
      xpos<-seq(from=0,to=(seqlen*xstep),by=xstep)
      ypos<-rep(3,times=seqlen)
      #n.ypos<-ypos-2
      
    } else {
      xpos<-rep(seq(from=0,to=(naas-1)*xstep,by=xstep),times=nline)
      ypos<-rep(seq(from=ymax,to=(ymax+(nline-1)*ystep),by=ystep),each=naas)
      xpos<-xpos[1:seqlen]
      ypos<-ypos[1:seqlen]
      #n.ypos<-ypos-2
    }
    
    # get the xpos and ypos for sm.res
    sm.res$xpos<-xpos[sm.res$position]
    sm.res$ypos<-ypos[sm.res$position]
    
    
    return(sm.res)
  }
  
  
}
  
