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

PlotSeq<-function(refseq,refstr,sm.res) {
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
  if (refstr=='S288C') {refstr<-'STC'}
  
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
    n.ypos<-ypos-2

  } else {
    xpos<-rep(seq(from=0,to=(naas-1)*xstep,by=xstep),times=nline)
    ypos<-rep(seq(from=ymax,to=(ymax+(nline-1)*ystep),by=ystep),each=naas)
    xpos<-xpos[1:seqlen]
    ypos<-ypos[1:seqlen]
    n.ypos<-ypos-2
  }
  
  
  seqplot<-as.data.frame(cbind(xpos,ypos))
  seqplot$refseq<-refseq
  
  # plot the numbering
  numplot<-as.data.frame(cbind(xpos,n.ypos))
  numplot$num<-NumSeq(refseq)
  
  # get the match and highlight in yellow
  # get the col of refstr match character
  if (length(sm.res)>0) {
    m.label<-sm.res[,which(colnames(sm.res)==refstr)]
    sm.uni.row<-!duplicated(sm.res$MSAposition)
    m.label<-m.label[sm.uni.row]
    m.xpos<-xpos[sm.res$MSAposition[sm.uni.row]]
    m.ypos<-ypos[sm.res$MSAposition[sm.uni.row]]
    hlplot<-as.data.frame(cbind(m.xpos,m.ypos))
    hlplot$hl<-m.label
    
    ggplot(data=seqplot)+
      coord_cartesian(xlim=c(0,xmax),ylim=c(0,ymax))+
      theme(panel.background=element_rect(fill='white'),
            plot.margin = margin(2, 2, 2, 2, "pt"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            legend.title=element_blank(),
            legend.position='none',
            axis.ticks.y = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())+
      geom_text(aes(x=xpos,y=ypos,label=refseq,fontface='bold'),size=5)+
      geom_label(data=hlplot,aes(x=m.xpos,y=m.ypos,label=hl),fontface='bold',
                 size=5,fill='yellow',label.padding = unit(0.15, "lines"),
                 label.size = 0)+
      geom_text(data=numplot,aes(x=xpos,y=n.ypos,label=num),size=2.5)
  } else {
    
    ggplot(data=seqplot)+
      coord_cartesian(xlim=c(0,xmax),ylim=c(0,ymax))+
      theme(panel.background=element_rect(fill='white'),
            plot.margin = margin(2, 2, 2, 2, "pt"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            legend.title=element_blank(),
            legend.position='none',
            axis.ticks.y = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())+
      geom_text(aes(x=xpos,y=ypos,label=refseq,fontface='bold'),size=5)+
      geom_text(data=numplot,aes(x=xpos,y=n.ypos,label=num),size=2.5)
    
  }
  
  
  
  
}
