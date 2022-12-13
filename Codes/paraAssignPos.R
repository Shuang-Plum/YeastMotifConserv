# assign position for each aas in seqs
paraAssignPos<-function(ndlali) {
  
  #ndlali<-as.data.frame(ndlali)
  
  naas<-50
  
  xstep=5
  ystep=(-12)
  yspace=3
  
  seqlen<-nrow(ndlali)/2
  nline<-ceiling(seqlen/naas)
  
  xmax<-xstep*(naas+1)
  ymax<-(-nline*ystep)
  
  if (nline<2) {
    xpos<-seq(from=0,to=((seqlen-1)*xstep),by=xstep)
    ypos.1<-rep(11,times=seqlen)
    ypos.2<-rep((11-yspace),times=seqlen)
    n1.ypos<-ypos.1+2
    n2.ypos<-ypos.2-2
    
  } else {
    xpos<-rep(seq(from=0,to=(naas-1)*xstep,by=xstep),times=nline)
    ypos.1<-rep(seq(from=ymax,to=(ymax+(nline-1)*ystep),by=ystep),each=naas)
    xpos<-xpos[1:seqlen]
    ypos.1<-ypos.1[1:seqlen]
    ypos.2<-ypos.1-yspace
    n1.ypos<-ypos.1+2
    n2.ypos<-ypos.2-2
  }
  
  
  ypos<-c(ypos.1,ypos.2)
  n.ypos<-c(n1.ypos,n2.ypos)
  
  ndlali$xpos<-c(xpos,xpos)
  ndlali$ypos<-ypos
  ndlali$n.ypos<-n.ypos
  
  return(ndlali)

    
}
