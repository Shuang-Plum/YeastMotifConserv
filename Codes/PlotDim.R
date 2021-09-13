# return plotseq dimensions
PlotDim<-function(refseq) {
  
  naas<-50
  
  xstep=3
  ystep=(-5)
  
  seqlen<-length(refseq)
  nline<-ceiling(seqlen/naas)
  
  ymax<-(-nline*ystep)
  
  return(ymax)
  
}