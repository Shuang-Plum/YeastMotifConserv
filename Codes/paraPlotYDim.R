# return plotseq dimensions
paraPlotYDim<-function(ndlali) {
  
  naas<-50
  ystep=(-12)
  
  seqlen<-nrow(ndlali)/2
  nline<-ceiling(seqlen/naas)
  
  ymax<-(-nline*ystep)
  
  return(ymax)
  
}