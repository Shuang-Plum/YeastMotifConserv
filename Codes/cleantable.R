# clean xpos and ypos from results
cleantable<-function(res) {
  
  res$xpos<-NULL
  res$ypos<-NULL
  
  return(res)
}