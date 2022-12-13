# return symbol freq
# aalist<-c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T',
#           'W','Y','V','B','Z','X','-')
# 
# ali<-sample(aalist,size=1012,replace=T)
# 
# ali<-rep('A',times=1012)
# 
# ali<-sample(c('A','-','X'),size=1012,replace=T,prob=c(0.5,0.2,0.3))

SymbolFreq<-function(ali) {
  
  aa<-c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T',
        'W','Y','V','B','Z','X','-')
  
  ali<-as.character(ali)
  temp<-table(ali)
  
  missaas<-setdiff(aa,names(temp))
  
  for (a in missaas){temp[a]<-0}
  
  # order temp in the same order for output
  temp<-temp[aa]
  
  temp<-round(temp/length(ali),digits=4)
  
  
  return(temp)
}
