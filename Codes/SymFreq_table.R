# generate symbol frequency table for shiny rendering
# takes in results from searchMotif
# need to source SymbolFreq.R

source('./Codes/SymbolFreq.R')
Symfreq_table<-function(refstr,sm.res) {
  
  
  if (length(sm.res)>0) {
    
    # get the properties
    #sm.res<-sm.res.list$sm.res
    ref.col<-which(colnames(sm.res)==refstr)
    symfreq<-sm.res[,c(1:2,ref.col,3)]
    symfreq$group<-paste0('Match_',symfreq$group)
    
    # add cols of each aas 
    aa<-c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T',
          'W','Y','V','B','Z','X','-')
    aas<-matrix(data=0,nrow=nrow(symfreq),ncol=length(aa))
    colnames(aas)<-aa
    symfreq<-cbind(symfreq,aas)
    rm(aas)
    
    # get the alignments
    ali.df<-sm.res[,4:(ncol(sm.res)-2)]
    
    # for each row (aligned position) calculate symbol frequency
    for (n in 1:nrow(symfreq)){
      
      symfreq[n,5:ncol(symfreq)]<-SymbolFreq(ali.df[n,])
      
    }
    
    # remove cols (aas) that does not appear at least once
    del.col<-vector(mode='numeric')
    for (i in 5:ncol(symfreq)) {
      if (sum(symfreq[,i])==0) {del.col<-c(del.col,i)}
    }
    
    symfreq[,del.col]<-NULL
    
    # add xpos and ypos for interactive clicking
    symfreq$xpos<-sm.res$xpos
    symfreq$ypos<-sm.res$ypos
    
    return(symfreq)
    
  }
  
}
