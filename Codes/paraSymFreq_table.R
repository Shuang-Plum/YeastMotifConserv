# generate symbol frequency table for shiny rendering
# takes in results from searchMotif
# need to source SymbolFreq.R

source('Codes/SymbolFreq.R')
paraSymfreq_table<-function(sm.res) {
  
  if (length(sm.res)>0) {
    
    # get the properties
    #sm.res<-sm.res.list$sm.res
    
    ref.col<-which(colnames(sm.res)=='S288C')
    symfreq<-sm.res[,c(1:3,ref.col,4)]
    
    
    # add cols of each aas 
    aa<-c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T',
          'W','Y','V','B','Z','X','-')
    aas<-matrix(data=0,nrow=nrow(symfreq),ncol=length(aa))
    colnames(aas)<-aa
    
    symfreq<-cbind(symfreq,aas)
    rm(aas)
    
    # get the alignments
    ali.df<-sm.res[,5:(ncol(sm.res)-2)]
    
    # for each row (aligned position) calculate symbol frequency
    for (n in 1:nrow(symfreq)){
      
      if(symfreq$proteinposition[n]!='-') {
        
        symfreq[n,6:ncol(symfreq)]<-SymbolFreq(ali.df[n,])
      } else {
        
        symfreq[n,6:ncol(symfreq)]<-NA
      }
    }
    
    # remove cols (aas) that does not appear at least once
    del.col<-vector(mode='numeric')
    for (i in 6:ncol(symfreq)) {
      if (sum(symfreq[,i],na.rm=T)==0) {del.col<-c(del.col,i)}
    }
    
    symfreq[,del.col]<-NULL
    
    # add xpos and ypos for interactive clicking
    symfreq$xpos<-sm.res$xpos
    symfreq$ypos<-sm.res$ypos
    
    # parse more info
    ORF<-strsplit(symfreq$group,'_')
    ORF<-sapply(ORF,'[[',1)
    
    paralog.pair<-rep(c(1:(nrow(symfreq)/2)),each=2)
    
    symfreq<-as.data.frame(cbind(paralog.pair,ORF,symfreq))
    
    #colnames(symfreq)[3]<-'detail'
    
    symfreq$group<-NULL
    
    return(symfreq)
    
  }
  
}
