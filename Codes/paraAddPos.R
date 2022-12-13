# add xpos and ypos positions to sm.res for yellow highlighted regions
paraAddPos<-function(ndlali, sm.res){
  
  if (length(sm.res)>0) {
    sm.res$xpos<-NA
    sm.res$ypos<-NA
    
    colnames(sm.res)[which(colnames(sm.res)=='STC')]<-'S288C'
    
    for (n in 1:nrow(sm.res)) {
      
      pt<-sm.res$proteinposition[n]
      # could match to '-' in needleman
      if (pt!='-') {
        gt<-unlist(strsplit(sm.res$group[n],'_'))[1]
        
        mt<-which(ndlali$ORF==gt & ndlali$ori.pos==pt)
        
        sm.res$xpos[n]<-ndlali$xpos[mt]
        sm.res$ypos[n]<-ndlali$ypos[mt]
      } else {
        # when matched to '-' refer xpos and ypos from above line
        sm.res$xpos[n]<-sm.res$xpos[n-1]
        # get the above line protein and ndl.pos
        # get the ypos from ndlali of the other gene same ndl.pos
        ag<-unlist(strsplit(sm.res$group[n-1],'_'))[1]
        a.pos<-as.numeric(sm.res$proteinposition[n-1])
        b.npos<-ndlali$ndl.pos[which(ndlali$ORF==ag & ndlali$ori.pos==a.pos)]
        bg<-setdiff(unique(ndlali$ORF),ag)
        sm.res$ypos[n]<-ndlali$ypos[which(ndlali$ORF==bg & ndlali$ndl.pos==b.npos)]
        
      }
    }
    
    return(sm.res)
    
  }
  
}
