# get motif alignments

paraSearchMotif<-function(paralogs,co_filep,paramotif,ndlali) {
  
  paralogs<-unlist(strsplit(paralogs,'_'))
  p1<-paralogs[1]
  p2<-paralogs[2]
  
  p1motif<-SearchMotif(p1,paramotif,'S288C',co_filep)
  p2motif<-SearchMotif(p2,paramotif,'S288C',co_filep)
  
  if(length(p1motif)>0) {
    p1motif$group<-paste0(p1,'_',p1motif$group)
  }
  
  if(length(p2motif)>0) {
    p2motif$group<-paste0(p2,'_',p2motif$group)
  }
  
  comb<-rbind(p1motif,p2motif)
  
  
  # add in needleman aligned position info in each paralog
  
  if (length(comb)>0) {
    
    comb.all<-comb[1,]
    comb.all<-comb.all[-1,]
    
    # create a NA temp for NA matches
    na.temp<-comb[1,]
    na.temp$group<-''
    na.temp[,2:1016]<-'-'
    
    
    comb.gene<-strsplit(comb$group,'_')
    comb.gene<-sapply(comb.gene,'[[',1)
    
    # mark whether this row has been moved up to match a previous paralog pos
    self.flag<-rep(F,times=nrow(comb))
    
    for (n in 1:nrow(comb)) {
      
      if (self.flag[n]==F) {
        
        ct.gene<-comb.gene[n]
        ct.pos<-as.numeric(comb$proteinposition[n])
        
        if(ct.gene==p1) {
          p.gene<-p2
        } else {
          p.gene<-p1
        }
        
        p.match<-ndlali$ndl.pos[which(ndlali$ORF==ct.gene & ndlali$ori.pos==ct.pos)]
        
        p.pos<-ndlali$ori.pos[which(ndlali$ORF==p.gene & ndlali$ndl.pos==p.match)]
        
        # the matched place could be NA for the paralog
        if (is.na(p.pos)) {
          
          na.temp$group<-paste0(p.gene,'_match_to_-')
          
          comb.all<-rbind(comb.all,comb[n,],na.temp)
          
        } else {
          
          # if the pos in paralog is already identified-move that row above
          # if not, add the row
          self.m<-which(comb.gene==p.gene & comb$proteinposition==p.pos)
          
          if (length(self.m)>0) {
            
            self.flag[self.m]<-T
            
            comb.all<-rbind(comb.all,comb[n,],comb[self.m,])
            
          } else {
            
            p.info<-SearchPosit(p.gene,as.character(p.pos),'1','S288C',co_filep)
            p.info$group<-paste0(p.gene,'_match_',ct.gene,'-',ct.pos)
            
            comb.all<-rbind(comb.all,comb[n,],p.info)
            
          }
          
        }
        
      }
      
    }
    return(comb.all)
  }
  
}
