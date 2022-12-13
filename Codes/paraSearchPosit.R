# get motif alignments


paraSearchPosit<-function(paralogs,co_filep,paraposit,paramlen,ndlali) {
  
  paralogs<-unlist(strsplit(paralogs,'_'))
  p1<-paralogs[1]
  p2<-paralogs[2]
  
  # only search for positions in first gene
  # and find corresponding positions with needleman alignment in the second gene
  p1posit<-SearchPosit(p1,paraposit,paramlen,'S288C',co_filep)
  
  
  if(length(p1posit)>0) {
    
    p1posit$group<-paste0(p1,'_pos_',p1posit$group)
    
    comb<-p1posit[1,]
    comb<-comb[-1,]
    
    # creat NA template for matched to - in needleman
    na.temp<-p1posit[1,]
    na.temp$group<-''
    na.temp[,2:1016]<-'-'
    
    # find correspoding match after needle alignment in second gene
    for (n in 1:nrow(p1posit)) {
      
      ct.pos<-p1posit$proteinposition[n]
      
      p.match<-ndlali$ndl.pos[which(ndlali$ORF==p1 & ndlali$ori.pos==ct.pos)]
      
      p.pos<-ndlali$ori.pos[which(ndlali$ORF==p2 & ndlali$ndl.pos==p.match)]
      
      # the matched place could be NA for the paralog
      if (is.na(p.pos)) {
        
        na.temp$group<-paste0(p2,'_match_to_-')
        comb<-rbind(comb,p1posit[n,],na.temp)
        
      } else {
        
        p.info<-SearchPosit(p2,as.character(p.pos),'1','S288C',co_filep)
        p.info$group<-paste0(p2,'_match_',p1,'-',ct.pos)
        
        comb<-rbind(comb,p1posit[n,],p.info)
        
      }
      
    }
    
    return(comb)
    
  }

  
}
