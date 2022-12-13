# recreate Mihalek04 pholygeny method
# input is the column of match 
# and the tree of phylogeny build by UPGMA method
# and cutree to 1-1011 nodes (1012 strains in total)

# aalist <- c("V","I","L","M",  "F","W","Y",  "S","T",
#         "N","Q",  "H","K","R",  "D","E",
#         "A","G",  "P",  "C",  "-","X")
# 
# aalist <- c("V","I","L","M",  "F","W","Y","S","T",
#         "N","Q",  "H","K","R",  "D","E",
#         "A","G",  "P",  "C")
# 
# ali<-sample(aalist,size=1012,replace=T)
# 
# ali<-rep('A',times=1012)
# 
# ali<-sample(c('A','-','X'),size=1012,replace=T,prob=c(0.05,0.5,0.45))
# 
# names(ali)<-row.names(tree)


gentropy<-function(gvec,gap.pen) {
  
  
  aa <- c("V","I","L","M",  "F","W","Y",  "S","T",
          "N","Q",  "H","K","R",  "D","E",
          "A","G",  "P",  "C",  "-","X")
  
  gvec[!gvec %in% aa]<-'X' # convert non-standard aas to X
  
  temp<-table(gvec)
  
  # remove gap and 'X' for calculating entropy
  temp<-temp[names(temp)!='-' & names(temp)!='X']
  
  # if only has 1 member, do not count in
  # if temp is empty after removing X and -, don't count in
  # (length(temp)==1 & temp[1]==1) | length(temp)==0
  if (length(temp)>1 | (length(temp)==1 & temp[1]>1 )) {
    
    ftemp<-temp/sum(temp)
    
    bmin<-min(length(gvec),20)
    lntemp<-log(ftemp,base=bmin)# make sure the results is between 0-1
    
    ge<-(-sum(ftemp*lntemp))
    # totally random =1  highly conserve =0
    
  } else {ge = 1} # only has 1 member or empty after removing gaps
  
  # reverse the score to be consistent with others
  ge = (1-ge)
  # totally random =0  highly conserve =1
  
  if (gap.pen) {
    # gap penalty decrease the score
    # ge = ge*(nongap aas percentage)
    gapaas<-sum(gvec=='X')+sum(gvec=='-')
    
    ge=ge*(1-(gapaas/length(gvec)))
    
  }
  
  return(ge)
  
}

# ali must be a named vec with the same row names as tree

# if the seq is random, not conserved, results should be close to 0
# if conserved, results should be close to 1

# group score, ge is reverted 
# gap penalty decrease the score, ge = ge*(nongap aas percentage)
# cause only after reversion, applying penalty this way makes sense

#tree<-read.csv('G:/My Drive/UNC/SeqMoid/Yeast1011/UPGMA_nodes.csv',row.names = 'X')

phylo_M<-function(ali,tree,gap.pen) {
  
  mscore<-0
  
  # check if ali has the same num of strains as in tree
  if (length(ali)==nrow(tree) & nrow(tree)==(ncol(tree)+1) & 
      length(names(ali))==nrow(tree)) {
    
    torder<-row.names(tree)
    
    ali<-ali[torder] # order the ali as the strains order in tree
    scount<-0
    for (n in 1:ncol(tree)) {
      
      # get the tree cut info
      treecut<-tree[,n]
      # n is also the same as subtrees
      gscore<-0
      gcount<-0 # count node that has member >1 and only normalized by this number
      for (g in 1:n) {
        gvec<-ali[treecut==g] # make sure get the right aas for each strain
        gvec<-toupper(gvec) # convert all char to upper case
        
        gscore<-gscore+gentropy(gvec,gap.pen)
        
        # if gap.pen applied
        # count all groups that has more member than the refstrain
        # even if it has refstr and all gaps -- ge=0
        # because this could be a ge=1
        # this helps gap penalty to work
        # if gap.pen is not turn on
        # only count groups that has a solid aas other than refstrain
        # for group with refstrain and all gaps would be treated as only has refstrain (1 member)
        # which should not be counted in
        if (gap.pen) {
          if (length(gvec)>1) {gcount<-gcount+1}
        } else {
          gtemp<-gvec[gvec!='X' & gvec!='-']
          if (length(gtemp)>1) {gcount<-gcount+1}
        }
        
      }
      
      # only normalize when gscore >0
      # when gscore=0, gcount could also be 0, then normalization does not make sense
      if(gscore>0) {
        gscore<-gscore/gcount # normalize gscore by the num of groups with memb>1
        # not the num of groups, cause group with 1 mem is not contributing to entropy
        
        # also count in this n node for normalization
        # it could be a highly conserved Aas with a lot of gaps (too sparse)
        # would result in a n node with 0 score
        scount<-scount+1
      }
      
      
      mscore<-mscore+gscore
      
      
    }
    
    mscore<-mscore/scount# normalized the results between 0-1
    return(mscore)
    
  } else {
    warning('num of matched strain is different from phylogeny tree 
            or tree is not 1:(n-1) 
            or input is not named with strains, 
            score not calculated')
    return(NA)
  }
  
}

