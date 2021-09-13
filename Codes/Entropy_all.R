# entropy based algorithm

# aalist <- c("V","I","L","M",  "F","W","Y",  "S","T",
#             "N","Q",  "H","K","R",  "D","E",
#             "A","G",  "P",  "C",  "-","X")
# 
# aalist <- c("V","I","L","M",  "F","W","Y","S","T",
#             "N","Q",  "H","K","R",  "D","E",
#             "A","G",  "P",  "C")
# 
# ali<-sample(aalist,size=1012,replace=T)
# 
# ali<-rep('A',times=1012)
# 
# ali<-rep('W',times=1012)
# 
# ali<-sample(c('A','-','X'),size=1012,replace=T,prob=c(0.91,0.045,0.045))
# 
# ali<-sample(c('A','-','X'),size=1012,replace=T,prob=c(0.2,0.4,0.4))
# 
# ali<-c(sample(aalist,size=500,replace=T),sample(c('X','-'),size=512,replace=T))

### shannon entropy ###

# ali can be just a vec without strain name

# if the seq is random, not conserved, result is close to 0
# if conserved, results should be close to 1

# gap penalty decrease the score

# shanE with gap.pen is shanE*(nongap aas percentage)

shannon_E<-function (ali,gap.pen) {
  
  aaSE <- c("V","I","L","M",  "F","W","Y",  "S","T",
          "N","Q",  "H","K","R",  "D","E",
          "A","G",  "P",  "C",  "-","X")
  ali<-toupper(ali) # covert char to upper case
  
  ali[!ali %in% aaSE]<-'X' # convert non-standard aas to X
  
  temp<-table(ali)
  
  # remove gap and 'X' for calculating entropy
  temp<-temp[names(temp)!='-' & names(temp)!='X']
  
  if (length(temp)==1) {
    shanE<-0
    
  } else {
    
    ftemp<-temp/sum(temp)
    
    bmin<-min(length(ali),20)
    lntemp<-log(ftemp,base=bmin)# make sure the results is between 0-1
    
    shanE<-(-sum(ftemp*lntemp))
    
    
  }
  
  # revert the score to be consistent with other scores
  shanE<-(1-shanE)
  
  if (gap.pen) {
    # gap penalty decrease the score
    # shanE with gap.pen is shanE*(nongap aas percentage)
    gapaas<-sum(ali=='X')+sum(ali=='-')
    
    shanE<-shanE*(1-(gapaas/length(ali)))
  }
  
  return(shanE)
  
}


# stereochemically sensitive entropy
# using shannon entropy
# Williamson95 9 category

# if the seq is random, not conserved, results should be close to 0
# if conserved, results should be close to 1

# gap penalty decrease the score

# stereoE with gap.pen is stereoE*(nongap aas percentage)

stereochem_E<-function (ali,gap.pen) {
  
  aaStE <- c("V","I","L","M",  "F","W","Y",  "S","T",
          "N","Q",  "H","K","R",  "D","E",
          "A","G",  "P",  "C",  "-","X")
  ali<-toupper(ali) # covert char to upper case
  
  ali[!ali %in% aaStE]<-'X' # convert non-standard aas to X
  
  temp<-table(ali)
  
  # remove gap and 'X' for calculating entropy
  
  temp<-temp[names(temp)!='-' & names(temp)!='X']
  
  
  if (length(temp)==1) {
    stereoE<-0
    
  } else {
    
    # convert table of aas into category table
    # Williamson '95
    #['V','L', 'I','M'], ['F','W','Y'], ['S','T'], ['N','Q'], 
    #['H','K','R'], ['D','E'], ['A','G'], ['P'], ['C']
    
    stemp<-vector(mode='numeric',length=9)
    
    # sum, na.rm can make aas not exists into 0 instead of NA
    stemp[1]<-sum(temp['V'],temp['L'],temp['I'],temp['M'],na.rm=T)
    stemp[2]<-sum(temp['F'],temp['W'],temp['Y'],na.rm=T)
    stemp[3]<-sum(temp['S'],temp['T'],na.rm=T)
    stemp[4]<-sum(temp['N'],temp['Q'],na.rm=T)
    stemp[5]<-sum(temp['H'],temp['K'],temp['R'],na.rm=T)
    stemp[6]<-sum(temp['D'],temp['E'],na.rm=T)
    stemp[7]<-sum(temp['A'],temp['G'],na.rm=T)
    stemp[8]<-sum(temp['P'],na.rm=T)
    stemp[9]<-sum(temp['C'],na.rm=T)
    
    stemp<-stemp[stemp!=0]
    
    fstemp<-stemp/sum(stemp)
    
    bmin<-min(length(ali),9)
    lnstemp<-log(fstemp,base=bmin)# make sure the results is between 0-1
    
    stereoE<-(-sum(fstemp*lnstemp))
    
    
  }
  # revert the score to be consistent with other scores
  stereoE<-(1-stereoE)
  
  if (gap.pen) {
    # gap penalty is decrease the score
    # stereoE with gap.pen is stereoE*(nongap aas percentage)
    gapaas<-sum(ali=='X')+sum(ali=='-')
    
    stereoE<-stereoE*(1-(gapaas/length(ali)))
  }
  
  return(stereoE)
  
}  


# if ali has the same freq as background dist -- totally random
# jsd returns 0
# if ali is highly preserved (only has 1 aas)
# jsd returns ~1

# gap.penalty decrease jsd -- jsd = jsd*(nongap aas percentage)

# jsd won't return 1 if all conserved as 'A' or any other aas
# cause it is comparing to the bkg.dist
# bkg.dist A.freq =0.078 W.freq=0.014
# all A would give jsd ~0.8, all W would give jsd~0.95
# the score is entropy but considering aas background frequency


#bkg.dist<-read.csv('BLOSUM62_dist.csv',header=T)

js_divergence<-function(ali,bkg.dist,gap.pen) {
  
  # get freq of aas
  aa<-c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T',
        'W','Y','V','B','Z','X','-')
  
  ali<-toupper(ali) # covert char to upper case
  
  ali[!ali %in% aa]<-'X' # convert non-standard aas to X
  
  temp<-table(ali)
  
  # remove non-common aas and gap
  temp<-temp[names(temp)!='B' & names(temp)!='Z' & names(temp)!='X' & names(temp)!='-']
  
  # get frequency
  ftemp<-temp/sum(temp)
  
  # make the missing aas freq as 0
  missa<-setdiff(names(bkg.dist),names(ftemp))
  
  for (i in missa) {ftemp[i]<-0}
  
  # sort ftemp as bkg.dist order
  ftemp<-ftemp[names(bkg.dist)]
  
  ftemp<-as.vector(ftemp)
  names(ftemp)<-names(bkg.dist)
  
  # make r distribution
  r.dist<-0.5*ftemp+0.5*bkg.dist
  
  jsd<-0
  
  
  
  for (n in 1:length(r.dist)) {
    # calculate distance with non-zero aas
    if (r.dist[n]!=0) {
      # avoid having 0*log2(0) which is NAN
      if (ftemp[n]==0) {
        jsd<-jsd+bkg.dist[n] * log2(bkg.dist[n]/r.dist[n])
      } else if (bkg.dist[n]==0) {
        jsd<-jsd+ftemp[n] * log2(ftemp[n]/r.dist[n])
      } else {
        jsd<-jsd+ftemp[n] * log2(ftemp[n]/r.dist[n])+bkg.dist[n] * log2(bkg.dist[n]/r.dist[n])
      }
      
    }
    
  }
  
  jsd<-jsd/2
  
  if (gap.pen) {
    # gap.penalty decrease jsd -- jsd = jsd*(nongap aas percentage)
    gapaas<-sum(ali=='X')+sum(ali=='-')+sum(ali=='B')+sum(ali=='Z')
    
    jsd<-jsd*(1-(gapaas/length(ali)))
  }
  
  jsd<-as.numeric(jsd)
  
  return(jsd)
  
}
