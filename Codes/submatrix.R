# substitution matrix

# BLOSUM62 after Karlin normalization, min =-1 max=1
# BLOSUM62 measures the logarithm for the ratio of the likelihood of two amino acids 
# appearing with a biological sense 
# and the likelihood of the same amino acids appearing by chance
# positiv score - likely substitution - conserved or makes biological sense
# negative score - less likely substitution - random or relaxed 

# normalized to 0-1 where conserved =1 random=0 to be consistent 
# all gap would lead to score =0.5
# score <0.5 meaning have substitute that has neg score

# gap.pen is inhereted in the submatrix 
# gap is penalized to score =0
# cannot remove gap and only focus on non-gap part
# gap score =0 but some aas substitute score is negative
# for some ali with neg aas sub scores, counting in gap > ignoring all gaps
# which does not make sense


#setwd("G:/My Drive/UNC/SeqMoid")
#submat<-read.csv('BLOSUM62_KarlinNorm.csv',header=T,row.names = 'AAS')

sub_matrix<-function(ali,submat) {
  
  aasub<-c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T',
               'W','Y','V','B','Z','X','-')
  
  colnames(submat)[24]="-"
  
  ali<-toupper(ali) # covert char to upper case
  
  ali[!ali %in% aasub]<-'X' # convert non-standard aas to X
  # this include '*' stop codon, may be improved to penalized more
  
  freq.aa   <- table(ali)
  
  unique.aa <- names(freq.aa)
  missing.aa <- unique.aa[!unique.aa %in% colnames(submat)]
  
  count <- 0
  score <- 0
  
  for(i in 1:length(unique.aa)) {
    aa.i <- unique.aa[i]
    freq.i <- freq.aa[i]
    for(j in i:length(unique.aa)) {
      aa.j <- unique.aa[j]
      freq.j <- freq.aa[j]
      ##sim <- mat[aa.i,aa.j]
      if(length(missing.aa)>0) {
        if(i==missing.aa || j==missing.aa) {
          sim <- 0
        } else {
          sim <- submat[aa.i,aa.j]
        }
      } else { sim <- submat[aa.i,aa.j] }
      
      ## number of comparisons
      if(aa.i == aa.j) {
        ncmp <- freq.i * (freq.i - 1)/2
      } else {
        ncmp <- freq.i * freq.j
      }
      count <- count + ncmp
      score <- score + (ncmp * sim)
    }
  }
  
  fscore<-score/count
  
  # normalize to 0-1 to be consistent with other scores
  smin<-min(submat) # the possible min of fscore 
  smax<-max(submat) # the possible max of fscore
  fscore<-(fscore-smin)/(smax-smin)
  
  fscore<-as.numeric(fscore)
  
  return(fscore)
  
}

