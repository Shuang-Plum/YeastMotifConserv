# calculate scores with different algos
# take results from searchMotif

source('Codes/submatrix.R')
source('Codes/Entropy_all.R')
source('Codes/phylo_ZOOM.R')
source('Codes/phylo_M.R')

# phyloFull to indicate whether to calculate full phylo score or just use ZOOM
# full phylo score takes significantly long time ~3-5min per site
ProConserve<-function(sm.res, refstr, gap.pen=T, phyloFull=F) {
  
  if (refstr=='S288C') {refstr<-'STC'}
  
  if (length(sm.res)>0) {
    # load predefined data
    # using BLOSUM62 dist and submat - can change to PAM30 if needed
    # BLOSUM62 submat is Karlin normalized to -1 to 1
    # using UPGMA tree built from SNPs distance of 1012 strains
    bkg.dist<-read.csv('Data/BLOSUM62_dist.csv',header=T,stringsAsFactors=F)
    submat<-read.csv('Data/BLOSUM62_KarlinNorm.csv',header=T,row.names = 'AAS',stringsAsFactors=F)
    tree<-read.csv('Data/UPGMA_nodes.csv',row.names = 'X',stringsAsFactors=F)
    
    # get the properties
    #sm.res<-sm.res.list$sm.res
    
    ref.col<-which(colnames(sm.res)==refstr)
    procon<-sm.res[,c(1:3,ref.col,4)]
    procon$group<-paste0('Match_',procon$group)
    
    if (colnames(procon)[4]=='STC') {
      colnames(procon)[4]<-'S288C'
    }
    
    # adding cols of algos
    procon$ShannonE<-NA
    procon$StereochemE<-NA
    procon$JSDivergence<-NA
    procon$PhyloZOOM<-NA
    procon$SubMatrix<-NA
    if (phyloFull) {procon$Phylo<-NA}
    
    # add the xpos and ypos for interactive clicking/brushing
    procon$xpos<-sm.res$xpos
    procon$ypos<-sm.res$ypos
    
    # get the data
    ali.df<-sm.res[,5:(ncol(sm.res)-2)]
    
    # loop thru each pos to calculate conserv scores
    for (n in 1:nrow(ali.df)) {
      
      ali<-ali.df[n,]
      
      procon$SubMatrix[n]<-round(sub_matrix(ali,submat),digits=4)
      
      procon$ShannonE[n]<-round(shannon_E(ali, gap.pen),digits=4)
      procon$StereochemE[n]<-round(stereochem_E(ali, gap.pen),digits=4)
      procon$JSDivergence[n]<-round(js_divergence(ali,bkg.dist,gap.pen),digits=4)
      
      procon$PhyloZOOM[n]<-round(phylo_ZOOM(ali,refstr,tree,gap.pen),digits=4)
      
      # this step takes a much longer time
      if(phyloFull) {procon$Phylo[n]<-round(phylo_M(ali,tree,gap.pen),digits=4)}
      
    }
    
    return(procon)
  }
  
  
}

