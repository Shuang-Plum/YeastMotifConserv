# get posit alignments

# library(seqinr)
# co_filep<-'E:/YeastMotifConserv_Shiny/CO_aligned'
# 
# gene<-'YOR267C'
# 
# posit<-'123,135,196,202'
# mlen<-2
# 
# refstr<-'STC'

SearchPosit<-function(gene,posit,mlen,refstr,co_filep) {
  
  co<-read.fasta(file=paste0(co_filep,'/',gene,'_aligned.fasta'),seqtype = 'AA',
                 as.string = F,strip.desc = T,)
  
  snames<-getName(co)
  snames<-strsplit(snames,split=paste0('_',gene,'_'),fixed=T)
  # get the strain name
  snames<-sapply(snames,'[[',1)
  # replace S288C to STC to be consistent
  snames<-gsub('S288C','STC',snames)
  
  if(refstr=='S288C') {refstr<-'STC'}
  
  refseq<-getSequence(co[which(snames==refstr)],as.string = F)
  
  allseq<-getSequence(co,as.string=F)
  
  refseq<-unlist(refseq)
  
  # name each aas according to its position
  names(refseq)<-c(1:length(refseq))
  
  # remove the insertion mark and end codon 
  # this is the refstr seq which position input corresponds to
  refseq.c<-refseq[which(refseq!='-' &  refseq!='*')]
  # save the names as the position in the original alignment
  refseq.pos<-names(refseq.c)
  # generate original pos as in protein
  proseq.pos<-c(1:length(refseq.c))
  
  # get position of input, which is the pos in refstr
  posit.n<-as.numeric(unlist(strsplit(posit,',',fixed=T)))
  mlen<-as.numeric(mlen)

  
  # only proceed if input is within refstr protein limit
  if (sum(posit.n<=0)==0 & (max(posit.n)+mlen-1)<=length(refseq.c)) {
    
    
    # initialize a dataframe to store the alignments at each position
    # each row is a position
    # each col is a strain
    match.num<-length(posit.n)
    ali.df<-data.frame(matrix(ncol=(length(snames)+4),nrow=(match.num*mlen)))
    colnames(ali.df)<-c('group','proteinposition','MSAposition','insertion',snames)
    ali.df$group<-rep(c(1:match.num),each=mlen)
    
    # when getting the alignments need to use the ori coordinates (MSA alignment) with '-' and '*"
    
    
    for (n in 1:match.num) {
      # get each match motif
      po<-posit.n[n]
      ori.coord<-refseq.pos[po:(po+mlen-1)]
      ori.coord<-as.numeric(ori.coord)
      
      pro.coord<-proseq.pos[po:(po+mlen-1)]
      
      temp.ins<-0
      # test if there is insertion that split the motif
      ori.coord.s<-ori.coord[1]
      ori.coord.exp.e<-ori.coord.s+mlen-1
      if(ori.coord[mlen]!=ori.coord.exp.e & mlen>1) {
        #temp.ins<-0.5
        # get insertion position
        ins.pos<-setdiff(c(ori.coord[1]:ori.coord[mlen]),ori.coord)
        # go thru all strains to see if ins.pos are all '-'
        for (s in 1:length(allseq)) {
          ins.seq<-allseq[[s]][ins.pos]
          if (sum(ins.seq!='-')>0) {temp.ins<-temp.ins+1}
        }
      }
      
      # get each aas inside a motif
      for (m in 1:mlen) {
        # get the ori position (MAS position) of the aas
        ori.co<-ori.coord[m]
        pro.co<-pro.coord[m]
        # get the aligned aas at that pos for all strains
        temp.ali<-sapply(allseq,'[[',ori.co)
        
        # corresponding row number in ali.df
        ali.rown<-(n-1)*mlen+m
        
        ali.df$proteinposition[ali.rown]<-pro.co
        ali.df$MSAposition[ali.rown]<-ori.co
        ali.df$insertion[ali.rown]<-round(temp.ins/length(snames),digits=4)
        ali.df[ali.rown,5:ncol(ali.df)]<-temp.ali
        
      }
      
    } 
    
    return(ali.df)
  } 
  
}
