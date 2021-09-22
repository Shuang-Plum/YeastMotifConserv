# get motif alignments
# 
# library(seqinr)
# co_filep<-'E:/YeastMotifConserv_Shiny/YeastMotifConserv/CO_aligned'
# 
# gene<-'YOR267C'
# 
# motif<-'QQQ'
# 
# gene<-'YAL008W'
# motif<-'KS'
# 
# refstr<-'S288C'

searchMotif<-function(gene,motif,refstr,co_filep) {
  
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
  refseq.c<-refseq[which(refseq!='-' &  refseq!='*')]
  # save the names as the position in the original alignment
  refseq.pos<-names(refseq.c)
  # concatenate into a string for searching motif
  refseq.c<-paste(refseq.c,collapse = '')
  
  
  motif.n<-unlist(strsplit(motif,',',fixed=T))
  # search for motif
  # use regular expression lookahead assertions to match motif in a repeats
  # (?=xx)
  # after match, it steps back and rematch so that CC can match twice in CCC
  motif.n<-toupper(motif.n)
  
  for (mn in 1:length(motif.n)) {
    
    pmotif<-paste0('(?=',motif.n[mn],')')
    motif.match<-gregexpr(pattern=pmotif,refseq.c,perl=T)
    
    # only get the match position
    motif.match<-as.vector(motif.match[[1]])
    
    
    # only proceed when there is a match
    if (motif.match[1]!=(-1)) {
      
      # get the motif length from match return
      # length would be different from word length if regrex is used
      motif.len<-gregexpr(pattern=motif.n[mn],refseq.c,perl=T)
      motif.len<-unique(attr(motif.len[[1]],'match.length'))
      
      # only support fixed length regexp search now
      # positive look-ahead assertion (?=...) match cannot return match length
      # therefore, for flexible length regrex, cannot determine end position properly
      # for fixed length regrex, use a regular match without the assertion to get the accurate word length
      # could be improved if there is an easy way to determine word length
      
      # only proceed if the regrex length is fixed and can be determined
      if (length(motif.len)==1) {
        
        # initialize a dataframe to store the alignments at each position
        # each row is a position
        # each col is a strain
        match.num<-length(motif.match)
        ali.df<-data.frame(matrix(ncol=(length(snames)+3),nrow=(match.num*motif.len)))
        colnames(ali.df)<-c('group','position','insertion',snames)
        ali.df$group<-paste0(motif.n[mn],'_',rep(c(1:match.num),each=motif.len))
        
        # when getting the alignments need to use the ori coordinates (MSA alignment) with '-' and '*"
        
        
        for (n in 1:match.num) {
          # get each match motif
          ma<-motif.match[n]
          ori.coord<-refseq.pos[ma:(ma+motif.len-1)]
          ori.coord<-as.numeric(ori.coord)
          
          temp.ins<-0
          # test if there is insertion that split the motif
          ori.coord.s<-ori.coord[1]
          ori.coord.exp.e<-ori.coord.s+motif.len-1
          if(ori.coord[motif.len]!=ori.coord.exp.e) {
            #temp.ins<-0.5
            # get insertion position
            ins.pos<-setdiff(c(ori.coord[1]:ori.coord[motif.len]),ori.coord)
            # go thru all strains to see if ins.pos are all '-'
            for (s in 1:length(allseq)) {
              ins.seq<-allseq[[s]][ins.pos]
              if (sum(ins.seq!='-')>0) {temp.ins<-temp.ins+1}
            }
          }
          
          # get each aas inside a motif
          for (m in 1:motif.len) {
            # get the ori position (MAS position) of the aas
            ori.co<-ori.coord[m]
            
            # get the aligned aas at that pos for all strains
            temp.ali<-sapply(allseq,'[[',ori.co)
            
            # corresponding row number in ali.df
            ali.rown<-(n-1)*motif.len+m
            
            ali.df$position[ali.rown]<-ori.co
            ali.df$insertion[ali.rown]<-round(temp.ins/length(snames),digits=4)
            ali.df[ali.rown,4:ncol(ali.df)]<-temp.ali
            
          }
          
        }
        
        if (!exists('ali.df.c')) {
          ali.df.c<-ali.df
        } else {
          ali.df.c<-rbind(ali.df.c,ali.df)
        }
        
      }
      
    }
    
  }
  if (exists('ali.df.c')) {return(ali.df.c)}
  
}
