# plot paralog pairs needleman alignment with matched places highlighted

PlotPara<-function(ndlali,sm.res,ymax) {
  
  ymax<-as.numeric(ymax)
  
  xmax<-5*(50+1)
  
  p1<-unique(ndlali$ORF)[1]
  p2<-unique(ndlali$ORF)[2]
  
  if (length(sm.res)>0) {
    
    hlplot<-sm.res[,c('S288C','xpos','ypos')]
    
    ggplot(data=ndlali)+
      coord_cartesian(xlim=c(0,xmax),ylim=c(0,ymax+5))+
      theme(panel.background=element_rect(fill='white'),
            plot.margin = margin(2, 2, 2, 2, "pt"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            legend.title=element_blank(),
            legend.position='none',
            axis.ticks.y = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())+
      geom_text(aes(x=xpos,y=ypos,label=needlealign,fontface='bold'),size=5)+
      geom_text(aes(x=xpos,y=n.ypos,label=label),size=2.5)+
      geom_label(data=hlplot,aes(x=xpos,y=ypos,label=S288C),fontface='bold',
                 size=5,fill='yellow',label.padding = unit(0.15, "lines"),
                 label.size = 0)+
      geom_text(aes(x=xpos[1],y=(n.ypos[1]+2),label=p1),size=3.5)+
      geom_text(aes(x=xpos[1],y=(n.ypos[1]-9),label=p2),size=3.5)
    
  } else {
    
    ggplot(data=ndlali)+
      coord_cartesian(xlim=c(0,xmax),ylim=c(0,ymax+5))+
      theme(panel.background=element_rect(fill='white'),
            plot.margin = margin(2, 2, 2, 2, "pt"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            legend.title=element_blank(),
            legend.position='none',
            axis.ticks.y = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())+
      geom_text(aes(x=xpos,y=ypos,label=needlealign,fontface='bold'),size=5)+
      geom_text(aes(x=xpos,y=n.ypos,label=label),size=2.5)+
      geom_text(aes(x=xpos[1],y=(n.ypos[1]+2),label=p1),size=3.5)+
      geom_text(aes(x=xpos[1],y=(n.ypos[1]-9),label=p2),size=3.5)
    
  }
  
  
}


