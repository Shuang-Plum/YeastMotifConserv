sftable_clean<-sftable()

if (length(sftable_clean)>0) {
  sftable_clean$xpos<-NULL
  sftable_clean$ypos<-NULL
  
  fstyle.level<-unique(sftable_clean$group)
  fstyle.color<-rep(c('#202020','#2c7da0'),times=round(length(fstyle.level)/2))
  fstyle.color<-fstyle.color[1:length(fstyle.level)]
  
  datatable(
    data = sftable_clean,
    options = list(pageLength = 10),
    rownames = FALSE
  ) %>% formatStyle('group',
                    target='row',
                    fontWeight='bold',
                    color=styleEqual(fstyle.level,fstyle.color))
  
} else {
  
  if (input$intype=='type.motif') {
    stop(safeError(
      'Input Motif does not exist in the target protein.'
    ))
  } else {
    stop(safeError(
      'Please input positive integer(s) within the protein length.
            Remember the max position allowed is (aas length)-(motif length).'
    ))
  }
  
}
