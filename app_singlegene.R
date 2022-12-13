# Load packages ----------------------------------------------------------------

library(seqinr)
library(DT)
library(shiny)
library(shinythemes)
library(shinycssloaders)
library(ggplot2)



# Load data and source code -----------------------------------------------------

source('Codes/SearchMotif.R')
source('Codes/SearchPosit.R')
source('Codes/GetSeq.R')
source('Codes/PlotSeq.R')
source('Codes/SeqPos.R')
source('Codes/SymFreq_table.R')
source('Codes/ProConserve.R')
source('Codes/PlotDim.R')
source('Codes/NumSeq.R')



# Define UI --------------------------------------------------------------------
ui <- fluidPage(
  
  #shinythemes::themeSelector(), # set theme selector
  theme = shinytheme("flatly"),
  
  titlePanel(title=span(img(src="icon1.jpg"),tags$p(HTML(paste0(strong("CoSMoS.c."), " - ", 
                                                                strong("Co"),"nserved ",strong("S"),"equence ",strong("Mo"),"tif in ",strong("S"),"accharomyces ", strong("c"),"erevisiae")))),
             windowTitle = "YeastMotifConserv"), # set title and tab title
  fluidRow(column(12,
                  h4('This algorithm can be used to search for and identify motifs conserved in the 1002 yeast genome project.')
                  )),
  sidebarLayout(
    sidebarPanel(width = 2,
                 
                 selectizeInput(
                   inputId = "gene",
                   label = "ORF Name:",
                   choices = NULL,
                   selected = NULL,
                   options=list(maxOptions = 3,maxItems = 1,
                                placeholder = 'select an ORF')
                 ),
                 
                 radioButtons("intype", "Search type:",
                              c("Motif" = "type.motif",
                                "Position" = "type.position"),
                              selected="type.motif"
                 ),
                 conditionalPanel(
                   condition = "input.intype == 'type.motif'",
                   
                   textInput(
                     inputId = "motif",
                     label = "Target Motif:",
                     placeholder = "LR,[SP]T,N.D"
                   )
                 ),
                 
                 conditionalPanel(
                   condition = "input.intype == 'type.position'",
                   
                   textInput(
                     inputId = "posit",
                     label = "Target Positions:",
                     placeholder = "9,25,198"
                   ),
                   numericInput(
                     inputId = "mlen",
                     label = "Motif length:",
                     value = 2,
                     min=1,
                     max=10
                   )
                 ),
                 
                 selectizeInput(
                   inputId = "refstr",
                   label = "Reference Strain Standard Name:",
                   choices = NULL,
                   selected = 'S288C',
                   options=list(maxOptions = 3,maxItems = 1,
                                placeholder = 'select a strain')
                 ),
                 
                 checkboxInput(
                   inputId = 'gap.pen',
                   label='apply gap penalty',
                   value=TRUE
                 ),
                 
                 actionButton('procal','Calculate',width='100%',class = "btn-success"),
                 
                 
                 br(),
                 hr(),
                 tags$p('For strain name details, refer to User Tips tab.'),
                 hr(),
                 tags$a("Strain details adapted from the 1002 Yeast Genome website.", 
                        href = "https://doi.org/10.1038/s41586-018-0030-5",target="_blank")
                 
                 
    ),
    
    mainPanel(width=10,
              wellPanel(id="genePanel",style = "height:50px; font-size: 22px; padding: 8px 15px",
                        uiOutput(outputId = 'genelink')),
              hr(),
              wellPanel(id = "seqPanel",style = "overflow-y:scroll; max-height: 200px",
                        #textOutput(outputId = 'trial'),
                        plotOutput(outputId = "geneseq", click = "plot_click")
              ),
              hr(),
              
              # set output tabs -- interactive, all symbol freq and all scores
              
              tabsetPanel(
                tabPanel("User Tips",
                         br(),
                         tags$p(strong('Search type - Motif:'), " input pattern is matched with regular expression using the same syntax and semantics as Perl (see examples below). The motif must be fixed in length."),
                         tags$p(HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),strong("Common examples:")),
                         tags$p(HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),strong("[ABC] - matches A or B or C;")),
                         tags$p(HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),strong("[^ABC] - matches anything except A, B, or C;")),
                         tags$p(HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),strong(". - matches anything once.")),
                         tags$p(HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),"Please refer to ", 
                                tags$a("this link",href="https://stat.ethz.ch/R-manual/R-devel/library/base/html/regex.html",target="_blank"), 
                                " and the detailed supporting information in Support Info tab for further details."),
                         br(),
                         tags$p(strong('Search type - Position:'), " input the position(s) of the first amino acids in the motif(s) and select the correct 'motif length'. Set 'motif length' to 1 if checking individual positions. In this case, insertion will be 0 for all positions."),
                         br(),
                         tags$p(strong('Standard Gene Name: '),'is displayed at the top panel on the right side. Gene name is hyperlinked to the corresponding SGD page.'),
                         br(),
                         tags$p(strong('Gap penalty: '),'when applied decreases the score if there are non-standard amino acids at the alignment position (target site).'),
                         br(),
                         tags$p(strong("Sequence Map "),"is plotted in the second panel on the right, with matched motif(s) or position(s) highlighted in yellow. '-' represents a gap in multi-sequence alignment (MSA) results, suggesting there is at least one other strain among the 1012 yeast strains that has an insertion in this position. MSA is performed with Clustal Omega 1.2.4. For further information about MSA, please refer to 'detailed supporting information' under Support Info."),
                         br(),
                         tags$p(strong('Selected Site: '),'display an interactive table showing stats for user selected amino acid position(s). Click on the yellow highlighted amino acids in the Sequence Plot above to check Symbol Frequency and Conservation Score for individual position. Click in between two amino acids to show both together for comparisons.'),
                         br(),
                         tags$p('Click on ',strong('Symbol Frequency '), 'and ', strong('Conservation Score '), 'for the full table of corresponding stats.'),
                         br(),
                         tags$p('All scores are normalized to: ', strong('conserved -> 1; relaxed -> 0.')),
                         br(),
                         hr(),
                         tags$p('Click here to download ',strong('Strain Name'),' details.'),
                         downloadButton("download_strname", "Download Strain Name sheet"),
                         br(),
                         br(),
                         tags$p('Amino acids will be displayed by their ', strong('one letter code '), '. Click here to download the code table.'),
                         downloadButton("download_aascode", "Download amino acid code table")
                ),
                tabPanel("Selected Site",
                         br(),
                         tags$p(strong('Symbol Frequency for selected site')), # bold
                         shinycssloaders::withSpinner(
                           dataTableOutput('slct_symfreq'),
                           type = 1, color = "#408000", size = 0.8,proxy.height=80
                         ),
                         hr(),
                         tags$p(strong('Conservation Scores for selected site')),
                         shinycssloaders::withSpinner(
                           dataTableOutput('slct_score'),
                           type = 1, color = "#408000", size = 0.8,proxy.height=80)
                ), 
                tabPanel("Symbol Frequency", 
                         br(),
                         shinycssloaders::withSpinner(
                           dataTableOutput("symbolfreqtable"),
                           type = 1, color = "#408000", size = 0.8,proxy.height=200
                         ),
                         hr(),
                         br(),
                         downloadButton("download_data_symfreq", "Download Symbol Frequency table")), 
                tabPanel("Conservation Score",
                         br(),
                         tags$i(strong('Refer to Support Info tab for details about each algorithm.')),
                         br(),
                         br(),
                         shinycssloaders::withSpinner(
                           dataTableOutput("scoretable"),
                           type = 1, color = "#408000", size = 0.8,proxy.height=200
                         ),
                         hr(),
                         br(),
                         downloadButton("download_data_score", "Download Conservation Scores Table")),
                tabPanel("Support Info",
                         br(),
                         tags$p(strong("Conservation Score interpretation:")),
                         br(),
                         tags$p(strong('Insertion: '),'calculate the percentage of amino acid insertion events within the motif among all strains in the database.'),
                         br(),
                         tags$p(strong('positions: '),'proteinposition is the position of the amino acid in the protein (not counting "-", gap in multi-sequence alignment). MSAposition is the position of amino acid in the multi-sequence alignment (counting the "-").'),
                         br(),
                         tags$p(tags$a(strong('SubMatrix: '),
                                       href = "https://doi.org/10.1128/jb.178.7.1881-1894.1996",target="_blank"), 
                                'likeliness of existing substitution based on BLOSUM62.'),
                         br(),
                         tags$p(tags$a(strong('ShannonE: '),
                                       href = "https://doi.org/10.1002/prot.340090107",target="_blank"),
                                'diversity of the target site.'),
                         br(),
                         tags$p(tags$a(strong('StereochemE: '),
                                       href = "https://doi.org/10.1006/jmbi.1999.2911",target="_blank"),
                                'stereochemical property of the diversity.'),
                         br(),
                         tags$p(tags$a(strong('JSDivergence: '),
                                       href = "https://doi.org/10.1093/bioinformatics/btm270",target="_blank"),
                                'the similarity of diversity with BLOSUM62 frequency (relaxed).'),
                         br(),
                         tags$p(tags$a(strong('PhyloZOOM: '),
                                       href = "https://doi.org/10.1016/j.jmb.2003.12.078",target="_blank"),
                                'diversity weighted by evolutionary distance to reference strain.'),
                         hr(),
                         tags$p('Click here to download detailed supporting information'),
                         downloadButton("download_supp", "Download SupportingInfo"),
                         br()
                ),
                tabPanel("Citation",
                         br(),
                         tags$p('Please cite Shuang Li and Henrik G. Dohlman, personal communication. shuang9@email.unc.edu; hdohlman@med.unc.edu.')
                  
                )
              )
              
    )
    
  )
  
)


# Define server ----------------------------------------------------------------

server <- function(input, output, session) {
  
  co_filep = 'CO_aligned'
  #gene.name.all<-dir(co_filep)
  gene.name.all<-as.character(read.csv('Data/ORFList.csv',header=F,stringsAsFactors=F))
  
  refstr.all<-as.character(read.csv('Data/RefStrainList.csv',header=F,stringsAsFactors=F))
  
  otg<-read.csv('Data/ORFtoGene.csv',header=T,stringsAsFactors=F)
  
  
  # updateRadioButtons(session,'intype')
  # updateTextInput(session, 'motif')
  # updateTextInput(session, 'posit')
  # updateNumericInput(session, 'mlen')
  # updateCheckboxInput(session,'gap.pen')
  
  updateSelectizeInput(session, 'gene', choices = gene.name.all, server = TRUE)
  
  updateSelectizeInput(session, 'refstr', choices = refstr.all, server = TRUE)
  
  sm.res.initial<-eventReactive(input$procal, {
    
    if (input$intype=='type.motif') {
      req(input$motif)
      
      SearchMotif(input$gene,input$motif,input$refstr,co_filep)
      
    } else {
      req(input$posit)
      req(input$mlen)
      
      SearchPosit(input$gene,input$posit,input$mlen,input$refstr,co_filep)
    }
  })
  
  
  refseq<-eventReactive(input$procal, {GetSeq(input$gene,input$refstr,co_filep)})
  
  
  sm.res<-eventReactive(input$procal, {SeqPos(refseq(),sm.res.initial())})
  
  sftable<-eventReactive(input$procal, {Symfreq_table(input$refstr,sm.res())})
  stable<-eventReactive(input$procal, {ProConserve(sm.res(), input$refstr, input$gap.pen, phyloFull=F)})
  
  ymax<-eventReactive(input$procal, {PlotDim(refseq())})
  
  # debugging module
  #output$trial<-renderText({
  #  otg[1:2,]
  #})
  
  
  observeEvent(input$procal, {
    
    output$geneseq <- renderPlot({
      
      PlotSeq(refseq(),input$refstr,sm.res()) 
    },
    width=1000,
    #height=1500)
    height=ymax()*10)
    
    output$genelink<-renderUI({
      orf.n<-which(otg$ORF==input$gene)
      tags$p(strong('Standard Gene Name: '), HTML('&nbsp;'),tags$a(strong(otg$Gene[orf.n]),
             href = otg$URL[orf.n],target="_blank"))
      
    })
    
  })
  

  output$symbolfreqtable<-renderDataTable({

    validate(
      if (input$intype=='type.motif') {
        need(length(sftable())>0, 'Input Motif does not exist in the target protein. Or Input Motif pattern is of flexible length. Please convert to fixed length pattern and redo search.')
      } else {
        need(length(sftable())>0, 'Please input positive integer(s) within the protein length. Remember the max position allowed is (total length)-(motif length).')
      }

    )


    sftable_clean<-sftable()
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

  })
  
  
  output$scoretable<-renderDataTable({
    
    validate(
      if (input$intype=='type.motif') {
        need(length(sftable())>0, 'Input Motif does not exist in the target protein. Or Input Motif pattern is of flexible length. Please convert to fixed length pattern and redo search.')
      } else {
        need(length(sftable())>0, 'Please input positive integer(s) within the protein length. Remember the max position allowed is (total length)-(motif length).')
      }
      
    )
    
    stable_clean<-stable()
    stable_clean$xpos<-NULL
    stable_clean$ypos<-NULL
    
    style.level<-unique(stable_clean$group)
    style.color<-rep(c('#202020','#2c7da0'),times=round(length(style.level)/2))
    style.color<-style.color[1:length(style.level)]
    
    datatable(
      data = stable_clean,
      options = list(pageLength = 10),
      rownames = FALSE
    ) %>% formatStyle('group',
                      target='row',
                      fontWeight='bold',
                      color=styleEqual(style.level,style.color))
    
  })
  
  output$slct_symfreq <- renderDataTable({
    slct_sf<-nearPoints(sftable(), input$plot_click,
                        xvar = "xpos", yvar = "ypos",
                        threshold=15, maxpoints=5) 
    slct_sf$xpos<-NULL
    slct_sf$ypos<-NULL
    datatable(
      data = slct_sf,
      rownames = FALSE
    )
  })

  
  
  output$slct_score <- renderDataTable({
    
    slct_s<-nearPoints(stable(), input$plot_click,
                       xvar = "xpos", yvar = "ypos",
                       threshold=15, maxpoints=5) 
    slct_s$xpos<-NULL
    slct_s$ypos<-NULL
    datatable(
      data = slct_s,
      rownames = FALSE
    )
  })
  
  
  output$download_data_symfreq <- downloadHandler(
    filename = function() {
      if (input$intype=='type.motif') {
        paste0(input$gene,'-',input$motif,'-',input$refstr,'_symbolfreq_table.csv')
      } else {
        paste0(input$gene,'-',input$posit,'-',input$refstr,'_symbolfreq_table.csv')
      }
      
    },
    content = function(file) { 
      write.csv(sftable(), file, row.names = F) 
    }
  )
  
  output$download_data_score <- downloadHandler(
    filename = function() {
      if (input$intype=='type.motif') {
        paste0(input$gene,'-',input$motif,'-',input$refstr,'_conservation_score_table.csv')
      } else {
        paste0(input$gene,'-',input$posit,'-',input$refstr,'_conservation_score_table.csv')
      }
      
    },
    content = function(file) { 
      write.csv(stable(), file, row.names = F) 
    }
  )
  
  
  output$download_strname <- downloadHandler(
    filename = 'Yeast1012StrainName.xlsx',
    content = function(file) { 
      file.copy("www/1012StrainName.xlsx", file) 
    }
  )
  
  output$download_supp <- downloadHandler(
    filename = 'YeastMotifConservationScore_SupportingInfo.pptx',
    content = function(file) { 
      file.copy("www/YeastMotifConservationScore_SupportingInfo.pptx", file) 
    }
  )
  
  output$download_aascode <- downloadHandler(
    filename = 'AminoAcidsCodetable.xlsx',
    content = function(file) { 
      file.copy("www/AAScode.xlsx", file) 
    }
  )

  
}

# Create a Shiny app object ----------------------------------------------------

shinyApp(ui = ui, server = server)
