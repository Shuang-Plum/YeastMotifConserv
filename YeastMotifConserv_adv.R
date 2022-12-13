# Load packages ----------------------------------------------------------------

library(seqinr)
library(DT)
library(shiny)
library(shinythemes)
library(shinycssloaders)
library(ggplot2)



# Load data and source code -----------------------------------------------------

source('./Codes/SearchMotif.R')
source('./Codes/SearchPosit.R')
source('./Codes/GetSeq.R')
source('./Codes/PlotSeq.R')
source('./Codes/SeqPos.R')
source('./Codes/SymFreq_table.R')
source('./Codes/ProConserve.R')
source('./Codes/PlotDim.R')



# Define UI --------------------------------------------------------------------
ui <- fluidPage(
  
  #shinythemes::themeSelector(), # set theme selector
  theme = shinytheme("flatly"),
  
  titlePanel("Yeast Motif Conservation Calculator",
             windowTitle = "YeastMotifConserv"), # set title and tab title
  
  sidebarLayout(
    sidebarPanel(width = 2,
                 
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
                     placeholder = "ST"
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
                   ),
                   tags$p("Set Motif length to 1 if checking individual positions. In this case, insertion is 0 for all positions."),
                 ),
                 
                 selectizeInput(
                   inputId = "gene",
                   label = "ORF Name:",
                   choices = NULL,
                   selected = 'YAL002W',
                   options=list(maxOptions = 5,maxItems = 1,
                                placeholder = 'select an ORF')
                 ),
                 
                 selectizeInput(
                   inputId = "refstr",
                   label = "Reference Strain Standard Name:",
                   choices = NULL,
                   selected = 'STC',
                   options=list(maxOptions = 3,maxItems = 1,
                                placeholder = 'select a strain')
                 ),
                 
                 checkboxInput(
                   inputId = 'gap.pen',
                   label='apply gap penalty',
                   value=TRUE
                 ),
                 
                 #br(),
                 tags$p("STC is S288C"),
                 
                 #br(),
                 tags$a("Refer here for Strain Name detail.", 
                        href = "1012StrainName.xlsx"),
                 hr(),
                 tags$a("Strain details adapted from the 1002 Yeast Genome website.", 
                        href = "https://doi.org/10.1038/s41586-018-0030-5")
                 
                 
    ),
    
    mainPanel(width=10,
              wellPanel(id = "seqPanel",style = "overflow-y:scroll; max-height: 200px",
                        #textOutput(outputId = 'trial'),
                        plotOutput(outputId = "geneseq", click = "plot_click"),
              ),
              hr(),
              
              # set output tabs -- interactive, all symbol freq and all scores
              
              tabsetPanel(
                tabPanel("Selected Site",
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
                tabPanel("SymbolFreqTable", 
                         shinycssloaders::withSpinner(
                           dataTableOutput("symbolfreqtable"),
                           type = 1, color = "#408000", size = 0.8,proxy.height=200
                         ),
                         hr(),
                         downloadButton("download_data_symfreq", "Download Symbol Frequency table")), 
                tabPanel("ConservScoreTable", 
                         shinycssloaders::withSpinner(
                           dataTableOutput("scoretable"),
                           type = 1, color = "#408000", size = 0.8,proxy.height=200
                         ),
                         hr(),
                         downloadButton("download_data_score", "Download Conservation Scores Table")),
                tabPanel("SupportInfo",
                         br(),
                         tags$p(strong("Conservation Score interpretation:")),
                         br(),
                         tags$p('All scores are normalized to have: ', strong('conserved -> 1; relaxed -> 0.')),
                         br(),
                         tags$p(strong('Gap penalty '),'when applied decreases the score if there are non-standard amino acids at the alignment position (target site).'),
                         br(),
                         tags$p(strong('Insertion '),'calculate the percentage of amino acid insertion event within the motif among all strains in the database.'),
                         br(),
                         tags$p(strong('SubMatrix: '),'likeliness of existing substitution based on BLOSUM62.'),
                         br(),
                         tags$p(strong('ShannonE: '),'diveristy of the target site.'),
                         br(),
                         tags$p(strong('StereochemE: '),'stereochemical property of the diversity.'),
                         br(),
                         tags$p(strong('JSDivergence: '),'the similarity of diversity with BLOSUM62 frequency (relaxed).'),
                         br(),
                         tags$p(strong('PhyloZOOM: '),'diversity weighted by evolutionary distance to refrence strain.'),
                         hr(),
                         tags$a("Refer here for detailed Supporting Information.", 
                                href = "YeastMotifConservationScore_SupportingInfo.pptx"),
                         br()
                )
              )
              
    )
    
  )
  
)


# Define server ----------------------------------------------------------------

server <- function(input, output, session) {
  
  co_filep = './CO_aligned'
  gene.name.all<-dir(co_filep)
  gene.name.all<-gsub('_aligned.fasta','',gene.name.all)
  
  refstr.all<-as.character(read.csv('./Data/RefStrainList.csv',header=F))
  
  updateRadioButtons(session,'intype')
  updateTextInput(session, 'motif')
  updateTextInput(session, 'posit')
  updateNumericInput(session, 'mlen')
  
  updateSelectizeInput(session, 'gene', choices = gene.name.all, server = TRUE)
  
  updateSelectizeInput(session, 'refstr', choices = refstr.all, server = TRUE)
  
  updateCheckboxInput(session,'gap.pen')

  
  gene.s<-reactive({input$gene})
  refstr.s<-reactive({input$refstr})
  
  gap.pen.s<-reactive({input$gap.pen})
  intype.s<-reactive({input$intype})
  
  motif.s<-reactive({input$motif})
  posit.s<-reactive({input$posit}) 
  mlen.s<-reactive({input$mlen})
  
  sm.res.initial<-reactive({
    if (intype.s()=='type.motif') {
      req(input$motif)
      
      searchMotif(gene.s(),motif.s(),refstr.s(),co_filep)
      
    } else {
      req(input$posit)
      req(input$mlen)
      
      searchPosit(gene.s(),posit.s(),mlen.s(),refstr.s(),co_filep)
    }
  })
  
  
  refseq<-reactive({
    req(input$gene)
    req(input$refstr)
    GetSeq(gene.s(),refstr.s(),co_filep)
           }) 
  
  sm.res<-reactive(SeqPos(refseq(),sm.res.initial())) 
  
  sftable<-reactive(Symfreq_table(refstr.s(),sm.res()))
  stable<-reactive(ProConserve(sm.res(), refstr.s(), gap.pen.s(), phyloFull=F))
  
  ymax<-reactive(PlotDim(refseq()))
  
  # output$trial<-renderText({
  #   req(input$posit)
  #   as.numeric(unlist(strsplit(input$posit,',',fixed=T)))
  # })
  
  observe({
    req(input$gene)
    req(input$refstr)
    req(input$intype)
    if (intype.s()=='type.motif') {
      req(input$motif)
    } else {
      req(input$posit)
      req(input$mlen)
      }
    output$geneseq <- renderPlot({
      req(input$gene)
      req(input$refstr)
      req(input$intype)
      if (intype.s()=='type.motif') {
        req(input$motif)
      } else {
        req(input$posit)
        req(input$mlen)
      }
      PlotSeq(refseq(),refstr.s(),sm.res()) 
    },
    width=1000,
    #height=1500)
    height=ymax()*10)
  })
  
  # output$geneseq <- renderPlot({
  #   req(input$gene)
  #   req(input$refstr)
  #   req(input$motif)
  #   PlotSeq(refseq(),sm.res())
  # },
  # width=1000,
  # height=1500)
  
  output$symbolfreqtable<-renderDataTable({
    req(input$gene)
    req(input$refstr)
    req(input$intype)
    if (intype.s()=='type.motif') {
      req(input$motif)
    } else {
      req(input$posit)
      req(input$mlen)
    }
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
    req(input$gene)
    req(input$refstr)
    req(input$intype)
    if (intype.s()=='type.motif') {
      req(input$motif)
    } else {
      req(input$posit)
      req(input$mlen)
    }
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
  
  # download data
  output$download_data_symfreq <- downloadHandler(
    filename = function() {
      paste0(gene.s(),'-',motif.s(),'-',refstr.s(),'_symbolfreq_table.csv')
    },
    content = function(file) { 
      write.csv(sftable(), file) 
    }
  )
  
  output$download_data_score <- downloadHandler(
    filename = function() {
      paste0(gene.s(),'-',motif.s(),'-',refstr.s(),'_conservation_score_table.csv')
    },
    content = function(file) { 
      write.csv(stable(), file) 
    }
  )
  
  
  output$slct_symfreq <- renderDataTable({
    req(input$gene)
    req(input$refstr)
    req(input$intype)
    if (intype.s()=='type.motif') {
      req(input$motif)
    } else {
      req(input$posit)
      req(input$mlen)
    }
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
    req(input$gene)
    req(input$refstr)
    req(input$intype)
    if (intype.s()=='type.motif') {
      req(input$motif)
    } else {
      req(input$posit)
      req(input$mlen)
    }
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
  
}

# Create a Shiny app object ----------------------------------------------------

shinyApp(ui = ui, server = server)