# Load packages ----------------------------------------------------------------
library(Biostrings)
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
source('Codes/paraGetSeq.R')
source('Codes/NdlPos.R')
source('Codes/paraSearchMotif.R')
source('Codes/paraSearchPosit.R')
source('Codes/paraAssignPos.R')
source('Codes/paraAddPos.R')
source('Codes/paraSymFreq_table.R')
source('Codes/paraProConserve.R')
source('Codes/paraPlotYDim.R')
source('Codes/PlotPara.R')
source('Codes/cleantable.R')


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
  tabsetPanel(
    
    tabPanel(title = HTML(paste(p(strong("Single Gene"),style = "font-size:18px"))),
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
                            
                            
                            hr(),
                            tags$p('Strain name details can be downloaded from User Tips tab.'),
                            hr(),
                            tags$a("Strains are adapted from the 1002 Yeast Genome paper.", 
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
                                    tags$p(strong('ORF Name: '),'All available ORF names are included in the dropdown list. Some ORFs are not included. Refer to Support Info for further details. The corresponding ',strong('Standard Gene Name '),'is displayed at the top panel on the right side. Gene name is hyperlinked to its SGD page.'),
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
                                    tags$p(strong('Search type - Position:'), " input the position(s) of the first amino acids in the motif(s) and select the correct 'motif length' (between 1 to 10). Position is defined as amino acid position in the protein. Set 'motif length' to 1 if checking individual positions. In this case, insertion will be 0 for all positions."),
                                    br(),
                                    tags$p(strong('Reference Strain: '),' the strain you are interested at (used as reference) for the algorithm PholyZOOM. Refer to Support Info for further details.'),
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
                                    hr(),
                                    tags$p('Click here to download ',strong('Strain Name'),' details.'),
                                    downloadButton("download_strname", "Download Strain Name sheet"),
                                    br(),
                                    br(),
                                    tags$p('Amino acids will be displayed by their ', strong('one letter code'), '. Click here to download the code table.'),
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
                                    tags$p(strong("positions: "),"proteinposition is the position of the amino acid in the protein (not counting '-', gap in multi-sequence alignment). MSAposition is the position of amino acid in the multi-sequence alignment (counting the '-')."),
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
                                    tags$p(strong('If you use CoSMoS.c., please cite '), tags$a(strong('this paper:'),href = 'https://doi.org/10.1101/2022.09.15.508132',target='_blank')),
                                    br(),
                                    tags$p('Shuang Li, Henrik Dohlman. Differential modification of protein paralogs reveals conserved sequence determinants of post-translational modification. (2022)'),
                                    hr(),
                                    tags$p('For source code, please refer to this ', 
                                           tags$a('github repository.', href = 'https://github.com/Shuang-Plum/YeastMotifConserv',target='_blank')),
                                    hr(),
                                    tags$p('For questions or to report problems, please email: shuang9@email.unc.edu; hdohlman@med.unc.edu.')
                                    
                           )
                         )
                         
               )
               
             )
             
             
    ),
    
    tabPanel(title = HTML(paste(p(strong("Paralogs"),style = "font-size:18px"))),
             sidebarLayout(
               sidebarPanel(width = 2,
                            
                            selectizeInput(
                              inputId = "paralogs",
                              label = "Paralog ORF Names:",
                              choices = NULL,
                              selected = NULL,
                              options=list(maxOptions = 3,maxItems = 1,
                                           placeholder = 'select paralogs')
                            ),
                            
                            tags$head(tags$style(HTML('#paralogs+ div>.selectize-input{font-size: 11px !important}
                                                      #paralogs+ div>.selectize-dropdown{font-size: 12px !important}'))),
                            
                            
                            radioButtons("paraintype", "Search type:",
                                         c("Motif" = "paratype.motif",
                                           "Position" = "paratype.position"),
                                         selected="paratype.motif"
                            ),
                            conditionalPanel(
                              condition = "input.paraintype == 'paratype.motif'",
                              
                              textInput(
                                inputId = "paramotif",
                                label = "Target Motif:",
                                placeholder = "LR,[SP]T,N.D"
                              )
                            ),
                            
                            conditionalPanel(
                              condition = "input.paraintype == 'paratype.position'",
                              
                              textInput(
                                inputId = "paraposit",
                                label = "Target Positions:",
                                placeholder = "9,25,198"
                              ),
                              numericInput(
                                inputId = "paramlen",
                                label = "Motif length:",
                                value = 2,
                                min=1,
                                max=10
                              )
                            ),
                            
                            checkboxInput(
                              inputId = 'paragap.pen',
                              label='apply gap penalty',
                              value=TRUE
                            ),
                            
                            br(),
                            tags$p('Reference strain is set to S288C only'),
                            br(),
                            
                            actionButton('paraprocal','Calculate',width='100%',class = "btn-success"),
                            
                            
                            hr()
               ),
               
               mainPanel(width=10,
                         wellPanel(id="paragenePanel",style = "height:50px; font-size: 22px; padding: 8px 15px",
                                   uiOutput(outputId = 'paragenelink')),
                         hr(),
                         wellPanel(id = "paraseqPanel",style = "overflow-y:scroll; max-height: 300px",
                                   #textOutput(outputId = 'trial'),
                                   plotOutput(outputId = "parageneseq", click = "paraplot_click")
                         ),
                         hr(),
                         
                         # set output tabs -- interactive, all symbol freq and all scores
                         
                         tabsetPanel(
                           tabPanel("User Tips",
                                    br(),
                                    tags$p(strong("Paralog ORF Names: "), "input pair of paralog ORF names with '_' in between. The sequence of the two ORFs matters for the following analysis. Refer to the paragraphs below about Search Types for details. All pairs of ORF names are included in the dropdown list. Some paralog pairs are not included. Refer to Support Info for further details. Their corresponding ",strong('Standard Gene Names '),'are displayed at the top panel on the right side. Gene names are hyperlinked to their SGD pages.'),
                                    br(),
                                    tags$p(strong("Search type - Motif:"), " input pattern is matched with regular expression using the same syntax and semantics as Perl (see examples below). The motif must be fixed in length. For Paralogs analysis, the motif(s) are matched for both ORFs. Additionally, for any motif match that exists in only one ORF, its corresponding sites in the other ORF are also included for analysis. Input ORF name sequence does not affect Motif search results."),
                                    tags$p(HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),strong("Common examples:")),
                                    tags$p(HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),strong("[ABC] - matches A or B or C;")),
                                    tags$p(HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),strong("[^ABC] - matches anything except A, B, or C;")),
                                    tags$p(HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),strong(". - matches anything once.")),
                                    tags$p(HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),"Please refer to ", 
                                           tags$a("this link",href="https://stat.ethz.ch/R-manual/R-devel/library/base/html/regex.html",target="_blank"), 
                                           " and the detailed supporting information in Support Info tab for further details."),
                                    br(),
                                    tags$p(strong("Search type - Position:"), " input the position(s) of the first amino acids in the motif(s) and select the correct 'motif length' (between 1 to 10). Position is defined as amino acid position in the protein. Set 'motif length' to 1 if checking individual positions. In this case, insertion will be 0 for all positions. For Paralogs analysis, the position(s) are searched against the first input ORF. Then the corresponding position(s) in the second ORF are also included for analysis. The order of ORF names affect Position search results. "),
                                    br(),
                                    tags$p(strong("Gap penalty: "),"when applied decreases the score if there are non-standard amino acids at the alignment position (target site)."),
                                    br(),
                                    tags$p(strong("Sequence Map "),"is plotted in the second panel on the right. Paralog sequences are plotted against each other with Needleman–Wunsch global alignment. '-' represents a gap in the alignment. ORF names are noted at the beginning of the sequences. Amino acid positions in the proteins (not counting '-') are noted above or below each sequence. Matched motif(s) or position(s) are highlighted in yellow. The Needleman–Wunsch alignment is performed with pairwiseAlignment from the Biostrings package in R."),
                                    br(),
                                    tags$p(strong("Selected Site: "),"display an interactive table showing stats for user selected amino acid position(s). Click on the yellow highlighted amino acids in the Sequence Plot above to check Symbol Frequency and Conservation Score for individual position. Click in between two amino acids to show both together for comparisons."),
                                    br(),
                                    tags$p('Click on ',strong('Symbol Frequency '), 'and ', strong('Conservation Score '), 'for the full table of corresponding stats. Results are arranged in a way that for each matched position, the two paralog sites are juxtaposed for better comparisions.'),
                                    br(),
                                    tags$p('All scores are normalized to: ', strong('conserved -> 1; relaxed -> 0.')),
                                    hr(),
                                    tags$p('Amino acids will be displayed by their ', strong('one letter code '), '. Click here to download the code table.'),
                                    downloadButton("paradownload_aascode", "Download amino acid code table")
                           ),
                           tabPanel("Selected Site",
                                    br(),
                                    tags$p(strong('Symbol Frequency for selected site')), # bold
                                    shinycssloaders::withSpinner(
                                      dataTableOutput('paraslct_symfreq'),
                                      type = 1, color = "#408000", size = 0.8,proxy.height=80
                                    ),
                                    hr(),
                                    tags$p(strong('Conservation Scores for selected site')),
                                    shinycssloaders::withSpinner(
                                      dataTableOutput('paraslct_score'),
                                      type = 1, color = "#408000", size = 0.8,proxy.height=80)
                           ), 
                           tabPanel("Symbol Frequency", 
                                    br(),
                                    shinycssloaders::withSpinner(
                                      dataTableOutput("parasymbolfreqtable"),
                                      type = 1, color = "#408000", size = 0.8,proxy.height=200
                                    ),
                                    hr(),
                                    br(),
                                    downloadButton("paradownload_data_symfreq", "Download Symbol Frequency table")), 
                           tabPanel("Conservation Score",
                                    br(),
                                    tags$i(strong('Refer to Support Info tab for details about each algorithm.')),
                                    br(),
                                    br(),
                                    shinycssloaders::withSpinner(
                                      dataTableOutput("parascoretable"),
                                      type = 1, color = "#408000", size = 0.8,proxy.height=200
                                    ),
                                    hr(),
                                    br(),
                                    downloadButton("paradownload_data_score", "Download Conservation Scores Table")),
                           tabPanel("Support Info",
                                    br(),
                                    tags$p(strong("Conservation Score interpretation:")),
                                    br(),
                                    tags$p(strong("positions: "),"proteinposition is the position of the amino acid in the protein (not counting '-', gap in multi-sequence alignment). MSAposition is the position of amino acid in the multi-sequence alignment (counting the '-')."),
                                    br(),
                                    tags$p(strong('Insertion: '),'calculate the percentage of amino acid insertion events within the motif among all strains in the database.'),
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
                                    downloadButton("paradownload_supp", "Download SupportingInfo"),
                                    br()
                           ),
                           tabPanel("Citation",
                                    br(),
                                    tags$p(strong('If you use CoSMoS.c., please cite '), tags$a(strong('this paper:'), href = 'https://doi.org/10.1101/2022.09.15.508132',target='_blank')),
                                    br(),
                                    tags$p('Shuang Li, Henrik Dohlman. Differential modification of protein paralogs reveals conserved sequence determinants of post-translational modification. (2022)'),
                                    hr(),
                                    tags$p('For source code, please refer to this ', 
                                           tags$a('github repository.', href = 'https://github.com/Shuang-Plum/YeastMotifConserv',target='_blank')),
                                    hr(),
                                    tags$p('For questions or to report problems, please email: shuang9@email.unc.edu; hdohlman@med.unc.edu.')
                                    
                           )
                         )
                         
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
  
  para.pair<-as.character(read.csv('Data/ParaList.csv',header=F,stringsAsFactors=F))
  
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
  
  sftable_out<-eventReactive(input$procal, {cleantable(sftable())})
  stable_out<-eventReactive(input$procal, {cleantable(stable())})
  
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


    sftable_clean<-sftable_out()

    fstyle.level<-unique(sftable_clean$group)
    fstyle.color<-rep(c('#202020','#2c7da0'),times=round(length(fstyle.level)/2))
    fstyle.color<-fstyle.color[1:length(fstyle.level)]

    datatable(
      data = sftable_clean,
      options = list(pageLength = 10, scrollX = T),
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
    
    stable_clean<-stable_out()
    
    style.level<-unique(stable_clean$group)
    style.color<-rep(c('#202020','#2c7da0'),times=round(length(style.level)/2))
    style.color<-style.color[1:length(style.level)]
    
    datatable(
      data = stable_clean,
      options = list(pageLength = 10, scrollX = T),
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
      write.csv(sftable_out(), file, row.names = F) 
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
      write.csv(stable_out(), file, row.names = F) 
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
  
  
  ##### paralogs tab
  
  para.pair<-as.character(read.csv('Data/ParaList.csv',header=F,stringsAsFactors=F))
  
  updateSelectizeInput(session, 'paralogs', choices = para.pair, server = TRUE)
  
  paraseq<-eventReactive(input$paraprocal, {paraGetSeq(input$paralogs,co_filep)})
  
  # needleman alignment
  para.ndlali<-eventReactive(input$paraprocal, {NdlPos(paraseq())})
  
  para.sm.res.initial<-eventReactive(input$paraprocal, {
    
    if (input$paraintype=='paratype.motif') {
      req(input$paramotif)
      
      paraSearchMotif(input$paralogs,co_filep,input$paramotif,para.ndlali())
      
    } else {
      req(input$paraposit)
      req(input$paramlen)
      
      paraSearchPosit(input$paralogs,co_filep,input$paraposit,input$paramlen,para.ndlali())
    }
  })
  
  para.ndlali.pos<-eventReactive(input$paraprocal, {paraAssignPos(para.ndlali())})
  
  para.sm.res<-eventReactive(input$paraprocal, {paraAddPos(para.ndlali.pos(),para.sm.res.initial())})
  
  parasftable<-eventReactive(input$paraprocal, {paraSymfreq_table(para.sm.res())})
  parastable<-eventReactive(input$paraprocal, {paraProConserve(para.sm.res(), input$paragap.pen, phyloFull=F)})
  
  parasftable_out<-eventReactive(input$paraprocal, {cleantable(parasftable())})
  parastable_out<-eventReactive(input$paraprocal, {cleantable(parastable())})
  
  paraymax<-eventReactive(input$paraprocal, {paraPlotYDim(para.ndlali.pos())})
  
  # debugging module
  #output$trial<-renderText({
  #  otg[1:2,]
  #})
  
  
  observeEvent(input$paraprocal, {
    
    output$parageneseq <- renderPlot({
      
      PlotPara(para.ndlali.pos(),para.sm.res(),paraymax()) 
    },
    width=1000,
    #height=1500)
    height=paraymax()*10)
    
    output$paragenelink<-renderUI({
      paras<-unlist(strsplit(input$paralogs,'_'))
      orf.1<-which(otg$ORF==paras[1])
      orf.2<-which(otg$ORF==paras[2])
      tags$p(strong('Standard Gene Name: '), HTML('&nbsp;'),
             tags$a(strong(otg$Gene[orf.1]),href = otg$URL[orf.1],target="_blank"),
             tags$a(strong('-')),
             tags$a(strong(otg$Gene[orf.2]),href = otg$URL[orf.2],target="_blank"))
      
    })
    
  })
  
  
  output$parasymbolfreqtable<-renderDataTable({
    
    validate(
      if (input$paraintype=='paratype.motif') {
        need(length(parasftable())>0, 'Input Motif does not exist in the paralogs. Or Input Motif pattern is of flexible length. Please convert to fixed length pattern and redo search.')
      } else {
        need(length(parasftable())>0, 'Please input positive integer(s) within the paralogs protein length. Remember the max position allowed is (total length)-(motif length).')
      }
      
    )
    
    
    parasftable_clean<-parasftable_out()
    
    parafstyle.level<-unique(parasftable_clean$paralog.pair)
    parafstyle.color<-rep(c('#202020','#2c7da0'),
                          times=round(length(parafstyle.level)/2))
    parafstyle.color<-parafstyle.color[1:length(parafstyle.level)]
    
    datatable(
      data = parasftable_clean,
      options = list(pageLength = 10, scrollX = T),
      rownames = FALSE
    ) %>% formatStyle('paralog.pair',
                      target='row',
                      fontWeight='bold',
                      color=styleEqual(parafstyle.level,parafstyle.color))
    
  })
  
  
  output$parascoretable<-renderDataTable({
    
    validate(
      if (input$paraintype=='paratype.motif') {
        need(length(parasftable())>0, 'Input Motif does not exist in the paralogs. Or Input Motif pattern is of flexible length. Please convert to fixed length pattern and redo search.')
      } else {
        need(length(parasftable())>0, 'Please input positive integer(s) within the paralogs protein length. Remember the max position allowed is (total length)-(motif length).')
      }
      
    )
    
    parastable_clean<-parastable_out()
    
    parastyle.level<-unique(parastable_clean$paralog.pair)
    parastyle.color<-rep(c('#202020','#2c7da0'),
                         times=round(length(parastyle.level)/2))
    parastyle.color<-parastyle.color[1:length(parastyle.level)]
    
    datatable(
      data = parastable_clean,
      options = list(pageLength = 10, scrollX = T),
      rownames = FALSE
    ) %>% formatStyle('paralog.pair',
                      target='row',
                      fontWeight='bold',
                      color=styleEqual(parastyle.level,parastyle.color))
    
  })
  
  output$paraslct_symfreq <- renderDataTable({
    paraslct_sf<-nearPoints(parasftable(), input$paraplot_click,
                            xvar = "xpos", yvar = "ypos",
                            threshold=15, maxpoints=6) 
    paraslct_sf$xpos<-NULL
    paraslct_sf$ypos<-NULL

    datatable(
      data = paraslct_sf,
      rownames = FALSE
    )
  })
  
  
  
  output$paraslct_score <- renderDataTable({
    
    paraslct_s<-nearPoints(parastable(), input$paraplot_click,
                           xvar = "xpos", yvar = "ypos",
                           threshold=15, maxpoints=6) 
    paraslct_s$xpos<-NULL
    paraslct_s$ypos<-NULL

    datatable(
      data = paraslct_s,
      rownames = FALSE
    )
  })
  
  
  output$paradownload_data_symfreq <- downloadHandler(
    filename = function() {
      if (input$paraintype=='paratype.motif') {
        paste0(input$paralogs,'-',input$paramotif,'_symbolfreq_table.csv')
      } else {
        paste0(input$paralogs,'-',input$paraposit,'_symbolfreq_table.csv')
      }
      
    },
    content = function(file) { 
      write.csv(parasftable_out(), file, row.names = F) 
    }
  )
  
  output$paradownload_data_score <- downloadHandler(
    filename = function() {
      if (input$paraintype=='paratype.motif') {
        paste0(input$paralogs,'-',input$paramotif,'_conservation_score_table.csv')
      } else{
        paste0(input$paralogs,'-',input$paraposit,'_conservation_score_table.csv')
      }
      
    },
    content = function(file) { 
      write.csv(parastable_out(), file, row.names = F) 
    }
  )
  
  
  output$paradownload_supp <- downloadHandler(
    filename = 'YeastMotifConservationScore_SupportingInfo.pptx',
    content = function(file) { 
      file.copy("www/YeastMotifConservationScore_SupportingInfo.pptx", file) 
    }
  )
  
  output$paradownload_aascode <- downloadHandler(
    filename = 'AminoAcidsCodetable.xlsx',
    content = function(file) { 
      file.copy("www/AAScode.xlsx", file) 
    }
  )

  
}

# Create a Shiny app object ----------------------------------------------------

shinyApp(ui = ui, server = server)
