# LIBRARIES
library(shiny)
library(openxlsx)
library(DT)
library(plotly)
library(shinyjs)
library(ggplot2)
library(reshape2)
library(igraph)
library(edgebundleR)


# SOURCES
source("data/processGC.R")
source("data/correlationGC.R")
source("data/processHILIC.R")
source("data/correlationHILIC.R")
source("data/processC18.R")
source("data/correlationC18.R")
options(shiny.sanitize.errors = TRUE)




ui <- shinyUI(
  navbarPage("AutismExposome Chemical Explorer",
    tabPanel("Home",
      absolutePanel( width = "90%", left = "5%", right = "5%",
        column(6,
            wellPanel(
              p(HTML(paste0(h4(a(target="_blank", href="https://autismexposome.org/", "AutismExposome"), " Project")))),
              hr(),
              p('The AutismExposome project aims to investigate what contributed a significant increase to the autism prevalence in the U.S. over the past decades, focusing on discovery of internal and external environmental exposures that are detectable in the blood from the patients with autism and their family members compared to individuals without autism.'),
              p('The project hypothesizes that environmental exposures contributing to autism create an exposure history reflected in transcriptomic and environmental metabolomic profiles. Therefore, this multi-faceted investigation provides an opportunity to search for novel factors in ASD that may play a role in ASD etiology and preliminary information for future research.'),
              p('The project is interrogating both child and parent exposomes, defined as the cumulative measure of internal and external environmental influences and associated biological responses (including exposures from environment, diet, behavior and endogenous process) throughout the lifespan.'),
              p(HTML(paste0('The AutismExposome project is funded by the National Institutes of Mental Health (NIMH ',
                            a(target="_blank", href="https://projectreporter.nih.gov/project_info_description.cfm?aid=8916991&icde=28673319", "R01MH107205"), ')'))),
              p(HTML(paste0('Check ', a(target="_blank", href="https://autismexposome.org/", "here"), ' AutismExposome project\'s web page for more information.')))
            )
        ),
        column(6,
            wellPanel(
              h4("Autism Exposome's Chemical Explorer"),
              hr(),
              p('The AutismExposome Chemical Explorer is an interactive web to explorer the different chemical measures obtained from a set of 93 samples including 28 Control, 22 Cases, 27 Father (cases related) and 16 Mother (cases related).'),
              p('This platform allows to explore the different levels of measured chemicals from both targeted (GC-MS) and un-targeted (LC-HRMS) metabolome data.'),
              p('Each tabs offers a specific set of data to explore.'),
              tags$ul(
              tags$li('The table tab allows to explore the amount of data recollected in each experiment as well as their annotation.'),
              tags$li('The boxplot tab allows to study the levels of each chemical in each set of samples.'),
              tags$li('The correlation tab allows to study the co-exposure to the chemicals within each samples group and can they be visualized as heatmaps.')
              )
            )
        )
      )
    ),
## GC-MS TAB ## ---------------------------------------------------------------
    navbarMenu("GC-MS",
      tabPanel("Table",
        sidebarLayout(position = "right",
          sidebarPanel(
            tags$h3("Measured Chemicals"),
            plotlyOutput("pieGCgroup")
          ),
          mainPanel(
            dataTableOutput('tableGC')
          )
        )
      ),
      tabPanel("Box-Plot",
        sidebarLayout(position = "left",
          sidebarPanel(
            selectInput("groupGC", "Select Group:", choices = ""),
            selectInput("metablGC", "Select Metabolite:", choices = "")
          ),
          mainPanel(
            plotlyOutput('boxplotGC')
          )
        )
      ),
      tabPanel("Correlation",
         sidebarPanel(
           selectInput("groupGCcircos", "Select Group:", 
             choices = c("All Samples", "Control Samples", "Case Samples",
                         "Mother Samples", "Father Samples"))
         ),
         mainPanel(
           edgebundleOutput("circosGC", width = "100%", height = "600px")
         )
      ),
      tabPanel("Detailed Correlation",
        plotlyOutput("corrGC", height = 600)
      )
    ),
## LC-HRMS/C18 TAB ## ---------------------------------------------------------
    navbarMenu("LC-HRMS/C18",
      tabPanel("Table",
        sidebarLayout(position = "right",
          sidebarPanel(
            tags$h3("Annotated Measured Chemicals"),
            plotlyOutput("pieC18group")
          ),
          mainPanel(
            dataTableOutput('tableC18')
          )
        )
      ),
      tabPanel("Box-Plot",
        sidebarLayout(position = "left",
          sidebarPanel(
            selectInput("groupC18", "Select Group:",
              choices = c("HMDB" = "HMDB", "KEGG" = "KEGG", "LIPID" = "LIPID"),
              selected = "HMDB"),
            selectInput("metabC18", "Select Metabolite:", choices = ""),
            selectInput("mztC18", "Select M/Z and Retention Time:", 
                        choices = "" )
          ),
          mainPanel(
            plotlyOutput('boxplotC18')
          )
        )
      ),
      tabPanel("Correlation",
        sidebarLayout(position = "left",
          sidebarPanel(
            selectInput("corC18set1", "Select Set 1:",
              choices = c("HMDB" = "HMDB", "KEGG" = "KEGG", "LIPID" = "LIPID"),
              selected = "HMDB"),
            selectInput("corC18set2", "Select Set 2:", choices = "")
          ),
          mainPanel(
            plotlyOutput('corrC18')
          )
        )
      ),
      tabPanel("KEGG Pathways",
        sidebarLayout(position = "left",
          sidebarPanel(
            selectInput("keggSelC18", "Select Pathway:", choices = ""),
            selectInput("keggGroupC18", "Select Group:", 
                        choices = c("All Samples", "Control Samples", "Case Samples",
                                    "Mother Samples", "Father Samples"))
          ),
          mainPanel(
            plotlyOutput("keggC18", height = "600px")
          )
        )
      )
    ),
## LC-HRMS/HILIC TAB ## -------------------------------------------------------
    navbarMenu("LC-HRMS/HILIC",
      tabPanel("Table",
        sidebarLayout(position = "right",
          sidebarPanel(
            tags$h3("Annotated Measured Chemicals"),
            plotlyOutput("pieHILICgroup")
          ),
          mainPanel(
            dataTableOutput('tableHILIC')
          )
        )
      ),
      tabPanel("Box-Plot",
        sidebarLayout(position = "left",
          sidebarPanel(
            selectInput("groupHILIC", "Select Group:",
              choices = c("HMDB" = "HMDB", "KEGG" = "KEGG", "LIPID" = "LIPID"),
              selected = "HMDB"),
            selectInput("metabHILIC", "Select Metabolite:", choices = ""),
            selectInput("mztHILIC", "Select M/Z and Retention Time:", 
              choices = "" )
          ),
          mainPanel(
            plotlyOutput('boxplotHILIC')
          )
        )
      ),
      tabPanel("Correlation",
        sidebarLayout(position = "left",
          sidebarPanel(
            selectInput("corHILICset1", "Select Set 1:",
            choices = c("HMDB" = "HMDB", "KEGG" = "KEGG", "LIPID" = "LIPID"),
            selected = "HMDB"),
            selectInput("corHILICset2", "Select Set 2:", choices = "")
          ),
          mainPanel(
            plotlyOutput('corrHILIC')
          )
        )
      ),
      tabPanel("KEGG Pathways",
        sidebarLayout(position = "left",
          sidebarPanel(
            selectInput("keggSelHILIC", "Select Pathway:", choices = ""),
            selectInput("keggGroupHILIC", "Select Group:", 
              choices = c("All Samples", "Control Samples", "Case Samples",
              "Mother Samples", "Father Samples"))
          ),
          mainPanel(
            plotlyOutput("keggHILIC", height = "600px")
          )
        )
      )
    )
  )
)

server <- function(session, input, output, event) {
  ## CREATE GC-MS PIE CHART ###################################################
  output$pieGCgroup <- renderPlotly({ pieGC() })
  ## CREATE C18 PIE CHART ## --------------------------------------------------
  output$pieC18group <- renderPlotly({ pieC18() })
  ## CREATE HILIC PIE CHART ## ------------------------------------------------
  output$pieHILICgroup <- renderPlotly({ pieHILIC() })
  
  
  ## CREATE GC-MS TABLE #######################################################
  dfGC <- reactiveValues(data = {
    rst <- gc_data$table[ , -5]
    rst
  })
  output$tableGC <- renderDataTable(datatable(dfGC$data,
    rownames = FALSE, escape = FALSE, filter = "top", 
    options = list(pageLength = 10), selection = "none")
  )
  ## CREATE C18 TABLE ## ------------------------------------------------------
  dfC18 <- reactiveValues(data = {
    rst <- c18_data$table[ , -6]
    rst
  })
  output$tableC18 <- renderDataTable(datatable(dfC18$data,
    rownames = FALSE, escape = FALSE, filter = "top",
    options = list(pageLength = 10), selection = "none")
  )
  ## CREATE HILIC TABLE ## ----------------------------------------------------
  dfHI <- reactiveValues(data = {
    rst <- hilic_data$table[ , -6]
    rst
  })
  output$tableHILIC <- renderDataTable(datatable(dfHI$data,
    rownames = FALSE, escape = FALSE, filter = "top",
    options = list(pageLength = 10), selection = "none")
  )
  
  ## FILL GC SELECTORS ########################################################
  updateSelectInput(session = session, inputId = "groupGC", 
    choices = familcyGC(), selected = familcyGC()[1])
  observeEvent(input$groupGC,
    updateSelectInput(session = session, inputId = "metablGC", 
      choices = namesGC(input$groupGC))
  )
  ## FILL C18 SELECTORS ## ----------------------------------------------------
  observeEvent(input$groupC18,
    updateSelectInput(session = session, inputId = "metabC18",
      choices = namesC18(input$groupC18))
  )
  observeEvent(input$metabC18, {
    mzt <- mztimeC18(input$metabC18)
    updateSelectInput(session = session, inputId = "mztC18",
      choices = mzt)
  })
  observeEvent(input$corC18set1, {
    opt <- list("HMDB" = c("KEGG", "LIPID"),
                "KEGG" = c("HMDB", "LIPID"),
                "LIPID" = c("HMDB", "LIPID"))
    updateSelectInput(session = session, inputId = "corC18set2",
                      choices = opt[[input$corC18set1]])
  })
  updateSelectInput(session = session, inputId = "keggSelC18",
                    choices = keggNameC18(), selected = keggNameC18()[1])
  ## FILL HILIC SELECTORS ## --------------------------------------------------
  observeEvent(input$groupHILIC,
    updateSelectInput(session = session, inputId = "metabHILIC",
      choices = namesHILIC(input$groupHILIC))
  )
  observeEvent(input$metabHILIC, {
    mzt <- mztimeHILIC(input$metabHILIC)
    updateSelectInput(session = session, inputId = "mztHILIC",
      choices = mzt)
  })
  observeEvent(input$corHILICset1, {
    opt <- list("HMDB" = c("KEGG", "LIPID"),
                "KEGG" = c("HMDB", "LIPID"),
                "LIPID" = c("HMDB", "LIPID"))
    updateSelectInput(session = session, inputId = "corHILICset2",
                      choices = opt[[input$corHILICset1]])
  })
  updateSelectInput(session = session, inputId = "keggSelHILIC",
                    choices = keggNameHILIC(), selected = keggNameHILIC()[1])
  
  
  ## CREATE GC BOXPLOT ########################################################
  observeEvent(input$metablGC,
    output$boxplotGC <- renderPlotly({ 
      boxplotGC(input$metablGC) 
    })
  )
  ## CREATE C18 BOXPLOT ## ----------------------------------------------------
  observeEvent(input$mztC18,
    output$boxplotC18 <- renderPlotly({
      boxplotC18(input$metabC18, input$mztC18)
    })
  )
  # ## CREATE HILIC BOXPLOT ## --------------------------------------------------
  observeEvent(input$mztHILIC,
    output$boxplotHILIC <- renderPlotly({
      boxplotHILIC(input$metabHILIC, input$mztHILIC)
    })
  )
  
  ## CREATE GC-MS CORRELATION HEATMAP #########################################
  output$corrGC <- renderPlotly({
    correlationGC()
  })
  observeEvent(input$groupGCcircos,
    output$circosGC <- renderEdgebundle({
      xx <- getCorGR(input$groupGCcircos)
      edgebundle(xx)
    })
  )
  ## CREATE C18 CORRELATION HEATMAP ## ----------------------------------------
  output$corrC18 <- renderPlotly({
    correlationC18(input$corC18set1, input$corC18set2)
  })
  output$keggC18 <- renderPlotly(
    getKEGGcorC18(keggNameC18()[1], input$keggGroupC18)
  )
  observeEvent(input$keggSelC18,
    output$keggC18 <- renderPlotly(
      getKEGGcorC18(input$keggSelC18, input$keggGroupC18)
    )
  )
  ## CREATE HILIC CORRELATION HEATMAP ## --------------------------------------
  output$corrHILIC <- renderPlotly({
    correlationHILIC(input$corHILICset1, input$corHILICset2)
  })
  output$keggHILIC <- renderPlotly(
    getKEGGcorC18(keggNameHILIC()[1], input$keggGroupHILIC)
  )
  observeEvent(input$keggSelHILIC,
     output$keggHILIC <- renderPlotly(
       getKEGGcorC18(input$keggSelHILIC, input$keggGroupHILIC)
     )
  )
}



shinyApp(ui = ui, server = server)

