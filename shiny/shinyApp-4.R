#### Install packages not yet installed ####

cran_packages <- c("rlang", "shiny", "DT", "ggplot2", "RColorBrewer", "heatmaply", "plotly", "reshape", "tidyverse","bslib", "scales", "ggrepel")
cran_installed_packages <- cran_packages %in% rownames(installed.packages())
bioconductor_packages <- c("limma", "edgeR")
bioconductor_installed_packages <- bioconductor_packages %in% rownames(installed.packages())
if (any(cran_installed_packages == FALSE)) {
  install.packages(cran_packages[!cran_installed_packages])
}

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager") 
  }

if (any(bioconductor_installed_packages == FALSE)) {
  
  suppressMessages({
    BiocManager::install(bioconductor_packages[!bioconductor_installed_packages], quiet = TRUE)
  })
}


#### Packages loading ####
invisible(lapply(cran_packages, library, character.only = TRUE))
invisible(lapply(bioconductor_packages, library, character.only = TRUE))

source("functions.R")

feature_style='border-left: 0.75px solid #f0f0f0;'
title_style='border-left: 0.75px solid #f0f0f0;border-bottom: 0.75px solid #f0f0f0;'
space1_style='height:40px'
space2_style='height:20px;border-left: 0.75px solid #f0f0f0;'
space3_style='height:20px;border-left: 0.75px solid #f0f0f0;border-top: 0.75px solid #f0f0f0;'

ui <- fluidPage(
  theme = bs_theme(bootswatch = "minty", base_font = font_google("Roboto")
                    ),
  tags$head(
    tags$style(
      HTML("
      @import url('https://fonts.googleapis.com/css2?family=Roboto');
        body {font-family: 'Roboto';},
        h1 {font-family: 'Roboto';}
        h2 {font-family: 'Roboto';}
        h5 {font-family: 'Roboto';}
        strong {font-family: 'Roboto';}",
        "#inTabset {
          position:fixed;
          width: 100%;
          height: 4.7%;
          background-color: white;
          margin-top: -30px;
          z-index:3;
          }", 
        ".tab-content {
        margin-top: 25px;
        z-index:1;
        }"
        )
      )
    ),
  
  div(
    titlePanel(
      h2("RNA interference and CRISPR/Cas9 screen results")
      ),
    style = "background-color:#ffffff;position:fixed; z-index: 1; top: -10px;width:100%;height:80px;"),
  br(), 
  br(),
  br(),
  sidebarLayout(
    sidebarPanel(
      style = "position:fixed;width:32%;background-color:#dfede7;",
      fileInput('files', strong('Upload workflow outputs'), accept = c(".rds"), multiple = T),
      selectInput("contrast", "Contrast", choices = ""),
      numericInput("FDR", "FDR threshold", value = 0.05, min = 0, max = 1, step = 0.01),
      numericInput("logFC", "LogFC threshold", value = 1.5, min = 0, max = 15, step = 0.5), 
      selectInput("export_gene", strong("Select gene to generate report on:"), choices = ""),
      downloadButton("report", "Generate report")
    ),
    mainPanel(
      tabsetPanel(id = "inTabset",
                  type = "tabs",
                  tabPanel(icon("house-user"),
                           fluidRow(column(12, style = 'height:20px')),
                           h5("Welcome!"),
                           h6("This site allows the exploration and visualization of results from RNA interference and CRISPR/Cas9 screens."),
                           h6("Once you have successfully run the snakemake workflow on your data, upload the workflow output where indicated on the left panel of this site. 
                             This workflow output can be found in the snakemake results directory, named “shiny.rds”."),
                           h6("You can set a threshold for false discovery rate (FDR) and log fold-change (FC) where indicated on the left, then explore the findings through the tabs above."),
                           h6("You can also generate a report on a gene of interest in your data.")
                           ),
                  tabPanel("Quality control",
                           fluidRow(column(12, style = 'height:20px')),
                           fluidRow(column(12, h5("This tab displays visualises quality control features"))),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, plotlyOutput('indexcounts'), style = feature_style)),
                           fluidRow(column(12, strong("Guide RNA counts per sample. Hover over each bar to get specific counts per sample"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, plotlyOutput('shcounts'), style = feature_style)),
                           fluidRow(column(12, strong("Density of guide RNA counts across all samples"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, plotlyOutput('bcv'), style = feature_style)),
                           fluidRow(column(12, strong("Variability in the screen as a functon of guide RNA abundance"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(6, plotlyOutput('mds'), style = feature_style),
                                    column(6, plotlyOutput('correctedmds'))),
                           fluidRow(column(12, strong("Multidimensional plot to visualise the relationships between samples, A) uncorrected and B) corrected for batch"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(6, plotlyOutput('pca'), style = feature_style),
                                    column(6, plotlyOutput('correctedpca'))),
                           fluidRow(column(12, strong("First two principal components to visualise the relationships between samples, A) uncorrected and B) corrected for batch"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(6, plotlyOutput('sampledist'), style = feature_style),
                                    column(6, plotlyOutput('correctedsampledist'))),
                           fluidRow(column(12, strong("Distances between samples, A) uncorrected and B) corrected for batch"), style = title_style)),
                           fluidRow(column(12, style = space1_style))
                           ),
                  tabPanel("Differential expression",
                           fluidRow(column(12, style = 'height:20px')),
                           fluidRow(column(12, h5("This tab displays differential expression across guide RNAs"))),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, plotlyOutput('hist'), style = feature_style)),
                           fluidRow(column(12, strong("Histogram of guide RNA p values, hover over bars to get further details"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                        
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(7, plotlyOutput('plotsmear'), style = feature_style),
                                    column(5, plotlyOutput('volcano'))),
                           fluidRow(column(7, strong("Log fold-change (FC) against log counts-per-million (CPM) of the guide RNAs. 
                                                 Hover over points to get specific information on each guide RNA"), style = title_style),
                                    column(5, strong("Volcano plot of guide RNA, hover over points to get specific information on each guide RNA"))),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(4, numericInput("topguides", "Select number of top guides", value = 5, min = 0, max = 100), style = feature_style)),
                           fluidRow(column(12, style = space2_style)),
                           fluidRow(column(6, plotlyOutput('de'), style = feature_style),
                                    column(6, plotlyOutput('correctedde'))),
                           fluidRow(column(12, strong("Differential expression (logCPM) across samples for the top guide RNAs, A) uncorrected and B) corrected for batch"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, selectInput("corrected", "Select uncorrected or corrected to view in the below table:", 
                                                           choices=c("Uncorrected", "Corrected")), style = feature_style)),
                           fluidRow(column(12, style = space2_style)),
                           fluidRow(column(12, DTOutput('detable'), style = feature_style)),
                           fluidRow(column(12, strong("Table of differential expression (logCPM) across samples for each guide RNA."), style = title_style)),
                           fluidRow(column(12, style = space1_style))
                           ),
                  tabPanel("Gene level",
                           fluidRow(column(12, style = 'height:20px')),
                           fluidRow(column(12, h5("This tab displays results at the gene level"))),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, DTOutput('camera'), style = feature_style)),
                           fluidRow(column(12, strong("camera"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, plotlyOutput('camerarank'), style = feature_style)),
                           fluidRow(column(12, strong("camera rank"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, DTOutput('genelevel'), style = feature_style)),
                           fluidRow(column(12, strong("genelevel"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, plotlyOutput('generank'), style = feature_style)),
                           fluidRow(column(12, strong("gene rank"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(6, selectInput("barcode_gene", "Select Gene", choices = ""), style = feature_style)),
                           fluidRow(column(12, style = space2_style)),
                           fluidRow(column(12, plotOutput('barcode'), style = feature_style)),
                           fluidRow(column(12, strong("barcode plot"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           ),
                  tabPanel("Enrichment",
                           fluidRow(column(12, style = 'height:20px')),
                           fluidRow(column(12, h5("This tab displays enrichment for gene ontology and KEGG pathways"))),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, DTOutput('gotable'), style = feature_style)),
                           fluidRow(column(12, strong("GO table"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(4, numericInput("topgos", "Select number of top terms to view"
                                                           , value = 10, min = 0, max =100), style = feature_style)),
                           fluidRow(column(12, style = space2_style)),
                           fluidRow(column(6, plotlyOutput("upgoplot"), style = 'height:600px;border-left: 0.75px solid #f0f0f0;'),
                                    column(6, plotlyOutput("downgoplot"), style = 'height:600px;')),
                           fluidRow(column(12, style = 'height:20px')),
                           fluidRow(column(12, strong("GO plots"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, DTOutput("keggtable"), style = feature_style)),
                           fluidRow(column(12, strong("KEGG table"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(4, numericInput("topkeggs", "Select number of top pathways to view"
                                                           , value = 10, min = 0, max = 100), style = feature_style)),
                           fluidRow(column(12, style = space2_style)),
                           fluidRow(column(6, plotlyOutput("upkeggplot"), style = 'height:600px;border-left: 0.75px solid #f0f0f0;'),
                                    column(6, plotlyOutput("downkeggplot"), style = 'height:600px;')),
                           fluidRow(column(12, style = 'height:20px')),
                           fluidRow(column(12, strong("KEGG plots"), style = title_style)),
                           fluidRow(column(12, style = space1_style))
                           ),
                  tabPanel("Comparing contrasts",
                           fluidRow(column(12, style = 'height:20px')),
                           fluidRow(column(12, h5("This tab compares results between contrasts"))),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(4, style = 'height:20px', style = feature_style),
                                    column(4, selectInput("contrast_gene", "Select Gene", choices=""), style = 'height:20px'),
                                    column(4, style = 'height:20px;')),
                           fluidRow(column(12, style = space2_style)),
                           fluidRow(column(12, plotlyOutput('guidecontrast'), style = feature_style)),
                           fluidRow(column(12, strong("Comparing guide RNAs between contrasts"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, DTOutput('FCcontrast'), style = feature_style)),
                           fluidRow(column(12, strong("Table of mean log fold-change of genes across contrasts"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           )
      )
    )
  )
)

server <- function(input, output, session) {
  data <- reactive({
    file <- input$files
    ext <- tools::file_ext(file$datapath[-(grep("-shiny.rds", input$files$name))])
    req(file)
    validate(need(ext == "rds", "Please upload a rds file"))
    readRDS(file$datapath[-(grep("-shiny.rds", input$files$name))])
  })
  contrast_data <- reactive({
    file <- input$files
    ext <- tools::file_ext(file$datapath[(grep("-shiny.rds", input$files$name))])
    req(file)
    validate(need(ext == "rds", "Please upload a rds file"))
    dat <- list()
    for (i in (grep("-shiny.rds", input$files$name))) {
      input <- readRDS(file$datapath[(i)])
      dat <- c(dat, input)
    }
    dat
  })
  contrasts <- reactive({
    name <- names(contrast_data())[which(str_detect(names(contrast_data()), "_lrt"))]
    split <- str_split(name, "_")
    as.vector(sapply(split,"[[",1))
  })
  observe({
    updateSelectInput(session,
                      "contrast",
                      choices = contrasts())
  })
  inputcontrast <- reactive({paste(input$contrast)})
  FCthres <- reactive({as.numeric(input$logFC)})
  FDRthres <- reactive({as.numeric(input$FDR)})
  topguides <- reactive({as.numeric(input$topguides)})
  topgos <- reactive({as.numeric(input$topgos)})
  topkeggs <- reactive({as.numeric(input$topkeggs)})
  
  ### Quality control tab ###
  output$indexcounts <- renderPlotly({
    indexcounts(data())
  })
  output$shcounts <- renderPlotly({
    guidecounts(data())
  })
  output$bcv <- renderPlotly({
    bcv(data())
  })
  output$mds <- renderPlotly({
    mds(data())
  })
  output$correctedmds <- renderPlotly({
    cor_mds(data())
  })
  output$pca <- renderPlotly({
    pca(data())
  })
  output$correctedpca <- renderPlotly({
    cor_pca(data())
  })
  output$sampledist <- renderPlotly({
    sampledist(data())
  })
  output$correctedsampledist <- renderPlotly({
    cor_sampledist(data())
  })
  
  ### Differential expression tab ###
  
  output$hist <- renderPlotly({
    hist(contrast_data(), inputcontrast())
  })
  output$plotsmear <- renderPlotly({
    plotsmear(contrast_data(), inputcontrast(), FCthres(), FDRthres())
  })
  output$volcano <- renderPlotly({
    volcano(contrast_data(), inputcontrast(), FDRthres())
  })
  output$de <- renderPlotly({
    de(data(), contrast_data(), inputcontrast(), topguides())
  })
  output$correctedde <- renderPlotly({
    cor_de(data(), contrast_data(), inputcontrast(), topguides())
  })
  output$detable <- renderDT({
    datatable(detable(data(), inputcontrast(), input$corrected))
  })
  ### Gene level tab ###
  output$camera <- renderDT({
    datatable(camera(contrast_data(), inputcontrast()))
  })
  output$camerarank <- renderPlotly({
    s <- input$camera_rows_selected
    camerarank(contrast_data(), inputcontrast(),FDRthres(), s)
  })
  output$genelevel <- renderDT({
    datatable(genelevel(contrast_data(), inputcontrast()))
  })
  output$generank <- renderPlotly({
    s <- input$genelevel_rows_selected
    generank(contrast_data(), inputcontrast(), FCthres(), s)
  })
  observe({
    updateSelectInput(session,
                      "barcode_gene",
                      choices = c(unique(data()$x$genes$Gene)))
  })
  barcodegene <- reactive({
    paste(input$barcode_gene)
  })
  output$barcode<- renderPlot({
    barcode(data(), contrast_data(), inputcontrast(), barcodegene())
  })
  
  ### Enrichment tab ###
  
  output$gotable <- renderDT({
    datatable(gotable(contrast_data(), inputcontrast()))
  })
  output$upgoplot <- renderPlotly({
    upgoplot(contrast_data(), inputcontrast(), topgos())
  })
  output$downgoplot <- renderPlotly({
    downgoplot(contrast_data(), inputcontrast(), topgos())
    
  })
  output$keggtable <- renderDT({
    datatable(keggtable(contrast_data(), inputcontrast()))
  })
  output$upkeggplot <- renderPlotly({
    upkeggplot(contrast_data(), inputcontrast(), topkeggs())
  })
  output$downkeggplot <- renderPlotly({
    downkeggplot(contrast_data(), inputcontrast(), topkeggs())
  })
  
  ### Comparing contrasts tab ###
  observe({
    updateSelectInput(session,
                      "contrast_gene",
                      choices = c(unique(data()$x$genes$Gene)))
  })
  contrastgene <- reactive({
    paste(input$contrast_gene)
  })
  output$guidecontrast <- renderPlotly({
    guidecontrast(contrast_data(), contrastgene())
  })

  output$FCcontrast <- renderDT({
    datatable(FCcontrast(data(), contrast_data(), inputcontrast()),
              selection = "single")
  })
  observe({
    updateSelectInput(session,
                      "export_gene",
                      choices = c(unique(data()$x$genes$Gene)))
  })
  exportgene <- reactive({
    paste(input$export_gene)
  })
  output$report <- downloadHandler(
    filename <- paste0(exportgene(), ".html"),
    content <- function(file) {
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      params <- base::list(gene = exportgene(),
                           data=data(),
                           contrast_data=contrast_data())
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
}

shinyApp(ui,server)
