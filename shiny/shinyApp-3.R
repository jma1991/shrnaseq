#### Install packages not yet installed ####

cran_packages <- c("rlang", "shiny", "DT", "ggplot2", "RColorBrewer", "heatmaply", "plotly", "reshape", "tidyverse","bslib")
cran_installed_packages <- cran_packages %in% rownames(installed.packages())
bioconductor_packages=c("limma", "edgeR")
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

ui <- fluidPage(
  tags$head(
    # Note the wrapping of the string in HTML()
    tags$style(HTML(
    "@import url('https://fonts.googleapis.com/css2?family=Roboto');
      h1 {
        font-family: 'Roboto';
      }")
    )
  ),
  
  div(titlePanel( h1("RNA interference and CRISPR/Cas9 screen results")), 
      style = "background-color:#ffffff;position:fixed; z-index: 1; top: -10px;width:100%;"),
  theme =  bs_theme(bootswatch = "minty",
                    base_font = font_google("Roboto")),
  br(), br(), br(),
  sidebarLayout(
    sidebarPanel(
      style = "position:fixed;width:32%;background-color:#dfede7;",
      selectInput("contrast", "Contrast", choices=""),
      numericInput("FDR", "FDR threshold", value=0.05, min=0, max=1),
      numericInput("logFC", "LogFC threshold", value=1.5, min=0, max=100)
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Upload",  
                           fluidRow(column(4, fileInput('files', 'Upload workflow outputs', accept=c(".rds"), multiple=T))
                                    )
                           ),
                  tabPanel("Quality control",
                           fluidRow(column(12,plotlyOutput('indexcounts'),style='padding:20px')
                                    ),
                           fluidRow(column(12,plotlyOutput('shcounts'),style='padding:20px')
                                    ),
                           fluidRow(column(12,plotlyOutput('bcv'),style='padding:20px')
                                    ),
                           fluidRow(column(6,plotlyOutput('mds'),style='padding:20px'),
                                    column(6,plotlyOutput('correctedmds'),style='padding:20px')
                                    ),
                           fluidRow(column(6,plotlyOutput('pca'),style='padding:20px'),
                                    column(6,plotlyOutput('correctedpca'),style='padding:20px')
                                    ),
                           fluidRow(column(6,plotlyOutput('sampledist'),style='padding:20px'),
                                    column(6,plotlyOutput('correctedsampledist'),style='padding:20px')
                                    )
                           ),
                  tabPanel("Differential expression",
                           fluidRow(column(12,plotlyOutput('hist'),style='padding:20px')
                                    ),
                           fluidRow(column(6,plotlyOutput('plotsmear'),style='padding:30px'),
                                    column(6,plotlyOutput('volcano'),style='padding:30px')
                                    ),
                           fluidRow(column(4, numericInput("topguides", "Select number of top guides", value=5, min=0, max=100), style='padding:40px;height:20px')
                                    ),
                           fluidRow(column(6,plotlyOutput('de'),style='padding:40px'),
                                    column(6,plotlyOutput('correctedde'),style='padding:40px')
                                    ),
                           
                           fluidRow(column(12, selectInput("corrected", "Select uncorrected or corrected to view in the below table:", choices=c("Uncorrected", "Corrected")))
                                    ),
                           fluidRow(column(12)
                                    ),
                           fluidRow(column(12, DTOutput('detable'))
                           )),
                  tabPanel("Gene level",
                           fluidRow(column(12, DTOutput('camera'),style='padding:20px')
                                    ),
                           fluidRow(column(12, plotlyOutput('camerarank'),style='padding:20px;height:500px')
                                    ),
                           fluidRow(column(12, DTOutput('genelevel'),style='padding:20px')
                                    ),
                           fluidRow(column(12, plotlyOutput('generank'),style='padding:20px;height:500px')
                                    ),
                           fluidRow(column(12, style='padding:20px;')
                                    ),
                           fluidRow(column(4, selectInput("barcode_gene", "Select Gene", choices = ""), style='padding:20px;height:20px')
                                    ),
                           fluidRow(column(12, plotOutput('barcode'),style='padding:20px')
                                    )
                           ),
                  tabPanel("Enrichment",
                           fluidRow(column(12, DTOutput('gotable'), style="padding:20px")
                                    ),
                           fluidRow(column(12, style='height:20px')
                                    ),
                           fluidRow(column(4, numericInput("topgos", "Select number of top terms to view", value=10, min=0, max=100), style='padding:40px;height:20px')
                                    ),
                           fluidRow(column(6, plotlyOutput("upgoplot"),style='padding:40px;height:600px'),
                                    column(6, plotlyOutput("downgoplot"), style='padding:40px;height:600px')
                                    ),
                           fluidRow(column(12, style='height:20px')
                                    ),
                           fluidRow(column(12, DTOutput("keggtable"),style='padding:40px')
                                    ),
                           fluidRow(column(12, style='height:20px')
                                    ),
                           fluidRow(column(4, numericInput("topkeggs", "Select number of top pathways to view", value=10, min=0, max=100), style='padding:40px;height:20px')
                                    ),
                           fluidRow(column(6, plotlyOutput("upkeggplot"),style='padding:40px;height:600px'),
                                    column(6, plotlyOutput("downkeggplot"), style="padding:40px;height:600px")
                                    ),                           
                           fluidRow(column(12, style='height:20px')
                                    )
                           ),
                  tabPanel("Comparing contrasts",
                           fluidRow(column(4, style='padding:30px;height:20px'),
                                    column(4, selectInput("contrast_gene", "Select Gene", choices=""), style='padding:30px;height:20px'),
                                    column(4, style='padding:30px;height:20px')
                                    ),
                           fluidRow(column(12, style='padding:30px;height:20px')),
                           fluidRow(column(12, plotlyOutput('guidecontrast'))
                                    ), 
                           fluidRow(column(12, h6("The plot above will not display genes with one guide RNA."))
                                    ),
                           fluidRow(column(12, style='padding:30px;height:20px')),
                           fluidRow(column(12, DTOutput('FCcontrast'), style='padding:20px')
                                    )
                           ),
                  tabPanel("Rmarkdown export",
                           fluidRow(column(6, selectInput("export_gene", "Select Gene", choices=""))
                           ),
                           fluidRow(column(6, downloadButton("report", "Generate report"))
                           )
                  )
      )
    )
  )
)

server <- function(input, output, session) {
  data=reactive({
    file <- input$files
    ext <- tools::file_ext(file$datapath[-(grep("-shiny.rds", input$files$name))])
    
    req(file)
    validate(need(ext == "rds", "Please upload a rds file"))
    
    readRDS(file$datapath[-(grep("-shiny.rds", input$files$name))])
  })
  contrast_data=reactive({
    file <- input$files
    ext <- tools::file_ext(file$datapath[(grep("-shiny.rds", input$files$name))])
    
    req(file)
    validate(need(ext == "rds", "Please upload a rds file"))
    dat=list()
    for (i in (grep("-shiny.rds", input$files$name))) {
      input=readRDS(file$datapath[(i)])
      dat=c(dat, input)
    }
    dat
  })
  contrasts=reactive({
    name=names(contrast_data())[which(str_detect(names(contrast_data()), "_lrt"))]
    split=str_split(name, "_")
    as.vector(sapply(split,"[[",1))
  })
  observe({
    updateSelectInput(session,
                      "contrast",
                      choices = contrasts())
  })
  inputcontrast=reactive({paste(input$contrast)})
  FCthres=reactive({as.numeric(input$logFC)})
  FDRthres=reactive({as.numeric(input$FDR)})
  topguides=reactive({as.numeric(input$topguides)})
  topgos=reactive({as.numeric(input$topgos)})
  topkeggs=reactive({as.numeric(input$topkeggs)})
  genes<-reactive({getgenes(data())})
  
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
    volcano(contrast_data(), inputcontrast())
  })
  output$de <- renderPlotly({
    de(data(), contrast_data(), inputcontrast(), topguides())
  })
  output$correctedde <- renderPlotly({
    cor_de(data(), contrast_data(), inputcontrast(), topguides())
  })
  output$detable <- renderDT({
    datatable(detable(data(), inputcontrast(), input$corrected), 
              
              caption = htmltools::tags$caption(style = 'caption-side: top; color:black; text-align: center; font-size:125% ;','Differential expression (logCPM)')
    )
  })
  ### Gene level tab ###
  output$camera <- renderDT({
    datatable(camera(contrast_data(), inputcontrast()), 
              caption = htmltools::tags$caption(style = 'caption-side: top; color:black; text-align: center; font-size:125% ;','Camera')
    )
  })
  output$camerarank <- renderPlotly({
    s = input$camera_rows_selected
    camerarank(contrast_data(), inputcontrast(), s)
  })
  output$genelevel <- renderDT({
    datatable(genelevel(contrast_data(), inputcontrast()),
              caption = htmltools::tags$caption( style = 'caption-side: top; color:black; text-align: center; font-size:125% ;','Gene level')
              
    )
  })
  output$generank <- renderPlotly({
    s = input$genelevel_rows_selected
    generank(contrast_data(), inputcontrast(), s)
  })
  observe({
    updateSelectInput(session,
                      "barcode_gene",
                      choices = c(genes()))
  })
  barcodegene=reactive({
    paste(input$barcode_gene)
  })
  output$barcode<- renderPlot({
    barcode(data(), contrast_data(), inputcontrast(), barcodegene())
  })
  
  ### Enrichment tab ###
  
  output$gotable <- renderDT({
    datatable(gotable(contrast_data(), inputcontrast()),
              
              caption = htmltools::tags$caption(style = 'caption-side: top; color:black; text-align: center; font-size:125% ;','Gene ontology')
              
    )
  })
  output$upgoplot <- renderPlotly({
    upgoplot(contrast_data(), inputcontrast(), topgos())
  })
  output$downgoplot <- renderPlotly({
    downgoplot(contrast_data(), inputcontrast(), topgos())
    
  })
  output$keggtable <- renderDT({
    datatable(keggtable(contrast_data(), inputcontrast()),
              
              caption = htmltools::tags$caption( style = 'caption-side: top; color:black; text-align: center; font-size:125% ;','KEGG pathway enrichment')
              
    )
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
                      choices = c(genes()))
  })
  contrastgene=reactive({
    paste(input$contrast_gene)
  })
  output$guidecontrast <- renderPlotly({
    guidecontrast(contrast_data(), contrastgene())
  })

  output$FCcontrast <- renderDT({
    datatable(FCcontrast(data(), contrast_data(), inputcontrast()),
              selection = "single",
              caption = htmltools::tags$caption( style = 'caption-side: top; color:black; text-align: center; font-size:125% ;','Mean logFC across contrasts')
              
    )
  })
  observe({
    updateSelectInput(session,
                      "export_gene",
                      choices = c(genes()))
  })
  exportgene=reactive({
    paste(input$export_gene)
  })
  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = paste0(exportgene(), ".html"),
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- base::list(gene = exportgene(),
                           data=data(),
                           contrast_data=contrast_data())
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
}

shinyApp(ui,server)
