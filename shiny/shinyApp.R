library(shiny)
library(DT)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(heatmaply)
library(plotly)
library(reshape)
library(limma)
library(tidyverse)

#### Functions ####
indexcounts=function(list) {
  
  df=melt(data.frame(t(colSums(list$x$counts))))
  ggplotly(
    ggplot(df, aes(x=variable, y=value)) +  
      geom_bar(stat = "identity", fill="cornflowerblue") +
      theme_minimal() +
      labs(title="Counts per sample", x="", y = "Counts") +
      theme(plot.title = element_text(hjust = 0.5))
  )
}
guidecounts=function(list) {
  df=melt(data.frame(t(rowSums(list$x$counts))))
  ggplotly(
    ggplot(df, aes(x=variable, y=value)) +  
      geom_bar(stat = "identity", fill="cornflowerblue") +
      theme_minimal() +
      labs(title="Counts per hairpin", x="", y = "Counts") +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.text.x = element_text(angle = 45))
  )
}
bcv=function(list) {
  df=data.frame(cbind(list$xglm$AveLogCPM, 
                      list$xglm$trended.dispersion, 
                      list$xglm$common.dispersion, 
                      list$xglm$tagwise.dispersion))
  colnames(df)=c("AveLogCPM", "trended.dispersion",
                 "commom.dispersion", "tagwise.dispersion")
  ggplotly(
    ggplot(df, aes(x=AveLogCPM)) +
      geom_line(aes(y=trended.dispersion), colour = "cornflowerblue") + 
      geom_line(aes(y=commom.dispersion), colour = "red") + 
      geom_point(aes(y=tagwise.dispersion), colour = "cornflowerblue") + 
      theme_classic() + 
      labs(title="BCV plot", x="Average log CPM", y = "Biological coefficient of variation") +
      theme(plot.title = element_text(hjust = 0.5))
  )
}
mds=function(list) {
  colour=rainbow(length(unique(list$x$samples$group)))
  plotMDS(list$x, labels=list$x$samples$group, col=colour,cex=0.8, main="MDS Plot")
  legend(par("usr")[2], par("usr")[4],legend=c(unique(list$x$samples$group)), col=colour, pch=16, box.lwd = 0,box.col = "white",bg = "white")
  
}
cor_mds=function(list) {
  colour=rainbow(length(unique(list$x$samples$group)))
  par(xpd = TRUE, mar = par()$mar + c(0, 0, 0, 5))
  plotMDS(list$corrected, labels=list$x$samples$group, col=colour,cex=0.8, main="Batch corrected MDS Plot")
  legend(par("usr")[2], par("usr")[4],legend=c(unique(list$x$samples$group)), 
         col=colour, pch=16, box.lwd = 0 ,box.col = "white",bg = "white")
}
pca=function(list) {
  list$x$samples$group=factor(list$x$samples$group, levels = c(unique(list$x$samples$group)))
  
  mat=cpm(list$x$counts, log=TRUE, prior.count = 1)
  var <- matrixStats::rowVars(mat)
  num <- min(500, length(var))
  ind <- order(var, decreasing = TRUE)[seq_len(num)]
  pca <- prcomp(t(mat[ind, ]))
  pct <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)
  dat <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], group = list$x$samples$group)
  
  if (length(unique(list$x$samples$group))<3) {
    n=3 
  } else {
    n=length(unique(list$x$samples$group))
  }
  col=brewer.pal(n = n, name = "Paired")
  ggplotly(
    ggplot(dat, aes_string(x = "PC1", y = "PC2", color = "group")) + 
      geom_point(size = 3) + 
      xlab(paste0("PC1: ", round(pct[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(pct[2] * 100), "% variance")) + 
      coord_fixed() +
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      labs(title="Principal component analysis") + 
      theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position="bottom") +
      scale_color_manual(values=col)
  )
}
cor_pca=function(list) {
  var <- matrixStats::rowVars(list$corrected)
  num <- min(500, length(var))
  ind <- order(var, decreasing = TRUE)[seq_len(num)]
  pca <- prcomp(t(list$corrected[ind, ]))
  pct <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)
  dat <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], group = list$x$samples$group)
  
  if (length(unique(list$x$samples$group))<3) {
    n=3 
  } else {
    n=length(unique(list$x$samples$group))
  }
  col=brewer.pal(n = n, name = "Paired")
  ggplotly(
    ggplot(dat, aes_string(x = "PC1", y = "PC2", color = "group")) + 
      geom_point(size = 3) + 
      xlab(paste0("PC1: ", round(pct[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(pct[2] * 100), "% variance")) + 
      coord_fixed() +
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      labs(title="Batch corrected principal component analysis") +
      theme(plot.title = element_text(hjust = 0.5)) +  theme(legend.position="bottom") +
      scale_color_manual(values=col)
  )
}
sampledist=function(list) {
  mat=cpm(list$x$counts, log=TRUE, prior.count = 1)
  sampleDists <- dist(t(mat))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(list$x$samples$group, list$x$samples$Replicate, sep = " - " )
  colnames(sampleDistMatrix) <- paste(list$x$samples$group, list$x$samples$Replicate, sep = " - " )
  col <- colorRampPalette( rev(brewer.pal(5, "Blues")) )(255)
  heatmaply(sampleDistMatrix, col=col, main = "Sample distances")
}
cor_sampledist=function(list) {
  sampleDists <- dist(t(list$corrected))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(list$x$samples$group, list$x$samples$Replicate, sep = " - " )
  colnames(sampleDistMatrix) <- paste(list$x$samples$group, list$x$samples$Replicate, sep = " - " )
  col <- colorRampPalette( rev(brewer.pal(5, "Blues")) )(255)
  heatmaply(sampleDistMatrix, col=col, main = "Batch corrected sample distances")
}
plotsmear=function(list, string, FCthres, FDRthres) {
  obj <- get('list')[[paste0(string, "_lrt")]]
  top2 = topTags(obj, n=Inf)
  top2ids = top2$table[top2$table$FDR<FDRthres,1]
  df=data.frame(obj$table)
  ggplotly(
    ggplot(df, aes(x=logCPM, y=logFC)) +geom_point() +
      geom_point(data = df %>% filter(row.names(df) %in% top2ids), color = "red") +
      geom_hline(yintercept=(-(FCthres)), linetype="dashed", color="green") +
      geom_hline(yintercept=0,  color="cornflowerblue") +
      geom_hline(yintercept=(FCthres), linetype="dashed", color="green") +
      theme_bw() + labs(title="Mean-difference plot") + theme(plot.title = element_text(hjust = 0.5))
  )
}
volcano=function(list, string) {
  obj <- get('list')[[paste0(string, "_lrt")]]
  df=data.frame(obj$table)
  ggplotly(
    ggplot(df, aes(x=logFC, y=-10*log10(PValue))) + geom_point() +
      geom_point() +
      theme_bw() + 
      labs(title="Volcano plot", x="M", y = "-10*log(P-value)") +
      theme(plot.title = element_text(hjust = 0.5)) 
  )
}
de=function(list, string, FCthres) {
  obj <- get('list')[[paste0(string, "_lrt")]]
  top2 = topTags(obj, n=Inf)
  selY <- rownames(top2$table)[abs(top2$table$logFC)>FCthres]
  
  mat <- cpm(list$x$counts, log=TRUE, prior.count = 1)
  colnames(mat) <- paste(list$x$samples$group, list$x$samples$Replicate, sep = " - " )
  mat = subset(mat, rownames(mat) %in% selY)
  col <- colorRampPalette( rev(brewer.pal(5, "Blues")) )(255)
  heatmaply(mat, col=col, main="Differential expression" )
}
cor_de=function(list, string, FCthres) {
  obj <- get('list')[[paste0(string, "_lrt")]]
  top2 = topTags(obj, n=Inf)
  selY <- rownames(top2$table)[abs(top2$table$logFC)>FCthres]
  
  colnames(list$corrected)= paste(list$x$samples$group, list$x$samples$Replicate, sep = " - " )
  corrected = subset(list$corrected, rownames(list$corrected) %in% selY)
  col <- colorRampPalette( rev(brewer.pal(5, "Blues")) )(255)
  heatmaply(corrected, col=col, main="Batch corrected differential expression")
}
hist=function(list, string) {
  obj <- get('list')[[paste0(string, "_lrt")]]
  df=data.frame(obj$table)
  ggplotly(
    ggplot(df, aes(x=PValue)) + 
      geom_histogram(bins=45,color="white",fill="cornflowerblue") +
      theme_classic() +
      labs(title="Histogram of hairpin P values", x="Hairpin p values", y = "Frequency") +
      theme(plot.title = element_text(hjust = 0.5)) 
  )
}
camera=function(list, string) {
  obj <- get('list')[[paste0(string, "_camera")]]
  colnames(obj)=c("nGuides", "Direction", "Pvalue", "FDR", "Gene")
  obj
}
genelevel=function(list, string) {
  obj <- get('list')[[paste0(string, "_genelevel")]]
  obj
}
guiderank=function(list, string, FCthres) {
  obj <- get('list')[[paste0(string, "_lrt")]]
  dat=data.frame(obj$table)
  dat$Rank = rank(dat$logFC)
  dat$guide=rownames(dat)
  ggplotly(
    ggplot(dat, aes(x=Rank, y=logFC, label=guide)) +
      geom_point(color="gray") +
      geom_point(data = dat %>% filter(dat$logFC>1), color = "red") +
      geom_point(data = dat %>% filter(dat$logFC<(-1)), color = "cornflowerblue") +
      geom_text(data = dat %>% filter(dat$logFC>FCthres), check_overlap = TRUE, nudge_x = max(dat$Rank)*0.05, size=2) +
      geom_text(data = dat %>% filter(dat$logFC<(-FCthres)), check_overlap = TRUE, nudge_x = max(dat$Rank)*0.05, size=2) +
      theme_classic() +
      labs(title="Guide rank plot") +
      theme(plot.title = element_text(hjust = 0.5)) 
    , height = 500
  )
}
getgenes=function(list) {
  genesymbollist = list()
  genesymbols=as.character(list$x$genes$Gene)
  unq = unique(genesymbols)
  unq = unq[!is.na(unq)]
  genes=NULL
  for (i in unq) {
    sel = genesymbols == i & !is.na(genesymbols)
    if (sum(sel) > 1) 
     genes=rbind(genes, i)
  }
  as.vector(genes)
}
barcode=function(list, string, gene) {
  obj <- get('list')[[paste0(string, "_lrt")]]
  genesymbollist = list()
  genesymbols=as.character(list$x$genes$Gene)
  unq = unique(genesymbols)
  unq = unq[!is.na(unq)]
  
  for (i in unq) {
    sel = genesymbols == i & !is.na(genesymbols)
    if (sum(sel) > 1) 
      genesymbollist[[i]] =which(sel)
  }
  barcodeplot(obj$table$logFC,index=genesymbollist[[gene]], main=paste0("Barcodeplot for ", gene),
              labels=c("Negative logFC", "Positive logFC"),
              quantile=c(-0.5,0.5))
}
gotable=function(list, string) {
  obj <- get('list')[[paste0(string, "_go")]]
  obj$GOID=rownames(obj)
  obj

}
upgoplot=function(list, string) {
  obj <- get('list')[[paste0(string, "_go")]]
  pos=obj[obj$P.Up<0.05,]
  pos=pos[order(pos$P.Up),]
  pos$logPvalue=-log(pos$P.Up)
  subset=pos[c(1:10),]
  subset$Term = str_wrap(subset$Term, width = 20)
  subset$Term <- factor(subset$Term, levels=rev(subset$Term))
  ggplotly(
    ggplot(subset, aes(x=logPvalue, y=Term, group=Up)) +
      geom_point(aes(size=Up),color="red") +
      theme_classic() +
      labs(title=paste0("GO enrichment of up regulated genes"), x="-log P value", y = "") +
      theme(plot.title = element_text(hjust = 0.5))
  , height = 600
  )
  
}
downgoplot=function(list, string) {
  obj <- get('list')[[paste0(string, "_go")]]
  neg=obj[obj$P.Down<0.05,]
  neg=neg[order(neg$P.Down),]
  neg$logPvalue=-log(neg$P.Down)
  subset=neg[c(1:10),]
  subset$Term = str_wrap(subset$Term, width = 20)
  subset$Term <- factor(subset$Term, levels=rev(subset$Term))
  ggplotly(
    ggplot(subset, aes(x=logPvalue, y=Term)) +
      geom_point(aes(size=Down),color="cornflowerblue") +
      theme_classic() +
      labs(title=paste0("GO enrichment of down regulated genes"), x="-log P value", y = "") +
      theme(plot.title = element_text(hjust = 0.5)) 
  , height = 600
  )
  
}
keggtable=function(list, string) {
  obj <- get('list')[[paste0(string, "_kegg")]]
  obj$KEGGID=rownames(obj)
  obj
  
}
upkeggplot=function(list, string) {
  obj <- get('list')[[paste0(string, "_kegg")]]
  pos=obj[obj$P.Up<0.05,]
  pos=pos[order(pos$P.Up),]
  pos$logPvalue=-log(pos$P.Up)
  subset=pos[c(1:10),]
  subset$Pathway = str_wrap(subset$Pathway, width = 20)
  subset$Pathway <- factor(subset$Pathway, levels=rev(subset$Pathway))
  ggplotly(
    ggplot(subset, aes(x=logPvalue, y=Pathway, group=Up)) +
      geom_point(aes(size=Up),color="red") +
      theme_classic() +
      labs(title=paste0("KEGG enrichment of up regulated genes"), x="-log P value", y = "") +
      theme(plot.title = element_text(hjust = 0.5)) 
    , height = 600 )
}
downkeggplot=function(list, string) {
  obj <- get('list')[[paste0(string, "_kegg")]]
  neg=obj[obj$P.Down<0.05,]
  neg=neg[order(neg$P.Down),]
  neg$logPvalue=-log(neg$P.Down)
  subset=neg[c(1:10),]
  subset$Pathway = str_wrap(subset$Pathway, width = 20)
  subset$Pathway <- factor(subset$Pathway, levels=rev(subset$Pathway))
  ggplotly(
    ggplot(subset, aes(x=logPvalue, y=Pathway)) +
      geom_point(aes(size=Down),color="cornflowerblue") +
      theme_classic() +
      labs(title=paste0("KEGG enrichment of down regulated genes"), x="-log P value", y = "") +
      theme(plot.title = element_text(hjust = 0.5)) 
    , height = 600
  )
}
guidecontrast=function(list, string, gene) {
  name=names(list)[which(str_detect(names(list), "_lrt"))]
  split=str_split(name, "_")
  contrasts=as.vector(sapply(split,"[[",1))
  compare=NULL
  for (i in contrasts) {
    obj <- get('list')[[paste0(i, "_lrt")]]
    compare=cbind(compare,obj$table$logFC)
  }
  colnames(compare)=contrasts
  rownames(compare)=rownames(obj)
  res=compare[which(str_detect(rownames(compare), gene)),]
  d=melt(t(res), id.vars=colnames(res))
  
  d %>% ggplot(aes(x=X1, y=value, fill=X1)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(colour = X2)) +
    theme_bw() +
    xlab("") +
    guides(fill = "none") +
    ggtitle(paste0(gene)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_discrete(name  ="Hairpin") +
    scale_y_continuous(name ="logFC") +
    scale_x_discrete(name ="Contrast")
}
meancontrast=function(list, string, gene) {
  
  name=names(list)[which(str_detect(names(list), "_lrt"))]
  split=str_split(name, "_")
  contrasts=as.vector(sapply(split,"[[",1))
  compare=NULL
  for (i in contrasts) {
    obj <- get('list')[[paste0(i, "_lrt")]]
    df=obj$table[which(str_detect(rownames(obj$table),  gene)),]
    compare=cbind(compare,mean(df$logFC))
  }
  colnames(compare)=contrasts
  d=melt(t(compare), id.vars=colnames(compare))
  colnames(d)=c("Contrast", "X2", "mean.logFC")
  ggplot(d, aes(x=Contrast, y=mean.logFC, fill=Contrast)) +
    geom_bar(stat="identity")  +
    theme_classic() +
    labs(title=paste0(gene), x="Contrast", y = "Mean logFC") +
    theme(plot.title = element_text(hjust = 0.5)) 
}
FCcontrast=function(list,string) {
  name=names(list)[which(str_detect(names(list), "_lrt"))]
  split=str_split(name, "_")
  contrasts=as.vector(sapply(split,"[[",1))
  vector=NULL
  dat=NULL
  for (gene in unique(list$x$genes$Gene)) {
    for (i in contrasts) {
      obj <- get('list')[[paste0(i, "_lrt")]]
      df=obj$table[which(str_detect(rownames(obj$table),  gene)),]
      vector=cbind(vector,mean(df$logFC))
    }
    vector=cbind(vector, gene)
    dat=rbind(dat,vector)
    vector=NULL
  }
  colnames(dat)=c(contrasts, "Gene")
  dat
}

ui <- fluidPage(
  titlePanel("RNA interference and CRISPR/Cas9 screen results"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("contrast", "Contrast", choices=""),
      numericInput("FDR", "FDR threshold", value=0.085, min=0, max=1),
      numericInput("logFC", "LogFC threshold", value=5, min=0, max=100)
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Upload",  fileInput('f1', 'Select workflow output', accept=c(".rds"), multiple=T)),
                  tabPanel("Quality control",
                           fluidRow(column(6,plotlyOutput('indexcounts'),style='padding:20px')
                                    ),
                           fluidRow(column(12,plotlyOutput('shcounts'),style='padding:20px')
                                    ),
                           fluidRow(column(6,plotlyOutput('bcv'),style='padding:20px')
                                    ),
                           fluidRow(column(6,plotOutput('mds'),style='padding:20px'),
                                    column(6,plotOutput('correctedmds'),style='padding:20px')
                                    ),
                           fluidRow(column(6,plotlyOutput('pca'),style='padding:20px'),
                                    column(6,plotlyOutput('correctedpca'),style='padding:20px')
                                    ),
                           fluidRow(column(6,plotlyOutput('sampledist'),style='padding:20px'),
                                    column(6,plotlyOutput('correctedsampledist'),style='padding:20px')
                                    )),
                  tabPanel("Differential expression",
                           fluidRow(column(6,plotlyOutput('plotsmear'),style='padding:20px'),
                                    column(6,plotlyOutput('volcano'),style='padding:20px')
                                    ),
                           fluidRow(column(6,plotlyOutput('de'),style='padding:20px'),
                                    column(6,plotlyOutput('correctedde'),style='padding:20px')
                                    ),
                           fluidRow(column(12,plotlyOutput('hist'),style='padding:20px')
                                    )),
                  tabPanel("Gene level",
                           fluidRow(column(4, DTOutput('camera'),style='padding:20px')
                                    ),
                           fluidRow(column(8, DTOutput('genelevel'),style='padding:20px')
                                    ),
                           fluidRow(column(12, plotlyOutput('guiderank'),style='padding:20px;height:500px')
                                    ),
                           fluidRow(column(4, selectInput("barcode_gene", "Select Gene", choices = ""), style='padding:20px;height:20px')
                                    ),
                           fluidRow(column(12, plotOutput('barcode'),style='padding:20px')
                                    )),
                  tabPanel("Enrichment",
                           fluidRow(column(12, DTOutput('gotable'), style="padding:20px")
                                    ),
                           fluidRow(column(6, plotlyOutput("upgoplot"),style='padding:20px;height:600px'),
                                    column(6, plotlyOutput("downgoplot"), style='padding:20px;height:600px')
                                    ),
                           fluidRow(column(12, DTOutput("keggtable"),style='padding:20px')
                                    ),
                           fluidRow(column(6, plotlyOutput("upkeggplot"),style='padding:20px;height:600px'),
                                    column(6, plotlyOutput("downkeggplot"), style="padding:20px;height:600px")
                                    )
                           ),
                  tabPanel("Comparing contrasts",
                           fluidRow(column(4, style='padding:30px;height:20px'),
                                    column(4, selectInput("contrast_gene", "Select Gene", choices=""), style='padding:30px;height:20px'),
                                    column(4, style='padding:30px;height:20px')
                                    ),
                           fluidRow(column(6, plotlyOutput('guidecontrast'), style='padding:30px'),
                                    column(6, plotlyOutput('meancontrast'), style='padding:30px')
                                    ),
                           fluidRow(column(12, DTOutput('FCcontrast'), style='padding:20px')
                                    )
                           ),
                  tabPanel("Rmarkdown export",
                           fluidRow(column(6, selectInput("export_gene", "Select Gene", choices=""), style='padding:30px;height:20px')
                           ),
                           fluidRow(column(6, downloadButton("report", "Generate report"), style='padding:30px;height:20px')
                           )
                  )
      )
    )
  )
)

server <- function(input, output, session) {
  list=reactive({
    file <- input$f1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "rds", "Please upload a rds file"))
    
    readRDS(file$datapath)
  })
  contrasts=reactive({
    name=names(list())[which(str_detect(names(list()), "_lrt"))]
    split=str_split(name, "_")
    as.vector(sapply(split,"[[",1))
  })
  observe({
    updateSelectInput(session,
                      "contrast",
                      choices = contrasts())
  })
  string=reactive({paste0(input$contrast)})
  FCthres=reactive({as.numeric(input$logFC)})
  FDRthres=reactive({as.numeric(input$FDR)})
  genes<-reactive({getgenes(list())})
  ### Quality control tab ###
  output$indexcounts <- renderPlotly({
    indexcounts(list())
  })
  output$shcounts <- renderPlotly({
      guidecounts(list())
  })
  output$bcv <- renderPlotly({
      bcv(list())
  })
  output$mds <- renderPlot({
      mds(list())
  })
  output$correctedmds <- renderPlot({
      cor_mds(list())
  })
  output$pca <- renderPlotly({
      pca(list())
  })
  output$correctedpca <- renderPlotly({
      cor_pca(list())
  })
  output$sampledist <- renderPlotly({
      sampledist(list())
  })
  output$correctedsampledist <- renderPlotly({
      cor_sampledist(list())
  })
  
  ### Differential expression tab ###
  
  output$plotsmear <- renderPlotly({
    plotsmear(list(), string(), FCthres(), FDRthres())
  })
  output$volcano <- renderPlotly({
    volcano(list(), string())
  })
  output$de <- renderPlotly({
    de(list(), string(), FCthres())
  })
  output$correctedde <- renderPlotly({
    cor_de(list(), string(), FCthres())
  })
  output$hist <- renderPlotly({
    hist(list(), string())
  })
  
  ### Gene level tab ###
  output$camera <- renderDT({
    datatable(camera(list(), string())
    )
  }, server=FALSE)
  output$genelevel <- renderDT({
    datatable(genelevel(list(), string())
    )
  })
  output$guiderank <- renderPlotly({
    guiderank(list(), string(), FCthres())
  })
  observe({
      updateSelectInput(session,
                       "barcode_gene",
                       choices = c(genes()))
  })
  barcodegene=reactive({
    paste0(input$barcode_gene)
    })
  output$barcode<- renderPlot({
    barcode(list(), string(), barcodegene())
  })
  
  ### Enrichment tab ###

  output$gotable <- renderDT({
    datatable(gotable(list(), string())
    )
  })
  output$upgoplot <- renderPlotly({
    upgoplot(list(), string())
  })
  output$downgoplot <- renderPlotly({
    downgoplot(list(), string())
    
  })
  output$keggtable <- renderDT({
    datatable(keggtable(list(), string())
    )
  })
  output$upkeggplot <- renderPlotly({
    upkeggplot(list(), string())
  })
  output$downkeggplot <- renderPlotly({
    downkeggplot(list(), string())
  })

  ### Comparing contrasts tab ###
  observe({
    updateSelectInput(session,
                      "contrast_gene",
                      choices = c(genes()))
  })
  contrastgene=reactive({
    paste0(input$contrast_gene)
  })
  output$guidecontrast <- renderPlotly({
    guidecontrast(list(), string(), contrastgene())
  })
  output$meancontrast <- renderPlotly({
    meancontrast(list(), string(), contrastgene())
  })
  output$FCcontrast <- renderDT({
    datatable(FCcontrast(list(), string())
    )
  })
  observe({
    updateSelectInput(session,
                      "export_gene",
                      choices = c(genes()))
  })
  exportgene=reactive({
    paste0(input$export_gene)
  })
  output$report <- downloadHandler(
    filename = "report.html",
    content = function(file) {
      #tempReport <-  file.path("report.Rmd")
      #file.copy("report.Rmd", tempReport, overwrite = TRUE)
      params <- list(g = exportgene())
      rmarkdown::render("~/report.Rmd", 
                        output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
}

shinyApp(ui,server)
