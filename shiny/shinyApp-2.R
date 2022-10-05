#### Install packages not yet installed ####

packages <- c("rlang", "shiny", "DT", "ggplot2", "RColorBrewer", "heatmaply", "plotly", "reshape", "limma", "tidyverse","bslib")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

#### Packages loading ####
invisible(lapply(packages, library, character.only = TRUE))

#### Functions ####
indexcounts=function(data) {
  
  df=melt(data.frame(t(colSums(data$x$counts))))
  ggplotly(
    ggplot(df, aes(x=variable, y=value)) +  
      geom_bar(stat = "identity", fill="cornflowerblue") +
      theme_minimal() +
      labs(title="Counts per sample", x="", y = "Counts") +
      theme(plot.title = element_text(hjust = 0.5))
  )
}
guidecounts=function(data) {
  df=melt(data.frame(t(rowSums(data$x$counts))))
  ggplotly(
    ggplot(df, aes(x=variable, y=value)) +  
      geom_bar(stat = "identity", fill="cornflowerblue") +
      theme_minimal() +
      labs(title="Counts per hairpin", x="", y = "Counts") +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.text.x = element_text(angle = 45))
  )
}
bcv=function(data) {
  df=data.frame(cbind(data$xglm$AveLogCPM, 
                      data$xglm$trended.dispersion, 
                      data$xglm$common.dispersion, 
                      data$xglm$tagwise.dispersion))
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
mds=function(data) {
  
  data$x$samples$group=factor(data$x$samples$group, levels = c(unique(data$x$samples$group)))
  
  mat <- cpm(data$x$counts, log = TRUE, prior.count = 1)
   
  var <- matrixStats::rowVars(mat)
   
  num <- min(500, length(var))
   
  ind <- order(var, decreasing = TRUE)[seq_len(num)]
   
  dst <- dist(t(mat[ind, ]))
  
  mds <- cmdscale(as.matrix(dst))
   
  dat <- data.frame(
     MD1 = mds[, 1], 
     MD2 = mds[, 2], 
     group = data$x$samples$group
  )

  ggplotly(
    ggplot(dat, aes(MD1, MD2, colour = group)) + 
      geom_point(size = 3) + 
      labs(x = "MDS 1", y = "MDS 2", colour = "Group") + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      labs(title="Multidimensional scaling") + 
      theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position="bottom") + 
      scale_color_brewer(palette = "Set3")
  )

}
cor_mds=function(data) {
  
  data$x$samples$group=factor(data$x$samples$group, levels = c(unique(data$x$samples$group)))
  
  mat <- data$corrected
  
  var <- matrixStats::rowVars(mat)
  
  num <- min(500, length(var))
  
  ind <- order(var, decreasing = TRUE)[seq_len(num)]
  
  dst <- dist(t(mat[ind, ]))
  
  mds <- cmdscale(as.matrix(dst))
  
  dat <- data.frame(
    MD1 = mds[, 1], 
    MD2 = mds[, 2], 
    group = data$x$samples$group
  )
  

  ggplotly(
    ggplot(dat, aes(MD1, MD2, colour = group)) + 
      geom_point(size = 3) + 
      labs(x = "MDS 1", y = "MDS 2", colour = "Group") +
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      labs(title="Batch corrected multidimensional scaling") + 
      theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position="bottom") + 
      scale_color_brewer(palette = "Set3")
  )

}
pca=function(data) {
  data$x$samples$group=factor(data$x$samples$group, levels = c(unique(data$x$samples$group)))
  
  mat=cpm(data$x$counts, log=TRUE, prior.count = 1)
  var <- matrixStats::rowVars(mat)
  num <- min(500, length(var))
  ind <- order(var, decreasing = TRUE)[seq_len(num)]
  pca <- prcomp(t(mat[ind, ]))
  pct <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)
  dat <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], group = data$x$samples$group)
  
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
      scale_color_brewer(palette="Set3")
  )
}
cor_pca=function(data) {
  var <- matrixStats::rowVars(data$corrected)
  num <- min(500, length(var))
  ind <- order(var, decreasing = TRUE)[seq_len(num)]
  pca <- prcomp(t(data$corrected[ind, ]))
  pct <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)
  dat <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], group = data$x$samples$group)
  
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
      scale_color_brewer(palette = "Set3")
  )
}
sampledist=function(data) {
  mat=cpm(data$x$counts, log=TRUE, prior.count = 1)
  sampleDists <- dist(t(mat))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  colnames(sampleDistMatrix) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  getPalette = colorRampPalette(brewer.pal(9, "Blues"))
  heatmaply(sampleDistMatrix, col=getPalette, main = "Sample distances") 
}
cor_sampledist=function(data) {
  sampleDists <- dist(t(data$corrected))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  colnames(sampleDistMatrix) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  getPalette = colorRampPalette(brewer.pal(9, "Blues"))
  heatmaply(sampleDistMatrix, col=getPalette, main = "Batch corrected sample distances")
}
hist=function(data, inputcontrast) {
  obj <- get('data')[[paste0(inputcontrast, "_lrt")]]
  df=data.frame(obj$table)
  ggplotly(
    ggplot(df, aes(x=PValue)) + 
      geom_histogram(bins=45,color="white",fill="cornflowerblue") +
      theme_classic() +
      labs(title="Histogram of hairpin P values", x="Hairpin p values", y = "Frequency") +
      theme(plot.title = element_text(hjust = 0.5)) 
  )
}
plotsmear=function(data, inputcontrast, FCthres, FDRthres) {
  obj <- get('data')[[paste0(inputcontrast, "_lrt")]]
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
volcano=function(data, inputcontrast) {
  obj <- get('data')[[paste0(inputcontrast, "_lrt")]]
  df=data.frame(obj$table)
  ggplotly(
    ggplot(df, aes(x=logFC, y=-10*log10(PValue))) + geom_point() +
      geom_point() +
      theme_bw() + 
      labs(title="Volcano plot", x="M", y = "-10*log(P-value)") +
      theme(plot.title = element_text(hjust = 0.5)) 
  )
}
de=function(data, inputcontrast, FCthres) {
  obj <- get('data')[[paste0(inputcontrast, "_lrt")]]
  top2 = topTags(obj, n=Inf)
  selY <- rownames(top2$table)[abs(top2$table$logFC)>FCthres]
  
  mat <- cpm(data$x$counts, log=TRUE, prior.count = 1)
  colnames(mat) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  mat = subset(mat, rownames(mat) %in% selY)
  getPalette = colorRampPalette(brewer.pal(9, "Blues"))
  heatmaply(mat, col=getPalette, main="Differential expression" )
}
cor_de=function(data, inputcontrast, FCthres) {
  obj <- get('data')[[paste0(inputcontrast, "_lrt")]]
  top2 = topTags(obj, n=Inf)
  selY <- rownames(top2$table)[abs(top2$table$logFC)>FCthres]
  
  colnames(data$corrected)= paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  corrected = subset(data$corrected, rownames(data$corrected) %in% selY)
  getPalette = colorRampPalette(brewer.pal(9, "Blues"))
  heatmaply(corrected, col=getPalette, main="Batch corrected differential expression")
}
detable=function(data, inputcontrast, FCthres, corrected) {
if (corrected=="Uncorrected") {
  mat <- cpm(data$x$counts, log=TRUE, prior.count = 1)
  colnames(mat) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  mat
  } else {
    colnames(data$corrected)= paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
    data$corrected
  }
}
camera=function(data, inputcontrast) {
  obj <- get('data')[[paste0(inputcontrast, "_camera")]]
  colnames(obj)=c("nGuides", "Direction", "Pvalue", "FDR", "Gene")
  obj[,c(5,1:4)]
}
genelevel=function(data, inputcontrast) {
  obj <- get('data')[[paste0(inputcontrast, "_genelevel")]]
  obj <- obj[order(obj$stouffers.pvalue),]
  colnames(obj)=c("Gene", "nGuides", 	"Mean logFC", "IQR logFC",
                  "Direction mean logFC",	"Direction smallest Pvalue",
                  "Stouffer's Pvalue",	"Stouffer's FDR")
  obj
}
guiderank=function(data, inputcontrast, FCthres, s) {
  obj <- get('data')[[paste0(inputcontrast, "_lrt")]]
  dat=data.frame(obj$table)
  dat$Rank = rank(dat$logFC)
  dat$guide=rownames(dat)
  ggplotly(
    ggplot(dat, aes(x=Rank, y=logFC, label=guide)) +
      geom_point(color="gray") +
      geom_point(data = dat %>% filter(dat$logFC>1), color = "red") +
      geom_point(data = dat %>% filter(dat$logFC<(-1)), color = "cornflowerblue") +
      {if (length(s)) geom_point(data = dat[s,], shape=1)} +
      geom_text(data = dat %>% filter(dat$logFC>FCthres), check_overlap = TRUE, nudge_x = max(dat$Rank)*0.05, size=2) +
      geom_text(data = dat %>% filter(dat$logFC<(-FCthres)), check_overlap = TRUE, nudge_x = max(dat$Rank)*0.05, size=2) +
      theme_classic() +
      labs(title="Guide rank plot") +
      theme(plot.title = element_text(hjust = 0.5)) 
    , height = 500
  )
}
getgenes=function(data) {
  genesymbollist = list()
  genesymbols=as.character(data$x$genes$Gene)
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
barcode=function(data, inputcontrast, gene) {
  obj <- get('data')[[paste0(inputcontrast, "_lrt")]]
  genesymbollist = list()
  genesymbols=as.character(data$x$genes$Gene)
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
gotable=function(data, inputcontrast) {
  obj <- get('data')[[paste0(inputcontrast, "_go")]]
  obj=obj[order(obj$P.Up),]
  obj$P.Up=signif(obj$P.Up,4)
  obj$P.Down=signif(obj$P.Down,4)
  colnames(obj)=c("Term", "Ont", "nGenes", "nGenes.Up", "nGenes.Down", "P.Up", "P.Down")
  obj

}
upgoplot=function(data, inputcontrast) {
  obj <- get('data')[[paste0(inputcontrast, "_go")]]
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
      labs(title="GO enrichment of up regulated genes", x="-log P value", y = "") +
      theme(plot.title = element_text(hjust = 0.5))
  , height = 600
  )
  
}
downgoplot=function(data, inputcontrast) {
  obj <- get('data')[[paste0(inputcontrast, "_go")]]
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
      labs(title="GO enrichment of down regulated genes", x="-log P value", y = "") +
      theme(plot.title = element_text(hjust = 0.5)) 
  , height = 600
  )
  
}
keggtable=function(data, inputcontrast) {
  obj <- get('data')[[paste0(inputcontrast, "_kegg")]]
  obj=obj[order(obj$P.Up),]
  obj$P.Up=signif(obj$P.Up,4)
  obj$P.Down=signif(obj$P.Down,4)
  colnames(obj)=c("Pathway", "nGenes", "nGenes.Up", "nGenes.Down", "P.Up", "P.Down")
  obj
}
upkeggplot=function(data, inputcontrast) {
  obj <- get('data')[[paste0(inputcontrast, "_kegg")]]
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
      labs(title="KEGG enrichment of up regulated genes", x="-log P value", y = "") +
      theme(plot.title = element_text(hjust = 0.5)) 
    , height = 600 )
}
downkeggplot=function(data, inputcontrast) {
  obj <- get('data')[[paste0(inputcontrast, "_kegg")]]
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
      labs(title="KEGG enrichment of down regulated genes", x="-log P value", y = "") +
      theme(plot.title = element_text(hjust = 0.5)) 
    , height = 600
  )
}
guidecontrast=function(data, inputcontrast, gene) {
  name=names(data)[which(str_detect(names(data), "_lrt"))]
  split=str_split(name, "_")
  contrasts=as.vector(sapply(split,"[[",1))
  compare=NULL
  for (i in contrasts) {
    obj <- get('data')[[paste0(i, "_lrt")]]
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
    ggtitle(paste(gene)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_discrete(name  ="Hairpin") +
    scale_y_continuous(name ="logFC") +
    scale_x_discrete(name ="Contrast")
}
meancontrast=function(data, inputcontrast, gene) {
  
  name=names(data)[which(str_detect(names(data), "_lrt"))]
  split=str_split(name, "_")
  contrasts=as.vector(sapply(split,"[[",1))
  compare=NULL
  for (i in contrasts) {
    obj <- get('data')[[paste0(i, "_lrt")]]
    df=obj$table[which(str_detect(rownames(obj$table),  gene)),]
    compare=cbind(compare,mean(df$logFC))
  }
  colnames(compare)=contrasts
  d=melt(t(compare), id.vars=colnames(compare))
  colnames(d)=c("Contrast", "X2", "mean.logFC")
  ggplot(d, aes(x=Contrast, y=mean.logFC, fill=Contrast)) +
    geom_bar(stat="identity")  +
    theme_classic() +
    labs(title=paste(gene), x="Contrast", y = "Mean logFC") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_brewer(palette="Set3")
}
FCcontrast=function(data,inputcontrast) {
  name=names(data)[which(str_detect(names(data), "_lrt"))]
  split=str_split(name, "_")
  contrasts=as.vector(sapply(split,"[[",1))
  vector=NULL
  dat=NULL
  for (gene in unique(data$x$genes$Gene)) {
    for (i in contrasts) {
      obj <- get('data')[[paste0(i, "_lrt")]]
      df=obj$table[which(str_detect(rownames(obj$table),  gene)),]
      vector=cbind(vector,mean(df$logFC))
    }
    vector=cbind(gene, vector)
    dat=rbind(dat,vector)
    vector=NULL
  }
  colnames(dat)=c(contrasts, "Gene")
  dat
}

ui <- fluidPage(
  tags$head(
    # Note the wrapping of the string in HTML()
    tags$style(HTML("
      @import url('https://fonts.googleapis.com/css2?family=Roboto');
      h1 {
        font-family: 'Roboto';
      }"))
  ),
  titlePanel(h1("RNA interference and CRISPR/Cas9 screen results")),
  theme =  bs_theme(bootswatch = "minty",
                    base_font = font_google("Roboto")),
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
                           fluidRow(column(6,plotlyOutput('mds'),style='padding:20px'),
                                    column(6,plotlyOutput('correctedmds'),style='padding:20px')
                                    ),
                           fluidRow(column(6,plotlyOutput('pca'),style='padding:20px'),
                                    column(6,plotlyOutput('correctedpca'),style='padding:20px')
                                    ),
                           fluidRow(column(6,plotlyOutput('sampledist'),style='padding:20px'),
                                    column(6,plotlyOutput('correctedsampledist'),style='padding:20px')
                                    )),
                  tabPanel("Differential expression",
                           fluidRow(column(12,plotlyOutput('hist'),style='padding:20px')
                                    ),
                           fluidRow(column(6,plotlyOutput('plotsmear'),style='padding:20px'),
                                    column(6,plotlyOutput('volcano'),style='padding:20px')
                                    ),
                           fluidRow(column(6,plotlyOutput('de'),style='padding:20px'),
                                    column(6,plotlyOutput('correctedde'),style='padding:20px')
                                    ),
                           fluidRow(column(4, style='padding:30px;height:20px'),
                                    column(4, selectInput("corrected", "Select uncorrected or corrected", choices=c("Uncorrected", "Corrected")), style='padding:30px;height:20px')
                                    ),
                           fluidRow(column(12, DTOutput('detable'),style='padding:30px')
                           )),
                  tabPanel("Gene level",
                           fluidRow(column(12, DTOutput('camera'),style='padding:20px')
                                    ),
                           fluidRow(column(12, DTOutput('genelevel'),style='padding:20px')
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
    file <- input$f1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "rds", "Please upload a rds file"))
    
    readRDS(file$datapath)
  })
  contrasts=reactive({
    name=names(data())[which(str_detect(names(data()), "_lrt"))]
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
    hist(data(), inputcontrast())
  })
  output$plotsmear <- renderPlotly({
    plotsmear(data(), inputcontrast(), FCthres(), FDRthres())
  })
  output$volcano <- renderPlotly({
    volcano(data(), inputcontrast())
  })
  output$de <- renderPlotly({
    de(data(), inputcontrast(), FCthres())
  })
  output$correctedde <- renderPlotly({
    cor_de(data(), inputcontrast(), FCthres())
  })
  output$detable <- renderDT({
    datatable(detable(data(), inputcontrast(), FCthres(), input$corrected), 
              
              caption = htmltools::tags$caption(style = 'caption-side: top; font-weight: bold; color:black; text-align: center; font-size:150% ;','Differential expression (logCPM)')
    )
  })
  ### Gene level tab ###
  output$camera <- renderDT({
    datatable(camera(data(), inputcontrast()), 
              
              caption = htmltools::tags$caption(style = 'caption-side: top; font-weight: bold; color:black; text-align: center; font-size:150% ;','Camera')
    )
  })
  output$genelevel <- renderDT({
    datatable(genelevel(data(), inputcontrast()),
              editable = "row", 
              caption = htmltools::tags$caption( style = 'caption-side: top; font-weight: bold; color:black; text-align: center; font-size:150% ;','Gene level')
              
    )
  })
  output$guiderank <- renderPlotly({
    s = input$genelevel_rows_selected
    guiderank(data(), inputcontrast(), FCthres(), s)
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
    barcode(data(), inputcontrast(), barcodegene())
  })
  
  ### Enrichment tab ###

  output$gotable <- renderDT({
    datatable(gotable(data(), inputcontrast()),
              
              caption = htmltools::tags$caption(style = 'caption-side: top; font-weight: bold; color:black; text-align: center; font-size:150% ;','Gene ontology')
              
    )
  })
  output$upgoplot <- renderPlotly({
    upgoplot(data(), inputcontrast())
  })
  output$downgoplot <- renderPlotly({
    downgoplot(data(), inputcontrast())
    
  })
  output$keggtable <- renderDT({
    datatable(keggtable(data(), inputcontrast()),
              
              caption = htmltools::tags$caption( style = 'caption-side: top; font-weight: bold; color:black; text-align: center; font-size:150% ;','KEGG pathway enrichment')
              
    )
  })
  output$upkeggplot <- renderPlotly({
    upkeggplot(data(), inputcontrast())
  })
  output$downkeggplot <- renderPlotly({
    downkeggplot(data(), inputcontrast())
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
    guidecontrast(data(), inputcontrast(), contrastgene())
  })
  output$meancontrast <- renderPlotly({
    meancontrast(data(), inputcontrast(), contrastgene())
  })
  output$FCcontrast <- renderDT({
    datatable(FCcontrast(data(), inputcontrast()),
              
              caption = htmltools::tags$caption( style = 'caption-side: top; font-weight: bold; color:black; text-align: center; font-size:150% ;','Mean logFC across contrasts')
              
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
    filename = "report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- base::list(gene = exportgene(), 
                           contrast=inputcontrast(),
                           data=data())
      
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
