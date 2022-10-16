#### Functions ####
indexcounts <- function(data) {
  
  df <- melt(data.frame(t(colSums(data$x$counts))))
  ggplotly(
    ggplot(df, aes(x=variable, y=value)) +  
      geom_bar(stat = "identity", fill="#b8dbcc") +
      theme_minimal() +
      labs(title="Counts per sample", x="", y = "Counts") +
      theme(plot.title = element_text(hjust = 0.5))
  )
}
guidecounts <- function(data) {
  df <- melt(data.frame(t(rowSums(data$x$counts))))
  ggplotly(
    ggplot(df, aes(x=variable, y=value)) +  
      geom_bar(stat = "identity", fill="#b8dbcc") +
      theme_minimal() +
      labs(title="Counts per guide RNA", x="", y = "Counts") +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.text.x = element_text(angle = 45))
  )
}
bcv <- function(data) {
  df <- data.frame(cbind(data$xglm$AveLogCPM, 
                      data$xglm$trended.dispersion, 
                      data$xglm$common.dispersion, 
                      data$xglm$tagwise.dispersion))
  colnames(df) <- c("AveLogCPM", "trended.dispersion",
                 "commom.dispersion", "tagwise.dispersion")
  ggplotly(
    ggplot(df, aes(x=AveLogCPM)) +
      geom_line(aes(y=trended.dispersion), colour = "#b8dbcc") + 
      geom_line(aes(y=commom.dispersion), colour = "red") + 
      geom_point(aes(y=tagwise.dispersion), colour = "#b8dbcc") + 
      theme_classic() + 
      labs(title="Biological coefficient of variation plot", x="Average log CPM", y = "Biological coefficient of variation") +
      theme(plot.title = element_text(hjust = 0.5))
  )
}
mds <- function(data) {
  
  data$x$samples$group <- factor(data$x$samples$group, levels = c(unique(data$x$samples$group)))
  
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
cor_mds <- function(data) {
  
  data$x$samples$group <- factor(data$x$samples$group, levels = c(unique(data$x$samples$group)))
  
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
pca <- function(data) {
  data$x$samples$group <- factor(data$x$samples$group, levels = c(unique(data$x$samples$group)))
  
  mat <- cpm(data$x$counts, log=TRUE, prior.count = 1)
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
cor_pca <- function(data) {
  data$x$samples$group <- factor(data$x$samples$group, levels = c(unique(data$x$samples$group)))
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
sampledist <- function(data) {
  mat <- cpm(data$x$counts, log=TRUE, prior.count = 1)
  sampleDists <- dist(t(mat))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  colnames(sampleDistMatrix) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  getPalette <- colorRampPalette(brewer.pal(9, "BuGn"))
  heatmaply(sampleDistMatrix, col=getPalette, main = "Sample distances") 
}
cor_sampledist <- function(data) {
  sampleDists <- dist(t(data$corrected))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  colnames(sampleDistMatrix) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  getPalette <- colorRampPalette(brewer.pal(9, "BuGn"))
  heatmaply(sampleDistMatrix, col=getPalette, main = "Batch corrected sample distances")
}
hist <- function(contrast_data, inputcontrast) {
  obj <- get('contrast_data')[[paste0(inputcontrast, "_lrt")]]
  df <- data.frame(obj$table)
  ggplotly(
    ggplot(df, aes(x=PValue)) + 
      geom_histogram(bins=45,color="white",fill="#b8dbcc") +
      theme_classic() +
      labs(title="Histogram of hairpin P values", x="Hairpin p values", y = "Frequency") +
      theme(plot.title = element_text(hjust = 0.5)) 
  )
}
plotsmear <- function(contrast_data, inputcontrast, FCthres, FDRthres) {
  obj <- get('contrast_data')[[paste0(inputcontrast, "_lrt")]]
  top2 <- topTags(obj, n=Inf)
  top2ids <- top2$table[top2$table$FDR<FDRthres,1]
  df <- data.frame(obj$table)
  df$Guide <- rownames(df)
  ggplotly(
    ggplot(df, aes(x=logCPM, y=logFC, text=Guide)) +geom_point() +
      geom_point(data = df %>% filter(row.names(df) %in% top2ids), color = "red") +
      geom_hline(yintercept=(-(FCthres)), linetype="dashed", color="#b8dbcc") +
      geom_hline(yintercept=0,  color="cornflowerblue") +
      geom_hline(yintercept=(FCthres), linetype="dashed", color="#b8dbcc") +
      theme_bw() + labs(title="Mean-difference plot") + theme(plot.title = element_text(hjust = 0.5))
  )
}
volcano <- function(contrast_data, inputcontrast) {
  obj <- get('contrast_data')[[paste0(inputcontrast, "_lrt")]]
  df <- data.frame(obj$table)
  df$Guide <- rownames(df)
  ggplotly(
    ggplot(df, aes(x=logFC, y=-10*log10(PValue), text=Guide)) + geom_point() +
      geom_point() +
      theme_bw() + 
      labs(title="Volcano plot", x="M", y = "-10*log(P-value)") +
      theme(plot.title = element_text(hjust = 0.5)) 
  )
}
de <- function(data, contrast_data, inputcontrast, topguides) {
  obj <- get('contrast_data')[[paste0(inputcontrast, "_lrt")]]
  top2 <- data.frame(topTags(obj, n=Inf))
  top2 <- top2[order(top2$logFC),]
  selY <- rownames(top2)[c(1:topguides)]
  
  mat <- cpm(data$x$counts, log=TRUE, prior.count = 1)
  colnames(mat) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  mat <- subset(mat, rownames(mat) %in% selY)
  getPalette <- colorRampPalette(brewer.pal(9, "BuGn"))
  heatmaply(mat, col=getPalette, main="Differential expression" )
}
cor_de <- function(data, contrast_data, inputcontrast, topguides) {
  obj <- get('contrast_data')[[paste0(inputcontrast, "_lrt")]]
  top2 <- data.frame(topTags(obj, n=Inf))
  top2=top2[order(top2$logFC),]
  selY <- rownames(top2)[c(1:topguides)]
  
  colnames(data$corrected) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  corrected <- subset(data$corrected, rownames(data$corrected) %in% selY)
  getPalette <- colorRampPalette(brewer.pal(9, "BuGn"))
  heatmaply(corrected, col=getPalette, main="Batch corrected differential expression")
}
detable <- function(data, inputcontrast, corrected) {
  if (corrected=="Uncorrected") {
    mat <- cpm(data$x$counts, log=TRUE, prior.count = 1)
    colnames(mat) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
    mat
  } else {
    colnames(data$corrected) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
    data$corrected
  }
}
camera <- function(contrast_data, inputcontrast) {
  obj <- get('contrast_data')[[paste0(inputcontrast, "_camera")]]
  colnames(obj) <- c("nGuides", "Direction", "Pvalue", "FDR", "Gene")
  obj[,c(5,1:4)]
}
camerarank <- function(contrast_data, inputcontrast, s) {
  obj <- get('contrast_data')[[paste0(inputcontrast, "_camera")]]
  colnames(obj) <- c("nGuides", "Direction", "Pvalue", "FDR", "Gene")
  
  # obj$Rank <- rank(-log(obj$Pvalue))
  # ggplotly(
  #   ggplot(obj, aes(x=Rank, y=(-log(Pvalue)), label=Gene)) +
  #     geom_point(aes(size=nGuides), color="#b8dbcc") +
  #     geom_point(aes(size=nGuides), color = "red",
  #                data = obj %>% filter(obj$Pvalue<0.05 & obj$Direction=="Up")) +
  #     geom_point(aes(size=nGuides), color = "cornflowerblue",
  #                data = obj %>% filter(obj$Pvalue<0.05 & obj$Direction=="Down")) +
  #     theme_classic() +
  #     labs(title="Camera rank plot", y="-log(P Value)") +
  #     theme(plot.title = element_text(hjust = 0.5)) +
  #     {if (length(s)) geom_point(aes(size=nGuides), data = obj[s,], shape=1) } +
  #     {if (length(s)) geom_text(data = obj[s,], check_overlap = TRUE, hjust = 0, nudge_x = (max(obj$Rank))*0.05, size=3)} 
  #   , height = 600
  # ) 

  res <- as.data.frame(obj)
  
  fdr <- 0.6

  res$Status <- factor("NS", levels = c("Up", "NS", "Down"))
  
  res$Status[res$Direction == "Up" & res$FDR < fdr] <- "Up"
  
  res$Status[res$Direction == "Down" & res$FDR < fdr] <- "Down"
    
  res$Pvalue <- -log10(res$Pvalue)
  
  res$Rank <- rank(res$Pvalue)
  
  col <- c(
    "Up"   = "#FF0000",
    "NS"   = "#B8DBCC",
    "Down" = "#6495ED"
  )

  res.s <- res[s, , drop = FALSE]
  
  plt <- ggplot(res, aes(x = Rank, y = Pvalue, colour = Status, size = nGuides, label = Gene)) + 
    geom_point() + 
    geom_point(data = res.s, shape = 1) + 
    geom_text(data = res.s, hjust = 0, nudge_x = 10) +  
    scale_colour_manual(values = col, breaks = names(col)) + 
    labs(
      title = "Camera rank plot",
      x = "Rank",
      y = "-log10(Pvalue)",
      colour = "Status"
    ) + 
    guides(size = "none") +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5))

  plt <- ggplotly(plt, height = 600)

  lab <- c(
    "Up"   = sprintf("Up (%s)", comma(sum(res$Status == "Up"))),
    "NS"   = sprintf("NS (%s)", comma(sum(res$Status == "NS"))),
    "Down" = sprintf("Down (%s)", comma(sum(res$Status == "Down")))
  )
  
  plt$x$layout$legend$title$text <- "Status"
  
  plt$x$data[[1]]$name <- lab["Up"]

  plt$x$data[[2]]$name <- lab["NS"]
  
  plt$x$data[[3]]$name <- lab["Down"]
  
  plt <- layout(p = plt, legend = list(itemsizing = "constant"))
  
  plt

}

genelevel <- function(contrast_data, inputcontrast) {
  obj <- get('contrast_data')[[paste0(inputcontrast, "_genelevel")]]
  colnames(obj) <- c("Gene", "nGuides", 	"Mean logFC", "IQR logFC",
                  "Direction mean logFC",	"Direction smallest Pvalue",
                  "Stouffer's Pvalue",	"Stouffer's FDR")
  row.names(obj) <- NULL
  obj
}
generank <- function(contrast_data, inputcontrast, s) {
  obj <- get('contrast_data')[[paste0(inputcontrast, "_genelevel")]]
  colnames(obj) <- c("Gene", "nGuides", 	"Mean logFC", "IQR logFC",
                  "Direction mean logFC",	"Direction smallest Pvalue",
                  "Stouffer's Pvalue",	"Stouffer's FDR")
  obj$Rank <- rank(obj$`Mean logFC`)
  ggplotly(
    ggplot(obj, aes(x=Rank, y=`Mean logFC`, label=Gene)) +
      geom_point(aes(size=nGuides), color="#b8dbcc") +
      geom_point(aes(size=nGuides), color = "red",
        data = obj %>% filter(obj$`Mean logFC`>1)) +
      geom_point(aes(size=nGuides), color = "cornflowerblue",
        data = obj %>% filter(obj$`Mean logFC`<(-1))) +
      theme_classic() +
      labs(title="Gene rank plot") +
      theme(plot.title = element_text(hjust = 0.5)) +
      {if (length(s)) geom_point(aes(size=nGuides), data = obj[s,], shape=1) } +
      {if (length(s)) geom_text(data = obj[s,], check_overlap = TRUE, hjust = 0, nudge_x = (max(obj$Rank))*0.05, size=3)} 
    , height = 600
  )
}
getgenes <- function(data) {
  genesymbollist <- list()
  genesymbols <- as.character(data$x$genes$Gene)
  unq <- unique(genesymbols)
  unq <- unq[!is.na(unq)]
  genes <- NULL
  for (i in unq) {
    sel <- genesymbols == i & !is.na(genesymbols)
    if (sum(sel) > 1) 
      genes=rbind(genes, i)
  }
  as.vector(genes)
}
barcode <- function(data, contrast_data, inputcontrast, gene) {
  obj <- get('contrast_data')[[paste0(inputcontrast, "_lrt")]]
  genesymbollist <- list()
  genesymbols <- as.character(data$x$genes$Gene)
  unq <- unique(genesymbols)
  unq <- unq[!is.na(unq)]
  
  for (i in unq) {
    sel <- genesymbols == i & !is.na(genesymbols)
    if (sum(sel) > 1) 
      genesymbollist[[i]]  <- which(sel)
  }
  barcodeplot(obj$table$logFC,index=genesymbollist[[gene]], main=paste0("Barcodeplot for ", gene),
              labels=c("Negative logFC", "Positive logFC"),
              quantile=c(-0.5,0.5))
}
gotable <- function(contrast_data, inputcontrast) {
  obj <- get('contrast_data')[[paste0(inputcontrast, "_go")]]
  obj <- obj[order(obj$P.Up),]
  obj$P.Up <- signif(obj$P.Up,4)
  obj$P.Down <- signif(obj$P.Down,4)
  colnames(obj) <- c("Term", "Ont", "nGenes", "nGenes.Up", "nGenes.Down", "P.Up", "P.Down")
  obj
  
}
upgoplot <- function(contrast_data, inputcontrast, topgos) {
  obj <- get('contrast_data')[[paste0(inputcontrast, "_go")]]
  pos <- obj[order(obj$P.Up),]
  pos$logPvalue <- -log(pos$P.Up)
  subset <- pos[c(1:topgos),]
  subset$Term <- str_wrap(subset$Term, width = 20)
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
downgoplot <- function(contrast_data, inputcontrast, topgos) {
  obj <- get('contrast_data')[[paste0(inputcontrast, "_go")]]
  neg <- obj[order(obj$P.Down),]
  neg$logPvalue <- -log(neg$P.Down)
  subset <- neg[c(1:topgos),]
  subset$Term <- str_wrap(subset$Term, width = 20)
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
keggtable <- function(contrast_data, inputcontrast) {
  obj <- get('contrast_data')[[paste0(inputcontrast, "_kegg")]]
  obj <- obj[order(obj$P.Up),]
  obj$P.Up <- signif(obj$P.Up,4)
  obj$P.Down <- signif(obj$P.Down,4)
  colnames(obj) <- c("Pathway", "nGenes", "nGenes.Up", "nGenes.Down", "P.Up", "P.Down")
  obj
}
upkeggplot <- function(contrast_data, inputcontrast, topkeggs) {
  obj <- get('contrast_data')[[paste0(inputcontrast, "_kegg")]]
  pos <- obj[order(obj$P.Up),]
  pos$logPvalue <- -log(pos$P.Up)
  subset <- pos[c(1:topkeggs),]
  subset$Pathway <- str_wrap(subset$Pathway, width = 20)
  subset$Pathway <- factor(subset$Pathway, levels=rev(subset$Pathway))
  ggplotly(
    ggplot(subset, aes(x=logPvalue, y=Pathway, group=Up)) +
      geom_point(aes(size=Up),color="red") +
      theme_classic() +
      labs(title="KEGG enrichment of up regulated genes", x="-log P value", y = "") +
      theme(plot.title = element_text(hjust = 0.5)) 
    , height = 600
    )
}
downkeggplot <- function(contrast_data, inputcontrast, topkeggs) {
  obj <- get('contrast_data')[[paste0(inputcontrast, "_kegg")]]
  neg <- obj[order(obj$P.Down),]
  neg$logPvalue <- -log(neg$P.Down)
  subset <- neg[c(1:topkeggs),]
  subset$Pathway <- str_wrap(subset$Pathway, width = 20)
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
guidecontrast <- function(contrast_data, gene) {
  name <- names(contrast_data)[which(str_detect(names(contrast_data), "_lrt"))]
  split <- str_split(name, "_")
  contrasts <- as.vector(sapply(split,"[[",1))
  compare <- NULL
  for (i in contrasts) {
    obj <- get('contrast_data')[[paste0(i, "_lrt")]]
    compare <- data.frame(cbind(compare,obj$table$logFC))
  }
  colnames(compare) <- contrasts
  compare$Guides <- rownames(obj$table)
  res <- data.frame(compare[which(str_detect(compare$Guides, gene)),])
  d <- melt(res)
  colnames(d) <- c("Guides","Contrast","logFC")
  
  ggplotly(
    d %>% ggplot(aes(x=Contrast, y=logFC)) +
      theme_bw() +
      ggtitle(paste(gene)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_boxplot(outlier.shape = NA, fill="#b8dbcc", show.legend=FALSE) +
      geom_jitter(aes(colour = Guides), show.legend = TRUE) +
      geom_hline(yintercept=0, linetype="dashed", color = "gray") +
      theme(legend.position="none")
  )
}
FCcontrast <- function(data, contrast_data, inputcontrast) {
  name <- names(contrast_data)[which(str_detect(names(contrast_data), "_lrt"))]
  split <- str_split(name, "_")
  contrasts <- as.vector(sapply(split,"[[",1))
  vector <- NULL
  dat <- NULL
  for (gene in unique(data$x$genes$Gene)) {
    for (i in contrasts) {
      obj <- get('contrast_data')[[paste0(i, "_lrt")]]
      df <- obj$table[which(str_detect(rownames(obj$table),  gene)),]
      vector <- cbind(vector,mean(df$logFC))
    }
    nguides <- nrow(data$x$genes[data$x$genes$Gene==gene,])
    vector <- cbind(gene, nguides, vector)
    dat <- rbind(dat,vector)
    vector <- NULL
  }
  colnames(dat) <- c("Gene", "nGuides", contrasts)
  dat
}
