pcaplt <- function(x, mat, pcatitle) {
    var <- matrixStats::rowVars(mat)
    num <- min(500, length(var))
    ind <- order(var, decreasing = TRUE)[seq_len(num)]
    pca <- prcomp(t(mat[ind, ]))
    pct <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)
    dat <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], group = x$samples$group)
    colors <- brewer.pal(length(unique(x$samples$group)),"Set3")

    ggplot(dat, aes_string(x = "PC1", y = "PC2", color = "group")) + 
    geom_point(size = 3) + 
    xlab(paste0("PC1: ", round(pct[1] * 100), "% variance")) + 
    ylab(paste0("PC2: ", round(pct[2] * 100), "% variance")) + 
    coord_fixed() +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(title=pcatitle) + 
    theme(plot.title = element_text(hjust = 0.5,  face="bold")) + theme(legend.position="bottom") +
    scale_color_manual(values=colors)

}

analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    library(matrixStats)
    library(ggplot2)
    library(RColorBrewer)

    x=readRDS(input$rds[1])
    x$samples$group=factor(x$samples$group, levels = c(unique(x$samples$group)))
    
    mat=cpm(x$counts, log=TRUE, prior.count = 1)
    png(output$plot[1], width=2500, height=2000, res=400)
    par(mar=c(5,5,6,5))
    print(pcaplt(x, mat, "Principal component analysis"))
    dev.off()

    #batch corrected
    mat=readRDS(input$rds[2])

    png(output$plot[2], width=2500, height=2000, res=400)
    par(mar=c(5,5,6,5))
    print(pcaplt(x, mat, "Batch corrected principal component analysis"))
    dev.off()

}
  
analysis(snakemake@input, snakemake@output, snakemake@log)