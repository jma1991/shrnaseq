analysis=function(input, output) {
    library(RColorBrewer)
    library(pheatmap)
    library(edgeR)
    x=readRDS(input$rds[1])
    lrt=readRDS(input$rds[2])
    top2 = topTags(lrt, n=Inf)
    y <- cpm(x$counts, log=TRUE, prior.count = 1)
  
    colnames(y) <- paste(x$samples$group, x$samples$Replicate, sep = " - " )
    selY <- rownames(top2$table)[abs(top2$table$logFC)>1.5]
    y = subset(y, rownames(y) %in% selY)

    colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)
    plt=pheatmap(t(y), col = colors, main="Differential expression across the groups (logCPM)",
         )

    png(output$plot, width=2000, height=2000, res=400)
    print(plt)
    dev.off()
}
  
analysis(snakemake@input, snakemake@output)