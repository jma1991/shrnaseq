analysis=function(input, output, params, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(RColorBrewer)
    library(pheatmap)
    library(edgeR)
    x=readRDS(input$rds[1])
    lrt=readRDS(input$rds[2])
    top2 = topTags(lrt, n=Inf)
    y <- cpm(x$counts, log=TRUE, prior.count = 1)
  
    colnames(y) <- paste(x$samples$group, x$samples$Replicate, sep = " - " )
    selY <- rownames(top2$table)[abs(top2$table$logFC)>params$FC]
    y = subset(y, rownames(y) %in% selY)
    colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)

    png(output$plot, width=2500, height=2000, res=400)
    pheatmap(t(y), col = colors, main="Differential expression across the groups (logCPM)")
    dev.off()
}
  
analysis(snakemake@input, snakemake@output, snakemake@params, snakemake@log)