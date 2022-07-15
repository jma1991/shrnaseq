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

    png(output$plot[1], width=2500, height=2000, res=400)
    pheatmap(t(y), col = colors, main="Differential expression across the groups (logCPM)")
    dev.off()

    #batch corrected
    corrected_lrt=readRDS(input$rds[3])
    corrected_top2 = topTags(corrected_lrt, n=Inf)
    corrected_y <- cpm(x$counts, log=TRUE, prior.count = 1)
  
    colnames(corrected_y) <- paste(x$samples$group, x$samples$Replicate, sep = " - " )
    corrected_selY <- rownames(corrected_top2$table)[abs(corrected_top2$table$logFC)>params$FC]
    corrected_y = subset(corrected_y, rownames(corrected_y) %in% corrected_selY)
    colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)

    png(output$plot[2], width=2500, height=2000, res=400)
    pheatmap(t(corrected_y), col = colors, main="Differential expression across the groups (logCPM)")
    dev.off()
}
  
analysis(snakemake@input, snakemake@output, snakemake@params, snakemake@log)