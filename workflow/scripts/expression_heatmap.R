exp_hm <- function(x, mat, selY, hmtitle) {
    colnames(mat) <- paste(x$samples$group, x$samples$Replicate, sep = " - " )  
    mat <- subset(mat, rownames(mat) %in% selY)
    colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)
    pheatmap(mat, col = colors,border_color=NA, main=hmtitle)
}

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
    
    lrt=readRDS(input$rds[1])
    top2 = topTags(lrt, n=Inf)
    selY <- rownames(top2$table)[abs(top2$table$logFC)>params$FC]

    x=readRDS(input$rds[2])
    mat <- cpm(x$counts, log=TRUE, prior.count = 1)
    
    png(output$plot[1], width=3000, height=3500, res=400)
    exp_hm(x, mat, selY, "Differential expression \n across the groups (logCPM)")
    dev.off()

    #batch corrected
    mat=readRDS(input$rds[3])
    
    png(output$plot[2], width=3000, height=3500, res=400)
    exp_hm(x, mat, selY, "Batch corrected differential expression \n across the groups (logCPM)")
    dev.off()
}
  
analysis(snakemake@input, snakemake@output, snakemake@params, snakemake@log)