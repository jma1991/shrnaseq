analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    library(DESeq2)
    library(DEFormats)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)

    x=readRDS(input$rds[1])
    mat=cpm(x$counts, log=TRUE, prior.count = 1)
    sampleDists <- dist(t(mat))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(x$samples$group, x$samples$Replicate, sep = " - " )
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

    png(output$plot[1], width=2000, height=1500, res=400)
    pheatmap(sampleDistMatrix,
            clustering_distance_rows = sampleDists,
            clustering_distance_cols = sampleDists,
            col = colors,
            main = "Sample distances")
    dev.off()

    #batch corrected
    mat=readRDS(input$rds[2])
    sampleDists <- dist(t(mat))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(x$samples$group, x$samples$Replicate, sep = " - " )
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

    png(output$plot[2], width=2000, height=1500, res=400)
    pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
            main = "Batch corrected sample distances")
    dev.off()
}
  
analysis(snakemake@input, snakemake@output, snakemake@log)