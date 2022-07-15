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
    dds = as.DESeqDataSet(x, design = ~ x$samples$group)
    vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$group, vsd$Replicate, sep = " - " )
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

    png(output$plot[1], width=2800, height=2000, res=400)
    pheatmap(sampleDistMatrix,
                clustering_distance_rows = sampleDists,
                clustering_distance_cols = sampleDists,
                col = colors)
    dev.off()

    #batch corrected
    corrected=readRDS(input$rds[2])
    dds = as.DESeqDataSet(corrected, design = ~ corrected$samples$group)
    vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$group, vsd$Replicate, sep = " - " )
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

    png(output$plot[2], width=2800, height=2000, res=400)
    pheatmap(sampleDistMatrix,
                clustering_distance_rows = sampleDists,
                clustering_distance_cols = sampleDists,
                col = colors)
    dev.off()
}
  
analysis(snakemake@input, snakemake@output, snakemake@log)