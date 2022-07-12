analysis=function(input, output) {
    library(edgeR)
    library(DESeq2)
    library(DEFormats)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    x=readRDS(input$rds)
    dds = as.DESeqDataSet(x, design = ~ x$samples$group)
    vsd <- vst(dds, blind = FALSE,  nsub=nrow(dds))
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$group, vsd$Replicate, sep = " - " )
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    plt=pheatmap(sampleDistMatrix,
                clustering_distance_rows = sampleDists,
                clustering_distance_cols = sampleDists,
                col = colors)

    png(output$plot, width=2800, height=2000, res=400)
    print(plt)
    dev.off()
}
  
analysis(snakemake@input, snakemake@output)