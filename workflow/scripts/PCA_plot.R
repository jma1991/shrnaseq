analysis=function(input, output) {
    library(edgeR)
    library(DESeq2)
    library(DEFormats)
    library(ggplot2)
    x=readRDS(input$rds)
    dds = as.DESeqDataSet(x, design = ~ x$samples$group)
    dds=rlog(dds)
    dds$group <- factor(dds$group, levels = c(unique(dds$group)))
    plt=plotPCA(dds, intgroup="group") + theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(title="Principal component analysis") 
    png(output$plot, width=2000, height=2000, res=400)
    print(plt)
    dev.off()
}
  
analysis(snakemake@input, snakemake@output)