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
    library(ggplot2)#

    x=readRDS(input$rds[1])
    dds = as.DESeqDataSet(x, design = ~ x$samples$group)
    dds=rlog(dds)
    dds$group <- factor(dds$group, levels = c(unique(dds$group)))
    plt=plotPCA(dds, intgroup="group") + theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(title="Principal component analysis") 
    png(output$plot[1], width=2000, height=2000, res=400)
    print(plt)
    dev.off()

    #batch corrected
    corrected=readRDS(input$rds[2])
    dds = as.DESeqDataSet(corrected, design = ~ corrected$samples$group)
    dds=rlog(dds)
    dds$group <- factor(dds$group, levels = c(unique(dds$group)))
    plt=plotPCA(dds, intgroup="group") + theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(title="Principal component analysis") 
    png(output$plot[2], width=2000, height=2000, res=400)
    print(plt)
    dev.off()

}
  
analysis(snakemake@input, snakemake@output, snakemake@log)