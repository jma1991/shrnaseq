analysis =function(input, output) {
    library(edgeR)
    x=readRDS(input$rds)
    sel=rowSums(cpm(x$counts)>0.5)>=3
    x=x[sel,]
    saveRDS(x, file = output$rds)
    png(output$plot,width=2000, height=2000, res=400)
        par(mfrow=c(2,1))
        barplot(colSums(x$counts), las=2, main="Counts per index", cex.names=0.5, cex.axis=0.8)
        barplot(rowSums(x$counts), las=2, main="Counts per hairpin", cex.names=0.5, cex.axis=0.8)
    dev.off()
}

analysis(snakemake@input, snakemake@output)