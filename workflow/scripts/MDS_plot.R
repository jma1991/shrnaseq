analysis=function(input, output) {
    library(edgeR)
    x=readRDS(input$rds)
    png(output$plot, width=2000, height=2000, res=400)
    plotMDS(x, labels=x$samples$group, col=rep(1:4, times=3), main="Another small screen: MDS Plot")
        legend("topright", legend=c(unique(x$samples$group)), col=1:4, pch=15)
    dev.off()
}

analysis(snakemake@input, snakemake@output)