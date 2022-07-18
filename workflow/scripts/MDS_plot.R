analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script 
    library(edgeR)
    x=readRDS(input$rds[1])
    png(output$plot[1], width=2000, height=2000, res=400)
    plotMDS(x, labels=x$samples$group, col=rep(1:4, times=3), main="MDS Plot")
        legend("topright", legend=c(unique(x$samples$group)), col=1:4, pch=15)
    dev.off()

    mat=readRDS(input$rds[2])
    png(output$plot[2], width=2000, height=2000, res=400)
    plotMDS(mat, labels=x$samples$group, col=rep(1:4, times=3), main="MDS Plot")
        legend("topright", legend=c(unique(x$samples$group)), col=1:4, pch=15)
    dev.off()
}

analysis(snakemake@input, snakemake@output, snakemake@log)