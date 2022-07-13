analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    xglm=readRDS(input$rds)
    png(output$plot, width=2000, height=2000, res=400)
        plotBCV(xglm, main="BCV Plot")
    dev.off()
}

analysis(snakemake@input, snakemake@output, snakemake@log)