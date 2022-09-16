analysis =function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    x=readRDS(input$rds)
    x=calcNormFactors(x)
    saveRDS(x,file=output$rds)
}

analysis(snakemake@input, snakemake@output, snakemake@log)