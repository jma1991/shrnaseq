analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script 
    library(limma)
    library(edgeR)
    x=readRDS(input$rds)

    mat=cpm(x$counts, log=TRUE, prior.count = 1)
    mat=removeBatchEffect(mat, batch=factor(x$samples$batch))

    saveRDS(mat, file=output$rds)
}
analysis(snakemake@input, snakemake@output, snakemake@log)