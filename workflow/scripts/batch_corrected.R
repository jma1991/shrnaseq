analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script 
    library(limma)
    library(edgeR)
    xglm=readRDS(input$rds)
    xglm$counts=round(exp(removeBatchEffect(log(xglm$counts), 
                                       batch=factor(xglm$samples$batch), 
                                       group=factor(xglm$samples$group))))
    saveRDS(xglm, file=output$rds)
}
analysis(snakemake@input, snakemake@output, snakemake@log)