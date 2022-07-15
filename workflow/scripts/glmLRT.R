analysis=function(input, output, params, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    matrix=readRDS(input$rds[1])
    fit=readRDS(input$rds[2])
    lrt = glmLRT(fit, contrast=matrix[, params$contrast])
    saveRDS(lrt,file=output$rds[1])

    #batch corrected
    corrected_fit=readRDS(input$rds[3])
    corrected_lrt = glmLRT(corrected_fit, contrast=matrix[, params$contrast])
    saveRDS(corrected_lrt,file=output$rds[2])

}

analysis(snakemake@input, snakemake@output, snakemake@params, snakemake@log)