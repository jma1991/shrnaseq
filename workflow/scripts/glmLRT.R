analysis=function(input, output, params) {
    library(edgeR)
    matrix=readRDS(input$rds[1])
    fit=readRDS(input$rds[2])
    lrt = glmLRT(fit, contrast=matrix[, params$contrast])
    saveRDS(lrt,file=output$rds)
}

analysis(snakemake@input, snakemake@output, snakemake@params)