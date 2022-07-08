analysis=function(input, output) {
    library(edgeR)
    fit=readRDS(input$rds)
    lrt = glmLRT(fit, contrast=c(0,0,-1,1))
    saveRDS(lrt,file=output$rds)
}

analysis(snakemake@input, snakemake@output)