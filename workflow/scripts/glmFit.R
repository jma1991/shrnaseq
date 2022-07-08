analysis=function(input, output) {
    library(edgeR)
    load(input$Rdata)
    fit = glmFit(xglm, des)
    saveRDS(fit,file=output$rds)
}

analysis(snakemake@input, snakemake@output)