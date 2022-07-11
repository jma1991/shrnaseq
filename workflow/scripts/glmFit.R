analysis=function(input, output) {
    library(edgeR)
    des=readRDS(input$rds[1])
    xglm=readRDS(input$rds[2])

    fit = glmFit(xglm, des)
    saveRDS(fit,file=output$rds)
}

analysis(snakemake@input, snakemake@output)