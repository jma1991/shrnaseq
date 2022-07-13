analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    des=readRDS(input$rds[1])
    xglm=readRDS(input$rds[2])

    fit = glmFit(xglm, des)
    saveRDS(fit,file=output$rds)
}

analysis(snakemake@input, snakemake@output, snakemake@log)