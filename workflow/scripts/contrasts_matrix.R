analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    des=readRDS(input$rds)
    conditions <- colnames(des)
    combinations <- expand.grid(conditions, conditions)
    contrasts <- apply(combinations, 1, paste, collapse = "-")
    matrix <- makeContrasts(contrasts = contrasts, levels = conditions)
    saveRDS(matrix,file=output$rds)
}

analysis(snakemake@input, snakemake@output, snakemake@log)