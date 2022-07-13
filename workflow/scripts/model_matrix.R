 analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    x=readRDS(input$rds)
    x$samples$group=factor(x$samples$group, levels = c(unique(x$samples$group)))
    des = model.matrix(~0+x$samples$group)
    colnames(des)= c(unique(x$samples$group))
    saveRDS(des,file=output$rds)
}

analysis(snakemake@input, snakemake@output, snakemake@log)