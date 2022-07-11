 analysis=function(input, output) {
    library(edgeR)
    x=readRDS(input$rds)
    des = model.matrix(~x$samples$group)
    saveRDS(des,file=output$rds)
}

analysis(snakemake@input, snakemake@output)