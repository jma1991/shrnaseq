analysis=function(input, output) {
    library(edgeR)
    lrt=readRDS(input$rds)
    thresh = 0.05
    top2 = topTags(lrt, n=Inf)
    top2ids = top2$table[top2$table$FDR<thresh,1]
    write.table(top2ids, output$txt, row.names=F, quote=F, col.names=F)
    saveRDS(top2ids, file=output$rds)
}

analysis(snakemake@input, snakemake@output)