analysis=function(input, output, params, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    thresh = params$threshold   
    
    lrt=readRDS(input$rds[1])
    top2 = topTags(lrt, n=Inf)
    top2ids = top2$table[top2$table$FDR<thresh,1]
    write.table(top2ids, output$tsv[1], row.names=F, quote=F, col.names=F)
    saveRDS(top2ids, file=output$rds[1])

    #batch_corrected
    lrt=readRDS(input$rds[2])
    top2 = topTags(lrt, n=Inf)
    top2ids = top2$table[top2$table$FDR<thresh,1]
    write.table(top2ids, output$tsv[2], row.names=F, quote=F, col.names=F)
    saveRDS(top2ids, file=output$rds[2])


}

analysis(snakemake@input, snakemake@output, snakemake@params, snakemake@log)