analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    library(metap)
    xglm=readRDS(input$rds[1])
    lrt=readRDS(input$rds[2])

    stouffers = NULL
    xglm$genes$Gene=as.character(xglm$genes$Gene)
    for (i in unique(xglm$genes$Gene)) {
    sel = xglm$genes$Gene == i & !is.na(xglm$genes$Gene)
    if (sum(sel) > 1) {
        p=sumz(lrt$table$PValue[which(sel)])
        pvals=cbind(i,p$p)
        stouffers=rbind(stouffers,pvals)
    }
    }
    colnames(stouffers)=c("Gene", "Stouffer's")

    write.table(stouffers, output$tsv, quote=F, row.names=F)
    saveRDS(stouffers,file=output$rds)

}

analysis(snakemake@input, snakemake@output, snakemake@log)