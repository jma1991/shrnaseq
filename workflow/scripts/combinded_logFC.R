analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    xglm=readRDS(input$rds[1])
    lrt=readRDS(input$rds[2])

    combinded_logFC=NULL
    xglm$genes$Gene=as.character(xglm$genes$Gene)
    for (i in unique(xglm$genes$Gene)) {
    sel = xglm$genes$Gene == i & !is.na(xglm$genes$Gene)
    if (sum(sel) > 1) {
        logFC=mean(lrt$table$logFC[which(sel)])
        logFC=cbind(i,logFC)
        combinded_logFC=rbind(combinded_logFC, logFC)
     }
    }

    colnames(combinded_logFC)=c("Gene", "combinded-logFC")

    write.table(combinded_logFC, output$tsv, quote=F, row.names=F)
    saveRDS(combinded_logFC,file=output$rds)

}

analysis(snakemake@input, snakemake@output, snakemake@log)