analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    library(metap)
    lrt=readRDS(input$rds)

    dat = NULL

    for (i in unique(lrt$genes$Gene)) {
    sel = lrt$genes$Gene == i & !is.na(lrt$genes$Gene)
    if (sum(sel) > 1) {
        stouffers=sumz(lrt$table$PValue[which(sel)])[2]
        nhairpins=length(which(sel))
        pvals=cbind(i,nhairpins, stouffers)
        dat=rbind(dat,pvals)
    }
    }
    dat=cbind(dat,p.adjust(dat[,"stouffers"], method="fdr"))
    colnames(dat)=c("Gene", "nhairpins", "stouffers", "FDR")

    write.table(dat, output$tsv, quote=F, row.names=F)
    saveRDS(dat,file=output$rds)

}

analysis(snakemake@input, snakemake@output, snakemake@log)