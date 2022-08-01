analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    lrt=readRDS(input$rds)

    dat=NULL
    for (i in unique(lrt$genes$Gene)) {
    sel = lrt$genes$Gene == i & !is.na(lrt$genes$Gene)
    logFC=mean(lrt$table$logFC[which(sel)])
    
    res=lrt$table[which(sel),]
    up=res[which(res$logFC>0),]
    down=res[which(res$logFC<0),]
    if ((min(up$Pvalue)>min(down$PValue)) | nrow(up)==0)
        dir="Down"
    if ((min(up$Pvalue)<min(down$PValue)) | nrow(down)==0)
        dir="Up"
    
    logFC=cbind(i,logFC, dir)
    dat=rbind(dat, logFC)
    
    }

    colnames(dat)=c("Gene", "mean-logFC", "Direction-(smallest-pval)")

    write.table(dat, output$tsv, quote=F, row.names=F)
    saveRDS(dat,file=output$rds)

}

analysis(snakemake@input, snakemake@output, snakemake@log)