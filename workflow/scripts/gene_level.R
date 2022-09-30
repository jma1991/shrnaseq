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
    # gene level
    dat=NULL

    for (i in unique(lrt$genes$Gene)) {

        sel = lrt$genes$Gene == i & !is.na(lrt$genes$Gene)
        
        #number of hairpins
        nhairpin=length(which(sel))
        
        if (nhairpin==1) {
            logFC=lrt$table$logFC[which(sel)]
            IQRFC=NA
            if (logFC>0) {
                dir="Up" 
            } else {
                dir="Down"
            }
            if (logFC>0) {
                dir_pval="Up" 
            } else {
                dir_pval="Down"
            }
            stouffers=NA
        }
        
        if (nhairpin>1) {
            #calculate average logFC
            logFC=mean(lrt$table$logFC[which(sel)])
            IQRFC=IQR(lrt$table$logFC[which(sel)])
            if (logFC >0) {
                dir="Up"
                } else {
                dir="Down"
                }
            
            #select up regulated hairpins
            res=lrt$table[which(sel),]
            up=res[which(res$logFC>0),]
            down=res[which(res$logFC<0),]
            options(warn=-1)
            if ((min(up$Pvalue)>min(down$PValue)) | nrow(up)==0)
                dir_pval="Down"
            if ((min(up$Pvalue)<min(down$PValue)) | nrow(down)==0)
                dir_pval="Up"
            options(warn=-0)
            pvals=lrt$table$PValue[which(sel)]
            pvals=pvals[!pvals==1]
            if (length(pvals)>1) {
            stouffers=as.numeric(sumz(lrt$table$PValue[which(sel)])[2])
            } else {
            stouffers=NA
            }
        }
        vector=cbind(i,nhairpin, logFC, IQRFC, dir, dir_pval, stouffers)
        dat=rbind(dat, vector)
    }
    dat=cbind(dat,p.adjust(dat[,"stouffers"], method="fdr"))
    colnames(dat)=c("gene", "nguides", "mean.logFC", "iqr.logFC", "direction.mean.logFC",
     "direction.smallest.pvalue", "stouffers.pvalue", "stouffers.FDR")

    write.table(dat, output$tsv, quote=F, row.names=F, sep=",")
    saveRDS(dat,file=output$rds)

}

analysis(snakemake@input, snakemake@output, snakemake@log)