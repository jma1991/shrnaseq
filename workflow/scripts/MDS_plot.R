analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script 
    library(edgeR)
    library(RColorBrewer)
    x=readRDS(input$rds[1])

    colors <- brewer.pal(length(unique(x$samples$group)),"Set3")

    png(output$plot[1], width=2500, height=1800, res=400)
    par(mar=c(5,5,3,7), xpd=TRUE)
    plotMDS(x, labels=x$samples$group, col=colors,cex=0.8, main="MDS Plot")
        legend("topright", inset=c(-0.35,0), legend=c(unique(x$samples$group)), 
        col=colors, pch=15, box.lwd = 0,box.col = "white",bg = "white")
    dev.off()
    
    ##batch corrected
    mat=readRDS(input$rds[2])
   
    colors <- brewer.pal(length(unique(x$samples$group)),"Set3")

    png(output$plot[2], width=2500, height=1800, res=400)
    par(mar=c(5,5,3,7), xpd=TRUE)
    plotMDS(mat, labels=x$samples$group, col=colors,cex=0.8, main="Batch corrected MDS Plot")
        legend("topright", inset=c(-0.35,0), legend=c(unique(x$samples$group)), 
        col=colors, pch=15,box.lwd = 0,box.col = "white",bg = "white")
    dev.off()
}

analysis(snakemake@input, snakemake@output, snakemake@log)