analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    lrt=readRDS(input$rds[1])
    png(output$plot[1], width=2000, height=2000, res=400)
    plot(lrt$table$logFC, -10*log10(lrt$table$PValue), main="Volcano plot", xlab="M", ylab="-10*log(P-val)")
    dev.off()

    #batch corrected 
    corrected_lrt=readRDS(input$rds[2])
    png(output$plot[2], width=2000, height=2000, res=400)
    plot(corrected_lrt$table$logFC, -10*log10(corrected_lrt$table$PValue), main="Volcano plot", xlab="M", ylab="-10*log(P-val)")
    dev.off()

}
  
analysis(snakemake@input, snakemake@output, snakemake@log)