analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    lrt=readRDS(input$rds)
    png(output$plot, width=3500, height=2000, res=400)
    hist(lrt$table$PValue, breaks = 40, main = "Histogram of hairpin P values", xlab="Hairpin p values")
    dev.off()
}
  
analysis(snakemake@input, snakemake@output, snakemake@log)