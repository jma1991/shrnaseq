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
    hist(lrt$table$PValue, breaks = 30, main = "Histogram of hairpin P values", xlab="Hairpin p values")
    dev.off()

    corrected_lrt=readRDS(input$rds[2])
    png(output$plot[2], width=2000, height=2000, res=400)
    hist(corrected_lrt$table$PValue, breaks = 30, main = "Histogram of hairpin P values", xlab="Hairpin p values")
    dev.off()
}
  
analysis(snakemake@input, snakemake@output, snakemake@log)