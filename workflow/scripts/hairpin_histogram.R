analysis=function(input, output) {
    library(edgeR)
    lrt=readRDS(input$rds)
    png(output$plot, width=2000, height=2000, res=400)
    hist(lrt$table$PValue,, breaks = 30, main = "Histogram of hairpin P values", xlab="Hairpin p values")
    dev.off()
}
  
analysis(snakemake@input, snakemake@output)