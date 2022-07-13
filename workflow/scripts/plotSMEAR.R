analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
  library(edgeR)
  lrt=readRDS(input$rds[1])
  top2ids=readRDS(input$rds[2])
  png(output$plot, width=2000, height=2000, res=400)
      plotSmear(lrt, de.tags=top2ids,
        pch=20, cex=0.6, main="Another small screen: logFC vs logCPM")
      abline(h = c(-1, 0, 1), col = c("dodgerblue", "yellow", "dodgerblue"), lty=2)
    dev.off()
}

analysis(snakemake@input, snakemake@output, snakemake@log)