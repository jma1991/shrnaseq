analysis=function(input, output, params, log) {
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
      pch=20, cex=0.6, main="logFC vs logCPM")
    abline(h = c(-(params$FC), 0, (params$FC)), col = c("dodgerblue", "yellow", "dodgerblue"), lty=2)
  dev.off()
}


analysis(snakemake@input, snakemake@output,snakemake@params, snakemake@log)