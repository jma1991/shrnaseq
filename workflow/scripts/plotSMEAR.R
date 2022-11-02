analysis=function(input, output, params, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
  library(edgeR)
  library(ggplot2)
  lrt=readRDS(input$rds[1])
  top2 <- topTags(lrt, n=Inf)
  top2ids <- top2$table[(top2$table$logFC>params$FC | top2$table$logFC<(-(params$FC))),1]
  selY=top2$table[top2$table$FDR<params$FDR,]
  df <- data.frame(lrt$table)
  df$Guide <- rownames(df)
  colors <- c("FDR sig." = "red")

  plt=ggplot(df, aes(x=logCPM, y=logFC, text=Guide)) +
    geom_point(color = "#b8dbcc") +
    geom_point(data = df[(row.names(df) %in% top2ids),], color = "#000000") +
    geom_point(data = df[(df$Guide %in% selY$ID),], aes(color = "FDR sig.")) +
    geom_hline(yintercept=(-(params$FC)), linetype="dashed", color="#b8dbcc") +
    geom_hline(yintercept=0,  color="cornflowerblue") +
    geom_hline(yintercept=(params$FC), linetype="dashed", color="#b8dbcc") +
    labs(title = "Mean difference plot",
    color="") +
    scale_color_manual(values = colors) +
    theme_classic()
  png(output$plot, width=2000, height=2000, res=400)
    print(plt)
  dev.off()
}


analysis(snakemake@input, snakemake@output,snakemake@params, snakemake@log)