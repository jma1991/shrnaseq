analysis=function(input, output, params, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    library(ggplot2)
    lrt=readRDS(input$rds)
    top2 <- topTags(lrt, n=Inf)
    selY=top2$table[top2$table$FDR<params$FDR,]
    df <- data.frame(lrt$table)
    df$Guide <- rownames(df)
    colors <- c("FDR sig." = "red")

    plt=ggplot(df, aes(x=logFC, y=-10*log10(PValue), text=Guide)) + geom_point() +
    geom_point() + 
    geom_point(data = df[(df$Guide %in% selY$ID),], aes(color = "FDR sig.")) +
    labs(title = "Volcano plot", x="M", y = "-10*log(P-value)", color="") +
    scale_color_manual(values = colors) +
    theme_classic()

    png(output$plot, width=2000, height=2000, res=400)
    print(plt)
    dev.off()

}
  
analysis(snakemake@input, snakemake@output,snakemake@params, snakemake@log)