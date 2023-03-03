analysis=function(input, output, params, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    library(ggplot2)
    camera <- read.table(input$tsv, header=T)
    colnames(camera) <- c("nGuides", "Direction", "Pvalue", "FDR", "Gene")
    res <- as.data.frame(camera)

    res$Status <- factor("NS", levels = c("Up", "NS", "Down"))
    res$Status[res$Direction == "Up" & res$FDR < params$FDR] <- "Up"
    res$Status[res$Direction == "Down" & res$FDR < params$FDR] <- "Down"
    res$Pvalue <- -log10(res$Pvalue)
    res$Rank <- rank(res$Pvalue)

    col <- c(
    "Up"   = "#FF0000",
    "NS"   = "#B8DBCC",
    "Down" = "#6495ED"
    )

    res.s <- res[res$FDR<params$FDR,]

    plt <- ggplot(res, aes(x = Rank, y = Pvalue, colour = Status, size = nGuides, label = Gene)) + 
    geom_point() + 
    geom_point(data = res.s, shape = 1, color="black") + 
    geom_text(data = res.s, hjust = 0, nudge_x = -20, color="black", size=3, check_overlap = T) +  
    scale_colour_manual(values = col, breaks = names(col)) + 
    theme_classic() +
    labs(title= "Camera rank",
        x = "Rank",
        y = "-log10(Pvalue)",
        colour = "Status"
    )
    png(output$plot, width=3000, height=3500, res=400)
    print(plt)
    dev.off()

}

analysis(snakemake@input, snakemake@output, snakemake@params, snakemake@log)