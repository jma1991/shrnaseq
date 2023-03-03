 analysis=function(input, output, params, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    library(ggplot2)
    genelevel <- read.table(input$tsv, header=T, sep=",")

    colnames(genelevel) <- c("Gene", "nGuides", "Mean logFC", "IQR logFC", "Direction mean logFC",	"Direction smallest Pvalue", "Stouffer's Pvalue", "Stouffer's FDR")
    
    res <- as.data.frame(genelevel)
    FC_index=c(paste0("> ", params$FC), paste0("< ", params$FC, " and > -", params$FC), paste0("< -", params$FC))

    res$LogFC <- factor(FC_index[2], levels = FC_index)
    res$LogFC[res$`Mean logFC`>params$FC] <- FC_index[1]
    res$LogFC[res$`Mean logFC`<(-(params$FC))] <- FC_index[3]
    res$Rank <- rank(res$`Mean logFC`)

    col <- c( "#FF0000","#B8DBCC", "#6495ED")
    names(col)=c(FC_index)
    res.up <-  res[order(-res$`Mean logFC`),]
    res.up <- res.up[res.up$`Mean logFC`>params$FC,][c(1:5),]
    res.down <- res[order(res$`Mean logFC`),]
    res.down <- res.down[res.down$`Mean logFC`<(-params$FC),][c(1:5),]

    plt <- ggplot(res, aes(x = Rank, y = `Mean logFC`, colour = LogFC, size = nGuides, label = Gene)) + 
      geom_point() + 
      geom_point(data = res.up, shape = 1, color="black") + 
      geom_point(data = res.down, shape = 1, color="black") + 
      geom_text(data = res.up, hjust = 0, nudge_x = -20, color="black", size=3, check_overlap = T) +  
      geom_text(data = res.down, hjust = 0, nudge_x = -20, color="black", size=3, check_overlap = T) + 
      theme_classic() +
      scale_colour_manual(values = col, breaks = names(col)) + 
      labs(title = "Gene rank",
        x = "Rank",
        y = "Mean log fold-change",
        colour = "LogFC"
      )  
    png(output$plot, width=3000, height=3500, res=400)
    print(plt)
    dev.off()

}

analysis(snakemake@input, snakemake@output, snakemake@params, snakemake@log)