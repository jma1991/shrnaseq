 analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)

    x=readRDS(input$rds)

    x$samples$group=factor(x$samples$group, levels = c(unique(x$samples$group)))
    
    # Set group factor

    if ("group" %in% names(x$samples)) {
        group <- factor(x$samples$group)
        n.group <- nlevels(group)
        is.group <- n.group > 1
        } else {
        is.group <- FALSE
        }

    # Set batch factor

    if ("batch" %in% names(x$samples)) {
        batch <- factor(x$samples$batch)
        n.batch <- nlevels(batch)
        is.batch <- n.batch > 1
        } else {
        is.batch <- FALSE
        }

    # Construct design matrix

    if (is.group & !is.batch) {
        des <- model.matrix(~ 0 + group)
        }
    if (is.group & is.batch) {
        des <- model.matrix(~ 0 + group + batch)    
    }

    # Rename group coefficients

    which.group <- seq_len(n.group)
    colnames(des)[which.group] <- levels(group)      
    saveRDS(des,file=output$rds)
}

analysis(snakemake@input, snakemake@output, snakemake@log)