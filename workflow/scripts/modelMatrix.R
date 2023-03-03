# Author: James Ashmore
# Copyright: Copyright 20203, Zifo Technologies Ltd.
# Email: james.ashmore@zifornd.com
# License: MIT

modelMatrix <- function(object) {

    # Get samples data

    data <- object$samples

    names <- colnames(data)

    # Set condition factor

    if ("Condition" %in% names) {

        condition <- factor(data$Condition)

        n.condition <- nlevels(condition)

        is.condition <- n.condition > 1

    } else {

        is.condition <- FALSE

    }

    # Set batch factor

    if ("Batch" %in% names) {

        batch <- factor(data$Batch)

        n.batch <- nlevels(batch)

        is.batch <- n.batch > 1

    } else {

        is.batch <- FALSE

    }

    # Set batch2 factor

    if ("Batch2" %in% names) {

        batch2 <- factor(data$Batch2)

        n.batch2 <- nlevels(batch2)

        is.batch2 <- n.batch2 > 1

    } else {

        is.batch2 <- FALSE

    }

    # Construct design matrix

    if (is.condition & !is.batch & !is.batch2) {

        design <- model.matrix(~ 0 + condition)

    }

    if (is.condition & is.batch & !is.batch2) {

        design <- model.matrix(~ 0 + condition + batch)

    }

    if (is.condition & !is.batch & is.batch2) {

        design <- model.matrix(~ 0 + condition + batch2)

    }

    if (is.condition & is.batch & is.batch2) {

        design <- model.matrix(~ 0 + condition + batch + batch2)

    }

    # Rename condition coefficients

    which.condition <- seq_len(n.condition)

    colnames(design)[which.condition] <- levels(condition)

    # Return design matrix

    design

}

main <- function(input, output, log) {

    # Log

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script

    library(edgeR)

    object <- readRDS(input$rds)

    design <- modelMatrix(object)

    saveRDS(design, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@log)
