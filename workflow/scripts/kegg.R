analysis=function(input, output, params, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    keg <- kegga(qlf, species=params$species)
    topkegg <- topKEGG(keg, sort="up")
}
