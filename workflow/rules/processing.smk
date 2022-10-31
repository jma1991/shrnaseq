rule processAmplicons:
    input: 
        index=config["samples"],
        guideRNAs=config["guideRNAs"],
        fastq=config["fastq"]      
    output: 
        rds="results/processAmplicons.rds"
    params:
        readfile2=config["readfile2"],
        hairpinBeforeBarcode=config["hairpinBeforeBarcode"]
    log:
        out = "logs/processAmplicons.out",
        err = "logs/processAmplicons.err"
    message:
        "Process raw read data from genetic sequencing screens"
    conda:
        "../envs/bioconductor-edger.yaml"
    script:
        "../scripts/processAmplicons.R"