rule processAmplicons:
    input: 
        index=config["samples"],
        hairpins=config["hairpins"],
        fastq=config["fastq"]
    output: 
        rds="results/processAmplicons.rds"
    log:
        out = "logs/processAmplicons.out",
        err = "logs/processAmplicons.err"
    message:
        "Process raw read data from genetic sequencing screens"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/processAmplicons.R"