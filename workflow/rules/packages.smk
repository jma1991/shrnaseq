rule install:
    output:
        "results/success.txt"
    conda:
        "../envs/bioconductor-edger.yaml"
    script:
        "../scripts/edgeR-package.R"