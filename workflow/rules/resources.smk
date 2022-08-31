rule organism:
    output:
        directory(expand("resources/bioconductor/organism/lib/R/library/{organism}", organism = config["organism"]))
    params:
        name = config["organism"],
        path = lambda wc, output: Path(output[0]).parents[3]
    log:
        out = "logs/organism.out",
        err = "logs/organism.err"
    message:
        "Install organism package: {params.name}"
    shell:
        "conda create --quiet --yes --prefix {params.path} --strict-channel-priority --override-channels --channel conda-forge --channel bioconda --channel defaults bioconductor-{params.name}=3.13 1> {log.out} 2> {log.err}"
