# Documentation

## Table of contents
* [Usage](#Usage)
* [Configuration](#Configuration)
* [Output](#Output)
* [Tests](#Tests)
* [References](#References)

## Usage 

The workflow can be excuted using the following command: 

```console
$ snakemake --use-conda --cores 1
```
For further details on Snakemake, see the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/).

## Configuration

The workflis configured by editing the following file:

- `config/config/yaml`

The config file must contain the following: 

| Name | Description | Example |
| --- | --- | --- |
| organism | Specify your organism | Mus musculus |
| samples   | Path to sample file | config/samples.tsv |
| hairpins | Path to hairpin sequence file  | config/hairpin.tsv |
| fastq | Path to fastq file | config/subset.fastq |
| contrasts | Specify conditions to be contrasted |   - Day2-Day10 |
| FDR | Specify FDR threshold | 0.05 |
| FC | Specify fold change threshold | 1.5 |

Contrasts for the differential expression analysis need to be defined like the following:

```
contrast:
  Day2_vs_Day10:
    - Day2-Day10
```

The sample tsv file need to be formatted with the following columns as below:

Columns:

| Column | Description |
| --- | --- | |
| ID | Sample ID |
| Sequences | Index sequences |
| group | Sample condition | 
| batch | Sample batch number | 

Example file format:

| ID | Sequences | group | batch |
| --- | --- | --- | --- |
| 1 | GAAA | Day2 | 1 |
| 2 | GAAC | Day10 | 2 |


## Output

After the workflow is successfully run the output is found in results and plots directory and will contain the following:

```console
results
├── contrasts_matrix.rds
├── corrected_counts.rds
├── camera.rds
├── camera.tsv
├── FDR_hairpins.rds
├── FDR-sig-hairpins.tsv
├── glmLRT.rds
├── top-ranked-hairpins.tsv
├── diff_rep_analysis.rds
├── filter_hairpins.rds
├── glmFit.rds
├── model_matrix.rds
└── processAmplicons.rds
```

```console
plots
├── BCV-plots.png
├── counts-index-hairpin.png
├── expression-heatmap.png
├── hairpin-histogram.png
├── MDS-plot.png
├── PCA-plot.png
├── plotSmear.png
├── sample-dist-heatmap.png
└── volcano-plot.png
```

See below for details of each output file. Any contrast-specific files include the contrast name. 

#### Data processing
| File | Format | Description |
| --- | --- | --- |
| `processAmplicons.rds` | RDS | DGEList object |
| `filter_hairpins.rds` | RDS | DGEList object |
| `counts-index-hairpins.png` | PNG | Counts of index and hairpin sequences |
| `model_matrix.rds` | RDS | Model design matrix |
| `corrected_counts.rds` | RDS | Batch corrected model design matrix |
| `contrasts_matrix.rds` | RDS | Contrasts matrix |

#### Quality control
| File | Format | Description |
| --- | --- | --- |
| `BCV-plots.png` | PNG | Biological coefficient of variation plot |
| `counts-index-hairpins.png` | PNG | Counts of index and hairpin sequences |
| `MDS-plot.png` | PNG | MDS plot |
| `corrected-MDS-plot.png` | PNG | Batch corrected MDS plot |
| `PCA-plot.png` | PNG | Principal component analysis plot |
| `corrected-PCA-plot.png` | PNG | Batch corrected PCA plot |
| `sample-dist-heatmap.png` | PNG | Heatmap of sample distances |
| `corrected-sample-dist-heatmap.png` | PNG | Batch corrected heatmap of sample distances |

#### Differential expression analysis 
| File | Format | Description |
| --- | --- | --- |
| `estimateDisp.rds` | RDS |  DGEList object |
| `glmFit.rds` | RDS | DGEGLM object |
| `{contrast}-glmLRT.rds` | RDS | DGELRT object |
| `{contrast}-expression-heatmap.png` | PNG | Heatmap of differential expression |
| `corrected-{contrast}-expression-heatmap.png` | PNG | Batch corrected heatmap of differential expression |
| `{contrast}-volcano-plot.png` | PNG | Volcano plot |
| `{contrast}-hairpin-histogram.png` | PNG | Histogram of hairpin p values |
| `{contrast}-top-ranked-hairpins.tsv` | TSV | Table of top ranked hairpins |
| `{contrast}-FDR-sig-hairpins.tsv` | TSV | Table of FDR significant hairpins |
| `{contrast}-FDR_hairpins.rds` | RDS | Vector of hairpins |
| `{contrast}-plotSmear.png` | PNG | Plots log-Fold Change versus log-CPM  |

#### Gene set analysis 
| File | Format | Description |
| --- | --- | --- |
| `{contrast}-camera.tsv` | TSV|  Table of gene results  |
| `{contrast}-camera.rds` | RDS | Dataframe |

## Tests

Test cases are in the `.test` directory. They are automatically executed via
continuous integration with GitHub Actions.
## References

- Snakemake - Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33

### Bioconductor

- edgeR
- DESeq2
- DEFormats
- mixOmics
- limma

### CRAN

- pheatmap
- RColorBrewer
- ggplot2
- matrixStats