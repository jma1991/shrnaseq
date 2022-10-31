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
$ snakemake --use-conda --cores all
```
For further details on Snakemake, see the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/).

## Configuration

The `config/` directory must contain the following:

- `config.yaml`
- `samples.tsv`
- `guideRNA.tsv`
- one or more `.fastq` files

The workflow can be configured by editing the `config/config.yaml` file which must contain the following:

| Name | Description | Example |
| --- | --- | --- |
| organism | Specify your bioconductor organism package | org.Mm.eg.db |
| samples   | Path to sample file | config/samples.tsv |
| guideRNAs | Path to guideRNA sequence file  | config/guideRNA.tsv |
| fastq | Path to fastq file | config/subset.fastq |
| contrasts | Specify conditions to be contrasted |   - Day10-Day2 |
| FDR | Specify FDR threshold | 0.05 |
| FC | Specify fold change threshold | 1.5 |

Contrasts for the differential expression analysis need to be defined like the following:

```
contrast:
  Day10_vs_Day2:
    - Day10-Day2
```

The sample tsv file need to be formatted with the following columns as below:

Columns:

| Column | Description |
| --- | --- |
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
├── estimateDisp.rds
├── filter_guideRNAs.rds
├── glmFit.rds
├── model_matrix.rds
├── norm.rds
├── processAmplicons.rds
├── shiny.rds
├── {contrast}-camera.rds
├── {contrast}-camera.tsv
├── {contrast}-FDR_guideRNAs.rds
├── {contrast}-FDR-sig-guideRNAs.tsv
├── {contrast}-gene-level.rds
├── {contrast}-gene-level.tsv
├── {contrast}-glmLRT.rds
├── {contrast}-goana.rds
├── {contrast}-goana.tsv
├── {contrast}-kegg.rds
├── {contrast}-kegg.tsv
├── {contrast}-shiny.rds
├── {contrast}-top-goana.tsv
├── {contrast}-top-kegg.tsv
└── {contrast}-top-ranked-guideRNAs.tsv
```

```console
plots
├── BCV-plots.png
├── corrected-MDS-plot.png
├── corrected-PCA-plot.png
├── corrected-sample-dist-heatmap.png
├── counts-index-guideRNA.png
├── MDS-plot.png
├── PCA-plot.png
├── sample-dist-heatmap.png
├── {contrast}-corrected-expression-heatmap.png
├── {contrast}-expression-heatmap.png
├── {contrast}-guideRNA-histogram.png
├── {contrast}-plotSmear.png
└── {contrast}-volcano-plot.png
```

See below for details of each output file. Any contrast-specific files include the contrast name. 

#### Data processing
| File | Format | Description |
| --- | --- | --- |
| `processAmplicons.rds` | RDS | DGEList object |
| `model_matrix.rds` | RDS | Model design matrix |
| `corrected_counts.rds` | RDS | Batch corrected model design matrix |
| `contrasts_matrix.rds` | RDS | Contrasts matrix |

#### Quality control
| File | Format | Description |
| --- | --- | --- |
| `filter_guideRNAs.rds` | RDS | DGEList object |
| `counts-index-guideRNAs.png` | PNG | Counts of index and guideRNA sequences |
| `norm.rds` | RDS | DGEList object |
| `BCV-plots.png` | PNG | Biological coefficient of variation plot |
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
| `{contrast}-guideRNA-histogram.png` | PNG | Histogram of guideRNA p values |
| `{contrast}-top-ranked-guideRNAs.tsv` | TSV | Table of top ranked guideRNAs |
| `{contrast}-FDR-sig-guideRNAs.tsv` | TSV | Table of FDR significant guideRNAs |
| `{contrast}-FDR_guideRNAs.rds` | RDS | Vector of guideRNAs |
| `{contrast}-plotSmear.png` | PNG | Plots log-Fold Change versus log-CPM  |

#### Gene level analysis 
| File | Format | Description |
| --- | --- | --- |
| `{contrast}-camera.tsv` | TSV |  Table of gene level results  |
| `{contrast}-camera.rds` | RDS | Dataframe |
| `{contrast}-gene-level.tsv` | TSV |  Table of gene level results  |
| `{contrast}-gene-level.rds` | RDS | Dataframe |

#### Gene enrichment analysis 
| File | Format | Description |
| --- | --- | --- |
| `{contrast}-goana.tsv` | TSV |  Table of gene ontology results  |
| `{contrast}-goana.rds` | RDS | Dataframe |
| `{contrast}-kegg.tsv` | TSV |  Table of KEGG pathway results  |
| `{contrast}-kegg.rds` | RDS | Dataframe |
| `{contrast}-top_goana.tsv` | TSV | Table of top gene ontology results |
| `{contrast}-top_kegg.tsv` | TSV | Table of top KEGG pathway results |

#### Files for shiny application input
| File | Format | Description |
| --- | --- | --- |
| `shiny.rds` | RDS |  List object  |
| `{contrast}-shiny.rds` | RDS |  List object  |

## Tests

Test cases are in the `.test/integration` directory. They are automatically executed via
continuous integration with GitHub Actions.

## References

- Snakemake - Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33

### Bioconductor

- edgeR
- DESeq2
- DEFormats
- mixOmics
- limma
- AnnotationDbi
- GO.db
- SummarizedExperiment

### CRAN

- pheatmap
- RColorBrewer
- ggplot2
- matrixStats
- metap