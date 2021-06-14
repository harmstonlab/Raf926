# Raf_DorsalVentralPatterning
## Data
* All data for this analysis can obtained from the server at `/home/data/projects/raf_study`. The aligned and quantified data can also be found in this repository under the folder `data/aligned`.
* Data from the [Koenecke et al. paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1057-2) can be found on GEO with the accession number [GSE68983](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68983). The data was downloaded using the `prefetch` command from the SRA Toolkit. The aligned and quantified data can also be found in this repository under the folder `data/koenecke_rnaseq`. 
* Annotation file used: `/home/shared/genomes/dm6/StarIndex/ensembl97/Drosophila_melanogaster.BDGP6.22.97.chr.gtf`
* Reference fasta file used: `/home/shared/genomes/dm6/StarIndex/ensembl97/dm6.fa`
* The dorsal target genes were obtained from [Stathopoulos et al.](https://www.cell.com/cell/fulltext/S0092-8674(02)01087-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867402010875%3Fshowall%3Dtrue)
* The TGF-beta marker genes were obtained from Table S2 of [Dominguez et al.](https://www-sciencedirect-com.libproxy1.nus.edu.sg/science/article/pii/S0378111916305418#ec0010)


## Analysis
The analysis was carried out using only the wild type samples (WT_4, WT_5, WT_6) and the mutant loss-of-function samples (S926glc_4, S926glc5, S926glc_6), and using the `QC.Rmd`, `DE.Rmd`, and `V_E.Rmd` scripts for quality control, differential expression, and visualiation and enrichment respectively. 

If running from scratch, these scripts can be found in the `analysis` folder and need to be run in the sequential order of `QC.Rmd` > `DE.Rmd` > `V_E.Rmd` as each script generates files required by the subsequent scripts. If not running from scratch, the outputs of each script can be found in their respective `output` subfolder in the `analysis` main folder. If the repository is cloned as is, and tje `.Rproj` file is used,  there should have no issues running the scripts.

The same folder structure applies to the subanalyses.

### Subanalysis
#### Analysing Mutant_GOF Samples 
For this sub-analysis, the same pairwise analysis as above was carried out between wild type samples and mutant gain-of-function samples (RafGOF_4, RafGOF_5, RafGOF_6). While the dorsal targets in the ectoderm do display the expected pattern (downregulated in GOF mutant), the same pattern could be seen in the dorsal targets of the neuroectoderm as well as mesoderm. To avoid cherry picking, the results of this sub-analysis were dropped from the manuscript.

The `.Rmd` scripts and `.html` files for this analysis can be found in the `analysis_GOF` folder.

#### Comparing dorsal targets to Koenecke et al's dataset
Here we looked at the RNA-seq data from [Koenecke et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1057-2). The data was aligned using STAR and RSEM (see Versions/Koenecke Analysis for versions) via a Snakefile, which can be found in the `analysis_koenecke` folder. Again, the same three scripts were used, and in the same order -- `QC.Rmd`, `DE.Rmd`, and `V_E.Rmd`.

## Versions
### Shell
```
4.2.46(2)-release x86_64-redhat-linux-gnu
```

### R
```
## R version 4.0.5 (2021-03-31)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
```

### Main Analysis
* **STAR version:** 2.7.1a
* **RSEM version:** v1.3.1

### Koenecke Analysis
* **STAR version:** 2.7.9a
* **RSEM version:** v1.3.1

The versions of R packages used in the analysis can be found in the **Session Info** section of each `.html` file.
