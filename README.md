# Sampling time-dependent artifacts in single-cell genomics studies

This repository contains all the scripts, notebooks and reports to reproduce the scRNA-seq analysis of our paper ["Sampling time-dependent artifacts in single-cell genomics studies"](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02032-0), published in Genome Biology in 2020. Here, we describe how to access the data, the most important packages and versions used, and how to navigate the directories and files in this repository.

## Data

All the raw data (fastqs) and expression matrices are available at the Gene Expression Omnibus (GEO) under [GSE132065](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132065). The data in this project can be broadly divided into 5 subprojects:

- [Smart-seq2](https://www.nature.com/articles/nprot.2014.006): includes a total of 4 96-well plates, with ids P2568, P2664, P2671 and P2672.
- [10X scRNA-seq](https://www.10xgenomics.com/products/single-cell/) data for Peripheral Blood Mononuclear Cells (PBMC): divided into two batches, which we named "JULIA_03" (cDNA libraries: AH9225 and AH9226, hashtag oligonucleotide (HTO) libraries: AH9223 and AH9224) and "JULIA_04" (cDNA libraries: AI0101 and AI0102, HTO libraries: AI0099 and AI0100).
- 10X scRNA-seq data for Chronic Lymphocytic Leukemia (CLL) cells: a total of 5 libraries, which are named after a combination of the donor id ("1220", "1472", "1892") and the temperature ("4ºC" or room temperature (RT)): 1220_RT, 1472_RT, 1892_RT, 1472_4C and 1892_RT.
- 10X scRNA-seq data for T-cell activation experiment (see methods): "Tcell_activation_day0_rep1", "Tcell_activation_day2_rep1", "Tcell_activation_day0_rep2" and "Tcell_activation_day1_rep2".
- [10X scATAC-seq](https://www.10xgenomics.com/products/single-cell-atac/) data for PBMC.
- 10X scATAC-seq data for CLL.


### Fastqs

As described in the paper, we multiplexed several sampling times into the same 10X Chip Channel using the [cell hashing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1) technology. To map the fastqs to the reference genome to obtain the single-cell gene expression matrices, we followed the ["Feature Barcoding Analysis"](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis) pipeline from [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger). This is an example of a cellranger run we used to map one of the libraries:

```{bash}
cellranger count --libraries libraries.csv --feature-ref feature_reference.csv --id 1472_RT --chemistry SC3Pv3 --expect-cells 5000 --localcores 12 --localmem 64 --transcriptome eference/human/refdata-cellranger-GRCh38-3.0.0/;
```

As you can see, a key input in this command is the feature_reference.csv which, according to 10X, "declares the set of Feature Barcoding reagents in use in the experiment. For each unique Feature Barcode used, this file declares a feature name and identifier, the unique Feature Barcode sequence associated with this reagent, and a pattern indicating how to extract the Feature Barcode sequence from the read sequence". This files can be easily created from the file "GSE132065_conditions_10X.tsv", available in both this GitHub repository and in GEO.


### Expression matrices
A total of 3 files per library are needed to reconstruct the full expression matrix:

1. barcodes*.tsv.gz: corresponds to the cell barcodes (column names).
2. features*.tsv.gz: corresponds to the gene/condition identifiers (row names). Moreover, it contains a columns that ideantifes genes ("Gene Expression") and experimental conditions ("Antibody Capture").
3. matrix*mtx.gz: expression matrix in sparse format.


To make our data as FAIR (findable, accessible, interoperable, reusable) as possible, we have deposited the gene expression matrices and the Seurat objects that are saved in each of the Rmarkdown notebooks in [this Zenodo respository](https://zenodo.org/record/7308457#.Y2veY-zMJAc). One can download it in 3 lines of code:

```{bash}
wget https://zenodo.org/record/7308457/files/MassoniBadosa2020_GenomeBiol_scRNAseq_data.zip
unzip MassoniBadosa2020_GenomeBiol_scRNAseq_data.zip
cd MassoniBadosa2020_GenomeBiol_scRNAseq_data
```

The next step after downloading it should be reading the README.md that is inside the MassoniBadosa2020_GenomeBiol_scRNAseq_data folder.


## Package versions

These are the versions of the most important packages used throughout all the analysis:

CRAN:

* [tidyverse 1.3.0](https://cran.r-project.org/web/packages/tidyverse/vignettes/paper.html)
* [Seurat 3.1.4](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419305598%3Fshowall%3Dtrue)

Bioconductor:

* [Scater 1.10.1](https://academic.oup.com/bioinformatics/article/33/8/1179/2907823)
* [Scran 1.10.2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0947-7)
* [GOstats 2.48.0](https://academic.oup.com/bioinformatics/article/23/2/257/204776)

**Note**: Two months before compiling the notebooks to release them together with the paper, we updated most Bioconductor packages. Thus, some versions reported in the sessionInfo() of the notebooks might be slightly different to the ones used to produce the figures of the article.


## File system and name scheme

This repository contains 4 different analysis directories (which correspond to the main blocks of the article) and 1 directory with the scripts to produce the figures of the article:

* 1-PBMC
* 2-CLL
* 3-T_cell_activation
* 4-Revision
* figures_scripts

The first 3 have a set of similar notebooks, which match the common pre-processing steps of any single-cell expression matrix:

1. Demultiplexing: classify each cell to its original condition based on the expression of HTO.
2. QC and normalization: filter out poor-quality cells and genes and normalize expression counts.
3. Dimensionality reduction, clustering and annotation of cell types.

Each notebook (\*.Rmd) has an associated report (\*.html). The reports are useful to visualize the results of each section as well as the diagnostic plots that we used to set the thresholds and parameters. For a quick inspection, one can copy the URL of the report in the [GitHub & BitBucket HTML Preview](https://htmlpreview.github.io/) (note that side bar and table will not be available with it).

Finally, the figures_scripts directory contains most of the scripts needed to produce the figures as they appear in the paper. The remaining supplementary figures are created either in the notebooks in 4-Revision, or were created by other coauthors.


## Other studies

Here is a list of other important benchmarking studies:

* [Single-cell sequencing reveals dissociation-induced gene expression in tissue subpopulations](https://www.nature.com/articles/nmeth.4437).
* [scRNA-seq assessment of the human lung, spleen, and esophagus tissue stability after cold preservation](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1906-x).
* [Dissociation of solid tumor tissues with cold active protease for single-cell RNA-seq minimizes conserved collagenase-associated stress responses](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1830-0).
* [Expression levels for many genes in human peripheral blood cells are highly sensitive to ex vivo incubation](https://www.nature.com/articles/6364098).
* [Sample processing obscures cancer-specific alterations in leukemic transcriptomes](https://www.pnas.org/content/111/47/16802).
* [A single-cell and single-nucleus RNA-Seq toolbox for fresh and frozen human tumors](https://www.nature.com/articles/s41591-020-0844-1).
* [Multi-omic Single-cell Atlas Reveals Broad Ex Vivo Effects on Human Immunobiology](https://www.biorxiv.org/content/10.1101/2020.10.18.344663v4).
* [ACME dissociation: a versatile cell fixation-dissociation method for single-cell transcriptomics](https://www.biorxiv.org/content/10.1101/2020.05.26.117234v2).
* [Systematic assessment of tissue dissociation and storage biases in single-cell and single-nucleus RNA-seq workflows](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02048-6).
* [Single Cell Sequencing Reveals Glial Specific Responses to Tissue Processing & Enzymatic Dissociation in Mice and Humans](https://www.biorxiv.org/content/10.1101/2020.12.03.408542v1).
* [Cryopreservation of human cancers conserves tumour heterogeneity for single-cell multi-omics analysis](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00885-z?elqTrackId=283abc61b2f44964b72a610cb61f2b11&elq=2137b432e40849a3a608a59ec01ee7d3&elqaid=31262&elqat=1&elqCampaignId=10610)
* [Single-cell transcriptome conservation in cryopreserved cells and tissues](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1171-9)
* [Stress relief: emerging methods to mitigate dissociation-induced artefacts](https://www.cell.com/trends/cell-biology/fulltext/S0962-8924(21)00096-9)
