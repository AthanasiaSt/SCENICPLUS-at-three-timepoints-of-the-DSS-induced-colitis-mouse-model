## Single-cell Gene Regulatory Network analysis with SCENIC+ at three timepoints of the DSS-induced colitis mouse model

This repository contains the code used to create all the plots associated with our manuscript: **Distinct enhancer-driven transcriptional networks shape intestinal fibroblast identities and regeneration-associated activation** 

In-house raw data from this work are deposited in GEO with the following accession numbers: 

**scRNA-seq**
  - Control colon (D0): GSE296873
  - Regenerating colon (D14): GSE304016
  - Acute inflammation (D7k) and control (D0k), from *Kinchen et al.*: GSE114374
  
**scATAC-seq**
  - Homeostasis (D0), acute inflammation (D7), regeneration (D14): GSE304192


Processed files are deposited in zenodo: **10.5281/zenodo.17258163**

## Methods Summary

### Preprocessing

-   Raw scRNA FASTQ files were processed using **10X Genomics CellRanger v7.1.0** (`cellranger count`, default parameters)
-   Raw scRNA FASTQ files from Kinchen dataset were downloaded and processed using **10X Genomics CellRanger v7.1.0** (`cellranger aggr`, default parameters including sequencing depth normalization)
-   Raw scATAC FASTQ files were processed using **10X Genomics CellRanger ATAC v2.1.0** (`cellranger-atac count`, default parameters)
-   Mouse reference transcriptome: `mm10` 

### Downstream Analysis

-   scRNA-seq QC, Harmony integration (In-house/Kinchen), clustering, plotting: **Seurat v5.1.0**
-   scATAC-seq QC, clustering, plotting, peaks, motifs: **ArchR v1.0.2**
-   Gene Regulatory Network inference: **SCENIC+ v1.0a1**

<sub>For more details on our methodology, see our manuscript ....
------------------------------------------------------------------------

## Citation

If you use this dataset, please cite the associated manuscript:.............
