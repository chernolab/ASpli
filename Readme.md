# ASpli

# An integrative R package for analysing alternative splicing using RNAseq

## Authors
Estefania Mancini, Andrés Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz

## Introduction
Alternative splicing (AS) is a common mechanism of post-transcriptional gene 
regulation in eukaryotic organisms that expands the functional and regulatory 
diversity of a single gene by generating multiple mRNA isoforms that encode 
structurally and functionally distinct proteins. 

Genome-wide analysis of AS has been a very active field of research since 
the early days of NGS (Next generation sequencing) technologies.  Since then, evergrowing data availability and the development of increasingly sophisticated analysis methods have uncovered the complexity of the general splicing repertoire.  

`ASpli` was specifically designed to integrate several independent signals in order to deal with the complexity that might arise in splicing patterns. Taking into account genome annotation information, `ASpli` considers bin-based signals along  with junction inclusion indexes in order to assess for statistically significant changes in read coverage. In addition, annotation-independent signals are estimated based on the complete set of experimentally detected splice junctions.  `ASpli` makes use of a generalized linear model framework (as implemented in `edgeR` R-package) to assess for the statistical  analysis of specific contrasts of interest. In this way, `ASpli` can provide a comprehensive description of genome-wide splicing alterations even for complex experimental designs. 

A typical `ASpli` workflow  involves: parsing the genome annotation into subgenic features called bins, overlapping read alignments against them, perform junction counting, fulfill inference tasks of differential bin and junction usage and, finally, report integrated splicing signals. At every step `ASpli` generates self-contained outcomes that, if required, can be easily exported and integrated into other processing pipelines. 

## Installation

    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
   
    BiocManager::install("ASpli")
   
    library(ASpli)
    


Note: **samtools** is also required for image creation when exporting integrated signals (reports can be generated without **samtools** if images are not required).

## Quick start
ASpli provides toy BAM and GTF files to introduce the working pipeline.
Here is an example for a pairwise comparison between 2 conditions (Control vs Treatment, 3 replicates each) using default parameters.

Extract features from genome, define *targets* data.frame with phenotype data, and *mBAMs* data.frame with phenotype data for merged BAMs:

```
library(ASpli)
library(GenomicFeatures)

# gtf preprocessing ----
gtfFileName <- aspliExampleGTF()
genomeTxDb  <- makeTxDbFromGFF( gtfFileName )

# feature extraction ----
features    <- binGenome( genomeTxDb )

#bams and target file ----
BAMFiles <- aspliExampleBamList()
targets  <- data.frame(row.names = paste0('Sample',c(1:6)),
                       bam = BAMFiles[1:6],
                       f1  = c( 'control','control','control','treatment','treatment','treatment'),
                       stringsAsFactors = FALSE)
mBAMs <- data.frame( bam = sub("_[012]","",targets$bam[c(1,4)]),
                     condition = c("control","treatment"))
```


Read counting against annotated features:
```
gbcounts <- gbCounts(features=features, targets=targets,
                     minReadLength = 100, maxISize = 50000)
gbcounts
```


Junction-based *de-novo* counting and splicing signal estimation:

```
asd <- jCounts(counts=gbcounts, features=features, minReadLength=100)
asd
```

Differential gene expression and bin usage signal estimation:
```
gb  <- gbDUreport(gbcounts, contrast = c(-1,1))
gb
```

Differential junction usage analysis:
```
jdur <- jDUreport(asd, contrast=c(-1,1))
jdur
```

Bin and junction signal integration:
```
sr <- splicingReport(gb, jdur, counts=gbcounts)
```

Summary of integration of splicing signals along genomic-regions. 
```
is <- integrateSignals(sr,asd)
```

Export results:
```
exportIntegratedSignals(is,sr=sr,
                          output.dir = "aspliExample",
                          counts=gbcounts,features=features,asd=asd,
                          mergedBams = mBAMs)
```

## Documentation and help
Entry point for ASpli documentation is ASpli vignette, available after installing ASpli from R:

    browseVignettes("ASpli")
    
If user has a question not answered in ASpli vignette, ASpli has an issue board with previous issues available in https://github.com/chernolab/ASpli.
If no previous issue answers the question, user can upload a new issue requesting help.

## Citation
*ASpli: Integrative analysis of splicing landscapes through RNA-seq assays*, Mancini E, Rabinovich A, Iserte J, Yanovsky M, Chernomoretz A, Bioinformatics, March 2021, 


