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
As of July 2020, ASpli is freely available at Bioconductor devel branch. It will be part of the next Bioconductor official release scheduled for October 2020

    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
   
    # The following line initializes usage of Bioc devel branch and 
    # should not be necessary after the next official release scheduled for October 2020
    if(Sys.Date()<"2020-11-01") BiocManager::install(version='devel')
   
    BiocManager::install("ASpli")
   
    library(ASpli)
    


Note: **samtools** is also required for image creation when exporting integrated signals (reports can be generated without **samtools** if images are not required).

## Quick start
Here is an example for a pairwise comparison between 2 conditions (Control vs Treatment, 3 replicates each) using default parameters.
Extract features from genome, define targets data.frame with phenotype data, and mBAMs data.frame with phenotype data for merged BAMs:

    genome   <- loadDb("txdb.sqlite")
    features <- binGenome(genome)
    targets  <- data.frame(bam = c("CT_1.BAM", "CT_2.BAM","CT_3.BAM", "TR_1.BAM", "TR_2.BAM", "TR_3.BAM"), 
                           genotype = c( "CT", "CT", "CT", "TR", "TR", "TR" ), stringsAsFactors = FALSE )
    mBAMs <- data.frame(bam=c("CT.BAM", "TR.BAM"),condition=c("CT","TR"))
    
Read counting against annotated features:

    counts <- gbCounts(features = features, targets = targets, minReadLength = 125L, maxISize = 50000)
    
Junction-based de-novo counting:

    asd <- jCounts(counts = counts, features = features, minReadLength =125L)
    
Differential signal estimation:

    gb   <- gbDUreport(counts, contrast = c(-1, 1))
    jdur <- jDUreport(asd, contrast = c(-1, 1))
    
Report and signal integration:

    sr <- splicingReport(gb, jdur, counts)
    is <- integrateSignals(sr,asd)
    
Export results:

    exportSplicingReports( sr, output.dir="results")
    exportIntegratedSignals(is,sr=sr, output.dir = "results", counts=counts,
                              features=features,asd=asd, mergedBams = mBAMs)

## Documentation and help
Entry point for ASpli documentation is ASpli vignette, available after installing ASpli from R:

    browseVignettes("ASpli")
    
If user has a question not answered in ASpli vignette, ASpli has an issue board with previous issues available in gitlab.
If no previous issue answers the question, user can upload a new issue requesting help.
