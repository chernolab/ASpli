\name{gbDUreport}

\alias{gbDUreport}

\title{
  Differential gene expression and differential bin usage estimation
}
\description{
  Estimate differential expression at gene level and differential usage at bin 
  level using diffSpliceDGE function from edgeR package.  
}

\usage{
  gbDUreport( counts, 
             minGenReads = 10, 
             minBinReads = 5,
             minRds = 0.05, 
             contrast = NULL, 
             ignoreExternal = TRUE, 
             ignoreIo = TRUE, 
             ignoreI = FALSE,
             filterWithContrasted = TRUE,
             verbose = TRUE,
             formula = NULL,
             coef = NULL)
}

\arguments{
  \item{ counts }{ 
    An object of class ASpliCounts
  }
  \item{ minGenReads }{
    Genes with at least an average of \code{minGenReads} reads for any condition 
    are included into the differential expression test.
    Bins from genes with at least an average of \code{minGenReads} reads for all
    conditions are included into the differential bin usage test.
    Default value is 10 reads.
  }
  \item{ minBinReads }{ 
    Bins with at least an average of \code{minGenReads} reads for any condition 
    are included into the differential bin usage test.
    Default value is 5 reads.
  }
  \item{ minRds }{
    Genes with at least an average of \code{minRds} read density for any 
    condition are included into the differential expression test.
    Bins from genes with at least an average of \code{minRds} read density for 
    all conditions are included into the differential bin usage test.
    Bins with at least an average of \code{minRds} read density for any 
    condition are included into the differential bin usage test.
    Default value is 0.05.
    }
  \item{ignoreExternal}{
    Ignore Exon Bins at the beginning or end of the transcript. 
    Default value is TRUE. }
  \item{ignoreIo}{ Ignore original introns.  Default TRUE }
  \item{ignoreI}{Ignore intron bins, test is performed only for exons. 
    Default FALSE }
  \item{contrast}{ 
    Either a formula or a contrast can be tested.
    If contrast is used, it defines the comparison between conditions to be tested.
    \code{contrast} should be a vector with length equal to the number of 
    experimental conditions defined by \code{targets}. The values of this vector
    are the coefficients that will be used to weight each condition, the order 
    of the values corresponds to the order given by \code{getConditions} 
    function. When \code{contrast} is NULL, defaults to a vector containing -1,
    as the first value, 1 as the second an zero for all the remaining values, 
    this corresponds to a pair comparison where the first condition is assumed 
    to be a control and the second condition is the treatment, all other 
    conditions are ignored. 
    Default = NULL}
  \item{filterWithContrasted}{
    A logical value specifying if bins, genes and junction will be filtered by
    read quantity and read density using data from those conditions that will 
    be used in the comparison, i.e. those which coefficients in contrast 
    argument are different from zero. The default value is TRUE, it is strongly
    recommended to do not change this value. }
  \item{verbose}{
    A logical value that indicates that detailed information about each step in 
    the analysis will be presented to the user. }
  \item{formula}{
    Either a formula or a contrast can be tested.
    If formula is used, complex tests can be run.
    \code{formula} should be a formula specifying which experimental conditions defined by \code{targets} to test.
    If coef is specified, then that coefficient will be tested. If not, it defaults to the last term in the formula.
    }
  \item{coef}{
    For formula only. The coefficient to be tested. If null the test defaults to the last term in the formula}    
}

\value{ 
  An ASpliDU object with results at \code{genes}, \code{bins} level.
  \item{genesDE}{
  symbol: gene symbol
  locus_overlap: genes overlaping the same locus
  gene_coordinates: gene coordinates
  start: gene start
  end: gene end
  length: gene length
  effective_length: gene effective length
  logFC: gene log2 fold change between conditions
  pvalue: p-value 
  gen.fdr: fdr corrected p-value for multiple testing
  }
  \item{binsDU}{
  feature: bin type
  event: type of event asigned by ASpli when bining.
  locus: gene locus
  locus_overlap: genes overlaping the same locus
  symbol: gene symbol
  gene_coordinates: gene coordinates
  start: bin start
  end: bin end
  length: bin length
  logFC: bin log2 fold change between conditions
  pvalue: p-value 
  bin.fdr: fdr corrected p-value for multiple testing
  }
  

}

\author{   Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz }

\examples{
  
  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  library(GenomicFeatures)
  genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
                                 package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )

  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:6)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'C','C','C','D','D','D'),
               subject = c(0, 1, 2, 0, 1, 2))
  
  # Read counts from bam files
   gbcounts  <- gbCounts( features = features, 
                        targets = targets, 
                        minReadLength = 100, 
                        maxISize = 50000,
                        libType="SE", 
                        strandMode=0)

# Test for factor1
  # Test for factor1 controlling for paired subject
  gbPaired   <- gbDUreport(gbcounts, formula = formula(~subject+factor1))
  
  # Show all genes and bins ordered by FDR
  genesDE(gbPaired)
  binsDU(gbPaired)
  
  # Test for factor1 without controlling for paired subject. 
  # Must change conditions inside gbcounts object to accommodate the contrast.
  gbcounts@targets$condition <- targets$factor1
  gbcounts@condition.order   <- c("C", "D") 
  gbContrast <- gbDUreport(gbcounts, contrast = c(1, -1))
  
  # Show all genes and bins ordered by FDR
  genesDE(gbContrast)
  binsDU(gbContrast)
  
  # Export results  
  writeDU( du = gbPaired, output.dir = paste0(tempdir(), "/gbPaired") )
  writeDU( du = gbContrast, output.dir = paste0(tempdir(), "/gbContrast") )
}

\seealso{
  \code{\link{edgeR}}, \code{\link{jDUreport}}
  Accessors: \code{\link{genesDE}}, \code{\link{binsDU}}
  Export: \code{ \link{writeDU} }
}


