\name{rds}
\alias{rds}

\title{
  Read density of gene and bins
}

\description{
  Read density of gene and bins is the quotient between the number of read 
  counts and the length of the feature. The results are appended into an 
  ASpliCounts object that must be given as argument. The explicit calculation 
  of read densities is usually not required because is automatically performed
  by \code{readCounts} function.
}

\usage{
  rds( counts, targets )
}

\arguments{
  \item{counts}{ An ASpliCounts object }
  \item{targets}{ A data frame containing sample, bam and experimental factors 
    columns}
}

\value{
  An ASpliCounts object containing read densities of genes and bins.
}

\examples{
    # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  #library(GenomicFeatures)
  #genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
  #                               package="ASpli") )
  #
  ## Create an ASpliFeatures object from TxDb
  #features <- binGenome( genomeTxDb )
#  
#  # Define bam files, sample names and experimental factors for targets.
#  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
#                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
#  targets <- data.frame( 
#               row.names = paste0('Sample_',c(1:6)),
#               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
#               factor1 = c( 'C','C','C','D','D','D') )
#  
#  # Load reads from bam files
#  bams <- loadBAM( targets )
#  
#  # Read counts from bam files
#  counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
#                          maxISize = 50000 )
#                          
#  # Calculates read densities
#  counts <- rds( counts, targets )
}

\author{ Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz  }
