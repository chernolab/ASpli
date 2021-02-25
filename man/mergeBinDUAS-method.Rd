\name{mergeBinDUAS}

\alias{mergeBinDUAS}

\title{Differential usage of bins and PSI/PIR.}

\description{ 
  This function merges the results of differential usage of bins, from an 
  ASpliDU object, with PSI/PIR and junction information, from an ASpliAS object.
  Also, a delta PSI/PIR value is calculated from a contrast.  
}

\usage{
  mergeBinDUAS( du, 
                as,
                targets, 
                contrast = NULL ) 
}

\arguments{
  
  \item{du}{ An object of class \code{ASpliDU} }
  
  \item{as}{ An object of class \code{ASpliAS} }
  
  \item{targets}{ A data frame containing sample, bam files and experimental 
    factor columns. 
  }
  \item{contrast}{ Define the comparison between conditions to be tested.
    \code{contrast} should be a vector with length equal to the number of 
    experimental conditions defined by \code{targets}. The values of this vector
    are the coefficients that will be used to weight each condition, the order 
    of the values corresponds to the order given by \code{getConditions} 
    function. When \code{contrast} is NULL, defaults to a vector containing -1,
    as the first value, 1 as the second an zero for all the remaining values, 
    this corresponds to a pair comparison where the first condition is assumed 
    to be a control and the second condition is the treatment, all other 
    conditions are ignored. 
    The default value is NULL.
  }
}

\author{ Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz  }


\value{
  A data frame containing feature, event, locus, locus_overlap, symbol, 
  gene coordinates, start of bin, end of bin, bin length, log-Fold-Change value,
  p-value, fdr corrected p-value, J1 inclusion junction, J1 junction counts for each 
  sample, J2 inclusion junction, J2 junction counts for each sample, J3
  exclusion junction, J3 junction counts for each sample, PSI or PIR value for 
  each bin, and delta PSI/PIR.
}

\seealso{
  \code{ \link{ASpliDU} },
  \code{ \link{ASpliAS} }
}

\examples{

  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  #library(GenomicFeatures)
  #genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
  #                               package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  #features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  #bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
  #                   "A_D_0.bam", "A_D_1.bam", "A_D_2.bam",
  #                   "B_C_0.bam", "B_C_1.bam", "B_C_2.bam", 
  #                   "B_D_0.bam", "B_D_1.bam", "B_D_2.bam" )
                     
  #targets <- data.frame( 
  #             row.names = paste0('Sample_',c(1:12)),
  #             bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
  #             factor1 = c( 'A','A','A','A','A','A','B','B','B','B','B','B'),
  #             factor2 = c( 'C','C','C','D','D','D','C','C','C','D','D','D') )
  
  # Load reads from bam files
  #bams <- loadBAM( targets )
  
  # Read counts from bam files
  #counts <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
  #                        maxISize = 50000 )
  
  # Calculate differential usage of genes and bins
  #du <- DUreport.norm( counts, targets , contrast = c(1,-1,-1,1))
  
  # Calculate PSI / PIR for bins and junction.
  #as <- AsDiscover( counts, targets, features, bams, readLength = 100, 
  #                        threshold = 5, cores = 1 )

  #mas <- mergeBinDUAS( du, as, targets, contrast =  c(1,-1,-1,1) )                     
  
}
