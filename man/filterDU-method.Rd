\name{filterDU}

\alias{filterDU}

\title{Filtering ASpliDU objects}

\description{
  ASpliDU object can be filtered to retain genes, bins or junction according to
  their fdr corrected p-value estimated and log-fold-change.
}

\usage{ 
  filterDU( 
    du , 
    what = c( 'genes','bins','junctions'),
    fdr = 1,
    logFC = 0,
    absLogFC = TRUE,
    logFCgreater = TRUE
  )
}

\arguments{
  \item{du}{ 
    An ASpliDU object }
  \item{what}{ 
    A character vector that specifies the kind of features that will be 
    filtered. Accepted values are 'genes', 'bins', 'junctions'. Multiple values
    can be passed at the same time. The default value is c( 'genes','bins', 
    'junctions')}
  \item{fdr}{
    A double value representing the maximum accepted value of fdr corrected 
    p-value to pass the filter. The default value is 1, the neutral value for 
    fdr filtering operation.  
  }
  \item{logFC}{
    A double value representing the cut-off for accepted values of 
    log-fold-change  to pass the filter. The default value is 0, the neutral 
    value for logFC filtering operation if logFCgreater and absLocFC arguments 
    are both TRUE. 
  }
  \item{absLogFC}{ 
    A logical value that specifies that the absolute value of log-fold-change 
    will be used in the filter operation. The default value is TRUE.}
  \item{logFCgreater}{ 
    A logical value that specifies that the log-fold-change value ( or 
    abs(log-fold-change) if absLogFC argument is TRUE) of features 
    must be greater than the cut-off value to pass the filter.
    The default value is TRUE.}

}

\value{
  A new ASpliDU object with the results of the filtering operations. 
  The elements of features that were not specified to be filtered are kept from 
  the input ASpliDU object.
}

\author{ Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz  }

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
  #                   "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  #targets <- data.frame( 
  #             row.names = paste0('Sample_',c(1:6)),
  #             bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
  #             factor1 = c( 'C','C','C','D','D','D') )
  
  # Load reads from bam files
  #bams <- loadBAM( targets )
  
  # Read counts from bam files
  #counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
  #                        maxISize = 50000 )
  
  # Calculate differential usage of junctions only 
  #du       <- DUreport.norm( counts, targets )
  
  # Filter by FDR
  #duFiltered1 <- filterDU( du, what=c('genes','bins'), 
  #   fdr = 0.01 )


  # Filter by logFC, only those that were up-regulated 
  #duFiltered2 <- filterDU( du, what=c('genes','bins'), 
  #   logFC = log( 1.5, 2 ),  absLogFC = FALSE )

}

\seealso{
  \code{\link{DUreport.norm}},
  \code{\link{DUreport.offset}},
  \code{\link{jDUreport}},
  \code{\link{gbDUreport}},
}
