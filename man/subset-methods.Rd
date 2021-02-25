\name{Subset ASpli objects}

\alias{subsetBams}
\alias{subsetTargets}
\alias{subset}

\title{ Subset ASpli objects }

\description{
  ASpli provides utility functions to easy subset ASpliCounts objects, 
  ASpliAS objects, targets data frame and lists GAlignments generated with 
  loadBAM function. The subset can be done selecting some of the experimental 
  conditions or samples names ( but not both ).
}

\usage{
  subset( x, ... )
  subsetBams( x, targets, select )
  subsetTargets( targets, select, removeRedundantExpFactors )
}

\arguments{
  \item{ x }{ 
    An ASpliCount or ASpliAS object for subset function, or list of GAlignments 
    for subsetrBams function.
  }
  \item{ targets }{
    A dataframe containing sample, bam and experimental factor columns.
  }
  \item{ select }{
    A character vector specifying the conditions or samples to be kept after 
    subset operation. It's assumed that condition names are different from 
    sample names. 
  }
  \item{ removeRedundantExpFactors }{
    When sub-setting the targets data frame, one or more experimental factors 
    can have only one value. If this argument is TRUE those experimental 
    factors are absent in the resulting target data frame. 
  }
  \item{ ... }{ Subsetting ASpliCounts and ASpliAS objects sub subset method 
    requires a targets argument and a select argument with the same 
    specifications that the arguments with the same name in subsetBams and 
    subsetTargets functions }
  
}
\value{ A data frame similar to x ( or targets for subsetTargets) with only the
 containing only the selected elements. }

\author{ Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz  }

\examples{

  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  #library(GenomicFeatures)
  #genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
  #                               package="ASpli") )
  #
  # Create an ASpliFeatures object from TxDb
  #features <- binGenome( genomeTxDb )
  #
  # Define bam files, sample names and experimental factors for targets.
  #bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
  #                   "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  #targets <- data.frame( 
  #             row.names = paste0('Sample_',c(1:6)),
  #             bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
  #             factor1 = c( 'C','C','C','D','D','D') )
  #
  # Load reads from bam files
  #bams <- loadBAM( targets )
  #
  # Read counts from bam files
  #counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
  #                        maxISize = 50000 )
#  
 # # Create ASpliAS object                        
  #as       <- AsDiscover( counts, targets, features, bams, readLength = 100, 
   #                       threshold = 5, cores = 1 )
    #                      
  # Define selection
  #select <- c('Sample_1', 'Sample_2', 'Sample_4', 'Sample_5')                       
  #
  # Subset target 
  #targets2 <- subsetTargets( targets, select )
  
  # Subset bams 
  #bams2 <- subsetBams( bams, targets, select )
  
  # Subset ASpliCounts object 
  #counts2 <- subset( counts, targets, select )
  
  # Subset ASpliAS object 
  #as2 <- subset( as, targets, select )
  
}

