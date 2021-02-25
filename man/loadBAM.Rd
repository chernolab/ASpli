\name{loadBAM}

\alias{loadBAM}

\title{ Load BAM files }

\description{
  Load BAM files into R session using a targets specification.
}

\usage{ loadBAM(targets, cores,  libType, strandMode) }

\arguments{
  \item{targets}{A data frame containing sample, bam and experimental factors 
    columns }
  \item{cores}{Number of processors to use}
  \item{libType}{Options are: "SE" or "PE"}
  \item{strandMode}{Options are: 0,1,2. See ?strandMode for more information}
}

\value{
  A list of GAlignments or GAlignmentPairs  . Each element of the list correspond to a samples BAM 
  file.
}

\author{ Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz  }

\examples{

  # Define bam files, sample names and experimental factors for targets.
  #bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
  #                   "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  
  #targets <- data.frame( 
  #             row.names = paste0('Sample_',c(1:6)),
  #             bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
  #             factor1 = c( 'C','C','C','D','D','D') )
               
  # Load reads from bam files
  #bams <- loadBAM( targets )
  
}
