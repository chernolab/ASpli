\name{getConditions}

\alias{getConditions}

\title{
  Retrieve condition names from a targets data frame.
}

\description{ Targets data frame contains experimental factors values for each 
  sample. This function generates a simple name for each unique condition 
  resulting from the combination of all experimental factors. The order of the 
  conditions given by \code{getConditions} is the same that in \code{contrast} 
  argument of \code{DUreport.norm} function. }

\usage{
  getConditions( targets )
}

\arguments{
  
  \item{ targets }{ A dataframe containing sample, bam and experimental factors 
    columns }

}

\value{
  A character vector with the names of the conditions derived from experimental
  factors.
}

\author{ Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz  }

\seealso{
  \code{\link{DUreport.norm}},
  \code{\link{DUreport.offset}} 
}

\examples{
  
  # Define bam files, sample names and experimental factors for targets.
  #bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
  #                   "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  #targets <- data.frame( 
  #             row.names = paste0('Sample_',c(1:6)),
  #             bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
  #             factor1 = c( 'C','C','C','D','D','D') )
  
  # Load reads from bam files.
  # Return value is c('C', 'D') in this example. 
  #conditions <- getConditions(targets)


  # Define bam files, sample names and experimental factors for targets.
  #bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", "A_C_3.bam", 
  #                   "A_D_0.bam", "A_D_1.bam", "A_D_2.bam", "A_D_3.bam"  )
  #targets <- data.frame( 
  #             row.names = paste0('Sample_',c(1:8)),
  #             bam = file.path( 'extdata', bamFileNames, package="ASpli" ),
  #             factor1 = c( 'C','C','C','C','D','D','D','D'),
  #             factor2 = c( 'E','E','F','F','E','E','F','F')  )
  
  # Load reads from bam files.
  # Return value is c("C_E", "C_F", "D_E", "D_F") in this example. 
  #conditions <- getConditions(targets)

}

