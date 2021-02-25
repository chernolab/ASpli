\name{splicingReport}

\alias{splicingReport}

\title{
  Splicing report
}

\description{
This function integrates bin and junction usage in a comprehensive report
}

\usage{
  splicingReport(bdu, jdu, counts)
}

\arguments{
  \item{bdu}{ 
    An object of class \code{ASpliDU}
}
  \item{jdu}{ 
    An object of class \code{ASpliJDU}
  }
  \item{counts}{ 
    An object of class \code{ASpliCounts}
  }
}

\value{ 
  An \code{ASpliSplicingReport} object with junction differential usage report. See vignette for more details
  
}

\author{   Andres Rabinovich, Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz }

\seealso{
  Accesors: \code{\link{binbased}}, 
            \code{\link{localebased}},
            \code{\link{anchorbased}},
  Export:   \code{\link{writeSplicingReport}}
  \code{\link{gbDUreport}}, \code{\link{jDUreport}}, \code{\link{ASpliSplicingReport}},  \code{\link{splicingReport}}
}


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
                           minReadLength = 100, maxISize = 50000,
                           libType="SE", 
                           strandMode=0)
jcounts   <- jCounts(counts = gbcounts, 
                     features = features, 
                     minReadLength = 100,
                     libType="SE", 
                     strandMode=0)
                     

# Test for factor1 controlling for paired subject
gbPaired <- gbDUreport(gbcounts, formula = formula(~subject+factor1))
jPaired  <- jDUreport(jcounts, formula = formula(~subject+factor1))
  
# Generate a splicing report merging bins and junctions DU
sr              <- splicingReport(gbPaired, jPaired, gbcounts)

# Access splicing report elements
sr  
localebased(sr)
anchorbased(sr)
binbased(sr)
}



