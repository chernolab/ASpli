\name{integrateSignals}

\alias{integrateSignals}

\title{
  Integrate signals
}
\description{
Integrates differential usage signals from different sources using overlaping regions. See vignette for more details
}

\usage{
  integrateSignals(sr = NULL, asd = NULL, bin.FC = 3, bin.fdr = 0.05,
                 nonunif = 1, usenonunif = FALSE, bin.inclussion = 0.2,
                 bjs.inclussion = 10.3, bjs.fdr = 0.01, a.inclussion =
                 0.3, a.fdr = 0.01, l.inclussion = 0.3, l.fdr = 0.01,
                 otherSources = NULL, overlapType = "any")
}

\arguments{
  \item{sr}{ 
    An object of class \code{ASpliSplicingReport}
  }
  \item{asd}{ 
    An object of class \code{ASpliDU}
  }  
  \item{bin.FC}{
    Filter bin signals by fold change. Actually, log2 fold change is return, so default would return only bin signlas with bin.fc > log2(3).
  }
  \item{bin.fdr}{
    Filter bin signals by fdr.
  }
  \item{nonunif}{
    Filter intronic bins with non uniform support (nonunif << 1 is uniform)
  }
  \item{usenonunif}{
    Use non uniformity as filter.
  }
  \item{bin.inclussion}{
    Filter bin signals by junction support with dPIR or dPSI accordingly.
  }
  \item{bjs.inclussion}{
    Filter annotated junction signals by junction inclussion with dPIR or dPSI accordingly.
  }
  \item{bjs.fdr}{
    Filter annotated junction signals by fdr.
  }
  \item{a.inclussion}{
    Filter anchor junction signals by junction inclussion with dPIR.
  }
  \item{a.fdr}{
    Filter anchor junction signals by fdr.
  }
  \item{l.inclussion}{
    Filter locale junction signals by junction inclussion with dPSI.
  }
  \item{l.fdr}{
    Filter locale junction signals by fdr.
  }
  \item{otherSources}{
    If user wants to compare ASpli results with results from other methods, otherSources must be a GenomicRange object with all the regions found with the other methods. It will be integrated with a new column next to signals information.  }
  
  \item{overlapType}{
  Type of regions overlap matching between the different signals. Defaults to "any" and can be any of the following: "any", "start", "end", "within", "equal".
  }    
}

\value{ 
  It returns A \code{ASpliIntegratedSignals} with all overlaping signals present in the region filtered by different parameters.
  }

\author{
  Andres Rabinovich, Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz
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
  sr       <- splicingReport(gbPaired, jPaired, gbcounts)
  is       <- integrateSignals(sr, jcounts)
  
  # Show integrate signals results and filters used
  signals(is)
  filters(is)

}

\seealso{
Accesors: \link{signals}, 
          \link{filters},
Export:   \link{exportIntegratedSignals}
  \code{\link{gbDUreport}}, \code{\link{jDUreport}}, \code{\link{ASpliSplicingReport}},  \code{\link{splicingReport}}, \
  code{\link{ASpliIntegratedSignals}}

}


