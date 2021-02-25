\name{exportIntegratedSignals}

\alias{exportIntegratedSignals}

\title{
  Export integrated signals.
}
\description{
  Export integrated signals in an easy to analyze HTML table.
}

\usage{
exportIntegratedSignals( is, output.dir="is", 
                         sr, counts, features, asd,
                         mergedBams, 
                         jCompletelyIncluded = FALSE, zoomRegion = 1.5, 
                         useLog = FALSE, tcex = 1, ntop = NULL, 
                         openInBrowser = FALSE, 
                         makeGraphs = TRUE, bforce=FALSE
                        )
}

\arguments{
  \item{is}{ 
    An object of class \code{ASpliIntegratedSignals}
  }
  \item{sr}{ 
    An object of class \code{ASpliSplicingReport}
  }
  \item{counts}{ 
    An object of class \code{ASpliCounts}
  }
  \item{features}{ 
    An object of class \code{ASpliFeatures}
  }
  \item{asd}{ 
    An object of class \code{ASpliAS}
  }  
  \item{output.dir}{HTML reports output directory}
  \item{mergedBams}{
    Dataframe with two columns, bams and conditions. Bams are paths to merged bams for each condition to be ploted. 
  }
  \item{jCompletelyIncluded}{If TRUE only plot junctions completely included in plot region. Else plot any overlapping junction in the region}
  \item{zoomRegion}{ Magnify plot region by this factor}
  \item{useLog}{Plot counts log}
  \item{tcex}{Text size}
  \item{ntop}{Only show n top signals}
  \item{openInBrowser}{Open reports in browser when done}
  \item{makeGraphs}{Generate graphs in reports}
  \item{bforce}{Force plot generation even if plot already exists}
}

\value{ 
  Produces html reports
}
\author{   Andres Rabinovich, Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz }

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
               factor1 = c( 'C','C','C','D','D','D'))
  
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
                     

  # Test for factor1
  gbPaired <- gbDUreport(gbcounts, contrast = c(1, -1))
  jPaired  <- jDUreport(jcounts, , contrast = c(1, -1))
  
  # Generate a splicing report merging bins and junctions DU
  sr       <- splicingReport(gbPaired, jPaired, gbcounts)
  is       <- integrateSignals(sr, jcounts)
  
  #Make merged bams dataframe
  mergedBamsFileNames <- c( "A_C.bam", "A_D.bam" )
  mergedBams <- data.frame(bams = system.file( 'extdata', mergedBamsFileNames, package="ASpli" ), 
  condition = c("C", "D"), stringsAsFactors = FALSE)
  
  # Export integrated signals
  exportIntegratedSignals(is, output.dir = paste0(tempdir(), "/is"), sr, gbcounts,
                        features, jcounts, mergedBams, makeGraphs = TRUE, bforce = TRUE )
  
}

\seealso{
  \code{\link{gbDUreport}}, \code{\link{jDUreport}}, \code{\link{ASpliSplicingReport}},  \code{\link{splicingReport}}, \
  code{\link{ASpliIntegratedSignals}}
}


