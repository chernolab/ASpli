\name{filterSignals}

\alias{filterSignals}

\title{
  Filter signals
}
\description{
  Filter signals
}

\usage{
    filterSignals( sr,
                   bin.FC = 3,
                   bin.fdr = 0.05,
                   nonunif = 1,
                   bin.inclussion = 0.1,
                   bjs.inclussion = 0.2,
                   bjs.fdr = 0.1,
                   a.inclussion = 0.3,
                   a.fdr = 0.05,
                   l.inclussion = 0.3,
                   l.fdr = 0.05, 
                   bDetectionSummary = FALSE)
}

\arguments{
  \item{sr}{ 
    An object of class \code{ASpliSplicingReport}
  }
  \item{bin.FC}{Description TODO}
  \item{bin.fdr}{Description TODO}
  \item{nonunif}{Description TODO}
  \item{bin.inclussion}{Description TODO}
  \item{bjs.inclussion}{Description TODO}
  \item{bjs.fdr}{Description TODO}
  \item{a.inclussion}{Description TODO}
  \item{a.fdr}{Description TODO}
  \item{l.inclussion}{Description TODO}
  \item{l.fdr}{Description TODO}
  \item{bDetectionSummary}{Description TODO}
}

\value{ TODO }

\author{
  Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz
}

\examples{
  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored.
  
  #as <- new("ASpliAS")
  #targets <- 1
  #a <- junctionDUreportExt(as, targets)
  
}

\seealso{
  \code{\link{ASpliSplicingReport}}
}


