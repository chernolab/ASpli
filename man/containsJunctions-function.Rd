\name{ Examine ASpliDU objects }

\alias{containsJunctions}
\alias{containsGenesAndBins}

\title{ Examine ASpliDU objects }

\description{
  AspliDU object may contain results of differential expression of genes, 
  differential usage of bins and junctions, however not everything is 
  calculated at the same or even present. Calculations for genes and bins can 
  be done independently from junctions. Functions \code{containsJunctions} and
  \code{containsGenesAndBins} allow to interrogate an ASpliDU object about the 
  kind of results it contain. 
}

\usage{
  containsJunctions( du )
  containsGenesAndBins( du )
}

\arguments{
  \item{ du }{ 
    An ASpliDU object.
  }
}
\value{ 
  A logical value that indicates that results for genes and bins, or results for
  junctions are available in the object. }

\author{ Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz  }

\examples{
# see ASpli package

}

