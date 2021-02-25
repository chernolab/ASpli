\name{write}

\alias{writeRds}
\alias{writeCounts}
\alias{writeAS}
\alias{writeDU}
\alias{writeJDU}
\alias{writeSplicingReport}
\alias{writeAll}

\title{ Write results }

\description{ Export tab delimited files in structured output }

\usage{
  writeCounts(counts, output.dir="counts")
  writeRds(counts, output.dir="rds")
  writeDU(du, output.dir="du")
  writeAS(as, output.dir="as")
  writeJDU(jdu, output.dir="jdu")
  writeSplicingReport(sr, output.dir="sr")
  writeAll(counts, du, as, output.dir="output")
}

\arguments{
  \item{counts}{An ASpliCounts object}
  \item{as}{An ASpliAS object}
  \item{du}{An ASpliDU object}
  \item{jdu}{An ASpliJDU object}
  \item{sr}{An ASpliSplicingReport object}
  \item{output.dir}{ Name of output folder (new or existing) }
}

\value{ Tab delimited files are exported in a tidy manner into output folder }

\author{ Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz  }

\seealso{ 
  \code{ \link{jCounts}}, 
  \code{ \link{binGenome}},
  \code{ \link{DUreport.norm}},
  \code{ \link{DUreport.offset}}
}
