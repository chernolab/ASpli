\name{ASpliCounts-class}

\docType{class}

\alias{ASpliCounts-class}
\alias{jCounts,ASpliCounts-method}
\alias{DUreport,ASpliCounts-method}
\alias{DUreportBinSplice,ASpliCounts-method}
\alias{junctionDUreport,ASpliCounts-method}
\alias{AsDiscover,ASpliCounts-method}
\alias{countsb,ASpliCounts-method}
\alias{countse1i,ASpliCounts-method}
\alias{countsg,ASpliCounts-method}
\alias{countsie2,ASpliCounts-method}
\alias{countsj,ASpliCounts-method}
\alias{targets,ASpliCounts-method}
\alias{countsb<-,ASpliCounts,data.frame-method}
\alias{countse1i<-,ASpliCounts,data.frame-method}
\alias{countsg<-,ASpliCounts,data.frame-method}
\alias{countsie2<-,ASpliCounts,data.frame-method}
\alias{countsj<-,ASpliCounts,data.frame-method}
\alias{targets<-,ASpliCounts,data.frame-method}
\alias{gbDUreport,ASpliCounts-method}
\alias{DUreport.norm,ASpliCounts-method}
\alias{DUreport.offset,ASpliCounts-method}
\alias{plotBins,ASpliCounts-method}
\alias{rdsb,ASpliCounts-method}
\alias{rdsg,ASpliCounts-method}
\alias{rdsb<-,ASpliCounts,data.frame-method}
\alias{rdsg<-,ASpliCounts,data.frame-method}
\alias{condition.order,ASpliCounts-method}
\alias{rds,ASpliCounts-method}
\alias{subset,ASpliCounts-method}
\alias{writeAll,ASpliCounts-method}
\alias{writeCounts,ASpliCounts-method}
\alias{writeRds,ASpliCounts-method}


\title{Class \code{"ASpliCounts"}}

\description{ Contains results of read overlaps against all feature levels 
  summarization }

\section{ Slots }{
  \describe{
    \item{\code{gene.counts}:}{Object of class \code{"data.frame"}  }
    \item{\code{exon.intron.counts}:}{Object of class \code{"data.frame"}  }
    \item{\code{junction.counts}:}{Object of class \code{"data.frame"} }
    \item{\code{e1i.counts}:}{Object of class \code{"data.frame"} }
    \item{\code{ie2.counts}:}{Object of class \code{"data.frame"}  }
    \item{\code{gene.rd}:}{Object of class \code{"data.frame"} }
    \item{\code{bin.rd}:}{Object of class \code{"data.frame"} }
    \item{\code{condition.order}:}{Object of class \code{"character"} }
    \item{\code{targets}:}{Object of class \code{"data.frame"} }
    \item{\code{.ASpliVersion}:}{
      ASpli version when this object was created. It should not be modified by the user.
    }     
  }
}

\section{Methods}{
  \describe{
    \item{AsDiscover}{ psi and pir metrics}
    \item{countsb}{bin counts accesor}
    \item{countse1i}{ e1i counts accesor}
    \item{countsg}{ gene counts accesor}
    \item{countsie2}{ie2 counts accesor}
    \item{countsj}{junction counts accesor}
    \item{DUreport_DEXSeq}{ differential expression and usage estimation 
      using DEXSeq }
    \item{DUreport}{differential expression and usage estimation using DEXSeq}
    \item{rdsb}{bin read densities accesor}
    \item{rdsg}{gen read densities acceesor}
    \item{rds}{compute read densities on genes and bins}
    \item{writeCounts}{Export count tables}
    \item{writeRds}{Export read density tables}
	 }
}

\author{ Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz}
