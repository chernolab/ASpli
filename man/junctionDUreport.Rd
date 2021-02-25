\name{junctionDUreport}

\alias{junctionDUreport}

\title{
  Differential junction usage estimation
}
\description{
  Estimate differential usage at junction level. When targets has only two 
  conditions, and contrast is not set, the estimation of differential expression
  and usage is done with an exact test, otherwise is estimated using a 
  generalized linear model.
}

\usage{
  junctionDUreport( counts, 
                    targets,
                    appendTo = NULL, 
                    minGenReads = 10, 
                    minRds = 0.05,
                    threshold = 5,
                    offset = FALSE, 
                    offsetUseFitGeneX = TRUE,
                    contrast = NULL, 
                    forceGLM = FALSE)
}

\arguments{
  \item{ counts }{ 
    An object of class \code{ASpliCounts}
  }
  \item{ targets }{
    A dataframe containing sample, bam and experimental factor columns.
  }
  \item{ appendTo }{ An object of class \code{ASpliDU} to which append the 
    results of junction differential usage. If \code{appendTo} is \code{NULL} a 
    new \code{ASpliDU} is created
  }  
  \item{ minGenReads }{
    Junctions within genes with at least an average of \code{minGenReads} reads
    for all conditions are included into the differential junction usage test.
    Default value is 10 reads.
  }
  \item{ minRds }{
    Junctions within genes with at least an average of \code{minRds} read 
    density for all conditions are included into the differential bin usage 
    test. Junctions with at least an average of \code{minRds} read density for
    any condition are included into the differential junction usage test.
    Default value is 0.05.
  }
  \item{threshold}{
    Junction with at least \code{threshold} counts are included into the 
    differential usage test. }
  \item{offset}{ Corrects junction counts using an offset matrix derived from 
    gene expression data. Default = FALSE}
  \item{offsetUseFitGeneX}{ Fit a GLM using gene counts to build the offset 
    matrix. This argument is used only when 'offset' argument is set to TRUE.
    The default value is TRUE } 
  \item{contrast}{ Define the comparison between conditions to be tested.
    \code{contrast} should be a vector with length equal to the number of 
    experimental conditions defined by \code{targets}. The values of this vector
    are the coefficients that will be used to weight each condition, the order 
    of the values corresponds to the order given by \code{getConditions} 
    function. When \code{contrast} is NULL, defaults to a vector containing -1,
    as the first value, 1 as the second an zero for all the remaining values, 
    this corresponds to a pair comparison where the first condition is assumed 
    to be a control and the second condition is the treatment, all other 
    conditions are ignored. 
    Default = NULL}
  \item{forceGLM}{ Force the use of a generalized linear model to estimate 
    differential expression and usage. Default = FALSE }
}

\value{ 
  An \code{ASpliDU} object with results of differential usage of junctions }

\author{
  Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz
}

\examples{
	#This function has been deprecated. Please see vignette for new pipeline.
}

\seealso{
  \code{\link{edgeR}}, \code{\link{DUreport}}
  Accessors: \code{\link{junctionsDU}}
  Export: \code{\link{writeDU}}
}


