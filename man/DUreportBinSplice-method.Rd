\name{DUreportBinSplice}

\alias{DUreportBinSplice}

\title{
  Differential gene expression and differential bin usage estimation
}
\description{
  Estimate differential expression at gene level and differential usage at bin 
  level using diffSpliceDGE function from edgeR package. This is an alternative
  approach to DUreport. The results at gene level are the same as the results
  from DUreport. The results at bin level are slightly different.    
}

\usage{
  DUreportBinSplice( counts, 
                     targets, 
                     minGenReads = 10, 
                     minBinReads = 5,
                     minRds = 0.05, 
                     contrast = NULL, 
                     forceGLM = FALSE, 
                     ignoreExternal = TRUE, 
                     ignoreIo = TRUE, 
                     ignoreI = FALSE,
                     filterWithContrasted = FALSE,
                     verbose = TRUE )
}

\arguments{
  \item{ counts }{ 
    An object of class ASpliCounts
  }
  \item{ targets }{
    A dataframe containing sample, bam and experimental factor columns.
  }
  \item{ minGenReads }{
    Genes with at least an average of \code{minGenReads} reads for any condition 
    are included into the differential expression test.
    Bins from genes with at least an average of \code{minGenReads} reads for all
    conditions are included into the differential bin usage test.
    Default value is 10 reads.
  }
  \item{ minBinReads }{ 
    Bins with at least an average of \code{minGenReads} reads for any condition 
    are included into the differential bin usage test.
    Default value is 5 reads.
  }
  \item{ minRds }{
    Genes with at least an average of \code{minRds} read density for any 
    condition are included into the differential expression test.
    Bins from genes with at least an average of \code{minRds} read density for 
    all conditions are included into the differential bin usage test.
    Bins with at least an average of \code{minRds} read density for any 
    condition are included into the differential bin usage test.
    Default value is 0.05.
    }
  \item{ignoreExternal}{
    Ignore Exon Bins at the beginning or end of the transcript. 
    Default value is TRUE. }
  \item{ignoreIo}{ Ignore original introns.  Default TRUE }
  \item{ignoreI}{Ignore intron bins, test is performed only for exons. 
    Default FALSE }
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
    differential expression. It is not used to differential usage of bins. 
    Default = FALSE }
  \item{filterWithContrasted}{
    A logical value specifying if bins, genes and junction will be filtered by
    read quantity and read density using data from those conditions that will 
    be used in the comparison, i.e. those which coefficients in contrast 
    argument are different from zero. The default value is FALSE, it is strongly
    recommended to do not change this value. }
  \item{verbose}{
    A logical value that indicates that detailed information about each step in 
    the analysis will be presented to the user. }
}

\value{ 
  An ASpliDU object with results at \code{genes}, \code{bins} level.}

\author{
  Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz
}

\examples{
	#This function has been deprecated. Please see vignette for new pipeline.
}

\seealso{
  \code{\link{edgeR}}, \code{\link{junctionDUreport}}
  Accessors: \code{\link{genesDE}}, \code{\link{binsDU}}
  Export: \code{ \link{writeDU} }
}


