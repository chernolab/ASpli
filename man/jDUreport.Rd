\name{jDUreport}

\alias{jDUreport}

\title{ Differential junction usage estimation }

\description{ 
This function estimates the differential usage of junctions combining different types of evidence

Differential junction usage is estimated using a combination of evidences


  }

\usage{
  jDUreport(asd, 
            minAvgCounts                       = 5, 
            contrast                           = NULL,
            filterWithContrasted               = TRUE,
            runUniformityTest                  = FALSE,
            mergedBams                         = NULL,
            maxPValForUniformityCheck          = 0.2,
            strongFilter                       = TRUE,
            maxConditionsForDispersionEstimate = 24,
            formula                            = NULL,
            coef                               = NULL,
            maxFDRForParticipation             = 0.05,
            useSubset                          = FALSE)
}

\details{
Estimation is made at junction level using diffSpliceDGE function from edgeR package. 
Junctions belonging to the same AS event comprises the event "set". 
Each junction is tested against this "set" in a similar fashion that bins are tested against their gene in diffSpliceDGE. Localec are clusters made of junctions that share an end with at least another junction in the cluster.
}
\arguments{
\item{asd}{ An object of class \code{ASpliAS} with results of PSI and PIR using experimental junctions}
\item{minAvgCounts}{ Minimum average counts for filtering}  
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
If NULL must provide a formula. }
    
\item{filterWithContrasted}{
    A logical value specifying if bins, genes and junction will be filtered by
    read quantity and read density using data from those conditions that will 
    be used in the comparison, i.e. those which coefficients in contrast 
    argument are different from zero. The default value is TRUE, it is strongly
    recommended to do not change this value.  }

\item{runUniformityTest}{ Run uniformity test on Intron Retention. Sometimes Mutually Exclusive Exons (MEX) events can be confused with Intron Retention events. This test compares the standard deviation of the inner intron region (11 bases from both ends) to the mean of both intron ends. Numbers closer to 0 mean the event is more probably an Intron Retention event than an MEX event. The test takes some time to run so it defaults to FALSE. }

\item{mergedBams}{     Path to merged bams for each testing condition. If no merged bams exist (for example, paired samples without replicates), use the same bams as targets.  }

\item{maxPValForUniformityCheck}{    To speed up uniformity test only check junctions with pval < maxPValForUniformityCheck  } 

\item{strongFilter}{  If strongFilter is TRUE, then we remove all events with at least one junction that doesn't pass the filter. }  

\item{maxConditionsForDispersionEstimate}{  In order to reduce resource usage, estimate dispersion for statistics tests with a reduced number of conditions.  }

\item{formula}{    Either a formula or a contrast can be tested.    If formula is used, complex tests can be run. 
\code{formula} should be a formula specifying which experimental conditions defined by \code{targets} to test.
    If coef is specified, then that coefficient will be tested. If not, it defaults to the last term in the formula.
    }

\item{coef}{
    For formula only. The coefficient to be tested. If null the test defaults to the last term in the formula}   
  \item{maxFDRForParticipation}{
   In order to calculate junctionPSI participation, only use significant junctions (ie junctions with FDR < maxFDRForParticipation).}   
  \item{useSubset}{    Experimental. It is strongly recomended to leave the default, FALSE.}

}

\value{ 
  An \code{ASpliJDU} object with results of differential usage at \code{junctions} level.
  \item{localec}{
  size: number of junctions belonging to the cluster.
  cluster.LR: likelihood ratio of cluster differential usage.
  pvalue: pvalue of cluster differential usage.
  FDR: fdr of cluster differential usage.
  range: cluster location.
  participation: participation of the significant junction (FDR < maxFDRForParticipation) presenting maximal participation value inside the cluster
  dParticipation: delta participation of the significant junction (FDR < maxFDRForParticipation) presenting maximal participation value inside the cluster 
  }
  \item{localej}{
  cluster: name of the cluster the junction belongs to
  log.mean: log of mean counts accross all conditions for this junction
  logFC: log fold change of junction accross conditions
  pvalue: pvalue of junction
  FDR: FDR of junction
  annotated: is junction annotated or new
  participation: the maximal participation value observed across contrasted condictions
  dParticipation: delta participation of the maximal participation value observed across contrasted condictions
  From dParticipation to the end, junction counts for all samples
  }  
  \item{anchorc}{
  cluster.LR: likelihood ratio of cluster differential usage.
  pvalue: pvalue of cluster differential usage.
  FDR: fdr of cluster differential usage.
  }   
  \item{anchorj}{
  log.mean: log of mean counts accross all conditions for this junction
  logFC: log fold change of junction accross conditions
  LR: likelihood ratio of junction differential usage.
  pvalue: pvalue of junction
  FDR: FDR of junction
  J1.pvalue: pvalue of J1 junction
  J2.pvalue: pvalue of J2 junction
  NonUniformity: if non uniformity test was performed, numbers closer to zero mean uniformity and closer to one mean non uniformity
  dPIR: junction delta PIR
  annotated: is junction annotated or new
  From annotated to the end, junction counts for all samples
  }   
  \item{jir}{
  J3: J3 junction/s
  logFC: log fold change of junction accross conditions
  log.mean: log of mean counts accross all conditions for this junction
  pvalue: pvalue of junction
  FDR: FDR of junction
  LR: likelihood ratio of junction differential usage.  
  NonUniformity: if non uniformity test was performed, numbers closer to zero mean uniformity and closer to one mean non uniformity
  dPIR: junction delta PIR
  multiplicity: do multiple junctions cross the region
  From multiplicity to the end, junction counts for all samples
  }      
  \item{jes}{
  event: type of event
  J3: J3 junction/s
  logFC: log fold change of junction accross conditions
  log.mean: log of mean counts accross all conditions for this junction
  pvalue: pvalue of junction
  FDR: FDR of junction
  LR: likelihood ratio of junction differential usage.  
  dPSI: junction delta PSI
  multiplicity: do multiple junctions cross the region
  From multiplicity to the end, junction counts for all samples
  }      
  \item{jalt}{
  event: type of event
  J3: J3 junction/s
  logFC: log fold change of junction accross conditions
  log.mean: log of mean counts accross all conditions for this junction
  pvalue: pvalue of junction
  FDR: FDR of junction
  LR: likelihood ratio of junction differential usage.  
  dPSI: junction delta PSI
  multiplicity: do multiple junctions cross the region
  From multiplicity to the end, junction counts for all samples   }      
  \item{contrast}{   Conditions contrasted by ASpli  }
}

\author{   Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz }

\seealso{
  Accesors: \code{\link{localec}}, 
            \code{\link{localej}},
            \code{\link{anchorc}},
            \code{\link{anchorj}},
            \code{\link{jir}},
            \code{\link{jes}},
            \code{\link{jalt}},
            \code{\link{junctionsDU}},
  Export:   \code{\link{writeJDU}},
  \code{\link{writeDU}},  \code{\link{edgeR}}, \code{\link{ASpliAS}}
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
  jPaired <- jDUreport(jcounts, formula = formula(~subject+factor1))
  
  # Show junctions information
  jPaired
  localej(jPaired)
  localec(jPaired)
  anchorj(jPaired)
  anchorc(jPaired)
  jir(jPaired)
  jes(jPaired)
  jalt(jPaired)
  
  # Export results  
  writeJDU( jPaired, output.dir = paste0(tempdir(), "/jPaired") )

}



