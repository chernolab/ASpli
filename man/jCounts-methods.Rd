\name{jCounts}
\alias{jCounts}
\alias{AsDiscover}

\title{ Report PSI, PJU and PIR using experimental junctions }

\description{ Summarize read overlaps against junctions.
              Report PSI, PJU, PIR and counts for experimental junctions.
              PSI or PIR metrics are calculated for each bin and experimental condition.
              The selection of which metric is used is based on the kind of splicing event 
              associated with each bin. 
}
 
\usage{ 
 jCounts( counts, features, minReadLength, threshold = 5, minAnchor = 10,libType="SE", strandMode=0 )
}
\arguments{
\item{counts}{
An object of class ASpliCounts.
}

\item{features}{ An object of class ASpliFeatures.}


\item{minReadLength}{ Minimum read length of sequenced library. It is used for computing 
    E1I and IE2 read summarization. Make sure this number is smaller than the maximum read length in every bam file, otherwise no junctions will be found.}

\item{threshold}{ Minimum  number of reads supporting junctions. Default=5 }
\item{minAnchor}{ An intronic junction must overlap completely and at least an minAnchor\% into the exon region and the intron region. The regions can be exon1-intron or intron-exon2.}
 \item{libType}{Defines how reads will be treated according their sequencing lybrary type (paired (PE, default) or single end (SE))}
 \item{strandMode}{controls the behavior of the strand getter. It indicates how the strand of a pair should be inferred from the strand of the first and last alignments in the pair. 
  0: strand of the pair is always *.
  1: strand of the pair is strand of its first alignment. This mode should be used when the paired-end data was generated using one of the following stranded protocols: Directional Illumina (Ligation), Standard SOLiD.
 2: strand of the pair is strand of its last alignment. This mode should be used when the paired-end data was generated using one of the following stranded protocols: dUTP, NSR, NNSR, Illumina stranded TruSeq PE protocol. For more information see ?strandMode}
}
\value{ 
An object of class ASpliAS.
Accesors: irPIR, esPSI, altPSI, junctionsPIR, junctionsPJU

\item{irPIR}{
event: Type of event asigned by ASpli when bining.
J1: Semicolon separated list of all the junctions with an end matching the start of the intron.
J2: Semicolon separated list of all the junctions with an end matching the end of the intron.
J3: Semicolon separated list of all the junctions overlaping the intron.
All the columns from J1 to J2 represent the J1 counts in the different samples for each bin. The counts are the sum of all the J1 junctions.
All the columns from J2 to J3 represent the J2 counts in the different samples for each bin. The counts are the sum of all the J2 junctions.
All the columns from J3 to the first condition represent the J3 counts in the different samples for each bin. The counts are the sum of all the J3 junctions.
The last columns are the PIR metrics calculated for each condition. The PIR metric is calculated as:  
  \deqn{ PIR = \frac{J1 + J2}{J1 + J2 + 2*J3}}{PIR = (J1 + J2)/(J1 + J2 + 2*J3) }
  Where the junctions are the sum by condition.
}

\item{altPSI}{ 
event: Type of event asigned by ASpli when bining.
J1(J2): Semicolon separated list of all the junctions with an end matching the end of alt5'SS(alt3'SS).
J3: Semicolon separated list of all the junctions with an end matching the start of alt5'SS or the start of alt3'SS.
All the columns from J1 to J2 represent the J1 counts in the different samples for each bin. The counts are the sum of all the J1 junctions.
All the columns from J2 to J3 represent the J2 counts in the different samples for each bin. The counts are the sum of all the J2 junctions.
All the columns from J3 to the first condition represent the J3 counts in the different samples for each bin. The counts are the sum of all the J3 junctions.
The last columns are the PSI metrics calculated for each condition. The PSI metric is calculated as:  
  \deqn{ PSI = \frac{J12}{J12 + J3}}{PSI = (J12)/(J12 + J3) }  
  Where J12 is J1 if it's an alt 5' event or J2 if it's an alt 3' event and the junctions are the sum by condition.
}

\item{esPSI}{ 
event: Type of event asigned by ASpli when bining
J1: Semicolon separated list of all the junctions with an end on the alternative exon.
J2: Semicolon separated list of all the junctions with an end on the alternative exon.
J3: Semicolon separated list of all the junctions overlaping the alternative exon.
All the columns from J1 to J2 represent the J1 counts in the different samples for each bin. The counts are the sum of all the J1 junctions.
All the columns from J2 to J3 represent the J2 counts in the different samples for each bin. The counts are the sum of all the J2 junctions.
All the columns from J3 to the first condition represent the J3 counts in the different samples for each bin. The counts are the sum of all the J3 junctions.
The PSI metric is calculated as:  
  \deqn{ PSI = \frac{J1 + J2}{J1 + J2 + 2*J3}}{PSI = (J1 + J2)/(J1 + J2 + 2*J3) }  
  Where the junctions are the sum by condition.  
}

\item{junctionsPIR}{
PIR metric for each experimental junction using e1i and ie2 counts. Exclusion 
junction is the junction itself. This output helps to discover new introns as 
well as new retention events.
hitIntron: If the junction matches a bin, the bin is shown here. 
hitIntronEvent: If the junction matches a bin, the type of event asigned by ASpli to this bin.
All the columns from hitIntronEvent up to the first repetition of the samples names in the columns, represent the J1 counts in the different samples for each region. 
From there to the next time the names of the columns repeat themselves, the J2 counts and from there to the first condition, the J3 counts. 
The last columns are the PIR metrics calculated for each condition. The PIR metric is calculated as:  
  \deqn{ PIR = \frac{J1 + J2}{J1 + J2 + 2*J3}}{PIR = (J1 + J2)/(J1 + J2 + 2*J3) }  
  Where the junctions are the sum by condition.   
}

\item{junctionsPJU}{
  Given a junction, it is possible to analyze if it shares start, end or both with another junction. If so, it is because there is alternative splicing. 
  Junction: name of the junction.
  gene: gene it belongs to.
  strand: gene strand.
  multipleHit: if other gene overlaps the gene the junction belongs to.
  symbol: gene symbol.
  gene_coordinates: gene coordinates.
  bin_spanned: semicolon separated list of all the bins spaned by this junction.
  j_within_bin: other junctions in the bins.
  StartHit: all the junctions sharing the start with this junction and $PJU_{J1}=J3/(J1+ J3)$ for each condition, 
  EndHit: all the junctions sharing the end with this junction and $PJU_{J2}=J3/(J2+ J3)$ for each condition. 
  All the columns between j_within_bin and StartHit are the counts for J3 in the different samples for each region. 
  From there to EndHit, the J1 counts and $PJU_{J1}=J3/(J1+ J3)$ for each condition. 
  Then after EndHit, the J2 counts and $PJU_{J2}=J3/(J2+ J3)$.
  Rownames are J3 range. StartHit is J1 range and EndHit is J2 range.
}

}
\author{ Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky and Ariel Chernomoretz }

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
  
  # Access summary and gene and bin counts and display them
  gbcounts
  countsg(gbcounts)
  countsb(gbcounts)
  
  # Access summary and junction counts and display them
  jcounts
  irPIR(jcounts)
  esPSI(jcounts)
  altPSI(jcounts)
  junctionsPIR(jcounts)
  junctionsPJU(jcounts)
  
  
  # Export data
  writeAS( as = jcounts, output.dir = paste0(tempdir(), "/only_as") )

}

\seealso{

Accesors: \code{\link{irPIR}}, \code{\link{altPSI}}, \code{\link{esPSI}}, 
\code{\link{junctionsPIR}}, \code{\link{junctionsPJU} }

Export: \code{ \link{writeAS} }

}
