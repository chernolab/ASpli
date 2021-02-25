\name{binGenome}

\alias{binGenome}

\title{ Feature coordinates extraction }

\description{ 
  Exons and introns are subdivided into new features called exon and intron bins
  and are then classified into exclusively exonic bins, exclusively intronic 
  bins or alternative splicing (AS) bins .
}

\usage{ binGenome(genome,  geneSymbols, logTo = "ASpli_binFeatures.log" ) }

\arguments{
  \item{genome}{ An object of class transcriptDb (TxDb) }
  \item{geneSymbols}{ A dataframe with symbol (common names) of TxDb genes. 
    If geneSymbols is NULL, gene name will be repeated }
  \item{ logTo }{ Filename where to print features extraction log } }

\details{ 
  Exon and intron coordinates are extracted from gene annotation, only those 
  from multi-exonic genes are saved for further evaluation. In case more than 
  one isoform exist, some exons and introns will overlap. Exons and introns are
  then disjoint into new features called exon and intron bins, and then they are
  classified into exclusively exonic bins, exclusively intronic bind or 
  alternative splicing bins (AS-bins), which are labeled according to which 
  alternative splicing event are assumed to came from: 
  
  \itemize{ 
     \item {ES: exon skipping} 
     \item {IR: intron retention}
     \item {Alt5|3'ss: alternative five/three prime splicing site}
     \item {"*" (ES*, IR*, AltSS*) means this AS bin/region is involved 
       simultaneously in more than one AS event type}
     \item {external:  from the beginning or the end of a transcript}
  }
  
  Subgenic features are labeled as follow (hypothetical GeneAAA):
  
  \itemize{
    \item {GeneAAA:E001: defines first exonic bin}
    \item {GeneAAA:I001: defines first intronic bin}
    \item {GeneAAA:Io001: defines first intron before disjoint into bins}
    \item {GeneAAA:J001: defines first junction}
  }

  Junctions are defined as the last position of five prime exon (donor position)
  and first position of three prime exon (aceptor position). Using TxDb object,
  it is possible to extract annotated/known junctions. This information will be
  useful for the analysis of "experimental" junctions (reads aligned with gaps). 
  Bins and junctions are labelled always in  5' to 3' sense. This notation is
  strand independent. It implies that bin / junction with lower numbering is 
  always at 5'.
}

\value{ An ASpliFeatures object. It is  a list of features using GRanges 
  format. }

\author{Estefania Mancini, Javier Iserte, Marcelo Yanovsky, 
  Ariel Chernomoretz}

\seealso{ 
  \code{\link{featuresg}}, 
  \code{\link{featuresb}},
  \code{\link{featuresj}} 
}


\examples{
  # Create a transcript DB from gff/gtf annotation file.
  library(GenomicFeatures)
  gtfFileName <- aspliExampleGTF()
  genomeTxDb <- makeTxDbFromGFF( gtfFileName )
  
  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )

  # Extract gene, bin and junctions features
  GeneCoord <- featuresg(features)
  BinCoord <- featuresb(features)
  JunctionCoord <- featuresj(features)
  
} 


