\name{gbCounts}

\alias{gbCounts}
\alias{readCounts}

\title{
  Summarize read overlaps
}

\description{ Summarize read overlaps against all feature levels }

\usage{
  gbCounts( features,  targets, minReadLength,
                 maxISize, minAnchor = 10, libType="SE", strandMode=0)
}

\arguments{
  
  \item{features}{ An object of class ASpliFeatures. It is a list of GRanges 
    at gene, bin and junction level
  }
  
  
  \item{targets}{ A dataframe containing sample, bam and experimental factors 
    columns }
  
  \item{minReadLength}{ Minimum read length of sequenced library. It is used for computing 
    E1I and IE2 read summarization. Make sure this number is smaller than the maximum read length in every bam file, otherwise no E1I or IE2 will be found.}
    
  \item{maxISize}{Maximum intron expected size. Junctions longer than this size
    will be dicarded }
  
  
  \item{minAnchor}{ Minimum percentage of read that should be aligned to 
    an exon-intron boundary.}
    
  \item{libType}{Defines how reads will be treated according their sequencing lybrary type (paired (PE, default) or single end (SE))}
  \item{strandMode}{controls the behavior of the strand getter. It indicates how the strand of a pair should be inferred from the strand of the first and last alignments in the pair. 
  0: strand of the pair is always *.
  1: strand of the pair is strand of its first alignment. This mode should be used when the paired-end data was generated using one of the following stranded protocols: Directional Illumina (Ligation), Standard SOLiD.
 2: strand of the pair is strand of its last alignment. This mode should be used when the paired-end data was generated using one of the following stranded protocols: dUTP, NSR, NNSR, Illumina stranded TruSeq PE protocol. For more information see ?strandMode}
}

\value{
  An object of class ASpliCounts. Each slot is a dataframe containing features 
  metadata and read counts. Summarization is reported at gene, bin, junction 
  and intron flanking regions (E1I, IE2).
  \item{countsg}{
  symbol: gene symbol
  locus_overlap: other genes overlaping this locus
  gene_coordinates: gene coordinates
  start: gene start
  end: gene end
  length: gene length
  effective_length: gene effective lengh
  From effective_length to the end, gene counts for all samples
  }
  \item{countsb}{
  feature: bin type
  event: type of event asigned by ASpli when bining.
  locus: gene locus
  locus_overlap: genes overlaping the same locus
  symbol: gene symbol
  gene_coordinates: gene coordinates
  start: bin start
  end: bin end
  length: bin length
  From length to the end, bin counts for all samples
  }  
  \item{countsj}{
  junction: annotated junction matching the current junction.
  gene: gene matching the current junction.
  strand: gene strand for the current junction in case a gene matches with the junction.
  multipleHit: semicolon separated list of junctions matching the current junction.
  symbol: gene symbol.
  gene_coordinates: gene coordinates.
  bin_spanned: semicolon separated list of all the bins spaned by this junction.
  j_within_bin: other junctions in the bins.
  From j_within_bin to the end, junction counts for all samples.
  }   
  \item{countse1i}{
  event: type of event asigned by ASpli when bining.
  locus: gene locus
  locus_overlap: genes overlaping the same locus
  symbol: gene symbol
  gene_coordinates: gene coordinates
  start: bin start
  end: bin end
  length: bin length
  From length to the end, bin counts for all samples
  }   
  \item{countsie2}{
  event: type of event asigned by ASpli when bining.
  locus: gene locus
  locus_overlap: genes overlaping the same locus
  symbol: gene symbol
  gene_coordinates: gene coordinates
  start: bin start
  end: bin end
  length: bin length
  From length to the end, bin counts for all samples
  }      
  \item{rdsg}{
  symbol: gene symbol
  locus_overlap: other genes overlaping this locus
  gene_coordinates: gene coordinates
  start: gene start
  end: gene end
  length: gene length
  effective_length: gene effective lengh
  From effective_length to the end, gene counts/effective_length for all samples
  }
  \item{countsb}{
  feature: bin type
  event: type of event asigned by ASpli when bining.
  locus: gene locus
  locus_overlap: genes overlaping the same locus
  symbol: gene symbol
  gene_coordinates: gene coordinates
  start: bin start
  end: bin end
  length: bin length
  From length to the end, bin counts/length for all samples
  }    
  \item{condition.order}{
    The order in which ASpli is reading the conditions. This is useful for contrast tests, in order to make sure which conditions are being contrasted.
  }
}

\author{   Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz  }

\seealso{
  Accesors: \code{\link{countsg}}, 
            \code{\link{countsb}},
            \code{\link{countsj}},
            \code{\link{countse1i}},
            \code{\link{countsie2}},
            \code{\link{rdsg}},
            \code{\link{rdsb}},
            \code{\link{condition.order}},
  Export:   \code{\link{writeCounts}}
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
  # Access summary and gene and bin counts and display them
  gbcounts
  countsg(gbcounts)
  countsb(gbcounts)

  # Export data
  writeCounts( gbcounts, output.dir = paste0(tempdir(), "/only_counts") )
  
}

