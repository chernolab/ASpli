\name{ASpli-package}
\alias{ASpli-package}
\alias{ASpli}
\docType{package}

\title{ Analysis of Alternative Splicing Using RNAseq }

\description{
  ASpli is an integrative and flexible package that facilitates the 
  characterization of genome-wide changes in AS under different experimental 
  conditions. ASpli analyzes the differential usage of introns, exons, and splice 
  junctions using read counts, and estimates the magnitude of changes in AS by 
  calculating differences in the percentage of exon inclusion or intron retention 
  using splice junctions. This integrative approach allows the identification of 
  changes in both annotated and novel AS events.
  
  ASpli allows users to produce self-explanatory intermediate outputs, based on 
  the aim of their analysis. A typical workflow involves parsing the genome 
  annotation into new features called bins, overlapping read alignments against 
  those bins, and inferring differential bin usage based on the number of reads 
  aligning to the bins and junctions.
}

\details{
  \tabular{ll}{
    Package: \tab ASpli\cr
    Type: \tab Package\cr
    Version: \tab 1.5.1\cr
    Date: \tab 2018-02-22\cr
    License: \tab  GPL\cr
    Depends: \tab methods, GenomicRanges, GenomicFeatures, edgeR, methods, 
    BiocGenerics, IRanges, GenomicAlignments, Gviz
    \cr
  }
}


\author{ Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz  }

\references{
\itemize{
  \item Acute effects of light on alternative splicing in light-grown plants. 
    Photochemistry and Photobiology. Mancini, E, Sanchez, S, Romanowsky, A, 
    Yanovsky, MJ. DOI: 10.1111/php.12550
  
  \item GEMIN2 attenuates the effects of temperature on alternative splicing and 
    circadian rhythms in Arabidopsis thaliana. Proceedings of the National Academy
    of Sciences. Schlaen, RG, Mancini, E, Sanchez, SE, Perez-Santangelo, S, 
    Rugnone, ML, Simpson, CG, Brown, JWS, Zhang, X, Chernomoretz, A, Yanovsky, 
    MJ. DOI:10.1073/pnas.1504541112
  
  \item Genome wide comparative analysis of the effects of PRMT5 and PRMT4/CARM1 
    arginine methyltransferases on the Arabidopsis thaliana transcriptome. BMC 
    Genomics. Hernando, E, Sanchez, S, Mancini, E, Yanovsky MJ. 
    DOI:10.1186/s12864-015-1399-2
  
  \item A role for LSM genes in the regulation of circadian rhythms. Proceedings 
    of the National Academy of Sciences. Perez Santangelo, S, Mancini, E, Francey, 
    LJ, Schlaen, RG, Chernomoretz, A, Hogenesch, JB, Yanovsky MJ. 
    DOI: 10.1073/pnas.1409791111 
  
  \item The dengue virus NS5 protein intrudes in the cellular spliceosome and 
    modulates splicing. PLOS Pathogens. De Maio, F., Risso, G., Iglesias, G., 
    Shah, P, Pozzi, B., Gebhard, L., Mammi, L., Mancini, E., Yanovsky, M., 
    Andino, R., Krogan, N., Srebrow, A. and Gamarnik, A. 
    DOI:10.1371/journal.ppat.1005841
  }
}

\examples{ 
library(GenomicFeatures)
gtfFileName <- aspliExampleGTF()
genomeTxDb <- makeTxDbFromGFF( gtfFileName )
features <- binGenome( genomeTxDb )
BAMFiles <- aspliExampleBamList()
targets <- data.frame( 
  row.names = paste0('Sample',c(1:12)),
  bam = BAMFiles,
  f1 = c( 'A','A','A','A','A','A',
          'B','B','B','B','B','B'),
  f2 = c( 'C','C','C','D','D','D',
          'C','C','C','D','D','D'),
  stringsAsFactors = FALSE)
  getConditions(targets)
  mBAMs <- data.frame(bam      = sub("_[02]","",targets$bam[c(1,4,7,10)]),
                      condition= c("A_C","A_D","B_C","B_D"))

gbcounts  <- gbCounts( features = features, 
                           targets = targets, 
                           minReadLength = 100, maxISize = 50000,
                           libType="SE", 
                           strandMode=0)
asd   <- jCounts(counts = gbcounts, 
                     features = features, 
                     minReadLength = 100,
                     libType="SE", 
                     strandMode=0)
                     
   
gb      <- gbDUreport(counts=gbcounts,
contrast = c( 1, -1, -1, 1 ) )
jdur    <- jDUreport(asd, 
                     contrast =  c( 1, -1, -1, 1 ) ,
                     mergedBams = mBAMs)
sr      <- splicingReport(gb, jdur, counts =gbcounts    )
is      <- integrateSignals(sr,asd)
}

\keyword{ RNA-seq alternative splicing analysis using bin coverage and junctions }
