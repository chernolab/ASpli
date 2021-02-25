aspliTargetsExample <- function( ) {
  
  bamfiles <- system.file( 'extdata', 
      c('A_C_0.bam', 'A_C_1.bam', 'A_C_2.bam',
          'A_D_0.bam', 'A_D_1.bam', 'A_D_2.bam',
          'B_C_0.bam', 'B_C_1.bam', 'B_C_2.bam',
          'B_D_0.bam', 'B_D_1.bam', 'B_D_2.bam' ),
      package = "ASpli") 
  
  targets <- data.frame(
      row.names = c( 'A_C_0', 'A_C_1', 'A_C_2',
          'A_D_0', 'A_D_1', 'A_D_2',
          'B_C_0', 'B_C_1', 'B_C_2',
          'B_D_0', 'B_D_1', 'B_D_2' ),
      bam = bamfiles,
      f1  = rep( c('A','B'), each=6 ),
      f2  = rep( c('C','D'), 2, each=3 ),
      stringsAsFactors = FALSE )
  
  return( targets )
  
}

aspliExampleGTF <- function( ) {
  
  system.file( 'extdata', 'genes.mini.gtf' , package = "ASpli" )
  
}

aspliExampleBamList <- function( ) {
  bamFileNames <- paste( 
      rep( c( 'A', 'B' ), each = 6 ), 
      rep( c( 'C', 'D' ), 2, each=3 ), 
      paste0( rep(0:2,2), '.bam'),
      sep='_' )
  
  system.file( 'extdata', bamFileNames , package = "ASpli")
}

aspliBamsExample <- function( ) {
  targets <- aspliTargetsExample()
  loadBAM( targets )
} 

aspliDUexample1 <- function( ) {
  .loadASpliObject( "ASpliDU", 'du', 
      c('genes', 'bins', 'junctions'),
      c( 'genesDE', 'binsDU', 'junctionsDU' ) )
}

aspliDUexample2 <- function( ) {
  .loadASpliObject( "ASpliDU", 'du', 
      c('genes', 'bins' ),
      c( 'genesDE', 'binsDU' ) )
}

aspliASexample <- function( ) {
  .loadASpliObject( "ASpliAS", 'as', 
      c('altPSI', 'esPSI', 'irPIR', 'join', 'junctionsPIR', 'junctionsPSI'),
      c( 'altPSI', 'esPSI', 'irPIR', 'joint', 'junctionsPIR', 'junctionsPSI') )
}

aspliJunctionDUexample <- function( ) {
  .loadASpliObject( "ASpliDU", 'du', 
      c( 'junctions'),
      c( 'junctionsDU' ) )
}

aspliCountsExample <- function( ) {
  .loadASpliObject( "ASpliCounts", 'counts', 
      c( 'bin.rd' , 'e1i.counts', 'exon.intron.counts', 'gene.counts', 
          'gene.rd', 'ie2.counts', 'junction.counts'),
      c( 'rdsb' , 'countse1i', 'countsb', 'countsg', 'rdsg', 'countsie2',
          'countsj') )
}

aspliFeaturesExample <- function( ) {
  
  # -------------------------------------------------------------------------- #
  # Create dummy gene features
  genesGRL <- GRangesList()
  genesGRL$GENE01 <- GRanges( 'reference', IRanges( start= c(1,301,501), 
          end=c( 300,500,700 ) ), '+')
  genesGRL$GENE02 <- GRanges( 'reference', IRanges( start= c(1001,1301,1601), 
          end=c( 1250,1400,1800 ) ), '+')
  genesGRL$GENE03 <- GRanges( 'reference', IRanges( start= c(2001,2251,2501), 
          end=c( 2250,2350,2800 ) ), '+')
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Create summy bin features
  binsGR <- GRanges( 'reference', IRanges( 
          start = c( 1, 301, 501, 1001, 1251, 1301, 1401, 1601, 2001, 2251, 
              2351, 2501 ),
          end =   c( 300,  500,  700, 1250, 1300, 1400, 1600, 1800, 2250, 2350,
              2500, 2800  ) ),
      strand = '+')
  
  names( binsGR ) <- c( 'GENE01:E001', "GENE01:E002", "GENE01:E003", 
      "GENE02:E001", "GENE02:I001", "GENE02:E002", "GENE02:I003", "GENE02:E003",
      "GENE03:E001", "GENE03:E002", "GENE03:I002", "GENE03:E003" )
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Create dummy junction features
  junctionsGR <- GRanges( 'reference', IRanges( 
              start = c( 300, 1250, 1250, 1400, 2250, 2350 ),
              end =   c( 501, 1301, 1601, 1601, 2501, 2501 ) ),
          strand = '+')
  
  names( junctionsGR ) <- c( "GENE01:J001", "GENE02:J001", "GENE02:J003", 
      "GENE02:J002", "GENE03:J002", "GENE03:J001" )
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Create dummy ASpliFeatures object
  features <- new("ASpliFeatures")
  featuresg( features ) <- genesGRL
  featuresb( features ) <- binsGR
  featuresj( features ) <- junctionsGR
  # -------------------------------------------------------------------------- #
  
  return( features )
  
}

.loadASpliObject <- function ( Class, prefix, slotNames, setters ) {
  obj <- new( Class )
  if (length (slotNames) == length( setters) ) {
    for ( i in 1:length( slotNames ) ) {
      fn <- paste0( prefix, '.',slotNames[i] )
      fn <- system.file("extdata", fn,  package="ASpli")
      
      slot.table <- read.table( fn, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE )
      obj <- do.call( paste0(setters[i],'<-'), list( obj , slot.table ) )
    }
  }
  return( obj )
}

