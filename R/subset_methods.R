subsetBams <- function( x, targets, select) {
  nSamples <- sum( select %in% rownames( targets ) )
  nCond    <- sum( select %in% getConditions( targets ))
  selectFilter <- NULL
  targets <- .condenseTargetsConditions( targets )
  
  # subset by condition
  if ( nCond >0 & nSamples == 0 )  { selectFilter <- targets$condition %in% select }
  
  # subset by sample
  if ( nSamples > 0 & nCond == 0 ) { selectFilter <- rownames( targets ) %in% select }
  
  return( x[selectFilter] )
    
}

.getExpFactors <- function( targets ) {
  expFactors <- colnames( targets )[ colnames( targets ) != 'condition' ]
  expFactors[2:length(expFactors)]
}

subsetTargets <- function( targets, select, removeRedundantExpFactors = FALSE ) {
  nSamples <- sum( select %in% rownames( targets ) )
  nCond    <- sum( select %in% getConditions( targets ))
  
  selectFilter <- NULL
  targets2 <- .condenseTargetsConditions( targets )
  
  # subset by condition
  if ( nCond >0 & nSamples == 0 )  { selectFilter <- targets2$condition %in% select }
  
  # subset by sample
  if ( nSamples > 0 & nCond == 0 ) { selectFilter <- rownames( targets ) %in% select }

  # Subset Targets
  result <- targets[selectFilter, ]
  
  # Remove conditions with a unique value if required by user
  if ( removeRedundantExpFactors ) {
    for ( ef in .getExpFactors( result ) ) {
      nValues = length( unique( result[, ef]) )
      if ( nValues  == 1) {
        result[, ef] <- NULL
      } 
    }
  }
  
  return( result )
}

.subset.ASpliCounts <- function(x, targets, select) {
  
  nac <- new(Class="ASpliCounts")
  nSamples <- sum( select %in% rownames( targets ) )
  nCond    <- sum( select %in% getConditions( targets ))
  selectFilter <- NULL
  targets <- .condenseTargetsConditions( targets )
  
  # subset by condition
  if ( nCond >0 & nSamples == 0 )  { selectFilter <- targets$condition %in% select }
  
  # subset by sample
  if ( nSamples > 0 & nCond == 0 ) { selectFilter <- rownames( targets ) %in% select }
  
  
  if ( is.null( selectFilter ) ) {
    stop( "subset ASpliCounts can be done with samples or conditions only" )
  }
  
  countsg( nac ) <- cbind( 
      .extractDataColumns(  countsg( x ), targets ), 
      .extractCountColumns( countsg( x ), targets )[ , selectFilter] )
  
  countsb( nac ) <- cbind( 
      .extractDataColumns(  countsb( x ), targets ), 
      .extractCountColumns( countsb( x ), targets )[ , selectFilter] )
  
  countsj( nac ) <- cbind( 
      .extractDataColumns(  countsj( x ), targets ), 
      .extractCountColumns( countsj( x ), targets )[ , selectFilter] )
  
  countse1i( nac ) <- cbind(
      .extractDataColumns(  countse1i( x ), targets ), 
      .extractCountColumns( countse1i( x ), targets )[ , selectFilter] ) 
  
  countsie2( nac ) <- cbind(
      .extractDataColumns(  countsie2( x ), targets ), 
      .extractCountColumns( countsie2( x ), targets )[ , selectFilter] )
  
  rdsg( nac ) <- cbind(
      .extractDataColumns(  rdsg( x ), targets ), 
      .extractCountColumns( rdsg( x ), targets )[ , selectFilter] )
  
  rdsb( nac ) <- cbind(
      .extractDataColumns(  rdsb( x ), targets ), 
      .extractCountColumns( rdsb( x ), targets )[ , selectFilter] )
  
  return( nac )
}

.subset.ASpliAS <- function( x, targets, select ) {
  nas <- new(Class="ASpliAS")
  nSamples <- sum( select %in% rownames( targets ) )
  nCond    <- sum( select %in% getConditions( targets ))
  targets <- .condenseTargetsConditions( targets )
  sampleNames <- rownames( targets )
  
  acceptedSamples <- NULL
  acceptedConditions <- NULL
  
  # subset by condition
  if ( nCond == length(select) & nSamples == 0 )  { 
    acceptedSamples <- sampleNames[ targets$condition %in% select ]
    acceptedConditions <- select
  }
  
  # subset by sample
  if ( nSamples  == length(select) & nCond == 0 ) { 
    acceptedSamples <- select
    acceptedConditions <- targets[ select , 'condition']
  }
  
  
  if ( is.null( acceptedSamples ) ) {
    stop( "subset ASpliAS can be done with samples or conditions only" )
  }
  
  createMask <- function( cNames, sampleNames, conditionNames, 
      acceptedSamples, acceptedConditions ) {
    
    isSampleName <- cNames %in% sampleNames 
    isConditionName <- cNames %in% conditionNames
    isAcceptedSample <- cNames %in% acceptedSamples
    isAcceptedCondition <- cNames %in% acceptedConditions 
    
    return( ! ( isSampleName | isConditionName ) | 
        isAcceptedSample | isAcceptedCondition )

  }
  
  # irPIR / esPSI / altPSI / join tienen la misma estructura de columnas, 
  # Genera un solo filtro para subsetear estos dataframes
  binColNames <- colnames( irPIR( x ) )

  binsFilter <- createMask( binColNames, sampleNames, getConditions( targets ),
      acceptedSamples, acceptedConditions )
  
  irPIR( nas ) <- irPIR( x )[ , binsFilter ]
  colnames( irPIR( nas ) ) <- binColNames[ binsFilter ]
  altPSI( nas ) <- altPSI( x )[ , binsFilter ]
  colnames( altPSI( nas ) ) <- binColNames[ binsFilter ]
  esPSI( nas ) <- esPSI( x )[ , binsFilter ]
  colnames( esPSI( nas ) ) <- binColNames[ binsFilter ]
  joint( nas ) <- joint( x )[ , binsFilter ]
  colnames( joint( nas ) ) <- binColNames[ binsFilter ]
  
  # Junction PIR filter
  jpirNames <- colnames( junctionsPIR( x ) )

  jpirFilter <- createMask( jpirNames, sampleNames, getConditions( targets ),
      acceptedSamples, acceptedConditions )
  
  junctionsPIR( nas ) <- junctionsPIR( x )[ , jpirFilter ]
  colnames( junctionsPIR( nas ) ) <- jpirNames[ jpirFilter ]
  
  # Junction PIR filter

  jpsiNames <- colnames( junctionsPJU( x ) )
  jpsiNames <- sub( "(.start$)|(.end$)", "", jpsiNames )
  
  jpsiFilter <- createMask( jpsiNames, sampleNames, getConditions( targets ),
      acceptedSamples, acceptedConditions )

  junctionsPJU( nas ) <- junctionsPJU( x )[ , jpsiFilter ] 
  colnames( junctionsPJU( nas ) ) <- jpsiNames[ jpsiFilter ]
  
  return( nas )
  
}