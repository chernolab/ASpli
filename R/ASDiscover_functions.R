

# Filter junctions that has in at least one condition a given (threshold) number
# of junction for all samples.
.filterJunctionBySample <- function( df0, targets, threshold) {

  # Simple function to compute the minimum of values of two vectors, by position
  # Example:
  # vecMin( c(1,2,3), c(2,0,0) ) -> c(1,0,0)
  vecMin <- function ( a, b ) {
    at <- a < b
    bt <- a >= b
    av <- rep( NA, length( a ) )
    bv <- rep( NA, length( b ) )
    av[ at ] <- a[ at ]
    av[ !at ] <- 0
    bv[ bt ] <- b[ bt ]
    bv[ !bt ] <- 0
    return( av + bv )
  }
  
  cropped <- .extractCountColumns( df0, targets ) 
  uniqueConditions <- unique( targets$condition )

  # Creates an matrix with Inf ( the neutral element for Min function)
  filter <- matrix( Inf ,
      ncol = length( uniqueConditions ) ,
      nrow = nrow( cropped ) )
  
  # Iterates over conditions and over samples of each condition
  # uses filter matrix to compute partial min operations
  for ( i in 1:length( uniqueConditions ) ) {
    byCond <- cropped[ , targets$condition == uniqueConditions[ i ] , drop = FALSE]
    for ( j in 1:ncol( byCond ) ) {
      filter[ , i ] <- vecMin( filter[ , i ] , byCond[ , j ]  )
    }
  }
  
  # Filter the initial dataframe
  filter <- rowSums( filter > threshold ) > 0
  return( df0[ filter,  ] )

}

# Recover junctions that corresponds exactly to introns.
.e1e2JPIR <- function( intronGRanges, jcounts, targets ) {
  
  # Creates a GRanges from junction names and modifies its to match the 
  # corresponding intron.
  jranges <- sort( .createGRangesExpJunctions( rownames( jcounts ) ) )
  start( jranges ) <- start( jranges ) + 1  
  end( jranges ) <- end( jranges ) - 1
  
  # Looks junctions that overlaps exactly with introns
  overlapped <- findOverlaps( jranges, intronGRanges, type="equal" ) 

  # reorder junction counts by query names  
  queryNames   <- names( jranges )[ queryHits( overlapped ) ]
  subjectNames <- names( intronGRanges )[ subjectHits( overlapped ) ]
  queryOrder   <- match( queryNames, rownames( jcounts ) )
  jc <- .extractCountColumns( jcounts, targets )[ queryOrder, ]
  
  result <- data.frame( J3 = queryNames, 
                        jbin = subjectNames,
                        jc,
                        stringsAsFactors = FALSE)
  return( result )
}

# Get junction PSI for junction overlapping different exon regions.
# The possible overlapping types are: start, end and within
.getJPSIByOverlap <- function ( jranges, exonGRanges, jcounts, targets, overlapType ) {
  
  jranges.edge <- jranges
  
  if ( overlapType == 'start' ) {
    start( jranges.edge ) <- end( jranges.edge )
  } else if ( overlapType == 'end' ){
    end( jranges.edge ) <- start( jranges.edge )
  }
  
  overlapped <- findOverlaps( exonGRanges, jranges.edge, type = overlapType )
  
  overlappedQueryNames <- names( exonGRanges[ queryHits( overlapped ) ] )
  overlappedSubjectNames <- names( jranges.edge[ subjectHits( overlapped ) ] )
  
  df1 <- data.frame( 
    names = overlappedQueryNames ,
    .extractCountColumns( jcounts[ subjectHits( overlapped ), ] ,targets ),
    row.names = NULL )
  
  result <- merge(
    x = data.frame( aggregate( . ~ names, data = df1, sum ) ),
    y = data.frame( aggregate( 
            overlappedSubjectNames ~ overlappedQueryNames, 
            FUN = paste, collapse=";")),
    by.x = 1,
    by.y = 1
  )
  rownames( result ) <- result$names
  result[ , 1 ] <-  NULL
  return( result )
  
}

.createGRangesExpJunctions <- function( jnames ) {
  split <- data.frame( matrix( unlist( strsplit( jnames,  "[.]" ) ),
          byrow = TRUE, ncol = 3 ),   
      row.names = jnames,
      stringsAsFactors = FALSE )  
  colnames(split)  <- c("chrom","start","end")
  split$start <- as.numeric(split$start)  
  split$end <- as.numeric(split$end)
  split$chrom <- as.character(split$chrom)
  jranges <- GRanges( seqnames = split$chrom, 
      ranges = IRanges(
          start= split$start, 
          end = split$end,
          names = jnames)    )
  # jranges <- sort(jranges)
  return( jranges )
}

.junctionsDiscover <- function( df,
                                minReadLength, 
                                targets, 
                                features, 
                                minAnchor,
                                bam) 
  
  {
  cores=1
  
  # This function get the counts of the junctions that overlaps the an 
  # intron/exon region of a bin. The junction must overlap completely and at 
  # least an minAnchor% into the exon region and the intron region.
  # The regions can be exon1-intron or intron-exon2. All junctions are assumed 
  # to correspond to a intron.
  getExonIntronCounts <- function( jranges, 
                                   targets, 
                                   bams, 
                                   minReadLength, 
                                   regionType, 
                                   minAnchor ) {
    
    cores = 1
    exonIntron <- jranges

    minAnchor <- round( minAnchor * minReadLength / 100 )
    
    start( exonIntron ) <- if ( regionType == 'e1i' ) start( jranges ) else end( jranges )
    start( exonIntron ) <- start( exonIntron ) - ( minReadLength - minAnchor ) - 1
    end( exonIntron )   <- if ( regionType == 'e1i' ) start( jranges ) else end( jranges )
    end( exonIntron )   <- end( exonIntron ) + ( minReadLength - minAnchor ) - 1
    
    
    hits <- lapply( bams, function( x ) { 
      x <- GRanges(x)    
    
countOverlaps( exonIntron, x, ignore.strand = TRUE, minoverlap = minReadLength )      
    } )    
    hits <- do.call( cbind.data.frame, hits )
    
    return( hits )
  }
  
  # Get the counts of junctions
  jcounts <- .extractCountColumns( df, targets )

  # Convert junctions names to GRanges (the GRanges object do noy have the same
  # order that the initial data ).
  jranges <- sort( .createGRangesExpJunctions( rownames( df ) ) )

  # Extract the junction that has no gaps
  ungappedBams <- lapply( bam, function( x ) { x[ njunc( x ) == 0 , ] } ) 
  # Get counts for junction overlapping the exon1-intron region and intron-exon2
  # region
  e1i <- getExonIntronCounts( jranges, 
                              targets, 
                              ungappedBams, 
                              minReadLength, 
                              'e1i',  
                              minAnchor)
  
  ie2 <- getExonIntronCounts( jranges, 
                              targets, 
                              ungappedBams, 
                              minReadLength , 
                              'ie2', 
                              minAnchor)  
  # Calculates the PIR value 
  j1 <- .sumByCond( e1i,     targets )
  j2 <- .sumByCond( ie2,     targets )
  j3 <- .sumByCond( jcounts, targets )
  pirValues <- ( j1 + j2 ) / ( j1 + j2 + 2 * j3 )
  
  
  # Search junctions that overlaps exacly with annotated introns
  intronBins <- featuresb( features )
  
  start( intronBins ) <- start( intronBins ) - 1
  end( intronBins ) <- end( intronBins ) + 1

  overlapped <- findOverlaps( jranges, intronBins, ignore.strand = TRUE, 
      type="equal" )
  
  queryNames <- names( jranges[ queryHits( overlapped ) ] )
  subjectNames <- names( intronBins[ subjectHits( overlapped ) ] )

  # sort the names of the junctions-matching bins and the event type of those 
  # bins by the order of the junction data in the GRanges object .
  orderIndex <- match( names( jranges ), queryNames )
  
  hitIntron <- subjectNames[ orderIndex ]
  
  hitIntronEvent <- mcols( intronBins )[ subjectHits( overlapped ) , 
      'event' ][ orderIndex ]

  # Creates string representing genomic coordinates of each junction
  genomicCoordinates <- paste( 
      seqnames( jranges ), 
      start( jranges ),
      end( jranges ),
      sep=".")
  
  # Creates result data.frame
  # construir con do.call( cbind )
  result <- do.call( cbind , 
      list( hitIntron = hitIntron, 
            hitIntronEvent = hitIntronEvent,
            e1i, 
            ie2,
            jcounts,
            pirValues ) )
  rownames( result ) <- genomicCoordinates

  return( result )

}

.junctionsPSI_SUM <- function( df, targets ) {
  
  # This function gets the complete set of all junctions and returns a data 
  # frame containing information of junctions that shares its start or end.
  # The result value is a data.frame that has :
  # 1. the names of the junctions shared by each junction (the 'Hit' column of 
  #    the data frame), 
  # 2. the sums ( by sample ) of the shared junctions 
  # 3. the junction ratio ( counts of given junction / sum of counts of shared junctions ).
  getBoundarySharedData <- function(
      allJunctionCounts,
      junctionsGRanges, 
      boundaryType = "start" ) {
    
    junctionNames <- rownames( allJunctionCounts )
    
    j.boundary <- findOverlaps( jranges, drop.self=TRUE, drop.redundant=FALSE, type=boundaryType)
    
    # Search for the names of the junctions that share the boundary
    queryNames   <- names( jranges[ queryHits( j.boundary ) ] )
    subjectNames <- names( jranges[ subjectHits( j.boundary ) ] )
    sharedNames  <- data.frame( 
        aggregate( subjectNames ~ queryNames, FUN = paste, collapse = ";" ) )
    
    # extract the counts of junctions 'indexed' by the name of the query junctions
    junctionCounts <- data.frame( 
        names = queryNames, 
        .extractCountColumns( allJunctionCounts[ subjectNames, ], targets ),
        row.names = NULL ) 
    

    
    # Sum the rows of the junctions that share the boundary
    sharedRowSum <- data.frame( aggregate( . ~ names, data = junctionCounts, sum ) )
    
    colnames( sharedRowSum )[-1] <- rownames(targets)
    rownames( sharedRowSum ) <- sharedRowSum$names
    sharedRowSum <- .extractCountColumns( sharedRowSum, targets  )
    
    
    # Creates an new empty data.frame to store the counts of summed counts 
    # ordered by the original junction data frame
    sharedRowSumOrdered <- data.frame( 
        row.names = junctionNames, 
        matrix( NA, 
            nrow = nrow( allJunctionCounts ), 
            ncol = ncol( sharedRowSum ) ) )
    colnames( sharedRowSumOrdered ) <- colnames( sharedRowSum )
    
    orderIndex <- match( row.names( sharedRowSum ), junctionNames ) 
    sharedRowSumOrdered[ orderIndex, ] <- sharedRowSum#OK
    sharedRowSumOrdered[ is.na( sharedRowSumOrdered ) ] <- 0
    
    # Reorder shared names
    orderIndex <- match( sharedNames$queryNames, junctionNames )
    sharedNamesOrdered <- as.character( rep("-", nrow( sharedRowSumOrdered ) ) )
    sharedNamesOrdered[ orderIndex ] <- sharedNames$subjectNames

    # Calculates jratio
    j1 <- .sumByCond( .extractCountColumns( allJunctionCounts, targets ), targets  )
    j2 <- .sumByCond( sharedRowSumOrdered , targets  )
    jratioResult <- j1 / ( j1 + j2 )
    colnames( jratioResult ) <- paste( colnames( jratioResult ), boundaryType, sep="." )
    
    #Fix for the columns starting with X (targets starting with a number)
    result <- data.frame( sharedNamesOrdered, sharedRowSumOrdered, jratioResult, check.names=FALSE  )
 
    # Sets the name of the 'Hit' column in the result dataframe 
    colnames( result ) [1] <- if ( boundaryType == 'start' ) 'StartHit' else 'EndHit'
    
    return( result )
    
  }
  
  # Get names of the junctions
  junctionNames <- rownames( df )
  
  # Create a Granges object from the names ( the names has the pattern:
  # 'chromosome.start.end' )
  jranges <- sort( .createGRangesExpJunctions( junctionNames ) )
  
  # Retrieve data from junctions that shares their starts
  sharedStartData <- getBoundarySharedData( df, jranges, "start" )

  # Retrieve data from junctions that shares their ends
  sharedEndData <- getBoundarySharedData( df, jranges, "end" )
  
  # Assign putative events 
#  pAS <- rep( "-", nrow( df ) )
#  pAS[ ( sharedStartData[ , 1 ] != "-" & df$strand == "-" ) | ( sharedEndData[ , 1 ] != "-" & df$strand == "+" ) ] <- "Alt5ss"
#  pAS[ ( sharedStartData[ , 1 ] != "-" & df$strand == "+" ) | ( sharedEndData[ , 1 ] != "-" & df$strand == "-" ) ] <- "Alt3ss"
#  pAS[ sharedStartData[ , 1 ] != "-" & sharedEndData[ , 1 ] != "-"] <- "ES" 
  
   #: construir con do.call( cbind, ... ) para evitar que aparezcan .1 y .2 
  # en los nombre de las muestras y puedan filtrarse
  result <- do.call( cbind, list (
          df,
          sharedStartData,
          sharedEndData
      #  , pAS )
          ) )
  return( result )
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
