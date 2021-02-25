# TODO: Add comment
# 
# Author: javier
###############################################################################


# ------------------------------------------------ #
# quantifier | grouping | type   | whom  | filter  #
# ------------------------------------------------ #
# min        | set      | count  | gene  | all     #
# avg        | cond.    | rd     | bin   | any     #
#            |          |        | junc. |         #
# -------------------------------------------------#

#bamFileNames <- c( 
#    "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
#    "A_D_0.bam", "A_D_1.bam", "A_D_2.bam",
#    "B_C_0.bam", "B_C_1.bam", "B_C_2.bam", 
#    "B_D_0.bam", "B_D_1.bam", "B_D_2.bam" )
#
#targets <- data.frame( 
#    row.names = sub( ".bam","",bamFileNames ),
#    bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
#    factor1 = c( 'C','C','C','D','D','D','C','C','C','D','D','D'),
#    factor2 = c( 'A','A','A','A','A','A','B','B','B','B','B','B'),
#    stringsAsFactors = FALSE )


.vecMin <- function ( a, b ) {
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

.rowMin <- function ( x ) {
  result <- x[,1]
  for ( i in 1:ncol(x) ) {
    result <- .vecMin(result, x[,i])
  }
  result
}

.extractValuesToFilter <- function ( counts, type, what, targets ) {
  
  typeAcceptedValues <- c('count','rd')
  whatAcceptedValues <- c('gene','bin','junction')
  
  if ( ! type %in% typeAcceptedValues || ! what %in% whatAcceptedValues ) {
    stop( simpleError( "Type or what argument has an unaccepted value.\nAccepted values are:\n\ttype: 'count','rd'\n\twhat: 'gene','bin','junction'") )
  }
  
  accFunction <- if ( type == 'count' && what == 'gene') { countsg } else 
                 if ( type == 'count' && what == 'bin' ) { countsb } else 
                 if ( type == 'count' && what == 'junction' ) { countsj } else
                 if ( type == 'rd' && what == 'gene' ) { rdsg } else
                 if ( type == 'rd' && what == 'bin' ) { rdsb } else
                 if ( type == 'rd' && what == 'junction' ) {
                   stop( simpleError( "junction read density is not defined"))
                 } else stop( "You should not be here" )

  targets <- .condenseTargetsConditions( targets )
  .extractCountColumns( accFunction(counts) , targets )  
  
}

.evaluateByGroup <- function( values, quantifier, grouping, targets ) {
  
  qFunc <- if ( quantifier == 'min' ) { .rowMin } else 
           if ( quantifier == 'avg' ) { rowMeans } else {
             stop( simpleError( "Quantifier not recognized. Accepted values are: 'min', 'avg'"))
           } 
  
  targets <- .condenseTargetsConditions(targets)
  
  colGroups <- if ( grouping == 'set' ) { list( rep( TRUE, ncol(values)) ) } else
               if ( grouping == 'cond') {
                 lapply( unique( targets$condition ), function( x )  targets$condition %in% x )
               } else
               { stop(simpleError( "Set argument has an unaccepted value.\nAccepted values are: 'set', 'cond'")) }
  
  result <- lapply( colGroups, function( g ) { qFunc( values[ ,g ] ) } )
  do.call( cbind, result )

}
.filterByGroup <- function( values, filter, cutoff ) {
  
  n <- if ( filter == 'any') { 1 } else 
       if ( filter == 'all') { ncol(values) } else { 
        stop( simpleError( "Filter argument has an unaccepted value. Accepted values are: 'any', 'all'" ) ) }
  values <- values >= cutoff
  
  rowSums( values ) >= n
  
}

#a <- .extractValuesToFilter( counts , type = 'count', 'gene', targets )
#b <- .evaluateByGrouping( a, 'min', 'set', targets ) 
#c1 <- .filterByGroup( b, 'any', 310 )

.filterReadCounts <- function( counts, targets, cutoffValue, type=NULL, 
    quantifier=NULL, what= NULL, grouping = NULL, filter = NULL ) {

  typeAcceptedValues <- c('count','rd')
  whatAcceptedValues <- c('gene','bin','junction')
  quantifierAcceptedValues <- c('min', 'avg')
  groupingAcceptedValues <- c( 'cond', 'set' )
  filteAcceptedValues <- c( 'all', 'any' )  
  
  what <- tolower( what )
  quantifier <- tolower( quantifier )
  grouping <- tolower( grouping )
  type <- tolower( type )
  filter <- tolower(filter)
  
  targets <-.condenseTargetsConditions( targets )
  values <- .extractValuesToFilter( counts , type , what, targets )
  values <- .evaluateByGroup( values, quantifier,grouping, targets ) 
  filterMask <- .filterByGroup( values, filter , cutoffValue)
  
  if ( what == 'gene') {
    countsg( counts ) <- countsg( counts )[ filterMask , ]
    rdsg( counts )    <- rdsg( counts )[ filterMask , ]
  } else 
  if ( what == 'bin') {
    countsb( counts ) <- countsb( counts )[ filterMask , ]
    rdsb( counts )    <- rdsb( counts )[ filterMask , ]
  } else 
  if ( what == 'junction') {
    countsj( counts ) <- countsj( counts )[ filterMask , ]
  }
  return( counts )
}

#a <- .filterReadCounts( counts, targets, 315, 'count', 'min', 'gene', 'cond', 'all' )
