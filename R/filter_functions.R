.filter.ASpliDU <- function( aspliDU, what = c( 'genes','bins','junctions'), fdr = 1, 
    logFC = 0, absLocFC = TRUE, logFCgreater = TRUE ) {
  
  whatValidValues <- c( 'genes','bins','junctions')
  
  if ( all( tolower( what ) %in% whatValidValues ) ) {

    for ( slotName in tolower(what) ) {
      
      currentSlot <- slot( aspliDU, slotName )
      
      fdrCol <- switch( match( slotName, whatValidValues ), 'gen.fdr', 'bin.fdr', 
          'fdr' )
      
      if ( nrow( currentSlot ) > 0 ) {

        logFCfun <- if ( logFCgreater ) function(a,b) a >= b else function(a,b) a <= b
        
        abFun <- if ( absLocFC ) abs else identity
        
        passFdrCutoff   <- currentSlot[ , fdrCol ] <= fdr
        passLogFcCutOff <- logFCfun( abFun( currentSlot[,'logFC'] ) ,  logFC )
        
        currentSlot <- currentSlot[ passFdrCutoff & passLogFcCutOff , ]
        
        slot( aspliDU, slotName ) <- currentSlot
                
      } else {
        warning( simpleWarning( paste(slotName, ' in ASpliDU object is empty.')  ) )
      }
      
    }
    
    return( aspliDU )
    
  } else {
    stop( simpleError( "Argument 'what' contains invalid values. Valid values are: c( 'genes','bins','junctions')" ) )        
  }
  
}
