.mergeBinDUAS <- function( du, as, targets, contrast = NULL ) { 
  
  targets <- .condenseTargetsConditions( targets )
  bas <- joint( as )
  bdu <- binsDU( du )
  bas$event <- NULL
  
  bas2 <- as.data.frame( matrix( NA, nrow = nrow( bdu ), ncol = ncol( bas ) ) )
  rownames( bas2 ) <- rownames( bdu )
  colnames( bas2 ) <- colnames( bas )
  
  bas2[ intersect( rownames( bas ), rownames( bdu ) ), ] <- bas[ intersect( rownames( bas ), rownames( bdu ) ),  ]
  bas2[,'J1'] <- levels( bas[,'J1'] )[ bas2[,'J1'] ]
  bas2[,'J2'] <- levels( bas[,'J2'] )[ bas2[,'J2'] ]
  bas2[,'J3'] <- levels( bas[,'J3'] )[ bas2[,'J3'] ]
  
  if ( is.null( contrast )) contrast <- .getDefaultContrasts( targets$condition )

  psir <- t( bas2[ , getConditions( targets ) ] )
  psir <- psir * contrast
  psir <- t( psir[ contrast != 0, ] )
  psir <- data.frame( delta = rowSums( psir ) )
  
  r <- cbind( bdu, bas2[ ,], psir )
  colnames( r ) <- c( colnames( bdu), colnames( bas2), 'delta' )
  return( r )
  
}