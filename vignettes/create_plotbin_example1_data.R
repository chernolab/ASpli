makeExample01Data <- function() {

  # Get example counts and as, just for object skeleton
  counts <- aspliCountsExample()
  as     <- aspliASexample()

  # Define synthetic targets
  targets <- data.frame(
      row.names = paste(rep( c("A","B","C","D"), each=12 ),"t",rep( c(1:12)*4, 4 ),sep = '.') , 
      bam = paste(rep( c("A","B","C","D"), each=12 ),"t",rep( c(1:12)*4, 4 ),sep =  '.'),
      treat = rep( c('A','B','C','D'), each=12 ) ,
      time  = rep( c('T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10',
              'T11', 'T12' ), 4),
      stringsAsFactors = FALSE )
  
  # Make synthetic gene expression data
  a <- countsg( counts )
  b <- matrix( runif( 48 * nrow(a), 0, 30 ), ncol = 48, nrow = nrow(a) )
  b <- b + 2000
  b <- apply( b, c(1,2), as.integer) 
  colnames( b ) <- paste( rep( c("A","B","C","D"), each=12 ), "t", rep( c(1:12)*4, 4 ), sep = '.' )
  rownames( b ) <- rownames( a )
  countsg( counts ) <- cbind( a[1:7], b )
  
  # Make synthetic bin counts data
  a <- countsb( counts )
  b <- matrix( runif( 48 * nrow(a), 15, 26 ), ncol = 48, nrow = nrow(a) )
  profile <- c( 0, 0.2, 0.4,  0.8,1,0.1,  0, 0.3, 0.4,  0.75, 0.95,0.1 )
  modprofile <- c( 1 , 0.7, 0.3, 0)
  profile <- rep ( profile, 4) * rep( modprofile, each = 12) 
  b <- b + matrix( rep( profile * 80, nrow(a) ), ncol = 48, byrow=TRUE)
  b <- apply( b, c(1,2), as.integer)
  
  colnames( b ) <- paste( rep( c("A","B","C","D"), each=12 ), "t", rep( c(1:12)*4, 4 ), sep = '.' )
  rownames( b ) <- rownames( a )

  countsb( counts ) <-cbind( a[1:9], b )
  

  # Make synthetic J1 junction counts
  a <- joint( as )
  
  bj1 <- matrix( runif( 48 * nrow(a), 0, 6 ), ncol = 48, nrow = nrow(a) )
  
  profile <- c( 0, 0.2, 0.4,  0.8,1,0.1,  0, 0.3, 0.4,  0.75, 0.95,0.1 )
  modprofile <- c( 1 , 0.7, 0.3, 0)
  profile <- rep ( profile, 4) * rep( modprofile, each = 12) 
  
  bj1 <- bj1 + matrix( rep( profile * 31, nrow(a) ), ncol = 48, byrow=TRUE)
  bj1[bj1 < 0] <- 0 
  bj1 <- apply( bj1, c(1,2), as.integer)
  
  
  # Make synthetic J2 junction counts
  bj2 <- matrix( runif( 48 * nrow(a), 0,6  ), ncol = 48, nrow = nrow(a) )
  
  profile <- c( 0, 0.2, 0.4,  0.8,1,0.1,  0, 0.3, 0.4,  0.75, 0.95,0.1 )
  modprofile <- c( 1 , 0.7, 0.3, 0)
  profile <- rep ( profile, 4) * rep( modprofile, each = 12) 
  
  bj2 <- bj2 + matrix( rep( profile * 28, nrow(a) ), ncol = 48, byrow=TRUE)
  bj2[bj2 < 0] <- 0 
  bj2 <- apply( bj2, c(1,2), as.integer)
  
  
  # Make synthetic J3 junction counts
  bj3 <- matrix( runif( 48 * nrow(a), 22, 27 ), ncol = 48, nrow = nrow(a) )
  
  profile <- c( 0, 0.2, 0.4,  0.8,1,0.1,  0, 0.3, 0.4,  0.75, 0.95,0.1 )
  modprofile <- c( 1 , 0.7, 0.3, 0)
  profile <- rep ( profile, 4) * rep( modprofile, each = 12) 
  
  bj3 <- bj3 - matrix( rep( profile * 22, nrow(a) ), ncol = 48, byrow=TRUE)
  bj3[bj3 < 0] <- 0 
  bj3 <- apply( bj3, c(1,2), as.integer)
  
  # Make synthetic PSI/PIR values
  psir <- (bj1 + bj2) / (bj1 + bj2 + 2 * bj3) 
  
  # update AS object
  asj <- cbind( a[,1:2], bj1, a[,15], bj2, a[,28], bj3, psir )
  
  colnames( asj ) <- c( 'event', 'J1', 
      paste(rep( c("A","B","C","D"), each=12 ),"t",rep( c(1:12)*4, 4 ),sep = '.'),
      c('J2'),
      paste(rep( c("A","B","C","D"), each=12 ),"t",rep( c(1:12)*4, 4 ),sep = '.'),
      c('J3'),
      paste(rep( c("A","B","C","D"), each=12 ),"t",rep( c(1:12)*4, 4 ),sep = '.'),
      getConditions( targets ) )
  
  joint(as) <- asj
  
  # Define factors and values to plot
  fv <- list( 
      time = c('T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11',
          'T12'  ), treat =  c('A','B','C', 'D') )
  
  # Make the plots
  plotBins( counts, as, 'GENE02:E002', factorsAndValues = fv, targets,
      legendAtSide = FALSE, outfolder = 'grplots', outfileType='pdf', deviceOpt=list(
          width = 7, height = 7 ) )
}