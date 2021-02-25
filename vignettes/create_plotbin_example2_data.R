makeExample02Data <- function() {

  # Get example counts and as, just for object skeleton
  counts <- aspliCountsExample()
  as     <- aspliASexample()

  # Define synthetic targets
  targets <- data.frame(
      row.names = paste( rep( c('A','B'),each=6), rep( c('C','D'),2,each=3) , rep( 0:2, 4), sep = '_') , 
      bam = "s" ,
      f1 = rep( c('A','B'),each=6) ,
      f2 = rep( c('C','D'),2,each=3),
      stringsAsFactors = FALSE )
  
  # Define factors and values to plot
  fv <- list( 
          f1 = c( 'A', 'B' ), 
          f2 = c( 'C' ) )
  
  # Make the plots
  plotBins( counts, as, 'GENE02:E002', factorsAndValues = fv, targets,
      legendAtSide = TRUE, outfolder = 'grplots'
      , outfileType='pdf', deviceOpt=list( width = 3, height = 7 ) ,
      innerMargins    = c( 2.1, 5.1, 1.1, 1.1 ), las.x=1
  )
}