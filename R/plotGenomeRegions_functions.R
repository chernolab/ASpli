# TODO agregar verbosity

.getBamMatrix <- function( conditionMatrix, mergedFiles ) { 
  
  apply( conditionMatrix, c(1,2), 
      function( x ) {
        a <- unlist( mergedFiles[ mergedFiles$cond == x , 'bams' ] )
        if (any( is.na(a) ) ) { return(NA) }
        return(a)
      }
  )
} 

.mergeBamsByCondition <- function ( targets , regions , tempFolder = './tmp', 
    keepOld = FALSE, verbose) {
  
  targets <- cbind( targets[ rep(1:nrow(targets), length( regions )),], 
      data.frame( region = rep( c(1:length(regions)), each = nrow(targets)) ) )
  targets$filt <- file.path( tempFolder, paste0( 'sample.', 1:nrow(targets) ,'.bam') )
  
  destBams <- data.frame( cond = getConditions(targets) )
  destBams$bams <- file.path( tempFolder, 
      paste0( 'cond.', getConditions(targets),'.bam') )
  
  file.exists( tempFolder) || dir.create( tempFolder , recursive = TRUE)
  if ( ! keepOld ) {
    if (verbose) message("Extracting and merging reads from bam files")
    for ( i in 1:nrow( targets ) ) {
      
      indexBam( targets$bam[i] )
      
      filterBam(
          file = targets$bam[i],
          destination = targets$filt[i],
          param = ScanBamParam( which = regions[targets$region[i]] ),
          indexDestination = TRUE )
      if (verbose) message(paste("Extraction completed:", as.character(i), "/", as.character(nrow( targets )) ) )
    }
    
    lapply( getConditions(targets), function( cond ) {
          sourceBams <- as.character( targets$filt[ targets$condition == cond ] )
          
          fn <- mergeBam( 
              files = sourceBams , 
              destination = file.path( tempFolder, paste0( 'cond.', cond,'.bam') ) , 
              indexDestination = TRUE, overwrite = TRUE )
        } )
    if (verbose) message("Merging completed" )
    
  } else {
    if (verbose) message("Re-using previously extracted and merged reads")
  }
  return( destBams )
}


.definePlottingRegions <- function( x, xIsBin, verbose, features, genomeTxDb ) {
  
  binFeatures  <- featuresb( features )
  geneFeatures <- featuresg( features )
#  geneCounts <- countsg( counts )
  
  if ( xIsBin ) {
    if (verbose) message( "Extracting gene regions of selected bins" )
    selectedGenes <- unique( mcols(featuresb( features )[ names( featuresb( features ) ) %in% x ])[,'locus'] )
    #selectedGenes <- unique( countsb( counts )[ rownames(  countsb( counts ) ) %in% x ,'locus'] )
  } else {
    if (verbose) message( "Extracting gene regions of selected genes" )
    selectedGenes <- x
  }
  if ( length( selectedGenes ) > 0 ) {
   
#    geneCoordinates <- ranges( geneFeatures[ selectedGenes ] )
#    
#    geneCoordinates <- geneCounts [ selectedGenes, 'gene_coordinates', drop = FALSE ]
#    
#    geneCoordinates <- as.character( geneCoordinates[ 
#            !duplicated( geneCoordinates), 1 ] )
#    geneCoordinates <- strsplit( geneCoordinates, '[-:]')
#    geneCoordinates <- do.call( rbind , geneCoordinates )
#    geneCoordinates <- data.frame(geneCoordinates, stringsAsFactors = FALSE)
#    
#    geneCoordinates[,2] <- as.integer( geneCoordinates[,2] )
#    geneCoordinates[,3] <- as.integer( geneCoordinates[,3] )
#    
#    regions <- GRanges( geneCoordinates[,1], 
#        IRanges( geneCoordinates[,2] , geneCoordinates[,3]) , strand='+')
#    
#    names( regions ) <- selectedGenes
    
    regions <- genes(genomeTxDb)[ selectedGenes ]
   
    return( regions )
    
  } else {
    if (verbose) message( "No valid gene regions recovered" )
    return( GRanges( NULL, IRanges( NULL, NULL ) ) )
  }
}

.arrangeLayout <- function( targets, layout, colors, plotTitles, verbose ) {
  
  .colorsFromConditionMatrix <- function( conditionMatrix , colors, verbose ) {
    if ( length( colors ) == 1 && tolower( colors ) == 'auto' ) {
      if (verbose) message( "Automatic selection of colors" )
      colors <- colors()[1:(length(conditionMatrix)) + 15]
      matrix( colors, ncol = ncol(conditionMatrix) )
    } else {
      if (verbose) message( "Using colors given by user" )        
      matrix( colors, 
          ncol = ncol( conditionMatrix ),
          nrow = nrow( conditionMatrix ),
          byrow = TRUE)
    }
  }
  
  # Check if auto layout is required
  if ( is.character( layout ) && tolower( layout ) == 'auto' ) {
    
    expFactors <- colnames( targets )[ 
        ! colnames( targets ) %in% c( 'bam', 'condition' ) ]
    
    nExpFactors <- length( expFactors )
    
    # ------------------------------------------------------------------------ #
    # Case 1: just one experimental factor
    # result matrix is 1 x n
    if ( verbose ) message( "Auto arrange for one experimental factor: 1 x n matrix")
    if ( nExpFactors == 1 ) {
      conditionMatrix <- matrix( getConditions(targets), 
          ncol = length( getConditions(targets) ) )
      plotTitles <- conditionMatrix
      colors <- .colorsFromConditionMatrix( conditionMatrix, colors, verbose )
    } 
    # ------------------------------------------------------------------------ #
    
    # ------------------------------------------------------------------------ #
    # Case 2: more than one experimental factor and n * m * ... conditions. 
    # result matrix is n * ( m * ... ), where n is the number of values of the 
    # first factor and m * ... are the number of combination of remaining 
    # factors.
    if ( nExpFactors > 1 ) {
      nComb <- prod( apply( targets[, expFactors], 2 , 
              function( x ) length( unique( x ) ) ) )
      
      if ( length( getConditions( targets ) == nComb ) ) {
        if ( verbose ) message( "Auto arrange for more than one experimental factor: m x n matrix")
        conditionMatrix <- matrix( getConditions( targets ), 
            nrow = length( unique( targets[expFactors][,1]) ),
            byrow = TRUE )
        plotTitles <- conditionMatrix
        colors <- .colorsFromConditionMatrix(conditionMatrix, colors, verbose)
        
      } else {
        if ( verbose ) message( "Auto arrange for more than one experimental factor: 1 x n matrix")
        conditionMatrix <- matrix( getConditions( targets ), 
            ncol = length( getConditions(targets) ) )
        plotTitles <- conditionMatrix     
        colors <- .colorsFromConditionMatrix(conditionMatrix, colors, verbose)
      }
    } 
    # ------------------------------------------------------------------------ #
    
  } else {
    
    # ------------------------------------------------------------------------ #
    # Case 3: No auto arrange required. 
    # Matching matrix dimensions of layout, colors and plotTitles matrixes is 
    # checked.
    
    if ( ! is.matrix( layout ) ) { 
      stop(simpleError("Layout must be a matrix or 'auto'" ))
    }
    if ( ! is.matrix( colors ) ) {
      colors <- .colorsFromConditionMatrix( layout, colors, verbose )
    }
    if ( is.character( plotTitles ) && plotTitles == 'auto' ) {
      plotTitles <- layout
    }

   if ( ( ncol( plotTitles ) != ncol(layout) ) | 
         ( nrow( plotTitles ) != nrow(layout) ) |
         ( ncol( colors )     != ncol(layout) ) | 
         ( nrow( colors )     != nrow(layout) ) ) {
      stop( simpleError(
              "The number of columns and rows of colors and plotTitles must match to layout" ))
    }
    if ( verbose ) message( "Using given layout matrix to arrange plots")
    conditionMatrix <- layout
    # ------------------------------------------------------------------------ #
  }
  result <- list()
  result$conditionMatrix <- conditionMatrix
  result$plotTitles <- plotTitles
  result$colors <- colors
  return( result )
}

.plotGenomicRegions <- function( x, genomeTxDb, features, targets, xIsBin = TRUE, 
    layout = 'auto', colors = 'auto', plotTitles = 'auto', sashimi = FALSE, 
    zoomOnBins = FALSE, deviceOpt = NULL, highLightBin = TRUE, outfolder = NULL, 
    outfileType = 'png', mainFontSize = 24, annotationHeight = 0.2,
    annotationCol = 'black', annotationFill = 'gray', annotationColTitle = 'black',
    preMergedBAMs = NULL, useTransparency = TRUE, tempFolder = 'tmp', 
    avoidReMergeBams= FALSE, verbose = TRUE  ) {
  
  targets <- .condenseTargetsConditions( targets )
  
  # -------------------------------------------------------------------------- #
  # Autoarrange 
  layoutPars <- .arrangeLayout( targets, layout, colors, plotTitles, verbose )
  conditionMatrix <- layoutPars$conditionMatrix
  colors <- layoutPars$colors
  plotTitles <- layoutPars$plotTitles
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Keep existing genes and bins in the data set
  if ( xIsBin ) {
    if (verbose) message( "Selecting Bins" )
    #x <- x[ x %in% rownames( countsb( counts ) ) ]
    x <- x[ x %in% names( featuresb( features ) ) ]
  } else {
    if (verbose) message( "Selecting Genes" )
    #x <- x[ x %in% rownames( countsg( counts ) ) ]
    x <- x[ x %in% names( featuresg( features ) ) ]
  }
  if ( length( x ) == 0 ) {
    stop( simpleError( "Bin and/or genes names are incorrect." ) )
  } 
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Collect and merge bam files
  regions <- .definePlottingRegions( x, xIsBin, verbose, features, genomeTxDb )
  if ( is.null( preMergedBAMs) ) {
    if (verbose) message("Using selected regions to extract and merge reads")
    mergedFiles <- .mergeBamsByCondition( targets, regions, 
        keepOld = avoidReMergeBams, tempFolder, verbose )
  } else {
    if (verbose) message("Using pre merged bams files")
    mergedFiles <- data.frame( 
        cond = rownames( preMergedBAMs ),
        bams = as.character( preMergedBAMs[,1] ), stringsAsFactors = FALSE)
  }
  bamfiles <- .getBamMatrix( conditionMatrix , mergedFiles  )
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Plot
  plotCounter = 0
  for ( xi in x ) {
    plotCounter <- plotCounter +1 
    #currentGene <- if ( xIsBin ) countsb( counts )[ xi, c('locus') ] else xi
    currentGene <- if ( xIsBin ) mcols(featuresb( features )[xi])[,'locus'] else xi
    highLightBin <- highLightBin & xIsBin

    zoomOnBins <- if ( ! zoomOnBins | ! xIsBin ) FALSE else  zoomOnBins 
    
    genLims <- c( start( regions[currentGene] ), end( regions[currentGene] ) )
#    binLims <- if ( xIsBin ) as.integer( countsb( counts )[ 
#                  xi, c( 'start', 'end' ) ] ) else NULL 
    binLims <- if ( xIsBin ) {
          df <- as.data.frame( featuresb( features )[xi] )
          rownames(df) <- df$names
          c( df[1, 'start'], df[1,'end']) 
        } else NULL
    
    if ( ! is.null( outfolder ) ) {
      dir.exists( outfolder ) || dir.create( outfolder )
      outfile <- file.path( outfolder, 
          .makeValidFileName( paste0( xi,'.gr.', outfileType ) ) )
    } else { 
      outfile <- NULL
    }
    if ( verbose ) message( paste( "Plotting", xi, "(",
              as.character(plotCounter),"/", as.character(length(x)) ,")" ) )
    
    .makeGenomeRegionPlot( 
      main = xi, 
      genLims = genLims ,
      binLims = binLims ,
      chromosome = as.character( seqnames( regions[currentGene] ) ),
      outfile = outfile,
      genome = genomeTxDb,
      bamFiles = bamfiles,
      plotNames = plotTitles,
      colors = colors,
      mainFontSize = mainFontSize,
      sashimi = sashimi,
      zoomOnBins = zoomOnBins,
      outfileType = outfileType,
      deviceOpt = deviceOpt,
      highLightBin = highLightBin,
      annotationHeight = annotationHeight,
      annotationCol = annotationCol,
      annotationFill = annotationFill,
      annotationColTitle = annotationColTitle,
      useTransparency = useTransparency )
      
   }
  # -------------------------------------------------------------------------- #
}

.makeGenomeRegionPlot <- function ( 
    main , 
    genLims , 
    binLims , 
    chromosome ,  
    genome , 
    outfile , 
    bamFiles = NULL ,    # a matrix of H x V 
    plotNames = NULL ,   # a matrix of H x V
    colors = NULL ,      # a matrix of H x V  
    mainFontSize = 24,
    sashimi = TRUE , 
    zoomOnBins = FALSE, 
    individualWidth = 1024,
    individualHeight = 800,
    highLightBin = TRUE,
    outfileType = 'png',
    deviceOpt = NULL,
    annotationHeight = 0.2,
    annotationCol = 'black',
    annotationFill = 'gray',
    annotationColTitle = 'black',
    useTransparency = FALSE ) {
  
  # -------------------------------------------------------------------------- #
  # Extrae la cantidad de plots horizontales y verticales
  hplots <- ncol ( bamFiles )
  vplots <- nrow ( bamFiles )
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Zoom on Bins
  if ( zoomOnBins  ) {
    minZoomFactor <- abs( binLims[2] - binLims[1] ) / abs(genLims[2]-genLims[1])
    zoomOnBins <- min( max( zoomOnBins, minZoomFactor ) , 1 )
    
    deltaB <- abs( binLims[2] - binLims[1] )
    deltaG <- deltaB / zoomOnBins 
    
    binMean    <- ( binLims[2] + binLims[1] ) / 2
    genLims[1] <- max( binMean - deltaG / 2, genLims[ 1 ] )
    genLims[2] <- min( binMean + deltaG / 2, genLims[ 2 ] )
    genLims[1] <- genLims[1] - as.integer( abs( genLims[ 2 ] - genLims[ 1 ] ) * 0.1 ) 
    genLims[2] <- genLims[2] + as.integer( abs( genLims[ 2 ] - genLims[ 1 ] ) * 0.1 ) 
    
    genLims <- as.integer( genLims )
  }
  # -------------------------------------------------------------------------- #

  # -------------------------------------------------------------------------- #
  # Functions to define tracks
  alnTrack <- function( rangeFile, name, colors, size  ) {
    AlignmentsTrack( 
        range = rangeFile,
        isPaired = FALSE, 
        name=name,
        col= colors,
        fill.coverage = colors,
        col.coverage = colors,
        background.title = if (useTransparency ) 'transparent' else { 'white' },
        type = c( 'coverage', 'sashimi' )[c( TRUE, sashimi ) ],
        sashimiHeight = 0.5,
        coverageHeight = 0.5,      
        lwd.sashimiMax = 3,
        lwd.sashimiMin = 0.5,        
        showTitle = TRUE,
        lwd = 0,
        col.axis = 'black',
        cex.axis = 0.6,
        cex.title = 0.6,
        col.title = 'black',
        size = size,
        margin= 20
    )
  }

  binHighLight <- function ( tracks ) {
    tracks <- tracks[ ! is.na (tracks )]
    colRGB <- if( useTransparency ) rgb(0.5,0.5,0.5, 0.3) else rgb(0.5,0.5,0.5) 
    fillRGB <- if( useTransparency ) rgb( 0.5,0.5,0.5,0.25) else rgb( 0.8,0.8,0.8 )   
    inBackGround <- ! useTransparency
    
    HighlightTrack(
        trackList = tracks, 
        start = binLims[1], 
        end = binLims[2],
        chromosome = chromosome,
#        col="transparent",
        col=colRGB,
        fill=fillRGB,
        lwd=1,
        frame = TRUE,
        inBackground=inBackGround
    )
  }
  # -------------------------------------------------------------------------- #

  # -------------------------------------------------------------------------- #
  # Define tracks
  alntracks <- lapply(1:hplots, function ( y ) { lapply(1:vplots, function( z ) { NA } ) } )
  for ( h in 1:hplots ) {
    for ( v in 1:vplots) {
      if ( ! is.na( bamFiles[ v, h ] ) ) {
        cTrack <- alnTrack ( 
            rangeFile = bamFiles [ v, h ],
            name      = plotNames[ v, h ],
            colors    = colors   [ v, h ],
            size      = 10000 * ( 1 - annotationHeight )  / sum( ! is.na( bamFiles[,h]))) 
        alntracks[[ h ]][[ v ]] <- cTrack 
      } 
    }
  }
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  geneAnnotationTrack <- GeneRegionTrack(genome, 
      chromosome = chromosome,
      start = genLims[1], 
      end = genLims[2],
      options(ucscChromosomeNames=FALSE),
      name = "", 
      transcriptAnnotation = "symbol",
      size = annotationHeight * 10000,  
      fill = annotationFill ,
      col = annotationCol,
      shape = "arrow",
      col.title = annotationColTitle )
  # -------------------------------------------------------------------------- #


  
  if ( highLightBin ) {
    columnTracks <- lapply( alntracks, function (x) {
          binHighLight ( x )
        } )
  } else {
    columnTracks <- alntracks
  }

  
  # -------------------------------------------------------------------------- #
  # creates the graphic device
  outputIsAValidFile <- ( ! is.null( outfile ) ) && 
      ( file.access( dirname( outfile ) , 2 ) == 0 ) 
  
  if ( outputIsAValidFile ) {
    
    devicesNames <- c( 'png', 'bmp', 'jpeg', 'tiff', 'pdf')
    devices      <- list( png, bmp, jpeg, tiff, pdf)
    
    deviceIndex <- match( tolower( outfileType ), devicesNames )
    
    if ( is.na(deviceIndex) ) { 
      dev.new()
      message( paste('Format',outfileType,'is not recognized.',
              'Select png, bmp, jpeg, tiff or pdf') )
      
    } else {
      device <- devices[[ deviceIndex ]]
      undefDeviceOpt <-  is.null( deviceOpt )
      
      if ( (! is.na( deviceIndex )) && deviceIndex %in%  c(1:4) ) {
        fileNameOpt <- list( filename = outfile)
        defaultDeviceOpt <- list( res = 200, height = individualHeight * vplots, 
            width = individualWidth  * hplots, 
            pointsize = 12, units = 'px')
      }
      
      if ( (! is.na( deviceIndex )) && deviceIndex ==  5 ) {
        fileNameOpt <- list( file = outfile)
        defaultDeviceOpt <- list( height = 8 * vplots, width = 8 * hplots, 
            pointsize = 12)        
      }
      
      deviceOpt <- if ( undefDeviceOpt ) 
            append( defaultDeviceOpt, fileNameOpt ) 
          else 
            append( deviceOpt, fileNameOpt )
      do.call( device, deviceOpt )    
      
    }
    
  } else {
    dev.new()
    if ( ! is.null ( outfile ) ) {
      message( paste( "File:",outfile,"cannot be created." ) )
    }
  }
  # -------------------------------------------------------------------------- #


  grid.newpage()
  
  pushViewport( 
      viewport( layout = grid.layout( 
              nrow = 2, 
              ncol = hplots , 
              heights = unit( c ( 1, 5 ) , "null") ) ) )
  
  grid.text(main, gp=gpar( fontsize = mainFontSize ) , 
      vp = viewport(layout.pos.row = 1, layout.pos.col = 1:hplots) )
  
  first <- TRUE
  for ( i in 1:hplots) {
    pushViewport(viewport(layout.pos.col=i, layout.pos.row=2))
    
    plotTracks( c( columnTracks [[ i ]] ,geneAnnotationTrack ), 
        chromosome        = chromosome, 
        from              = genLims[1], 
        to                = genLims[2], 
        min.height        = 0, 
        minCoverageHeight = 0, 
        min.width         = 2, 
        min.distance      = 5,
        background.title  = 'white',
        margin            = 10,
        innerMargin       = 5,
        add               = TRUE
    ) 
    first <- FALSE
    popViewport(1)
  }
  # -------------------------------------------------------------------------- #
  # Close graphic device 
  if ( outputIsAValidFile ) {
    output <- capture.output( dev.off() )
  } 
  # -------------------------------------------------------------------------- #
}

