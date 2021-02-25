# Make plot for a set of bins
.plotBins <- function( counts, as, bin, factorsAndValues, targets, main ,
    colors ,  panelTitleColors, panelTitleCex, innerMargins , 
    outerMargins, useBarplots, barWidth, barSpacer, las.x, useHCColors,
    legendAtSide, outfolder, outfileType, deviceOpt) {
  
  for ( cBin in bin) {
    print(cBin)
    filename <- if ( ! is.null( outfolder )) {
          file.path( outfolder , .makeValidFileName( paste0(cBin,'.pr.',outfileType ) ) )
        } else {
          outfolder
        }
    
    .plotSingleBin( counts, as, cBin, factorsAndValues, targets, main, colors, 
        panelTitleColors, panelTitleCex, innerMargins, outerMargins, 
        useBarplots, barWidth, barSpacer, las.x, useHCColors, legendAtSide,
        filename, outfileType, deviceOpt )
  }
}

# ---------------------------------------------------------------------------- #
# .getHCColor return a color that should have high contrast with a another given
# color and also with the white background 
.getHCColor <- function ( color ) {
  rgbCol <- col2rgb ( color , alpha = TRUE)
  brightness <- mean( rgbCol [ 1:3,1 ] )
  if ( brightness < 95 ) {
    if ( rgbCol [ 1,1 ] > rgbCol [ 2,1 ] & rgbCol [ 1,1 ] > rgbCol [ 3,1 ] ) {
      return( rgb( 0.3,0.5,1,1 ) )
    }
    if ( rgbCol [ 2,1 ] > rgbCol [ 1,1 ] & rgbCol [ 2,1 ] > rgbCol [ 3,1 ] ) {
      return( rgb( 1, 0.1,0.1,1 ) )
    }
    return( rgb( 0.3,1, 0.5,1 ) )
  } else {
    return( rgb( 0, 0, 0 ,1 ) )
  }
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# This function get a vector data for a bin property, with many values as
# samples. Values are averaged by condition a then returned as a list of
# vectors grouped by the main factor values.
.prepareData <- function( samplesProfileData, gridPanels, targets, factorsAndValues ) {
  
  mainFactorIndex <- length( factorsAndValues )
  targets$samplesNames <- rownames( targets)
  
  data <- lapply ( 1:nrow( gridPanels ) , function ( i ) {
        
        samplesData <- merge( targets, gridPanels[ i, , drop = FALSE ], 
            by = colnames( gridPanels[ i, , drop = FALSE ] ) )
        
        mainFactorValues <- factorsAndValues[[ mainFactorIndex ]]
        
        plotData <- sapply ( mainFactorValues, 
            function( x ) { 
              samples <- samplesData[ 
                  samplesData[ , names( factorsAndValues[ mainFactorIndex ] ) ] == x  , 
                  'samplesNames' ]
              rowMeans ( samplesProfileData[ rownames( targets ) %in% samples ] )
            } )
        
        plotData[ is.na( plotData )] <- 0
        return( plotData )
      } )
  return(data)
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# .makeLegends creates a top-left legend over a single plot.
.makeLegends <- function( main , useHCColors, col, panelTitleColor, panelTitleCex ) {
  legend( 
      x = "topleft", 
      legend = main,
      adj = c( 0 , 0 ),
      yjust = 0.5,
      bty = "n",
      cex = panelTitleCex,
      text.col = if ( useHCColors ) .getHCColor( col ) else { panelTitleColor } )
}
# ---------------------------------------------------------------------------- #


# -------------------------------------------------------------------------- #
# Function to make a single xy plot 
.plotLines <- function ( data, dataCons, gridPanels, nPoints, xTicksLabels, 
    factorsAndValues, mainFactorIndex, minValue, maxValue, main, col, 
    legendAtSide, las.x , panelTitleCex ) {
  
  plot( x = 1:length( dataCons ) , 
      y = dataCons,
      pch = 20,
      type = 'p',
      col = 'gray',
      ylim = c( minValue, maxValue),
      xlab='',
      xaxt = "n",
      ylab = if ( legendAtSide ) main else '',
      yaxt = 'n',
      cex.lab = panelTitleCex )
  
  secFactorTable <- table( gridPanels[ , ncol(gridPanels)] )
  st = 0
  for ( i in 1:nrow( secFactorTable )) {
    from <-  (st + 1)
    to <- st + secFactorTable[i] * nPoints
    at <- c(from:to)
    axis( side = 1, 
        at = at, 
        labels= xTicksLabels[from:to], 
        las=las.x , 
        cex.axis=0.7, tck=0,
        mgp=c(3,0.3,0),
        col= rgb(0.9,0.9,0.9) )
    st <- to
  }
  
  st = 0
  for( i in 1:nrow( gridPanels )) { 
    lines( x = st + 1:length( factorsAndValues[[ mainFactorIndex ]] ) , 
        y = data[[i]],
        col = col)
    st <- st + length( factorsAndValues[[mainFactorIndex]] )
  }
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Function to make a single bar plot 
.plotBars <- function( dataCons, nPoints, barWidth, barSpacer, maxValue, main, 
    useHCColors, panelTitleColor, xTicksLabels, col, las.x, legendAtSide, panelTitleCex) {
  spacers <- c( (1 - barWidth) / barWidth, barSpacer)
  dataCons <- matrix( dataCons, ncol = nPoints, byrow = TRUE ) 
  barplot ( height = t( dataCons ),
      col = col,
      ylim = c( 0, maxValue),
      xlab = '',
      ylab = if ( legendAtSide ) main else '',
      yaxt = 'n',
      width = barWidth,
      space = spacers,
      border = '#DDDDDDFF',
      beside = TRUE,
      las = 2,
      cex.lab = panelTitleCex)
  
  ticks <- matrix( 1, 
      ncol = ncol(dataCons),
      nrow = nrow(dataCons))
  spacerTick <- barWidth * barSpacer - ( 1 - barWidth )
  ticks <- cbind( ticks, spacerTick )
  
  ticks <- as.vector( t( ticks ) )
  
  ac <- 0
  tickSum <- c()
  for ( i in 1 : (length( ticks ) -1) ) {
    tickSum <- append( tickSum , ac + ticks[i])
    ac <- tickSum[i]
  }
  
  tickSum <- tickSum - ( 1 - ( spacerTick + (1-barWidth) ) ) + barWidth / 2
  
  
  tlblmat <- matrix ( xTicksLabels, byrow=TRUE, ncol = ncol(dataCons) ) 
  tlblmat <- cbind( tlblmat, '' )
  xTicksLabels <- as.vector( t(tlblmat) )[1:( length( tickSum ) )]
  
  axis( 
      side = 1, 
      at = tickSum , 
      labels = xTicksLabels, 
      las = las.x , 
      cex.axis = 0.7, 
      tck=0,
      mgp = c(3,0.3,0),
      col = rgb( 0.9,0.9,0.9, 0 ),
      padj = 0.5 )
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Create the labels for factor corresponding to the x axis of plots. 
.makeXticksLabels <- function( justOneFactor, gridPanels, factorsAndValues, 
    mainFactorIndex, data ) {
  
  if ( ! justOneFactor ) {
    
    xTicksLabels <- rep( apply( gridPanels, 1, 
            function( x ) paste0(x,collapse=  ".") ), sapply( data, length  ) )
    xTicksLabels <- apply( gridPanels, 1, function( x ) paste0(x,collapse=  ".") )
    xTicksLabels <- expand.grid( factorsAndValues[[ mainFactorIndex]], xTicksLabels )
    xTicksLabels <- paste0( xTicksLabels[,2] ,".",  xTicksLabels[,1] )
    
  } else {
    xTicksLabels <- factorsAndValues[[mainFactorIndex]]
  }
  return( xTicksLabels )
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Main plotting function
.makePlot <- function ( data, useBarPlot = FALSE, main, col, factorsAndValues, 
    justOneFactor, gridPanels, nPoints, legendAtSide, las.x, panelTitleCex,
    panelTitleColor, useHCColors, barWidth, barSpacer ) {
  
  mainFactorIndex <- length( factorsAndValues )

  # Get the labels of ticks in the x axis.
  xTicksLabels <- .makeXticksLabels( justOneFactor , gridPanels, 
      factorsAndValues, mainFactorIndex, data )
  
  # Get the y axis minimum and maximum values. Used to scale the plots correctly.
  maxValue <- max ( sapply ( data, max ))
  minValue <- min ( sapply ( data, min ))
  
  # Make a unique vector will all data
  dataCons <- Reduce( function( a, b) { a <- append( a, b )}, data )

  # Remove box around plots
  par( bty = 'n')
  
  # Make the main plot
  if ( ! useBarPlot ) {
    .plotLines (data, dataCons , gridPanels , nPoints, xTicksLabels, 
        factorsAndValues, mainFactorIndex, minValue, maxValue, main, 
        col = col, legendAtSide = legendAtSide , las.x, panelTitleCex )
  } else {
    .plotBars(dataCons, nPoints, barWidth, barSpacer, maxValue, main, 
        useHCColors, panelTitleColor, xTicksLabels, col, las.x, legendAtSide, 
        panelTitleCex )
  }
  
  # Add legend inside plot is required
  if ( ! legendAtSide ) {
    .makeLegends( main, useHCColors, col, panelTitleColor, panelTitleCex )
  }
  
  # Draw y axis ticks
  axis( 2, las=1, cex.axis=0.7 )
  
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# .plotbins draw all plots for a given bin
.plotSingleBin <- function( 
    counts, 
    as,
    bin, 
    factorsAndValues, 
    targets,
    main            = NULL,
    colors          = c( '#2F7955', '#79552F', '#465579', '#A04935', '#752020', 
        '#A07C35') ,
    panelTitleColors = '#000000',
    panelTitleCex   = 1,
    innerMargins    = c( 2.1, 3.1, 1.1, 1.1 ),
    outerMargins    = c( 0, 0, 2.4, 0 ), 
    useBarplots     = NULL,
    barWidth        = 0.9,
    barSpacer       = 0.4,
    las.x           = 2,
    useHCColors     = FALSE,
    legendAtSide    = TRUE,
    outfile         = NULL,
    outfileType     = 'png',
    deviceOpt       = NULL
    ) {
  
  # -------------------------------------------------------------------------- #
  # Manipulate factors
  # If just one factor, add another fictional factor with only one value to
  # make computation easier.
  justOneFactor <- FALSE
  if ( length ( factorsAndValues) == 1 ) {
    factorsAndValues <- append( factorsAndValues, list( single = 1 ) )
    colnames( joint(as) ) [ 
        match( unique( targets[,names(factorsAndValues)[[1]]] ), colnames( joint(as) ))
    ] <- paste( unique( targets[,names(factorsAndValues)[[1]]] ), '1',sep='_')
    targets <- cbind( targets, single = 1 )
    justOneFactor <- TRUE
  }
  # Create condition names
  targets <- .condenseTargetsConditions(targets)
  # factorsAndValues is reversed, this is required later in expand.grid.
  factorsAndValues <- factorsAndValues[ rev( names( factorsAndValues ) ) ]
  # Get the index of the first factor before factor list is reversed.
  mainFactorIndex  <- length( factorsAndValues )
  # Get the number of points of main factor
  nPoints          <- length( factorsAndValues[[ mainFactorIndex ]] )
  # -------------------------------------------------------------------------- #

  # -------------------------------------------------------------------------- #
  # Get the combination of all factors other than the main 
  gridPanels <- expand.grid( factorsAndValues[-mainFactorIndex] )
  
  # Subset from all factor combinations those that are present in the targets
  condInTargets <-  unique( apply( targets[ , names( gridPanels ) , drop= FALSE], 1,
          function(x) paste0(x,collapse = "_") ) )
  
  condInGrid <- apply( gridPanels, 1,
      function(x) paste0(x,collapse = "_") )
  
  gridPanels <- gridPanels[ condInGrid %in% condInTargets , , drop=FALSE ]
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Extract data from bin
  binCounts <- .extractCountColumns( countsb( counts )[ bin,] , targets )
  geneLocus <- countsb( counts )[ bin,'locus']
  
  geneCounts <- .extractCountColumns( countsg( counts )[ geneLocus,] , targets )
  
  psir <- joint( as )[ bin, as.character( unique( targets$condition ) ) ]
  psir [ is.na(psir) ] <- 0
  
  J123 <- joint( as )[ bin,  colnames( joint( as ) ) %in%  rownames(targets)]
  ac <- 0
  J1   <- J123[ ( ac + 1 ) : ( ac + nrow(targets) ) ]
  ac <- ac + nrow(targets)
  J2   <- J123[ ( ac + 1 ) : ( ac + nrow(targets) ) ]
  ac <- ac + nrow(targets)
  J3   <- J123[ ( ac + 1 ) : ( ac + nrow(targets) ) ]
  # -------------------------------------------------------------------------- #

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
        defaultDeviceOpt <- list( res = 200, height = 1500, width = 600, 
            pointsize = 12, units = 'px')
      }
      
      if ( (! is.na( deviceIndex )) && deviceIndex ==  5 ) {
        fileNameOpt <- list( file = outfile)
        defaultDeviceOpt <- list( height = 7.5, width = 3, 
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
  
  # -------------------------------------------------------------------------- #
  # Draw the plot

  # Define colors. Repeat colors if required
  colors <- rep_len( colors, 6 )
  panelTitleColors <- rep_len( panelTitleColors, 6)

  # Sets plot layout  
  par( mfrow= c( 6 , 1 ) )
  
  # Set inner and outer margins
  par( oma = outerMargins,
       mar = innerMargins )
   
  # If useBarplot is not defined then set it as true if main factor has two or
  # less points to show
  useBarplots <- ( ( ! is.null ( useBarplots ) )  &&  useBarplots ) |
                 ( is.null( useBarplots ) ) && ( ( nPoints <= 2 ) ) 
  
  # Add Bin count data
  data <- .prepareData( binCounts, gridPanels, targets, factorsAndValues )
  .makePlot( data, useBarplots, main = "Bin counts", col = colors[1], 
      factorsAndValues, justOneFactor, gridPanels, nPoints, legendAtSide,
      las.x, panelTitleCex, panelTitleColors[1], useHCColors, barWidth, barSpacer )
  
  # Add Gene count data
  data <- .prepareData( geneCounts, gridPanels, targets, factorsAndValues )
  .makePlot( data, useBarplots, main = "Gene counts", col = colors[2], 
      factorsAndValues, justOneFactor, gridPanels, nPoints, legendAtSide,
      las.x, panelTitleCex, panelTitleColors[2], useHCColors, barWidth, barSpacer )
  
  # Add PSI/PIR count data
  # There is only one psir value for each condition ( instead for each sample as
  # other profiles being plotted. The value for condition is repeated for each
  # sample in that condition and then processed as the other data.
  psir <- psir[, rep( colnames(psir), as.vector( table(targets$condition)[ unique(targets$condition) ] ))]
  colnames(psir) <- rownames( targets )
  data <- .prepareData( psir , gridPanels, targets, factorsAndValues )
  .makePlot( data, useBarplots, main = "PSI / PIR", col = colors[3], 
      factorsAndValues, justOneFactor, gridPanels, nPoints, legendAtSide,
      las.x, panelTitleCex, panelTitleColors[3], useHCColors, barWidth, barSpacer )
  
  # Add J1 Exc count data
  data <- .prepareData( J1, gridPanels, targets, factorsAndValues )
  .makePlot( data, useBarplots, main = "Incl. junc 1", col = colors[4], 
      factorsAndValues, justOneFactor, gridPanels, nPoints, legendAtSide,
      las.x, panelTitleCex, panelTitleColors[4], useHCColors, barWidth, barSpacer )
  
  # Add J2 Exc count data
  data <- .prepareData( J2, gridPanels, targets, factorsAndValues )
  .makePlot( data, useBarplots, main = "Incl. junc 2", col = colors[5], 
      factorsAndValues, justOneFactor, gridPanels, nPoints, legendAtSide,
      las.x, panelTitleCex, panelTitleColors[5], useHCColors, barWidth, barSpacer )
  
  # Add J3 Exc count data
  data <- .prepareData( J3, gridPanels, targets, factorsAndValues )
  .makePlot( data, useBarplots, main = "Exclu. junc", col = colors[6], 
      factorsAndValues, justOneFactor, gridPanels, nPoints, legendAtSide,
      las.x, panelTitleCex, panelTitleColors[6], useHCColors, barWidth, barSpacer )
  
  # Draw the main Title of the plot
  if ( is.null(main) ) { main = bin }
  title( main = bin, outer=TRUE)
  # -------------------------------------------------------------------------- #

  # -------------------------------------------------------------------------- #
  # Close graphic device if required
  if ( outputIsAValidFile ) {
    output <- capture.output( dev.off() )
  } 
  # -------------------------------------------------------------------------- #
  
}


#-------------------------------------------------
#Plot splicing pattern
#-------------------------------------------------
# colnames(reportes$col_16_pcp_16@binbased)[21]   <- "junction.nonuniformity"
# colnames(reportes$col_16_pcp_16@anchorbased)[7] <- "junction.nonuniformity"
# colnames(reportes$col_23_col_16@binbased)[21]   <- "junction.nonuniformity"
# colnames(reportes$col_23_col_16@anchorbased)[7] <- "junction.nonuniformity"
# colnames(reportes$pcp_23_pcp_16@binbased)[21]   <- "junction.nonuniformity"
# colnames(reportes$pcp_23_pcp_16@anchorbased)[7] <- "junction.nonuniformity"
# colnames(reportes$col_23_pcp_23@binbased)[21]   <- "junction.nonuniformity"
# colnames(reportes$col_23_pcp_23@anchorbased)[7] <- "junction.nonuniformity"

#region=NULL;exones=NULL;genePlot=TRUE;zoomRegion=1.5;chrMap=NULL;hCov=0.7;hJun=0.3;useLog=FALSE
#bamFiles <- c("col0_16.star.bam","pcp_16.star.bam")
#mergedBAMs <- data.frame(bam=bamFiles,condition=c("col_16","pcp_16"))
#region <- iss$region[1]
#iss <- integrateSignals(sr, asd)
#chrMap <- as.character(1:5)
#names(chrMap) <- paste0("Chr",1:5)
#.plotSplicingPattern(region, iss, counts, f, mergedBAMs, sr, exones, chrMap = chrMap)
.countJbyCondition<-function(jcoords,counts){
  x<-countsj(counts)[rownames(jcoords),]
  mjsum <-c()
  for(ifila in 1:nrow(x)){
    jsum<-c()
    for(ccond in counts@condition.order){
      cnames <- rownames(targets(counts))[which(targets(counts)$condition%in%ccond)]
      jsum <- c(jsum,sum(x[ifila,cnames]))
    }
    mjsum<-rbind(mjsum,jsum)
  }
  colnames(mjsum)<-counts@condition.order
  rownames(mjsum)<-rownames(x)
  return(mjsum)
} 


# region: coordenadas genomicas de una region de interes (e.g. Chr1:1054200-1054285). 
#         Debe ser elegida a partid del objeto iss (ver abajo)
# iss   : objeto obtenido a partir de integrateSplicingSignals
# counts: ASpliCounts object
# f     : ASpliFeatures object
# mergedBAMs: data frame con dos columnas: full_path_file_name_bam_file condition (compatible con nomenclatura de aspli)
# sr    : ASpliSplicingReport object
# exones: si GRangesList object que contiene GRanges de bines para cada transcripto del genoma. 
#         si es NULL no se plotea estructura de variantes
# genePlot: si TRUE se grafica el esquema de splicing de la region en el contexto del gen al cual pertenece
#           si FALSE se centra la figura en la region
# jCompletelyIncluded: si TRUE solo se grafican junturas no diferenciales que esten completamente incluidas en la region de ploteo
#                      si FALSE se grafican junturas no diferenciales que overlapeen con la region de ploteo
# zoomRegion: factor para magnificar el area ploteada centrada en la region de interes
# chrMap    : vector nombrado que permite mapear nombres de cromosomas *utilizados en la denofinicion de region) en aquellos 
#             utilizados en bams
# 
# chrMap <- 1:5
# names(chrMap) <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
# mergedBams <- c("/xdata1/porcupine_data/STAR/mergedBAMS/col0_16.star.bam",
#                 "/xdata1/porcupine_data/STAR/mergedBAMS/col0_23.star.bam",
#                 "/xdata1/porcupine_data/STAR/mergedBAMS/pcp_16.star.bam",
#                 "/xdata1/porcupine_data/STAR/mergedBAMS/pcp_23.star.bam")
# mergedBAMs <- data.frame(full_path_file_name_bam_file = mergedBams, condition = c("col_16", "col_23", "pcp_16", "pcp_23"))
#.plotSplicingPattern(region, iss, counts, f, mergedBAMs, sr, chrMap = chrMap)

.plotSplicingPattern<-function(region=NULL,iss,counts,f,mergedBAMs,sr,asd,genePlot=TRUE,jCompletelyIncluded=TRUE,
                              zoomRegion=1.5,useLog=FALSE,tcex=2){
  #region <- r
  #iss <- is
  #counts
  #f <- features
  #mergedBAMs <- mergedBams
  
  if(length(transcriptExons(f))==1){
    warning("transcriptExons(f) has only one element. This could imply that something is wrong with the ASpliFeature object.\n")
  }
  
  variants = transcriptExons(f)
  
  
  #alturas relativas de paneles de coverage y junturas
  hCov=0.7
  hJun=0.3
  
  greg <- region
  nConditions <- nrow(mergedBAMs)
  
  iiss  <- iss[iss$region %in% greg,]
  if(nrow(iiss)==0){
    warning("No valid region provided. Out.\n")
    return()
  }
  
  geneName<-as.character(iiss[,"locus"])
  
  #encuentra x: GRanges para graficar
  if(!genePlot){
    
    #nombre de cromosoma en aspli
    aspli.chr <- strsplit2(strsplit2(iiss$region,":")[1], "Chr")
    if(ncol(aspli.chr) == 2){
      aspli.chr <- aspli.chr[, 2]
    }else{
      aspli.chr <- aspli.chr[, 1]
    }
    
    #nombre de cromosomas en features
    features.chr<-levels(seqnames(featuresb(f))@values)
    
    if(is.na(match(aspli.chr,features.chr))){
      #if(is.null(chrMap)){ #chrMap no se define mas asi que deberia entrar siempre aca
        warning(paste("No se pudo mapear nombres de cromosomas.",
                      "\n aspli.chr=",aspli.chr,
                      "\n features.chr=",paste(features.chr,collapse="/")))
      #}else{
       # chr <- features.chr[match(aspli.chr,features.chr)]
      #}
    }else{
      chr <- aspli.chr
    }
    
    roi   <- as.numeric(strsplit2(strsplit2(iiss$region,":")[2],"-"))
    #si hay una J3 != NA la uso para definir el rango
    if(!is.na(iiss$J3)){
      #si J3 se movio en un cluster, uso el cluster para definir el rango
      #  if(iiss$jl==1){
      #    #TODO
      #  }else{
      
      #roiJ3   <- as.numeric(strsplit2(as.character(iiss$J3),".",fixed=TRUE)[2:3])
      
      #chequeo si tengo mas de una J3
      j3aux <- strsplit2(iiss$J3,";")
      roiJ3<-c()
      for(ij3 in seq_along(j3aux)){
        roiJ3   <- c(roiJ3,as.numeric(strsplit2(j3aux[ij3],".",fixed=TRUE)[2:3]))
      }
      roiJ3 <- range(roiJ3)
      
      #  }
      
      roi[1] <- min(roi,roiJ3)
      roi[2] <- max(roi,roiJ3)
    }
    
    delta <- roi[2]-roi[1]
    
    zroi  <- c(roi[1]-delta*zoomRegion/2 , roi[2]+delta*zoomRegion/2)
    zdelta<-zroi[2]-zroi[1]
    
    gr    <- as(paste0(chr,":",zroi[1],"-",zroi[2]),"GRanges")
    hits  <- findOverlaps(gr,featuresb(f))
    bins  <-featuresb(f)[subjectHits(hits)]
    
    bins  <- bins[mcols(bins)$feature!="Io",]
    
    #redefino la zoomedroi en base a los bines que overlapeaban la zroi original
    delta <- max(end(bins))-min(start(bins))
    #zroi  <- c(start(bins)[1]-delta*zoomRegion/2 , end(bins)[length(bins)]+delta*zoomRegion/2)
    zroi  <- c(start(bins)[1], end(bins)[length(bins)])
    
  }else{
    
    #identifico coordenadas de biones ebinsonicos
    #iE <- which(mcols(featuresb(f))$locus%in%geneName & mcols(featuresb(f))$feature=="E"
    iE <- which(mcols(featuresb(f))$locus%in%geneName)
    if(length(iE)==0){
      #Is this a monoexonic gene?
      if(length(featuresg(f)[[geneName]])==1){
        bins <- featuresg(f)[geneName]
        mcols(bins)<-data.frame(locus=geneName,bin="E001",feature="mono E",
                                symbol=mcols(bins)$symbol,
                                locus_overlap=mcols(bins)$locus_overlap,
                                class="E",
                                event="novel",
                                eventJ="-")
        zroi <- roi <- as.numeric(c(start(bins)[1],end(bins)[length(bins)]))
      }else{
       warning("Something terribly wrong happened...No annotation data
              for",geneName," could be found!\nProbably a monoexonic gene...")
       return()
      }
    }else{
     bins <- featuresb(f)[iE,]
     bins <- bins[mcols(bins)$feature!="Io",]
     zroi <- roi <- c(start(bins)[1],end(bins)[length(bins)])
    }
    
    
    
  }
  
  #Armado del layout
  if(is.null(variants)){
    hh <- c(rep(c(hCov,hJun)/nConditions,nConditions),0.1)
    hh <- hh/sum(hh)
  }else{
    hh <- c(rep(c(hCov,hJun)/nConditions,nConditions),(hCov+hJun)/nConditions)
    hh <- hh/sum(hh)
  }
  #layout(matrix(c((2*nConditions+1):1,rep(1+2*nConditions+1:nConditions,each=2),(3*nConditions+2)),ncol=2),width=c(0.8,.2),height=hh)
  layout(matrix((2*nConditions+1):1),widths=c(0.8,.2),heights=hh)
  
  par(mar=c(.5, 1.1, .5, 1.1))
  
  rownames(mergedBAMs)<-mergedBAMs[,1]
  
  
  # Si hay variants...dibujo variantes
  if(!is.null(variants)){
    #transcriptGene<-strsplit2(names(variants),".",fixed=TRUE)[,1]
    transcriptGene <- mcols(variants)[,"gene"]
    
    ig<-which(transcriptGene%in%geneName)
    tex <- variants[ig]
    tex <- tex[order(names(tex))]
    
    limvariants<-unlist(lapply(tex,function(x){
      return(c(start(x),end(x)))
    }))
    
    plot(0,typ="n",xlim=range(limvariants),ylim=c(0,length(tex)+3),axes=FALSE,xlab="",ylab="")
    
    if(length(tex)>1 | !genePlot){
      for(itex in 1:length(tex)){
        bbins <- tex[[itex]]
        y<-length(tex)-(itex-1)
        lines(c(start(bbins[1]),end(bbins[length(bbins)])),c(y,y),col="gray")
        rect(start(bbins),y-.3,end(bbins),y+0.3,col=rgb(.7,.7,.7),border="black")
        
      }
    }
    ycollapsed <- length(tex)+2
  }else{
    ycollapsed <- 0
  }
  
  # Bines del genoma anotado que overlapean con el roi
  iE   <- mcols(bins)$feature%in%c("E","mono E")
  iroi <- (start(bins)>=roi[1] & start(bins)<roi[2])
  
  nbines <- length(bins)
  delta <- end(bins[nbines])-start(bins[1])
  
  if(is.null(variants)){
    plot(0,typ="n",xlim=c(start(bins[1]),end(bins[nbines])),ylim=c(-.5,.5),axes=FALSE,xlab="",ylab="")
  }
  lines(c(start(bins[1]),end(bins[nbines])),c(ycollapsed,ycollapsed),col="gray")
  rect(start(bins)[iE],ycollapsed-.45,end(bins)[iE],ycollapsed+.45,col="white")
  
  #remarco la roi original
    xroi   <- as.numeric(strsplit2(strsplit2(iiss$region,":")[2],"-"))
    lines(xroi,rep(ycollapsed-.6,2),col="black",lwd=3)
    
  #grafico el zoom
  if(!genePlot){
    lines(c(start(bins[1]),par("usr")[1]),c(ycollapsed+.45,par("usr")[4]),lty=1,col="gray")
    lines(c(end(bins[nbines]),par("usr")[2]),c(ycollapsed+.45,par("usr")[4]),lty=1,col="gray")
  }
  
  #bines diferenciales
  
  #que bines diferenciales hay en la region?
  iE<-which(names(bins)%in%iiss$bin)
  if(length(iE)>0){
    for(iie in seq_along(iE)){
      if(mcols(bins)$feature[iE[iie]]%in%c("Io","I")){
        lines(c(start(bins[iE[iie]]),end(bins[iE[iie]])),
              c(ycollapsed,ycollapsed),col="orange",lwd=3)
      }else{
        if(iiss$b!=0 | iiss$bjs!=0){
         if(iiss$b!=0 & iiss$bjs==0) cc <- "orange"
         if(iiss$b!=0 & iiss$bjs!=0) cc <- "red"
         if(iiss$b==0 & iiss$bjs!=0) cc <- "yellow"
         rect(start(bins[iE[iie]]),
              ycollapsed-.45,
              end(bins[iE[iie]]),
              ycollapsed+.45,col=cc,border=NA)
        } 
      }
    }
  }
  
  
  lines(par("usr")[1:2],rep(ycollapsed-1,each=2),lty=2,col="gray")
  
  # Junturas
  jcount1<-jcount2<-jcount0<-nj0<-nj1<-nj2<-0
  # aspli.chr <- strsplit2(strsplit2(iiss$region,":")[1], "[cC]hr")
  # if(ncol(aspli.chr) == 2){
  #   aspli.chr <- aspli.chr[, 2]
  # }else{
  #   aspli.chr <- aspli.chr[, 1]
  # }
  aspli.chr <- strsplit2(iiss$region,":")[1]
  
  
  
  
  #que junturas estan dentro de la zona?
  jsplit <- strsplit2(rownames(junctionsPJU(asd)),".",fixed=TRUE)
  if(jCompletelyIncluded){
    ijsplit <- which(jsplit[,1]==aspli.chr & as.numeric(jsplit[,2])>=zroi[1] & as.numeric(jsplit[,3])<=zroi[2])
  }else{
    ijsplit <- which(jsplit[,1]==aspli.chr & 
                       ((as.numeric(jsplit[,2])>=zroi[1] & as.numeric(jsplit[,2])<=zroi[2]) |
                          (as.numeric(jsplit[,3])>=zroi[1] & as.numeric(jsplit[,3])<=zroi[2])) )
  }
  

  # busco anchor junctions
  js <- unique(anchorbased(sr)$junction)
  aux <- strsplit2(js,".",fixed=TRUE)
  # ijs <-which(aux[,1]==aspli.chr &
  #               as.numeric(aux[,2])>=zroi[1] &
  #               as.numeric(aux[,2])<=zroi[2])  #ACA HAbia un error!
  if(jCompletelyIncluded){
    ijs <- which(aux[,1]==aspli.chr & as.numeric(aux[,2])>=zroi[1] & as.numeric(aux[,3])<=zroi[2])
  }else{
    ijs <- which(aux[,1]==aspli.chr & 
                   ((as.numeric(aux[,2])>=zroi[1] & as.numeric(aux[,2])<=zroi[2]) |
                      (as.numeric(aux[,3])>=zroi[1] & as.numeric(aux[,3])<=zroi[2])) )
  }
  jcoords1<-c()
  if(length(ijs)>0){
    #jcoords1          <- matrix(as.numeric(aux[ijs,]),ncol=3)
    jcoords1          <- data.frame(aux[ijs,,drop=FALSE],stringsAsFactors = FALSE)
    jcoords1[,2] <- as.numeric(jcoords1[,2])
    jcoords1[,3] <- as.numeric(jcoords1[,3])
    rownames(jcoords1)<-js[ijs]
    nj1               <- nrow(jcoords1)
  }
  
  #analizo locale 
  js <- unique(localebased(sr)$junction)
  aux <- strsplit2(js,".",fixed=TRUE)
  # ijs <-which(aux[,1]==aspli.chr &
  #               as.numeric(aux[,2])>=zroi[1] &
  #               as.numeric(aux[,2])<=zroi[2]) #ACA HAbia un error!
  if(jCompletelyIncluded){
    ijs <- which(aux[,1]==aspli.chr & as.numeric(aux[,2])>=zroi[1] & as.numeric(aux[,3])<=zroi[2])
  }else{
    ijs <- which(aux[,1]==aspli.chr & 
                   ((as.numeric(aux[,2])>=zroi[1] & as.numeric(aux[,2])<=zroi[2]) |
                      (as.numeric(aux[,3])>=zroi[1] & as.numeric(aux[,3])<=zroi[2])) )
  }
  jcoords2<-c()
  if(length(ijs)>0){
    #jcoords2          <- matrix(as.numeric(aux[ijs,]),ncol=3)
    jcoords2          <- aux[ijs,]
    jcoords2          <- data.frame(aux[ijs,,drop=FALSE],stringsAsFactors = FALSE)
    jcoords2[,2] <- as.numeric(jcoords2[,2])
    jcoords2[,3] <- as.numeric(jcoords2[,3])
    
    rownames(jcoords2)<-js[ijs]
    nj2               <- nrow(jcoords2)
  }
  
  # junturas en la region que no son anchor ni locale
  jcoords0<-c()
  if(length(ijsplit)>0){ 
    jcoords0 <- data.frame(chr=jsplit[ijsplit,1],#chr=as.numeric(jsplit[ijsplit,1]),
                           start=as.numeric(jsplit[ijsplit,2]),
                           end=as.numeric(jsplit[ijsplit,3]))
    rownames(jcoords0)<-rownames(junctionsPJU(asd))[ijsplit]
    
    jcoords0 <- jcoords0[!rownames(jcoords0)%in%unique(c(rownames(jcoords1),rownames(jcoords2))),]
    nj0               <- nrow(jcoords0)                       
  }
  
  
  #si detecto locale o anchor lo marco  
  if(nj1>0) abline(v=unique(c(jcoords1[,2],jcoords1[,3])),col="lightblue",lty=3)
  if(nj2>0) abline(v=unique(c(jcoords2[,2],jcoords2[,3])),col="lightgreen",lty=3)
  
  #dibujo paneles por condicion: junturas y coverage en la region de interes
  for(icond in 1:nConditions){
    # print(paste0("samtools depth -r ",
    #              paste0(aspli.chr,":",zroi[1],"-",zroi[2]," "),
    #              mergedBAMs[icond, 1]))
    ad <- system(paste0("samtools depth -r ",
                        paste0(aspli.chr,":",zroi[1],"-",zroi[2]," "),
                        mergedBAMs[icond, 1]), intern = TRUE)
    #ad <- matrix(as.numeric(strsplit2(ad,"\t")),ncol=3)
    ad <- strsplit2(ad,"\t")
    ad <- data.frame(ad,stringsAsFactors = FALSE)
    ad[,2]<-as.numeric(ad[,2])
    ad[,3]<-as.numeric(ad[,3])
    yylim <- range(ad[,3])

    nrep <- table(counts@targets$condition)[mergedBAMs[icond,2]]
    #me quedo con las que pasan un filtro de minima   
    jjcoords0<-jcoords0
    if(nj0>0){
      jcount0 <- .countJbyCondition(jcoords0,counts)[,mergedBAMs[icond,2],drop=FALSE]
      jok     <- rownames(jcount0)[jcount0>5*nrep]
      nj0     <- length(jok)
      if(nj0>0){
        jjcoords0  <- jcoords0[jok,,drop=FALSE]
        jcount0    <- jcount0[jok,,drop=FALSE]
      }  
    }
    
    jjcoords1<-jcoords1
    if(nj1>0){
      jcount1 <- .countJbyCondition(jcoords1,counts)[,mergedBAMs[icond,2],drop=FALSE]
      #   jok     <- rownames(jcount1)[jcount1>5*nrep]
      #   nj1     <- length(jok)
      #   if(nj1>0)  jjcoords1  <- jcoords1[jok,,drop=FALSE]
    }
    
    jjcoords2<-jcoords2
    if(nj2>0){
      jcount2 <- .countJbyCondition(jcoords2,counts)[,mergedBAMs[icond,2],drop=FALSE]
      #   jok     <- rownames(jcount2)[jcount2>5*nrep]
      #   nj2     <- length(jok)
      #   if(nj2>0)  jjcoords2  <- jcoords2[jok,,drop=FALSE]
    }
    

    #panel junturas
    if(FALSE){
      plot(0,typ="n",xlim=c(start(bins[1]),end(bins[nbines])),ylim=c(1,nj0+nj1+nj2+3),axes=FALSE,xlab="",ylab="")
      jmaxcount <- max(c(jcount0,jcount1,jcount2))
      if(nj0>0){
        jcounts <- jcount0
        jcoords <- jjcoords0
        ww <- jcounts/jmaxcount*3
        for(ij in 1:nj0){
          lines(jcoords[ij,2:3],rep(ij,2),lwd=ww[rownames(jcoords)[ij],1],col="lightgray")
          points(jcoords[ij,2:3],rep(ij,2),pch=18,cex=0.5,col="lightgray")
        }
      }
      if(nj1>0){
        jcounts <- jcount1
        jcoords <- jjcoords1
        ww <- jcounts/jmaxcount*3
        
        abline(v=unique(c(jcoords[,2],jcoords[,3])),col="lightblue",lty=3)
        for(ij in 1:nj1){
          lines(jcoords[ij,2:3],rep(nj0+ij,2),lwd=ww[rownames(jcoords)[ij],1],col="lightblue")
          points(jcoords[ij,2:3],rep(nj0+ij,2),pch=18,cex=0.5,col="lightblue")
          text(mean(as.numeric(jcoords[ij,2:3])),nj0+ij,jcounts[ij],pos=3,cex=tcex)
        }
      }
      if(nj2>0){
        jcounts <- jcount2
        jcoords <- jjcoords2
        ww <- jcounts/jmaxcount*3
        
        abline(v=unique(c(jcoords[,2],jcoords[,3])),col="lightgreen",lty=3)
        for(ij in 1:nj2){
          lines(jcoords[ij,2:3],rep(nj0+nj1+ij,2),lwd=ww[rownames(jcoords)[ij],1],col="lightgreen")
          points(jcoords[ij,2:3],rep(nj0+nj1+ij,2),pch=18,cex=0.5,col="lightgreen")
          text(mean(as.numeric(jcoords[ij,2:3])),nj0+nj1+ij,jcounts[ij],pos=3,cex=tcex)
        }
      }
    }
    
    #plot(0,typ="n",xlim=c(start(bins[1]),end(bins[nbines])),ylim=c(1,nj0+nj1+nj2+3),axes=FALSE,xlab="",ylab="")
    
    extra <- 5
    njlevels <- nj0+nj1+nj2+extra
    
    j12      <- intersect(rownames(jcoords1),rownames(jcoords2))
    njlevels <- nj0+nj1+nj2-length(j12)+extra
    
    xx<-as.numeric(c(start(bins[1]),end(bins[nbines])))
    plot(0,typ="n",xlim=xx,ylim=c(1,njlevels),axes=FALSE,xlab="",ylab="")
    jmaxcount <- max(c(jcount0,jcount1,jcount2))
    if(nj0>0){
      jcounts <- jcount0
      jcoords <- jjcoords0
      ww <- jcounts/jmaxcount*3
      for(ij in 1:nj0){
        lines(jcoords[ij,2:3],rep(ij,2),lwd=ww[rownames(jcoords)[ij],1],col="lightgray")
        points(jcoords[ij,2:3],rep(ij,2),pch=18,cex=0.5,col="lightgray")
      }
    }
    
    if(njlevels-nj0>extra){
      jcounts<-c()
      if(nj1>0) jcounts<-jcount1
      if(nj2>0){
          if(length(j12)>0){
            jcounts<-rbind(jcounts,jcount2[!rownames(jcount2)%in%j12,,drop=FALSE])
          }else{  
           jcounts<-rbind(jcounts,jcount2)
          }
      }
      #jcoords <- t(apply(cbind(rownames(jcounts,jcounts)),1,function(x){as.numeric(strsplit2(x,"[.]"))}))
      jcoords <- t(apply(cbind(rownames(jcounts,jcounts)),1,function(x){strsplit2(x,"[.]")}))
      rownames(jcoords)<-rownames(jcounts)
    
      ccolor  <- rep("lightgreen",nrow(jcounts))
      names(ccolor)<-rownames(jcounts)
      if(nj2>0) ccolor[rownames(jcount2)]<-"lightblue"
      if(length(j12)>0) ccolor[j12]<-"orange"
      
      ww <- jcounts/jmaxcount*3
      
      abline(v=unique(c(jcoords[,2],jcoords[,3])),col="lightblue",lty=3)
      for(ij in 1:length(ccolor)){
        #yij <- nj0+ij
        yij  <- ij * njlevels*0.8/(nj1+nj2-length(j12))
        lines(jcoords[ij,2:3],rep(yij,2),lwd=max(1,ww[rownames(jcoords)[ij],1]),col=ccolor[ij])
        points(jcoords[ij,2:3],rep(yij,2),pch=18,cex=2,col=ccolor[ij])
        text(mean(as.numeric(jcoords[ij,2:3])),yij,jcounts[ij],pos=2,cex=tcex)
      }
      
    }

    
    #coverage
    xx <- as.numeric(c(start(bins[1]),end(bins[nbines])))
    if(useLog){
      yylim[yylim==0]<-0.1
      ad[ad[,3]==0,3]<-0.1
      plot(0.01,typ="n",xlim=xx,ylim=yylim,axes=FALSE,xlab="",ylab="",log="y")
    }else{
      plot(0,typ="n",xlim=xx,ylim=yylim,axes=FALSE,xlab="",ylab="")
    }
    polygon(c(ad[1,2],ad[,2],ad[nrow(ad),2]),c(0.01,ad[,3],0.01)
            ,border=NA,col=topo.colors(nConditions,1)[icond],fillOddEven = TRUE)
    lines(c(ad[1,2],ad[,2],ad[nrow(ad),2]),  c(0.01,ad[,3],0.01),
          col=topo.colors(nConditions,1)[icond])
    text(ad[1,2],0.01,paste0("[0-",max(ad[,3]),"]"),cex=tcex,adj=c(-.25,-.5))
    
    if(nj1>0)abline(v=unique(c(jcoords1[,2],jcoords1[,3])),col="lightgreen",lty=3)
    if(nj2>0)abline(v=unique(c(jcoords2[,2],jcoords2[,3])),col="lightblue",lty=3)
    
    #remarco la roi original
    if(icond==1){
      xroi   <- as.numeric(strsplit2(strsplit2(iiss$region,":")[2],"-"))
      lines(xroi,rep(par("usr")[3],2),col="black",lwd=3)
      # rect(xroi[1],par("usr")[3],xroi[2],par("usr")[3]+diff(par("usr")[4:3])*.1,col="black",density=2)
    }
    legend("topleft", as.character(mergedBAMs$condition[icond]),bty="n")
  }
  mtext(paste(iiss$locus,iiss$region),line=-.5,cex=0.8)
  
  
  #vamos por el texto:
  if(FALSE){
    # Junturas Anchor
    plot(1,axes=FALSE,type="n",xlab="",ylab="")
    if(iiss$ja==1){
      col1 <- c( "junction","junction.annotated","junction.fdr","junction.nonuniformity","junction.participation","cluster.fdr","dPIR")
      a1 <- sr@anchorbased[sr@anchorbased$junction%in%iiss$J3,col1]
      if(nrow(a1)>0){
        a1[col1[3:7]]<-signif(a1[col1[3:7]],2)
        a1<-cbind(names(a1),t(a1))
        a1<-apply(a1,1,function(x){return(paste(x,collapse=": "))})
        names(a1)<- c("junction","annotated","j.fdr","non-unif","participation","cluster.fdr","D_PIR")
        a1<-paste("  ",a1,sep = "")
        a1<-c("Junction Anchorage:",a1)
        legend("topleft",a1,cex=tcex,bty="n",inset=0.02) 
        
        # a2 <- sr@anchorbased[sr@anchorbased$junction%in%iiss$J3,8+c(6,4,2,5,3,1)]*nrep
        # a2 <- cbind(names(a2),t(a2))
        # a2 <- apply(a2,1,function(x){return(paste(x,collapse=": "))})
        # a2<-paste("  ",a2,sep = "")
        
        # legend("bottomleft",a2,cex=tcex,bty="n",inset=0.02)
      }else{
        a1<-c("Junction Anchorage:","  no J3 passed the filter.")
        legend("topleft",a1,cex=tcex,bty="n",inset=0.02)   
      }
    }
    
    # Junturas locale
    plot(1,axes=FALSE,type="n",xlab="",ylab="")
    if(iiss$jl==1){
      col1 <- c( "junction","junction.annotated","junction.fdr","junction.participation","cluster.fdr")
      a1 <- sr@localebased[sr@localebased$junction%in%iiss$J3,col1]
      a1[col1[3:5]]<-signif(a1[col1[3:5]],2)
      a1<-cbind(names(a1),t(a1))
      a1<-apply(a1,1,function(x){return(paste(x,collapse=": "))})
      names(a1)<- c("junction","annotated","j.fdr","participation","cluster.fdr")
      a1<-paste("  ",a1,sep = "")
      a1<-c("Junction Locale:",a1)
      
      legend("topleft",a1,cex=tcex,bty="n",inset=0.02) 
      
      
    }
    
    # bin info
    plot(1,axes=FALSE,type="n",xlab="",ylab="")
    if(iiss$b==1 | iiss$bjs==1){
      col1 <- c( "bin","bin.fdr","junction.dPIR","junction.dPSI")
      a1 <- sr@binbased[sr@binbased$bin%in%iiss$bin,col1]
      a1[col1[2:4]]<-signif(a1[col1[2:4]],2)
      names(a1)<-c("bin","bin.fdr","D_PIR","D_PSI")
      a1<-cbind(names(a1),t(a1))
      a1<-apply(a1,1,function(x){return(paste(x,collapse=": "))})
      a1<-paste("  ",a1,sep = "")
      a1<-c("Bin:",a1)
      legend("topleft",a1,cex=tcex,bty="n",inset=0.02) 
      
    }
  }
  
}


.plotSplicingPattern.orig<-function(region=NULL,iss,counts,f,mergedBAMs,sr,asd,genePlot=TRUE,jCompletelyIncluded=TRUE,
                               zoomRegion=1.5,useLog=FALSE,tcex=1){
  #region <- r
  #iss <- is
  #counts
  #f <- features
  #mergedBAMs <- mergedBams
  
  if(length(transcriptExons(f))==1){
    warning("transcriptExons(f) has only one element. This could imply that something is wrong with the ASpliFeature object.\n")
  }
  
  exones = transcriptExons(f)
  
  
  #alturas relativas de paneles de coverage y junturas
  hCov=0.7
  hJun=0.3
  
  greg <- region
  nConditions <- nrow(mergedBAMs)
  
  iiss  <- iss[iss$region %in% greg,]
  if(nrow(iiss)==0){
    warning("No valid region provided. Out.\n")
    return()
  }
  
  geneName<-as.character(iiss[,"locus"])
  
  #encuentra x: GRanges para graficar
  if(!genePlot){
    
    #nombre de cromosoma en aspli
    aspli.chr <- strsplit2(strsplit2(iiss$region,":")[1], "Chr")
    if(ncol(aspli.chr) == 2){
      aspli.chr <- aspli.chr[, 2]
    }else{
      aspli.chr <- aspli.chr[, 1]
    }
    
    #nombre de cromosomas en features
    features.chr<-levels(seqnames(featuresb(f))@values)
    
    if(is.na(match(aspli.chr,features.chr))){
      #if(is.null(chrMap)){ #chrMap no se define mas asi que deberia entrar siempre aca
      warning(paste("No se pudo mapear nombres de cromosomas.",
                    "\n aspli.chr=",aspli.chr,
                    "\n features.chr=",paste(features.chr,collapse="/")))
      #}else{
      # chr <- features.chr[match(aspli.chr,features.chr)]
      #}
    }else{
      chr <- aspli.chr
    }
    
    roi   <- as.numeric(strsplit2(strsplit2(iiss$region,":")[2],"-"))
    #si hay una J3 != NA la uso para definir el rango
    if(!is.na(iiss$J3)){
      #si J3 se movio en un cluster, uso el cluster para definir el rango
      #  if(iiss$jl==1){
      #    #TODO
      #  }else{
      
      #roiJ3   <- as.numeric(strsplit2(as.character(iiss$J3),".",fixed=TRUE)[2:3])
      
      #chequeo si tengo mas de una J3
      j3aux <- strsplit2(iiss$J3,";")
      roiJ3<-c()
      for(ij3 in seq_along(j3aux)){
        roiJ3   <- c(roiJ3,as.numeric(strsplit2(j3aux[ij3],".",fixed=TRUE)[2:3]))
      }
      roiJ3 <- range(roiJ3)
      
      #  }
      
      roi[1] <- min(roi,roiJ3)
      roi[2] <- max(roi,roiJ3)
    }
    
    delta <- roi[2]-roi[1]
    
    zroi  <- c(roi[1]-delta*zoomRegion/2 , roi[2]+delta*zoomRegion/2)
    zdelta<-zroi[2]-zroi[1]
    
    gr    <- as(paste0(chr,":",zroi[1],"-",zroi[2]),"GRanges")
    hits  <- findOverlaps(gr,featuresb(f))
    bins  <-featuresb(f)[subjectHits(hits)]
    
    bins  <- bins[mcols(bins)$feature!="Io",]
    
    #redefino la zoomedroi en base a los bines que overlapeaban la zroi original
    delta <- max(end(bins))-min(start(bins))
    #zroi  <- c(start(bins)[1]-delta*zoomRegion/2 , end(bins)[length(bins)]+delta*zoomRegion/2)
    zroi  <- c(start(bins)[1], end(bins)[length(bins)])
    
  }else{
    
    #identifico coordenadas de biones ebinsonicos
    #iE <- which(mcols(featuresb(f))$locus%in%geneName & mcols(featuresb(f))$feature=="E"
    iE <- which(mcols(featuresb(f))$locus%in%geneName)
    if(length(iE)==0){
      #Is this a monoexonic gene?
      if(length(featuresg(f)[[geneName]])==1){
        bins <- featuresg(f)[geneName]
        mcols(bins)<-data.frame(locus=geneName,bin="E001",feature="mono E",
                                symbol=mcols(bins)$symbol,
                                locus_overlap=mcols(bins)$locus_overlap,
                                class="E",
                                event="novel",
                                eventJ="-")
        zroi <- roi <- as.numeric(c(start(bins)[1],end(bins)[length(bins)]))
      }else{
        warning("Something terribly wrong happened...No annotation data
                for",geneName," could be found!\nProbably a monoexonic gene...")
        return()
      }
    }else{
      bins <- featuresb(f)[iE,]
      bins <- bins[mcols(bins)$feature!="Io",]
      zroi <- roi <- c(start(bins)[1],end(bins)[length(bins)])
    }
    
    
    
    }
  
  #Armado del layout
  if(is.null(exones)){
    hh <- c(rep(c(hCov,hJun)/nConditions,nConditions),0.1)
    hh <- hh/sum(hh)
  }else{
    hh <- c(rep(c(hCov,hJun)/nConditions,nConditions),(hCov+hJun)/nConditions)
    hh <- hh/sum(hh)
  }
  #layout(matrix(c((2*nConditions+1):1,rep(1+2*nConditions+1:nConditions,each=2),(3*nConditions+2)),ncol=2),width=c(0.8,.2),height=hh)
  layout(matrix((2*nConditions+1):1),widths=c(0.8,.2),heights=hh)
  
  par(mar=c(.5, 1.1, .5, 1.1))
  
  rownames(mergedBAMs)<-mergedBAMs[,1]
  
  
  # Si hay exones...dibujo variantes
  if(!is.null(exones)){
    transcriptGene<-strsplit2(names(exones),".",fixed=TRUE)[,1]
    
    ig<-which(transcriptGene%in%geneName)
    tex <- exones[ig]
    tex <- tex[order(names(tex))]
    
    limExones<-unlist(lapply(tex,function(x){
      return(c(start(x),end(x)))
    }))
    
    plot(0,typ="n",xlim=range(limExones),ylim=c(0,length(tex)+3),axes=FALSE,xlab="",ylab="")
    
    if(length(tex)>1 | !genePlot){
      for(itex in 1:length(tex)){
        bbins <- tex[[itex]]
        y<-length(tex)-(itex-1)
        lines(c(start(bbins[1]),end(bbins[length(bbins)])),c(y,y),col="gray")
        rect(start(bbins),y-.3,end(bbins),y+0.3,col=rgb(.7,.7,.7),border="black")
        
      }
    }
    ycollapsed <- length(tex)+2
  }else{
    ycollapsed <- 0
  }
  
  # Bines del genoma anotado que overlapean con el roi
  iE   <- mcols(bins)$feature%in%c("E","mono E")
  iroi <- (start(bins)>=roi[1] & start(bins)<roi[2])
  
  nbines <- length(bins)
  delta <- end(bins[nbines])-start(bins[1])
  
  if(is.null(exones)){
    plot(0,typ="n",xlim=c(start(bins[1]),end(bins[nbines])),ylim=c(-.5,.5),axes=FALSE,xlab="",ylab="")
  }
  lines(c(start(bins[1]),end(bins[nbines])),c(ycollapsed,ycollapsed),col="gray")
  rect(start(bins)[iE],ycollapsed-.45,end(bins)[iE],ycollapsed+.45,col="white")
  
  #remarco la roi original
  xroi   <- as.numeric(strsplit2(strsplit2(iiss$region,":")[2],"-"))
  lines(xroi,rep(ycollapsed-.6,2),col="black",lwd=3)
  
  #grafico el zoom
  if(!genePlot){
    lines(c(start(bins[1]),par("usr")[1]),c(ycollapsed+.45,par("usr")[4]),lty=1,col="gray")
    lines(c(end(bins[nbines]),par("usr")[2]),c(ycollapsed+.45,par("usr")[4]),lty=1,col="gray")
  }
  
  #bines diferenciales
  
  #que bines diferenciales hay en la region?
  iE<-which(names(bins)%in%iiss$bin)
  if(length(iE)>0){
    for(iie in seq_along(iE)){
      if(mcols(bins)$feature[iE[iie]]%in%c("Io","I")){
        lines(c(start(bins[iE[iie]]),end(bins[iE[iie]])),c(ycollapsed,ycollapsed),col="orange",lwd=3)
      }else{
        if(iiss$b!=0 | iiss$bjs!=0){
          if(iiss$b!=0 & iiss$bjs==0) cc <- "orange"
          if(iiss$b!=0 & iiss$bjs!=0) cc <- "red"
          if(iiss$b==0 & iiss$bjs!=0) cc <- "yellow"
          rect(start(bins[iE[iie]]),ycollapsed-.45,end(bins[iE[iie]]),ycollapsed+.45,col=cc,border=NA)
        } 
      }
    }
  }
  
  
  lines(par("usr")[1:2],rep(ycollapsed-1,each=2),lty=2,col="gray")
  
  # Junturas
  jcount1<-jcount2<-jcount0<-nj0<-nj1<-nj2<-0
  aspli.chr <- strsplit2(strsplit2(iiss$region,":")[1], "[cC]hr")
  if(ncol(aspli.chr) == 2){
    aspli.chr <- aspli.chr[, 2]
  }else{
    aspli.chr <- aspli.chr[, 1]
  }
  
  #que junturas estan dentro de la zona?
  jsplit <- strsplit2(rownames(junctionsPJU(asd)),".",fixed=TRUE)
  if(jCompletelyIncluded){
    ijsplit <- which(jsplit[,1]==aspli.chr & as.numeric(jsplit[,2])>=zroi[1] & as.numeric(jsplit[,3])<=zroi[2])
  }else{
    ijsplit <- which(jsplit[,1]==aspli.chr & 
                       ((as.numeric(jsplit[,2])>=zroi[1] & as.numeric(jsplit[,2])<=zroi[2]) |
                          (as.numeric(jsplit[,3])>=zroi[1] & as.numeric(jsplit[,3])<=zroi[2])) )
  }
  
  
  # busco anchor junctions
  js <- unique(anchorbased(sr)$junction)
  aux <- strsplit2(js,".",fixed=TRUE)
  # ijs <-which(aux[,1]==aspli.chr &
  #               as.numeric(aux[,2])>=zroi[1] &
  #               as.numeric(aux[,2])<=zroi[2])  #ACA HAbia un error!
  if(jCompletelyIncluded){
    ijs <- which(aux[,1]==aspli.chr & as.numeric(aux[,2])>=zroi[1] & as.numeric(aux[,3])<=zroi[2])
  }else{
    ijs <- which(aux[,1]==aspli.chr & 
                   ((as.numeric(aux[,2])>=zroi[1] & as.numeric(aux[,2])<=zroi[2]) |
                      (as.numeric(aux[,3])>=zroi[1] & as.numeric(aux[,3])<=zroi[2])) )
  }
  jcoords1<-c()
  if(length(ijs)>0){
    jcoords1          <- matrix(as.numeric(aux[ijs,]),ncol=3)
    rownames(jcoords1)<-js[ijs]
    nj1               <- nrow(jcoords1)
  }
  
  #analizo locale 
  js <- unique(localebased(sr)$junction)
  aux <- strsplit2(js,".",fixed=TRUE)
  # ijs <-which(aux[,1]==aspli.chr &
  #               as.numeric(aux[,2])>=zroi[1] &
  #               as.numeric(aux[,2])<=zroi[2]) #ACA HAbia un error!
  if(jCompletelyIncluded){
    ijs <- which(aux[,1]==aspli.chr & as.numeric(aux[,2])>=zroi[1] & as.numeric(aux[,3])<=zroi[2])
  }else{
    ijs <- which(aux[,1]==aspli.chr & 
                   ((as.numeric(aux[,2])>=zroi[1] & as.numeric(aux[,2])<=zroi[2]) |
                      (as.numeric(aux[,3])>=zroi[1] & as.numeric(aux[,3])<=zroi[2])) )
  }
  jcoords2<-c()
  if(length(ijs)>0){
    jcoords2          <- matrix(as.numeric(aux[ijs,]),ncol=3)
    rownames(jcoords2)<-js[ijs]
    nj2               <- nrow(jcoords2)
  }
  
  # junturas en la region que no son anchor ni locale
  jcoords0<-c()
  if(length(ijsplit)>0){ 
    jcoords0 <- data.frame(chr=as.numeric(jsplit[ijsplit,1]),
                           start=as.numeric(jsplit[ijsplit,2]),
                           end=as.numeric(jsplit[ijsplit,3]))
    rownames(jcoords0)<-rownames(junctionsPJU(asd))[ijsplit]
    
    jcoords0 <- jcoords0[!rownames(jcoords0)%in%unique(c(rownames(jcoords1),rownames(jcoords2))),]
    nj0               <- nrow(jcoords0)                       
  }
  
  
  #si detecto locale o anchor lo marco  
  if(nj1>0) abline(v=unique(c(jcoords1[,2],jcoords1[,3])),col="lightblue",lty=3)
  if(nj2>0) abline(v=unique(c(jcoords2[,2],jcoords2[,3])),col="lightgreen",lty=3)
  
  #dibujo paneles por condicion: junturas y coverage en la region de interes
  for(icond in 1:nConditions){
    # print(paste0("samtools depth -r ",
    #              paste0(aspli.chr,":",zroi[1],"-",zroi[2]," "),
    #              mergedBAMs[icond, 1]))
    ad <- system(paste0("samtools depth -r ",
                        paste0(aspli.chr,":",zroi[1],"-",zroi[2]," "),
                        mergedBAMs[icond, 1]), intern = TRUE)
    ad <- matrix(as.numeric(strsplit2(ad,"\t")),ncol=3)
    yylim <- range(ad[,3])
    
    nrep <- table(counts@targets$condition)[mergedBAMs[icond,2]]
    #me quedo con las que pasan un filtro de minima   
    jjcoords0<-jcoords0
    if(nj0>0){
      jcount0 <- .countJbyCondition(jcoords0,counts)[,mergedBAMs[icond,2],drop=FALSE]
      jok     <- rownames(jcount0)[jcount0>5*nrep]
      nj0     <- length(jok)
      if(nj0>0)  jjcoords0  <- jcoords0[jok,,drop=FALSE]
    }
    
    jjcoords1<-jcoords1
    if(nj1>0){
      jcount1 <- .countJbyCondition(jcoords1,counts)[,mergedBAMs[icond,2],drop=FALSE]
      #   jok     <- rownames(jcount1)[jcount1>5*nrep]
      #   nj1     <- length(jok)
      #   if(nj1>0)  jjcoords1  <- jcoords1[jok,,drop=FALSE]
    }
    
    jjcoords2<-jcoords2
    if(nj2>0){
      jcount2 <- .countJbyCondition(jcoords2,counts)[,mergedBAMs[icond,2],drop=FALSE]
      #   jok     <- rownames(jcount2)[jcount2>5*nrep]
      #   nj2     <- length(jok)
      #   if(nj2>0)  jjcoords2  <- jcoords2[jok,,drop=FALSE]
    }
    
    
    #panel junturas
    if(FALSE){
      plot(0,typ="n",xlim=c(start(bins[1]),end(bins[nbines])),ylim=c(1,nj0+nj1+nj2+3),axes=FALSE,xlab="",ylab="")
      jmaxcount <- max(c(jcount0,jcount1,jcount2))
      if(nj0>0){
        jcounts <- jcount0
        jcoords <- jjcoords0
        ww <- jcounts/jmaxcount*3
        for(ij in 1:nj0){
          lines(jcoords[ij,2:3],rep(ij,2),lwd=ww[rownames(jcoords)[ij],1],col="lightgray")
          points(jcoords[ij,2:3],rep(ij,2),pch=18,cex=0.5,col="lightgray")
        }
      }
      if(nj1>0){
        jcounts <- jcount1
        jcoords <- jjcoords1
        ww <- jcounts/jmaxcount*3
        
        abline(v=unique(c(jcoords[,2],jcoords[,3])),col="lightblue",lty=3)
        for(ij in 1:nj1){
          lines(jcoords[ij,2:3],rep(nj0+ij,2),lwd=ww[rownames(jcoords)[ij],1],col="lightblue")
          points(jcoords[ij,2:3],rep(nj0+ij,2),pch=18,cex=0.5,col="lightblue")
          text(mean(as.numeric(jcoords[ij,2:3])),nj0+ij,jcounts[ij],pos=3,cex=tcex)
        }
      }
      if(nj2>0){
        jcounts <- jcount2
        jcoords <- jjcoords2
        ww <- jcounts/jmaxcount*3
        
        abline(v=unique(c(jcoords[,2],jcoords[,3])),col="lightgreen",lty=3)
        for(ij in 1:nj2){
          lines(jcoords[ij,2:3],rep(nj0+nj1+ij,2),lwd=ww[rownames(jcoords)[ij],1],col="lightgreen")
          points(jcoords[ij,2:3],rep(nj0+nj1+ij,2),pch=18,cex=0.5,col="lightgreen")
          text(mean(as.numeric(jcoords[ij,2:3])),nj0+nj1+ij,jcounts[ij],pos=3,cex=tcex)
        }
      }
    }
    
    #plot(0,typ="n",xlim=c(start(bins[1]),end(bins[nbines])),ylim=c(1,nj0+nj1+nj2+3),axes=FALSE,xlab="",ylab="")
    
    extra <- 5
    njlevels <- nj0+nj1+nj2+extra
    
    j12      <- intersect(rownames(jcoords1),rownames(jcoords2))
    njlevels <- nj0+nj1+nj2-length(j12)+extra
    
    xx<-as.numeric(c(start(bins[1]),end(bins[nbines])))
    plot(0,typ="n",xlim=xx,ylim=c(1,njlevels),axes=FALSE,xlab="",ylab="")
    jmaxcount <- max(c(jcount0,jcount1,jcount2))
    if(nj0>0){
      jcounts <- jcount0
      jcoords <- jjcoords0
      ww <- jcounts/jmaxcount*3
      for(ij in 1:nj0){
        lines(jcoords[ij,2:3],rep(ij,2),lwd=ww[rownames(jcoords)[ij],1],col="lightgray")
        points(jcoords[ij,2:3],rep(ij,2),pch=18,cex=0.5,col="lightgray")
      }
    }
    
    if(njlevels-nj0>extra){
      jcounts<-c()
      if(nj1>0) jcounts<-jcount1
      if(nj2>0){
        if(length(j12)>0){
          jcounts<-rbind(jcounts,jcount2[!rownames(jcount2)%in%j12,,drop=FALSE])
        }else{  
          jcounts<-rbind(jcounts,jcount2)
        }
      }
      jcoords <- t(apply(cbind(rownames(jcounts,jcounts)),1,function(x){as.numeric(strsplit2(x,"[.]"))}))
      rownames(jcoords)<-rownames(jcounts)
      
      ccolor  <- rep("lightgreen",nrow(jcounts))
      names(ccolor)<-rownames(jcounts)
      if(nj2>0) ccolor[rownames(jcount2)]<-"lightblue"
      if(length(j12)>0) ccolor[j12]<-"orange"
      
      ww <- jcounts/jmaxcount*3
      
      abline(v=unique(c(jcoords[,2],jcoords[,3])),col="lightblue",lty=3)
      for(ij in 1:length(ccolor)){
        #yij <- nj0+ij
        yij  <- ij * njlevels*0.8/(nj1+nj2-length(j12))
        lines(jcoords[ij,2:3],rep(yij,2),lwd=max(1,ww[rownames(jcoords)[ij],1]),col=ccolor[ij])
        points(jcoords[ij,2:3],rep(yij,2),pch=18,cex=0.5,col=ccolor[ij])
        text(mean(as.numeric(jcoords[ij,2:3])),yij,jcounts[ij],pos=2,cex=tcex)
      }
      
    }
    
    
    #coverage
    xx <- as.numeric(c(start(bins[1]),end(bins[nbines])))
    if(useLog){
      plot(0.01,typ="n",xlim=xx,ylim=yylim,axes=FALSE,xlab="",ylab="",log="y")
    }else{
      plot(0,typ="n",xlim=xx,ylim=yylim,axes=FALSE,xlab="",ylab="")
    }
    polygon(c(ad[1,2],ad[,2],ad[nrow(ad),2]), c(0.01,ad[,3],0.01)
            ,border=NA,col=topo.colors(nConditions,1)[icond],fillOddEven = TRUE)
    lines(c(ad[1,2],ad[,2],ad[nrow(ad),2]),  c(0.01,ad[,3],0.01),
          col=topo.colors(nConditions,1)[icond])
    text(ad[1,2],0.01,paste0("[0-",max(ad[,3]),"]"),cex=tcex,adj=c(-.25,-.5))
    
    if(nj1>0)abline(v=unique(c(jcoords1[,2],jcoords1[,3])),col="lightblue",lty=3)
    if(nj2>0)abline(v=unique(c(jcoords2[,2],jcoords2[,3])),col="lightgreen",lty=3)
    
    #remarco la roi original
    if(icond==1){
      xroi   <- as.numeric(strsplit2(strsplit2(iiss$region,":")[2],"-"))
      lines(xroi,rep(par("usr")[3],2),col="black",lwd=3)
      # rect(xroi[1],par("usr")[3],xroi[2],par("usr")[3]+diff(par("usr")[4:3])*.1,col="black",density=2)
    }
    legend("topleft", as.character(mergedBAMs$condition[icond]),bty="n")
  }
  mtext(paste(iiss$locus,iiss$region),line=-.5,cex=0.8)
  
  
  #vamos por el texto:
  if(FALSE){
    # Junturas Anchor
    plot(1,axes=FALSE,type="n",xlab="",ylab="")
    if(iiss$ja==1){
      col1 <- c( "junction","junction.annotated","junction.fdr","junction.nonuniformity","junction.participation","cluster.fdr","dPIR")
      a1 <- sr@anchorbased[sr@anchorbased$junction%in%iiss$J3,col1]
      if(nrow(a1)>0){
        a1[col1[3:7]]<-signif(a1[col1[3:7]],2)
        a1<-cbind(names(a1),t(a1))
        a1<-apply(a1,1,function(x){return(paste(x,collapse=": "))})
        names(a1)<- c("junction","annotated","j.fdr","non-unif","participation","cluster.fdr","D_PIR")
        a1<-paste("  ",a1,sep = "")
        a1<-c("Junction Anchorage:",a1)
        legend("topleft",a1,cex=tcex,bty="n",inset=0.02) 
        
        # a2 <- sr@anchorbased[sr@anchorbased$junction%in%iiss$J3,8+c(6,4,2,5,3,1)]*nrep
        # a2 <- cbind(names(a2),t(a2))
        # a2 <- apply(a2,1,function(x){return(paste(x,collapse=": "))})
        # a2<-paste("  ",a2,sep = "")
        
        # legend("bottomleft",a2,cex=tcex,bty="n",inset=0.02)
      }else{
        a1<-c("Junction Anchorage:","  no J3 passed the filter.")
        legend("topleft",a1,cex=tcex,bty="n",inset=0.02)   
      }
    }
    
    # Junturas locale
    plot(1,axes=FALSE,type="n",xlab="",ylab="")
    if(iiss$jl==1){
      col1 <- c( "junction","junction.annotated","junction.fdr","junction.participation","cluster.fdr")
      a1 <- sr@localebased[sr@localebased$junction%in%iiss$J3,col1]
      a1[col1[3:5]]<-signif(a1[col1[3:5]],2)
      a1<-cbind(names(a1),t(a1))
      a1<-apply(a1,1,function(x){return(paste(x,collapse=": "))})
      names(a1)<- c("junction","annotated","j.fdr","participation","cluster.fdr")
      a1<-paste("  ",a1,sep = "")
      a1<-c("Junction Locale:",a1)
      
      legend("topleft",a1,cex=tcex,bty="n",inset=0.02) 
      
      
    }
    
    # bin info
    plot(1,axes=FALSE,type="n",xlab="",ylab="")
    if(iiss$b==1 | iiss$bjs==1){
      col1 <- c( "bin","bin.fdr","junction.dPIR","junction.dPSI")
      a1 <- sr@binbased[sr@binbased$bin%in%iiss$bin,col1]
      a1[col1[2:4]]<-signif(a1[col1[2:4]],2)
      names(a1)<-c("bin","bin.fdr","D_PIR","D_PSI")
      a1<-cbind(names(a1),t(a1))
      a1<-apply(a1,1,function(x){return(paste(x,collapse=": "))})
      a1<-paste("  ",a1,sep = "")
      a1<-c("Bin:",a1)
      legend("topleft",a1,cex=tcex,bty="n",inset=0.02) 
      
    }
  }
  
}



