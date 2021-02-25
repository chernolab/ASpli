.DUreport <- function( 
    counts, 
    targets, 
    minGenReads  = 10,
    minBinReads  = 5,
    minRds = 0.05,
    offset = FALSE,
    offsetAggregateMode = c( "geneMode", "binMode" )[1],
    offsetUseFitGeneX = TRUE,
    contrast = NULL,
    ignoreExternal = TRUE,
    ignoreIo = TRUE, 
    ignoreI = FALSE,
    filterWithContrasted = TRUE,
    verbose = FALSE,
    threshold = 5
    # ---------------------------------------------------------------------- #
    # Comment to disable priorcounts usage in bin normalization 
    # , priorCounts = 0 
    # ---------------------------------------------------------------------- #
  ) {
  
  # Create result object                   
  du <- new( Class="ASpliDU" )
  
  # Generate conditions combining experimental factors
  targets <- .condenseTargetsConditions( targets ) 
  
  # Filter genes and calculates differential usage of genes
  du <- .DUreportGenes( du, counts, targets, minGenReads, minRds, contrast, 
      filterWithContrasted, verbose )

  # Filter bins and calculates differential usage of bins 
  du <- .DUReportBins( du, 
                       counts, 
                       targets, 
                       minGenReads, 
                       minBinReads, 
                       minRds,
                       offsetAggregateMode, 
                       offsetUseFitGeneX, 
                       offset, 
                       contrast, 
                       ignoreExternal, 
                       ignoreIo, 
                       ignoreI,
                       filterWithContrasted,
                       verbose ) 
  # ------------------------------------------------------------------------ #
  
  #Adds junctionDUreport if offset is FALSE
  if(offset == FALSE) du <- .junctionDUreport(counts = counts, targets = targets, appendTo = du, 
                                              minGenReads = minGenReads, minRds = minRds, threshold = threshold,
                                              contrast = contrast)

  return( du )
}

.junctionDUreport <- function ( 
    counts, 
    targets, 
    appendTo = NULL, 
    minGenReads = 10,
    minRds = 0.05,
    threshold = 5,
    offset   = FALSE,
    offsetUseFitGeneX = TRUE,
    contrast = NULL,
    forceGLM = FALSE 
    # ------------------------------------------------------------------------ #
    # Comment to disable priorcounts usage in bin normalization 
    # , priorCounts = 0 
    # ------------------------------------------------------------------------ #
) {
  
  du <- if ( is.null( appendTo ) ) new( Class = "ASpliDU" ) else appendTo
  
  targets <- .condenseTargetsConditions( targets )
  
  df0 <- countsg(counts)
  
  dfG0 <- .filterByReads(
      df0 = df0,
      targets = targets,
      min = minGenReads,
      type = "all")
  
  dfGen <- .filterByRdGen(
      df0 = dfG0,
      targets = targets,
      min = minRds,
      type = "all" )
  
  df0 <- countsj(counts)[countsj(counts)[,"gene"]%in%rownames(dfGen),]
  
  df <- .filterJunctionBySample( df0 = df0, 
      targets = targets, 
      threshold = threshold)  #mean > one of the condition
  
  df <- df[ df$multipleHit == "-",]
  
  if( offset ){
    
    warning( simpleWarning( "Junctions DU with offsets is not fully tested. Use results with caution") )
    mOffset <- .getOffsetMatrix( 
        df, 
        dfGen,
        targets,
        offsetAggregateMode = 'geneMode',
        offsetUseFitGeneX= offsetUseFitGeneX )
  } else {
    mOffset <- NULL
  }

  junctionsdeSUM <- .junctionsDU_SUM(
      df = df,
      dfGen = dfGen,
      targets = targets,
      mOffset = mOffset,
      contrast = contrast,
      priorCounts = 0  
      # ------------------------------------------------------------------ # 
      # Comment to disable priorcounts usage in normalizefeatuebygen.         
#          ,priorCounts = priorCounts
      
  # ------------------------------------------------------------------ # 
  )
  
  
  du@junctions <- junctionsdeSUM
  message( "Junctions DU completed" )
  
  return( du )
}


.DUreportBinSplice <- function (  
    counts, 
    targets,
    minGenReads  = 10,
    minBinReads  = 5,
    minRds = 0.05,
    contrast = NULL,
    forceGLM = FALSE,
    ignoreExternal = TRUE, 
    ignoreIo = TRUE, 
    ignoreI = FALSE,
    filterWithContrasted = TRUE,
	  verbose = TRUE,
    formula = NULL,
    coef = NULL) {

  if(is.null(contrast) & is.null(formula)){
    stop("Must provide either contrast or formula.")
  }
  if(is.null(targets)){
    targets <- counts@targets
  }
  
  # Create result object                   
  du <- new( Class="ASpliDU" )
  
  #Contrast comes before formula
  if(!is.null(contrast)){
    du@contrast        <- contrast
    formula            <- NULL
  }else{
    design             <- model.matrix( formula, data = targets )
    contrast           <- ginv(design)
    if(is.null(coef)){
      coef <- nrow(contrast)
    }
    contrast           <- contrast[coef, ] #We add one because the first element is the intercept
    contrast[abs(contrast) < 1e-5] <- 0  
    contrastAggregated <- aggregate(contrast~targets$condition,FUN=sum)
    contrast <- setNames(contrastAggregated[,2],contrastAggregated[,1])[counts@condition.order]
    du@contrast        <- contrast    
  }
  if(!.hasSlot(counts, ".ASpliVersion")){
    counts@.ASpliVersion = "1" #Last version before 2.0.0 was 1.14.0. 
  }
  if(counts@.ASpliVersion == "1"){
    targets <- .condenseTargetsConditions( targets ) 
  }else{
    names(du@contrast) <- counts@condition.order
  }

  # Generate conditions combining experimental factors
  #targets <- .condenseTargetsConditions( targets ) 
  
  # Filter genes and calculates differential usage of genes
  du <- .DUreportGenes( du, counts, targets, minGenReads, minRds, contrast, 
       filterWithContrasted, verbose, formula, coef )
  message("Genes DE completed")

  # Filter bins and calculates differential usage of bins 
  du <- .DUReportBinsWithDiffSplice( counts, targets, contrast, du, minGenReads, 
      minBinReads, minRds, ignoreExternal, ignoreIo, ignoreI, filterWithContrasted, formula, coef )
  message("Bins DE completed")
  
  return( du )
  
}

.DUreportGenes <- function (
    du, 
    counts, 
    targets, 
    minGenReads  = 10,
    minRds = 0.05,
    contrast = NULL,
    filterWithContrasted = TRUE,
    verbose = FALSE,
    formula = NULL,
    coef    = NULL) {
  
  dfGen <- .filterGenes( counts, targets, minGenReads, minRds, contrast, 
      filterWithContrasted, verbose )
  
  genesde <- .genesDE( df=dfGen, targets = targets,
      contrast = contrast, verbose, formula, coef  )
  
  genesDE( du ) <- genesde
  
  return( du )
  
}

.DUReportBins <- function( 
    du, 
    counts, 
    targets, 
    minGenReads, 
    minBinReads, 
    minRds,
    offsetAggregateMode, 
    offsetUseFitGeneX, 
    offset, 
    contrast, 
    ignoreExternal, 
    ignoreIo, 
    ignoreI,
    filterWithContrasted = TRUE,
    verbose = FALSE) {
  
  # Filter bins
  filteringResult <- .filterBins(counts, targets, minGenReads, minBinReads, 
      minRds, ignoreIo, contrast, filterWithContrasted, verbose )
  dfGen <- filteringResult$genes
  df2   <- filteringResult$bins
  
  if ( verbose ) message( "Differential usage of bins:")
  
  # Set offset matrix if required
  if( offset ) {
    if ( verbose ) message( "  Using an offset matrix")
    mOffset <- .getOffsetMatrix(
        df2,
        dfGen,
        targets,
        offsetAggregateMode = offsetAggregateMode,
        offsetUseFitGeneX   = offsetUseFitGeneX,
        verbose)
    mOffset <- mOffset[ rownames( df2 ), ]
  } else {
    mOffset <- NULL
  }
  
  # Calculate bins DU
  binsdu <- .binsDU(
      df = df2,
      dfGen = dfGen,
      targets = targets,
      mOffset = mOffset,
      contrast = contrast,
      ignoreExternal = ignoreExternal,
      ignoreIo = ignoreIo, 
      ignoreI = ignoreI
      # -------------------------------------------------------------------- # 
      # Comment to disable priorcounts usage 
      , priorCounts = 0, 
      # , priorCounts = priorCounts 
      # -------------------------------------------------------------------- #
      verbose
  ) 
  
  binsDU( du ) <- binsdu
  if ( verbose ) message( "Differential usage of bins done")
  return( du )
}

.DUReportBinsWithDiffSplice <- function( counts, targets, contrast, du, 
    minGenReads, minBinReads, minRds, ignoreExternal, ignoreIo, ignoreI, 
	filterWithContrasted, formula, coef ) {
  
  # Filter bins
  countData <- .filterBins( counts, targets, minGenReads, minBinReads, minRds, 
      ignoreIo, contrast, filterWithContrasted )

  binsdu <- .binsDUWithDiffSplice( countData[[2]], targets, contrast, 
		  ignoreExternal, ignoreIo, ignoreI, formula, coef )[[1]]  
  #colnames(binsdu)[1]<-"locus"
  
  binsDU( du ) <- binsdu
  
  return( du )
}



.filterGenes <- function( counts, targets, minGenReads, minRds, contrast= NULL, 
    filterWithContrasted = TRUE, verbose = FALSE ) {
  
  if ( verbose ) { message( "Filtering genes:") }
  
  dfG0 <- .filterByReads( df0 = countsg( counts ), targets = targets,
      min = minGenReads, type = "any", contrast, 
      onlyContrast = filterWithContrasted, verbose )
  
  dfGen <- .filterByRdGen( df0 = dfG0, targets = targets,
      min = minRds, type = "any", contrast, 
      onlyContrast = filterWithContrasted, verbose ) 
  
  if ( verbose ) { message( "Filtering genes done") }
  
  return( dfGen )
  
} 

.filterBins <- function( counts, targets, minGenReads, minBinReads, minRds, 
    ignoreIo, contrast, filterWithContrasted, verbose = FALSE ) {
  
  if ( verbose ) { message( "Filtering bins") }

  if ( verbose ) { message( "  Selection of expressed genes") }
  
  dfG0 <- .filterByReads( df0 = countsg( counts ), targets = targets,
      min = minGenReads, type = "all" , contrast, 
      filterWithContrasted , verbose )
  
  dfGen <- .filterByRdGen( df0 = dfG0, targets = targets, 
      min = minRds, type = "all" , contrast, 
      filterWithContrasted, verbose )
  
  if ( verbose ) { message( "  Selection of expressed genes done") }
  
  dfBin <- countsb(counts)[countsb(counts)[,"locus"]%in%row.names(dfGen),]
  
  if( ignoreIo ) dfBin <- dfBin[dfBin[,"feature"]!="Io",]
  
  if ( verbose ) { message( "  Selection of bins of expressed genes") }
  
  df1 <- .filterByReads( df0=dfBin, targets=targets, 
      min=minBinReads, type="any", contrast, 
      filterWithContrasted, verbose )
  
  df2 <- .filterByRdBinRATIO( dfBin=df1, dfGen=dfGen,
      targets=targets, min=minRds, type="any", contrast , 
      filterWithContrasted, verbose )
  
  if ( verbose ) { message( "Filtering bins done") }
  
  return( list( 'genes'=dfGen, 'bins' = df2 ) )
  
}


.filterByReads <- function( df0, targets, min, type, contrast = NULL, 
    onlyContrast = FALSE, verbose = FALSE ) {

  # subset a working data frame
  cropped <- .extractCountColumns( df0, targets )
  if ( verbose ) { message( "  Filtering by reads.") }
  
  return( .filterByGeneric( cropped, df0, targets, min, type, contrast,
          onlyContrast, verbose  ) )

}


.filterByRdGen <- function( df0, targets, min, type, contrast = NULL, 
    onlyContrast = FALSE, verbose = FALSE ) {
  # keeps those genes which avg rd > min in any condition
  cropped <- .extractCountColumns( df0, targets ) / df0$effective_length
  if ( verbose ) { message( "  Filtering by read density.") }
  
  return( .filterByGeneric( cropped, df0, targets, min, type, contrast,
          onlyContrast, verbose ) )
}

.filterByGeneric <- function( cropped, df0, targets, min, type, contrast = NULL, 
    onlyContrast = FALSE, verbose = FALSE ) {
  
  # Get contrasted conditions
  contrastedConditions <- .getFilteringConditions( targets, contrast, 
      onlyContrast )
  if (verbose) message( paste("  Filtering using", paste( contrastedConditions, 
                collapse = ",") ,"conditions") )
  
  # TODO: vectorize this code
  list <- matrix( unlist(
          lapply( contrastedConditions , 
              function( x ) {
                rowMeans( cropped[ , targets$condition == x , drop = FALSE] ) >= min } )),  
      nrow = nrow( cropped ), 
      byrow = FALSE )
  
  
  if (type=="any")  {
    if (verbose) message( paste("  Filtering any condition with mean minimum value",min) )    
    filteredDF <- df0[ rowSums( list ) > 0             , ]
  } else  {
    if (verbose) message( paste("  Filtering all conditions with mean minimum value",min) )    
    filteredDF <- df0[ rowSums( list ) == length( contrastedConditions ) , ]
  }
  
  return (filteredDF)
}


.filterByRdBinRATIO <- function( 
    dfGen,
    dfBin,
    targets, 
    min, 
    type,
    contrast = NULL,
    onlyContrast = FALSE,
    verbose = FALSE) {

  # -------------------------------------------------------------------------- #
  # Get densities for genes and bins
  genes.rd <- .extractCountColumns( dfGen, targets ) / dfGen$effective_length 
  bins.rd  <- .extractCountColumns( dfBin, targets ) / dfBin$length
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Get contrasted conditions
  contrastedConditions <- .getFilteringConditions( targets, contrast, 
      onlyContrast )
  if (verbose) message( paste("  Filtering using", paste( contrastedConditions, 
                collapse = ",") ,"conditions") )
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Calculates avg. bin density by condition
  # TODO: vectorize this code
  avRdBin <- matrix( unlist(
          lapply( contrastedConditions , 
              function( x ) { 
                rowMeans( bins.rd[ , targets$condition == x, drop = FALSE] ) 
              } ) ),  
      nrow = nrow( bins.rd ), 
      byrow = FALSE) 
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Calculates avg. gene density by condition
  # TODO: vectorize this code
  avRdGen <- matrix( unlist(
          lapply( contrastedConditions, 
              function(x) rowMeans(genes.rd[ , targets$condition == x, drop = FALSE]))),  
      nrow = nrow(genes.rd), 
      byrow = FALSE)#Ok
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Calculates bin to gene read density ratio
  gen.rdb <- avRdGen[ match( dfBin$locus, rownames( dfGen ) )  , ]
  bin.gen.rd <-  avRdBin / gen.rdb
  bin.gen.rd[ is.na( bin.gen.rd ) ] <- 0 
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Apply filter
  if ( type == "any" )  { 
    if (verbose) message( paste("  Filtering any condition with mean minimum value",min) )    
    dfBin <- dfBin[ rowSums( bin.gen.rd >= min ) > 0 ,]
  } else{ 
    if (verbose) message( paste("  Filtering all conditions with mean minimum value",min) )    
    dfBin <- dfBin[ rowSums( bin.gen.rd >= min ) == length( contrastedConditions ) ,]
  }
  # -------------------------------------------------------------------------- #
  
return (dfBin)
}

.getFilteringConditions <- function ( targets, contrast, onlyContrast ) {

  allConditions <- getConditions( targets )
  
  constrast <- if ( is.null( contrast) ) .getDefaultContrasts( targets$condition )
  
  if ( onlyContrast && ( ! is.null( contrast  ) ) ) {
    return( allConditions[ (contrast != 0) ] ) 
  } else {
    return( allConditions )
  }   
}

# ---------------------------------------------------------------------------- #
.normalizeByGenFeature <- function( feature, gene, targets, priorCounts = 0,
    verbose = FALSE) {

  if ( verbose ) message( "  Using bin/gene normalization")
  
  # Search the column with the gene name in the features 
  colLocus <- match( c( "gene", "locus" ), colnames( feature ) )
  colLocus <- colLocus[ ! is.na( colLocus ) ][1]
  
  # Repeat gene rows by the number of bins in that gene 
  f.index <- match( feature[ , colLocus ], rownames(gene) ) 
  counts.genes.f <- gene[f.index,]

  # Extract count data
  counts.genes.f  <- .extractCountColumns( counts.genes.f, targets)
  counts.features <- .extractCountColumns( feature, targets) 
  
  # Calculate mean of gene counts 
  gen.mean.f <- rowMeans( counts.genes.f )   
  
  # Apply normalization
  normalizedCounts <- round( ( as.matrix(  counts.features ) + priorCounts ) / 
                       ( as.matrix( counts.genes.f ) + priorCounts ) * 
                       ( gen.mean.f + priorCounts ) )
  
  # Set invalid results to zero
  normalizedCounts[ is.na( normalizedCounts ) | is.infinite( normalizedCounts) ] <- 0 
  
  # Create result dataframe
  normalizedFull <- cbind( .extractDataColumns( feature , targets ) , normalizedCounts )

  return( normalizedFull )
}
# ---------------------------------------------------------------------------- #

.getDefaultContrasts <- function ( conditions ) {
  contrast <- rep( 0, length( unique( conditions ) ) )
  contrast[1:2] <- c(-1,1)
  return( contrast )
}

.genesDE <- function( df, targets, contrast = NULL,
    verbose = FALSE, formula = NULL, coef = NULL ) { 
  
  if (verbose) message("Genes differential expression:")
  
  if( is.null( contrast ) ) contrast <- .getDefaultContrasts(targets$condition)
  
  if ( verbose ) {
    msg <- paste(mapply( paste0,contrast[contrast!=0],
            getConditions(targets)[contrast!=0 ]),collpase = " ") 
    message("  Contrast:",msg)
  }
  
  cols <- match( rownames( targets ), colnames( df ) )
  
  group <- targets$condition
  
  groupFactor <- factor( group, unique( group ), ordered = TRUE )
  
  er <- DGEList( counts = df[ , cols ], samples=targets, group=groupFactor)
  er <- calcNormFactors( er )
  
  justTwoConditions <- sum( contrast != 0 ) == 2
  
  if( justTwoConditions & !  is.null(formula)){
    if (verbose) message("  Running exact test")
    er   <- estimateDisp( er , robust=TRUE)
    capture.output( er   <- estimateDisp( er, robust=TRUE ) )
    pair <- which( contrast != 0 )
    
    # at this point pair must have just two elements.
    # the first element of pair is assumed to be the control. 
    # Therefore, pair must be adjusted according to contrast 
    pair <- if (contrast[pair][1] == -1) pair else rev(pair)
    et   <- exactTest(er, pair=pair)
  } else {
    if (verbose) message("  Running GLM LRT")
    if(is.null(formula)){
      design   <- model.matrix( ~0 + groupFactor, data = er$samples )
      captured <- capture.output(
          er     <- estimateDisp( er, design = design , robust=TRUE)
      )
      glf    <- glmFit( er, design = design)  
      et     <- glmLRT( glf, contrast = contrast)
    }else{
      if(verbose) message("  Running with formula")
      design   <- model.matrix( formula, data = er$samples )
      captured <- capture.output(
        er     <- estimateDisp( er, design = design , robust=TRUE)
      )
      glf    <- glmFit( er, design = design)  
      et     <- glmLRT( glf, coef = coef )
    }
  } 
  
  fdr.gen <- p.adjust( et$table$PValue, method="BH" )
  
  cols <- match(rownames( targets ), colnames( df ) )
  geneData <- .extractDataColumns( df, targets )
  genesFull <- data.frame( geneData ,
      logFC = as.numeric( et$table$logFC ), 
      pvalue = as.numeric( et$table$PValue ), 
      gen.fdr = as.numeric(fdr.gen), 
      stringsAsFactors = FALSE)
  rownames( genesFull ) <- rownames( df )
  
  if (verbose) message("Genes differential expression done")
  
  return( genesFull )
}


.binsDU <- function (
      df,
      dfGen,
      targets,                   
      ignoreExternal= TRUE, 
      ignoreIo = TRUE, 
      ignoreI = FALSE,
      contrast = NULL,

      mOffset  = NULL,
      priorCounts = 0,
      verbose = TRUE) {

  # -------------------------------------------------------------------------- #
  # Filtrar bins de acuerdo de diferentes criterios
  df = df[ ! ignoreExternal | df$event != "external" ,] 
  df = df[ ! ignoreIo | df$feature != "Io" ,] 
  df = df[ ! ignoreI | df$feature != "I" ,] 
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Normalize bins by gen or filter offsets
  if( is.null( mOffset ) ){
    df <- .normalizeByGenFeature( feature=df, gene=dfGen, targets = targets, 
        priorCounts = priorCounts, verbose )
  } else {
    mOffset <- mOffset[ rownames( df ), ]
  }
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # perform edgeR extact or glm test 
  et <- .edgeRtest( df, dfGen, targets, mOffset, contrast, verbose )
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Build result dataframe
  bin.fdr <- p.adjust( et$table$PValue, method="BH" )
  logFC   <- et$table$logFC
  pvalue  <- et$table$PValue
  
  cols <- match( rownames( targets ), colnames( df ) )

  splicing_full <- data.frame ( df[ , -cols ],
      logFC, 
      pvalue,
      bin.fdr,
      stringsAsFactors = FALSE )
  
  rownames( splicing_full ) <- rownames( df )
  
  return( splicing_full )
  # -------------------------------------------------------------------------- #
  
}

.binsDUWithDiffSplice <- function( countData, targets, contrast, 
    ignoreExternal = FALSE, ignoreIo = TRUE, ignoreI = FALSE, formula = NULL, coef = NULL){
  
  # Filter bins
  countData = countData[ ! ignoreExternal | countData$event != "external" ,] 
  countData = countData[ ! ignoreIo | countData$feature != "Io" ,] 
  countData = countData[ ! ignoreI | countData$feature != "I" ,] 
  
  # Define group and contrast
  #targets$condition <- rep(c("adjacent_non_tumor_tissue", "tumor_tissue"), each=3)
  group <- targets$condition
  if( is.null( contrast ) ) contrast <- .getDefaultContrasts( group )
  
  # make DU analysis
  y <- DGEList( counts = .extractCountColumns( countData, targets ),
      group = factor( group, levels = getConditions(targets), ordered = TRUE) ,
      genes = .extractDataColumns(countData, targets) )       

  # TODO: Este filtro es muy resctrictivos
  #  keep <- rowSums( cpm( y ) > 1) >= 2
  #  y <- y[ keep, , keep.lib.sizes = FALSE ]
	y <- calcNormFactors( y )
  
  
  # model.matrix sort columns alphabetically if formula has characters instead 
  # of factors. Therefore to preserve order group is converted to ordered factors 
  groupFactor <- factor( group, unique( group ), ordered = TRUE )
  
  if(is.null(formula)){
    design <- model.matrix( ~0 + groupFactor, data = y$samples )
    
    y   <- estimateDisp( y, design, robust=TRUE )
    fit <- glmFit( y, design )
    ds  <- diffSpliceDGE( fit, contrast = contrast, geneid = "locus", exonid = NULL, verbose = FALSE )
  }else{
    message("Running with formula")
    y$samples <- targets
    design <- model.matrix( formula, data = y$samples )
    
    y   <- estimateDisp( y, design , robust=TRUE)
    fit <- glmFit( y, design )
    ds  <- diffSpliceDGE( fit, coef = coef, geneid = "locus", exonid = NULL, verbose = FALSE )    
  }
  tsp <- topSpliceDGE( ds, test = "exon", FDR = 2, number = Inf )
  tspg<- topSpliceDGE( ds, test = "gene", FDR = 2, number = Inf )
  
  # make column names equal to the results of DUReport method
  colnames( tsp )[ match( 'FDR', colnames( tsp )) ] <- 'bin.fdr'
  colnames( tsp )[ match( 'P.Value', colnames( tsp )) ] <- 'pvalue'
  tsp$exon.LR <- NULL
  
  #tsp <- tsp[,c("locus","logFC","pvalue","bin.fdr")]
  #colnames(tsp)[1]<-"J3"
  
  rownames(tspg) <- tspg$locus
  tspg <- tspg[,c("gene.LR","P.Value","FDR")]
  colnames(tspg)[1:3] <- c("cluster.LR","cluster.pvalue","cluster.fdr")
  
  return( list(junction=tsp,cluster=tspg) )
  
  
 # return( tsp )
  
} 

.junctionsDU_SUM <- function( df, 
                              dfGen,
                              targets, 
                              mOffset = NULL,
                              contrast = NULL,
                              priorCounts = 0 ) {

  # -------------------------------------------------------------------------- #
  # Inner function to compute junction ratio
  jratio <- function( junctions ) {
    junctions[ is.na( junctions ) ] <- 0 
    return( junctions[1] / ( junctions[1] + junctions[2] ) )
  }
  # -------------------------------------------------------------------------- #
  
  if( is.null( mOffset ) ) {
    df <- .normalizeByGenFeature( feature=df, gene=dfGen, targets, priorCounts )
  }
  
  et <- .edgeRtest( df, dfGen, targets, mOffset, contrast)
  
  fdr <- p.adjust( et$table$PValue, method="BH" )
  logFC <- et$table$logFC
  pvalue <- et$table$PValue
  
  group<- targets$condition
  
  jranges <- .createGRangesExpJunctions( rownames( df ) )

# TODO:  Es Necesario seguir preguntando por una version vieja de IRanges ?
# Ademas, las dos llamadas a la version vieja llaman diferente al parametro
# ignore.redundant ( en la otra llamada es ignoreRedundant )
#  if ( packageVersion("IRanges") < 2.6) {
#    j.start <- findOverlaps( jranges, ignoreSelf=TRUE, ignore.redundant=FALSE,
#        type="start")  
#  } else {
    j.start <- findOverlaps( jranges, drop.self=TRUE, drop.redundant=FALSE,
        type="start")
#  }
  
  jjstart <- as.data.frame( j.start )
  
  jjstart$queryHits <- names(jranges[jjstart$queryHits])
  jjstart$subjectHits <- names(jranges[jjstart$subjectHits])
  
  
  # BUG FIX: aggregate fails when no junction overlaps 
  if ( length( j.start ) > 0) {
    shareStart <- data.frame(aggregate(subjectHits ~ queryHits, 
            data = jjstart, paste, collapse=";"))
  } else {
    shareStart <- data.frame( 
        subjectHits = character(0), 
        queryHits = character(0))
  }
  
  start <- ncol( df ) - nrow(targets) + 1 
  end   <- ncol( df ) 

  dfCountsStart <- data.frame( names=jjstart$queryHits, 
      df[ jjstart$subjectHits, start:end],
      row.names=NULL) 
  
  # BUG FIX: aggregate fails with 0-rows dfCountsStart. 
  if ( nrow( dfCountsStart  ) > 0 ) {
    dfSumStart <- data.frame(aggregate(. ~ names, data = dfCountsStart, sum))
  } else {
    dfSumStart <- data.frame( names = character(0) )
    for ( i in 1:ncol( dfCountsStart ) ) {
      dfSumStart[, i+1] <- integer(0)
    }
    colnames( dfSumStart )[2:ncol(dfSumStart)] <- colnames( dfCountsStart )
  }
  
  
  sumJ <- paste(colnames(dfSumStart), "jsum", sep=".")
  colnames(dfSumStart) <-  sumJ
  rownames(dfSumStart) <- dfSumStart$names.jsum
  
  dfSumStart$names.jsum <- NULL

  dffStart <- data.frame(matrix(NA, nrow =  nrow(df), ncol = ncol(dfSumStart)) )
  rownames(dffStart) <- rownames(df)
  colnames(dffStart) <- colnames(dfSumStart)  
  mSumStart <- match( row.names(dfSumStart), row.names(dffStart)) 
  
  dffStart[mSumStart,] <-  dfSumStart
  dffStart[is.na(dffStart)] <- 0
  
  mbin_start_hit <- match(shareStart$queryHits, row.names(dffStart))
  #aca reacomodo el bin_start_hit con el indexado de dffStart
  bin_start_hit <- as.character(rep("-", nrow(dffStart)) )
  bin_start_hit[mbin_start_hit] <- shareStart$subjectHits
  ################################################
  
  cols <- match( rownames( targets ),colnames(df))
  
  
  ratioStart <- data.frame( df[,cols],dffStart)
  colnames(ratioStart) <- rep(rownames( targets ),2)
  #aca hay que armar un df itnermedio con la suma por condicion:

 

  ff <- rep(group,2)
  colnames(ratioStart) <- paste(ff, rep(1:2,each=length(group)))
  dfSum <- t(apply(ratioStart, 1, function(x){tapply(as.numeric(x), 
                INDEX=colnames(ratioStart), 
                sum)}))
  colnames(dfSum) <- rep(unique(group), each=2) # old version
  
  jratioStartRes <- t(apply(dfSum, 1, function(x){tapply(as.numeric(x), 
                INDEX=colnames(dfSum), 
                jratio )}))

#  if (packageVersion("IRanges")<2.6) { 
#    j.end <- findOverlaps(jranges, 
#        ignoreSelf=TRUE,
#        ignoreRedundant=FALSE,
#        type="end")
#    
#  } else {
    j.end <- findOverlaps( jranges, 
        drop.self = TRUE,
        drop.redundant = FALSE,
        type = "end" )
#  }

  jjend <- as.data.frame(j.end)
  jjend$queryHits <- names(jranges[jjend$queryHits])
  jjend$subjectHits <- names(jranges[jjend$subjectHits])
  
  
  # BUG FIX: aggregate fails when no junction overlaps 
  if ( length( j.end ) > 0) {
    shareEnd <- data.frame(aggregate(subjectHits ~ queryHits, 
            data = jjend, paste, collapse=";"))
  } else {
    shareEnd <- data.frame( 
        subjectHits = character(0), 
        queryHits = character(0))
  }

  
  
  dfCountsEnd <- data.frame( names=jjend$queryHits, 
      df[jjend$subjectHits,start:end],
      row.names=NULL) #recover counts
  

  # BUG FIX: aggregate fails with 0-rows dfCountsStart. 
  if ( nrow( dfCountsEnd  ) > 0 ) {
    dfSumEnd <- data.frame(aggregate(. ~ names, data = dfCountsEnd, sum))
  } else {
    dfSumEnd <- data.frame( names = character(0) )
    for ( i in 1:ncol( dfCountsEnd ) ) {
      dfSumEnd[, i+1] <- integer(0)
    }
    colnames( dfSumEnd )[2:ncol(dfSumEnd)] <- colnames( dfCountsEnd )
  }
  
  sumJ <- paste(colnames(dfSumEnd), "jsum", sep=".")
  colnames(dfSumEnd) <- sumJ
  rownames(dfSumEnd) <- dfSumEnd$names.jsum
  dfSumEnd$names.jsum <- NULL
  dffEnd =data.frame(matrix(NA, nrow =  nrow(df), 
          ncol = ncol(dfSumEnd)) )
  rownames(dffEnd) <- rownames(df)
  colnames(dffEnd) <- colnames(dfSumEnd)
  ########################################################################
  mSumEnd <- match(row.names(dfSumEnd), row.names(dffEnd))
  dffEnd[mSumEnd,] <- dfSumEnd
  dffEnd[is.na(dffEnd)] <- 0
  ########################################################################
  mbin_end_hit <- match(shareEnd$queryHits, row.names(dffEnd))
  bin_end_hit <- rep("-", nrow(dffEnd))
  bin_end_hit[mbin_end_hit] <- shareEnd$subjectHits
  ########################################################################
  
  ratioEnd <- data.frame(df[,cols],dffEnd)
  ff <- rep(group,2)
  colnames(ratioEnd) <- paste(ff,rep(1:2, each = length(group)))
  dfSum <- t(apply(ratioEnd, 1, function(x){tapply(as.numeric(x), 
                INDEX=colnames(ratioEnd), sum  )}))

  colnames(dfSum) <- rep(unique(group), each = 2) # new version
  
  jratioEndRes <- t(apply(dfSum, 1, function(x){tapply(as.numeric(x), 
                INDEX=colnames(dfSum), jratio )}))

  
  et_merge <- data.frame(df,                       
      logFC,
      pvalue, 
      fdr,
      bin_start_hit,
      dffStart,
      jratioStartRes,
      bin_end_hit,
      dffEnd,
      jratioEndRes )
  return(et_merge)
  
}

# ---------------------------------------------------------------------------- #
# Calculates the matrix of offset values for normalization
# TODO: Que es el valor 10^-4 que aparece, se puede pasar como parametro?
.getOffsetMatrix <- function( df, dfGen, targets, 
    offsetAggregateMode = c( "geneMode","binMode" )[1], 
    offsetUseFitGeneX = TRUE, verbose = FALSE) {
  
  locus  <- df[, na.omit( match( c("locus","gene"), colnames(df) ) ) ]
  
  if( offsetAggregateMode=="geneMode" ) {
    if(offsetUseFitGeneX){
      if ( verbose ) message( "  Using geneMode for offset matrix with GLM fit")
      countData = dfGen[,rownames(targets)]
      
      yg <- DGEList( counts = countData, 
                     group  = targets$condition, 
                     genes  = data.frame( locus = rownames(countData) ) )
      
      #filter lowcount genes
      # TODO: Los datos que llegan aca ya estan filtrados por low counts usando
      # los filtros de ASpli. Es necesario incluir este filtro ?
#      keep <- rowSums(cpm(yg)>1) >= 2
#      yg <- yg[keep, , keep.lib.sizes=FALSE]
      
      yg     <- calcNormFactors(yg)
      fc     <- targets$condition
      groupFactor <- factor( fc, unique( fc ), ordered = TRUE )
      design <- model.matrix(~0+groupFactor)
      yg     <- estimateDisp(yg,design, robust=TRUE)
      fitg   <- glmFit(yg,design)
      maux   <- fitg$fitted.values + 10^-4
      rownames( maux ) <- rownames( dfGen)
      
    } else {
      if ( verbose ) message( "  Using geneMode for offset matrix with gene counts")     
      maux   <- as.matrix( dfGen[,rownames(targets)] ) + 10^-4
      rownames( maux ) <- rownames( dfGen )
    }
  } else {
    # TODO:  Es necesario vectorizar esto ?
    if ( verbose ) message( "  Using binMode for offset matrix")     
    
    if ( "junction" %in% colnames( df ) ) {
      stop( simpleError( "Differential usage with offset and 'binMode' is not available." ))
    }
    
    a <- by( data = df[,c("feature",rownames(targets))],
             INDICES = as.factor( locus ),
             FUN = function(x){ 
               apply( x [x[,1] == "E", 2:ncol(x) ], 2 , sum )
             })
    maux<- matrix(unlist(a),byrow=TRUE,ncol=nrow(targets))+10^-4
    rownames(maux)<-names(a)
    colnames(maux)<-names(a[[1]])
  }
  #mOffset <- do.call( rbind, list( A = maux[ locus ,] ) )
  mOffset <- maux[ locus ,]
  rownames( mOffset )<-rownames(df)
  return(mOffset)
}

.setDefaultOffsets <- function ( aDGEList , mOffset) {

   aDGEList$offset <- log(mOffset) 
   return ( aDGEList )
}

.edgeRtest <- function( 
    df,
    dfGen,
    targets,
    mOffset = NULL,
    contrast = NULL,
    verbose = FALSE) {
  
  cols <- match( rownames( targets ), colnames( df ) )
  
  group <- targets$condition
  
  if( is.null( contrast ) ) contrast <- .getDefaultContrasts(group)
  
  er <- DGEList( counts  = df[,cols],
                 samples = targets,
                 group   = factor( group, levels = unique(group), ordered = TRUE ) )
  
  er <- calcNormFactors(er)
  
  justTwoConditions <- sum( contrast != 0 ) == 2
  
  # TODO: Forzar GLM no tiene efecto si se pasa un offset. 
  if( is.null( mOffset ) & justTwoConditions ){

    captured <- capture.output(
      er   <- estimateDisp( er, robust=TRUE )
    )
    pair <- which( contrast != 0 )
    testResult   <- exactTest( er, pair = pair )
    
#    message( "ExactTest... done" )
    
  } else {
    if (verbose) message("Running GLM LRT")

    if( ! is.null( mOffset ) ) er <- .setDefaultOffsets( er, mOffset ) 
    
    groupFactor <- factor( group, unique( group ), ordered = TRUE )
    design     <- model.matrix( ~0 + groupFactor, data = er$samples )
    captured <- capture.output(
        er         <- estimateDisp( er, design = design , robust=TRUE) 
    )      
    glf        <- glmFit( er, design = design )  
    testResult <- glmLRT( glf, contrast = contrast )
    
#    message( "glmLRT... done" )

  } 
  return( testResult )
}


