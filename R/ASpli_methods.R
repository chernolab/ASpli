# ---------------------------------------------------------------------------- #
# Class definitions
setClass( Class = "ASpliFeatures",
          representation = representation(
            genes = "GRangesList",
            bins = "GRanges",
            junctions = "GRanges",
            transcriptExons = "GRangesList"))

setClass( Class = "ASpliCounts",
          representation = representation(
            gene.counts = "data.frame", 
            exon.intron.counts = "data.frame",
            junction.counts = "data.frame",
            e1i.counts = "data.frame", 
            ie2.counts = "data.frame",
            gene.rd = "data.frame",
            bin.rd = "data.frame", 
            targets = "data.frame",
            condition.order = "character",
            .ASpliVersion = "character"))

setClass( Class="ASpliAS",
          representation = representation(
            irPIR = "data.frame",
            altPSI = "data.frame",
            esPSI = "data.frame",
            junctionsPIR = "data.frame",
            junctionsPJU = "data.frame",
            join = "data.frame", 
            targets = "data.frame",
            .ASpliVersion = "character") )

setClass( Class = "ASpliDU",
          representation = representation(
            genes = "data.frame",
            bins = "data.frame",
            junctions = "data.frame",
            contrast = "numeric",
            .ASpliVersion = "character" ))

setClass( Class = "ASpliJDU",
          representation = representation(
            localec = "data.frame",
            localej = "data.frame",
            anchorc = "data.frame",
            anchorj = "data.frame",
            jir     = "data.frame",
            jes     = "data.frame",
            jalt    = "data.frame",
            contrast= "numeric",
            .ASpliVersion = "character"))

setClass( Class = "ASpliSplicingReport",
          representation = representation(
            binbased    = "data.frame",
            localebased = "data.frame",
            anchorbased = "data.frame",
            contrast    = "numeric",
            .ASpliVersion = "character"))

setClass( Class = "ASpliIntegratedSignals",
          representation = representation(
            signals      = "data.frame",
            filters      = "data.frame",
            .ASpliVersion = "character"))
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Set methods
setGeneric ( name = "binGenome", 
             def = function( genome, geneSymbols = NULL, 
                             logTo = "ASpli_binFeatures.log" ) standardGeneric( "binGenome" ) )

setMethod(
  f = "binGenome",
  signature = "TxDb",
  definition = function ( genome, geneSymbols = NULL, logTo = "ASpli_binFeatures.log") {
    
    #Normalize seqnames. If . present in name, changes it to _ and warns the user
    if(length(grep("[.]", seqlevels(genome)) > 0)){
      seqlevels(genome) <- gsub("[.]", "_", seqlevels(genome))
      warning("Some seqnames had a '.' present in their names. ASpli had to normalize them using '_'.")
    }
    
    features <- new( Class = "ASpliFeatures" )
    
    if ( is.null( geneSymbols ) ) {
      # Recupera los nombres de los genes
      geneSymbols <- data.frame( names( transcriptsBy(genome) ), stringsAsFactors = FALSE)
      row.names(geneSymbols) <- names( transcriptsBy(genome) )
      colnames(geneSymbols)  <- "symbol" 
    }
    if ( ! is.null ( logTo ) ) {
      sink( file = logTo )
    }
    genes.by.exons <- .createGRangesGenes( genome, geneSymbols ) 
    
    msg <- paste( "* Number of extracted Genes =" , length( genes.by.exons ) )
    message( msg )
    if ( sink.number() > 0 ) cat( paste( msg, "\n" ) )
    
    
    exon.bins <- .createGRangesExons(genome, geneSymbols)
    #add locus_overlap
    index <- match(exon.bins@elementMetadata$locus, names(genes.by.exons))
    locus_overlap <- rep("-", length(exon.bins))
    locus_overlap <- genes.by.exons@elementMetadata$locus_overlap[index]
    mcols(exon.bins) <- append(mcols(exon.bins), DataFrame(locus_overlap=locus_overlap))
    
    msg <- paste( "* Number of extracted Exon Bins =", length( exon.bins ) )
    message( msg )
    if ( sink.number() > 0 ) cat( msg, "\n" ) 
    
    
    intron.tot <- .createGRangesIntrons(genome, geneSymbols)
    #add locus_overlap
    index <- match(intron.tot@elementMetadata$locus, names(genes.by.exons))
    locus_overlap <- rep("-", length(intron.tot))
    locus_overlap <- genes.by.exons@elementMetadata$locus_overlap[index]
    mcols(intron.tot) <- append( mcols(intron.tot), 
                                 DataFrame(locus_overlap=locus_overlap) )
    
    msg <- paste( "* Number of extracted intron bins =", length( intron.tot ) )
    message( msg )
    if ( sink.number() > 0 ) cat( paste( msg, "\n" ) )
    
    transcripts <- .createGRangesTranscripts(genome)
    
    msg <- paste( "* Number of extracted trascripts =" , length( unlist( transcripts ) ) )
    message( msg )
    if ( sink.number() > 0 ) cat( paste( msg,"\n") )
    
    junctions <- .createGRangesJunctions( genome ) 
    
    #add locus_overlap
    index <- match( junctions@elementMetadata$locus, names( genes.by.exons ) )
    locus_overlap <- rep( "-", length( junctions ) )
    locus_overlap <- genes.by.exons@elementMetadata$locus_overlap[ index ]
    mcols(junctions) <- append( mcols( junctions ), 
                                DataFrame( locus_overlap = locus_overlap ) )
    
    msg <- paste("* Number of extracted junctions =", length(junctions) )
    message( msg)
    if ( sink.number() > 0 ) cat( paste( msg, "\n" ) )
    
    intron.bins <- intron.tot[ intron.tot@elementMetadata$feature == "I"  ]
    intron.orig <- intron.tot[ intron.tot@elementMetadata$feature == "Io" ]
    
    class  <- rep("fullI", length( intron.orig ))
    mcols( intron.orig ) <- append( mcols( intron.orig ), DataFrame( class=class ))
    event  <- rep("-", length( intron.orig ))
    eventJ <- rep("-", length( intron.orig ))
    mcols( intron.orig ) <- append( mcols( intron.orig ), DataFrame( event=event ))
    mcols( intron.orig ) <- append( mcols( intron.orig ), DataFrame( eventJ=eventJ ))
    
    exons <- exons( genome )
    
    exons.introns <- .findAsBin( exons, exon.bins, intron.bins, transcripts, 
                                 junctions, logTo )
    fullT <- c(exons.introns,intron.orig)
    
    transcriptExons <- exonsBy(genome, use.names=TRUE)

    tByGene    <- transcriptsBy(genome, by = "gene")
    
    a         <- unlist(tByGene)
    #geneNames <- DataFrame(gene=names(a),tx=mcols(a)$tx_name)
    geneNames <- DataFrame(gene=names(a))
    rownames(geneNames)<-mcols(a)$tx_name
    mcols( transcriptExons ) <- append ( mcols( transcriptExons ), geneNames[names(transcriptExons),,drop=FALSE] )

    features@genes <- genes.by.exons
    features@bins <- fullT
    features@junctions <- junctions
    features@transcriptExons <- transcriptExons
    
    
    return( features ) 
  })


setGeneric( name = "rds",
            def = function( counts, targets ) standardGeneric("rds") )

# TODO: Las densidades de reads de genes y bins se calculan dos veces. Una 
# vez aca y otra vez cuando se hace el filtrado para hacer DU.
setMethod(
  f= "rds",
  signature = "ASpliCounts",
  definition = function(counts, targets) {
    geneStart <- ncol(countsg(counts))-nrow(targets)+1
    gene.rd <- cbind( countsg(counts)[,1:geneStart-1], 
                      countsg(counts)[,geneStart:ncol(countsg(counts))] / 
                        countsg(counts)$effective_length
    )
    binStart <- ncol(countsb(counts))-nrow(targets)+1
    bin.rd <- cbind(countsb(counts)[, 1:binStart-1], 
                    countsb(counts)[,binStart:ncol(countsb(counts))]
                    /countsb(counts)$length)
    
    tb <- match(bin.rd$locus, rownames(gene.rd))
    rdfinalb=cbind(bin.rd, bin.rd[,binStart:ncol(bin.rd)]
                   /gene.rd[tb,geneStart:ncol(countsg(counts))])
    
    counts@gene.rd <- gene.rd
    counts@bin.rd <- rdfinalb
    # counts@targets <- .condenseTargetsConditions(targets)
    # group                  <- count@targets$condition
    # counts@condition.order <- levels(factor( group, unique( group ), ordered = TRUE ))
    return(counts)
  }
)

# gbCounts is a wrapper around readCounts for improved legibility
setGeneric (
  name = "gbCounts",
  def = function( features, 
                  targets, minReadLength, maxISize, 
                  minAnchor = 10,
                  libType="SE",
                  strandMode=0)
    standardGeneric("gbCounts") )

setMethod(
  f = "gbCounts",
  signature = "ASpliFeatures",
  definition = function( features, targets,  minReadLength,  
                         maxISize, minAnchor = 10,
                         libType="SE",
                         strandMode=0) {
    counts <- readCounts( features = features,
                          bam = NULL, 
                          targets = targets, 
                          readLength = minReadLength, 
                          maxISize = maxISize, 
                          minAnchor = minAnchor,
                          libType=libType,
                          strandMode=strandMode)
    counts@.ASpliVersion = "2" #Marks ASpliCounts object with the ASpli update 2.0.0
    return(counts)
  }
)

# readCounts
setGeneric (
  name = "readCounts",
  def = function( features, 
                  bam,
                  targets, 
                  cores = 1, 
                  readLength, 
                  maxISize, 
                  minAnchor = 10,
                  libType=libType,
                  strandMode=strandMode
                  )
    standardGeneric("readCounts") )

setMethod(
  f = "readCounts",
  signature = "ASpliFeatures",
  definition = function( features, bam, targets, cores = 1, readLength,  
                         maxISize, 
                         minAnchor = 10,
                         libType=libType,
                         strandMode=strandMode) {

    if(!is.null(bam)){
      .Deprecated("gbCounts")
    }
    minReadLength <- readLength
    cores <- 1 #Allways use 1 core.
 
    #Create result object
    counts <- new(Class="ASpliCounts")
    counts@.ASpliVersion = "1" #Last version before 2.0.0 was 1.14.0.
    
    #Generates sample names in case there arent any
    targets <- .generateSamplesNames(targets)
    counts@targets <- .condenseTargetsConditions(targets) #ACH
    group                  <- counts@targets$condition
    counts@condition.order <- levels(factor( group, unique( group ), ordered = TRUE ))
    
    #Minimal anchors
    minAnchor <- if ( ! is.null(minAnchor) ) minAnchor else 10
    minA <- round( minAnchor * minReadLength / 100 )
    ptm <- proc.time()
    if(is.null(bam)) {
      ntargets <- nrow(targets)
    }else{
      ntargets <- 1
    }
    
    for(target in 1:ntargets){
      
      if(ntargets > 1 | is.null(bam)){
        #Verbose
        message(paste("Summarizing", rownames(targets)[target]))
        
        #Load bam from current target#
        #aca hay que pasarle el parametro de SE o PE, y strandMode
        bam <- loadBAM(targets[target, ], cores = NULL,
                       libType=libType, 
                       strandMode=strandMode) #With cores = NULL wont print deprecated message
      }      
      
      # Count Genes
      gene.hits <- .counterGenes( bam, featuresg( features ))
      if(ncol(counts@gene.counts) == 0){
        counts@gene.counts <- gene.hits
      }else{
        counts@gene.counts <- cbind(counts@gene.counts, .extractCountColumns(gene.hits, targets[target, ]))
        colnames(counts@gene.counts)[ncol(counts@gene.counts)] <- rownames(targets)[target]
      }
      if(ntargets == 1) message("Read summarization by gene completed")
      
      # Count exons 
      bins <- featuresb( features )
      exons.hits <- .counterBin( bam, bins, gene.hits)
      if(ncol(counts@exon.intron.counts) == 0){
        counts@exon.intron.counts <- exons.hits
      }else{
        counts@exon.intron.counts <- cbind(counts@exon.intron.counts, .extractCountColumns(exons.hits, targets[target, ]))
        colnames(counts@exon.intron.counts)[ncol(counts@exon.intron.counts)] <- rownames(targets)[target]      
      }
      if(ntargets == 1) message( "Read summarization by bin completed" )
      
      # Count introns
      introns <- c( bins[ mcols(bins)$feature == "I" ], 
                    bins[ mcols(bins)$feature == "Io"],
                    bins[ mcols(bins)$eventJ  == "IR"])
      
      # Count exon1 - intron regions
      e1i <- introns
      start( e1i ) <- start( introns ) - ( minReadLength - minA )
      end( e1i )   <- start( introns ) + ( minReadLength - minA )
      e1i.hits     <- .counterJbin(bam, e1i, gene.hits, minReadLength)
      if(ncol(counts@e1i.counts) == 0){
        counts@e1i.counts <- e1i.hits
      }else{
        counts@e1i.counts <- cbind(counts@e1i.counts, .extractCountColumns(e1i.hits, targets[target, ]))
        colnames(counts@e1i.counts)[ncol(counts@e1i.counts)] <- rownames(targets)[target]         
      }
      if(ntargets == 1) message("Read summarization by ei1 region completed")
      
      # Count intron - exon2 regions
      ie2 <- introns
      start( ie2 ) <- end( introns ) - ( minReadLength - minA )
      end( ie2 )   <- end( introns ) + ( minReadLength - minA )
      ie2.hits     <- .counterJbin( bam, ie2, gene.hits, minReadLength )
      if(ncol(counts@ie2.counts) == 0){
        counts@ie2.counts <- ie2.hits
      }else{
        counts@ie2.counts <- cbind(counts@ie2.counts, .extractCountColumns(ie2.hits, targets[target, ]))
        colnames(counts@ie2.counts)[ncol(counts@ie2.counts)] <- rownames(targets)[target]   
      }
      if(ntargets == 1) message("Read summarization by ie2 region completed")

      # Count junctions
      junction.hits    <- .counterJunctions( features, bam, maxISize )
      if(ncol(counts@junction.counts) == 0){
        counts@junction.counts <- junction.hits
      }else{
        dt1                    <- data.table(counts@junction.counts, keep.rownames = TRUE)
        dt2                    <- data.table(.extractCountColumns(junction.hits, targets[target, ]), keep.rownames = T)
        dt3                    <- data.frame(merge(dt1, dt2, by="rn", all.x=T, all.y=TRUE))
        for(s in c("junction", "gene", "strand", "multipleHit", "symbol", "gene_coordinates", "bin_spanned", "j_within_bin")){
          dt3[, s]           <- as.character(dt3[, s])
          junction.hits[, s] <- as.character(junction.hits[, s])
        }
        rownames(dt3)          <- dt3[, "rn"]
        dt3                    <- dt3[, -1]
        dt3[dt2$rn, 1:8]       <- .extractDataColumns(junction.hits, targets[target, ])
        counts@junction.counts <- dt3
        counts@junction.counts[is.na(counts@junction.counts)] <- 0
      }
      if(ntargets == 1) message("Junction summarization completed")
      if(length(grep("NA", rownames(counts@junction.counts))) > 0){
        print(target)
        break
      }
      if(length(grep("NA", rownames(junction.hits ))) > 0){
        print(target)
      }
      gc()
      if(ntargets > 1){     
        sptm <- (proc.time() - ptm)[3]/target/60        
        message(paste("ETA:", round(sptm*(nrow(targets) - target)), "min"))      
      }
    }
    
    for(s in c("junction", "gene", "strand", "multipleHit", "symbol", "gene_coordinates", "bin_spanned", "j_within_bin")){
      counts@junction.counts[, s] <- as.factor(counts@junction.counts[, s])
    }
    colnames(counts@junction.counts)[9:ncol(counts@junction.counts)] <- rownames(targets)
    junctions.order <- sort(rownames(counts@junction.counts))
    junctions.order <- strsplit2(junctions.order, "[.]")
    junctions.order <- GRanges(seqnames=junctions.order[, 1], IRanges(start=as.numeric(junctions.order[, 2]), end=as.numeric(junctions.order[, 3])))
    junctions.order <- sort(junctions.order)
    junctions.order <- paste(junctions.order@seqnames, junctions.order@ranges@start, (junctions.order@ranges@start+junctions.order@ranges@width-1), sep=".")
    counts@junction.counts <- counts@junction.counts[junctions.order, ]
    
    # Create result object
    counts <- rds( counts, targets )
    gc()
    return(counts)
    
  }
)

# jCounts is a wrapper around AsDiscover for improved legibility
setGeneric (
  name= "jCounts",
  def = function( counts, 
                  features, 
                  minReadLength, 
                  threshold = 5, 
                  minAnchor = 10,
                  libType="SE",
                  strandMode=0) standardGeneric("jCounts") )

setMethod(
  f = "jCounts",
  signature = "ASpliCounts",
  definition = function( counts, 
                         features, 
                         minReadLength, 
                         threshold = 5, 
                         minAnchor = 10,
                         libType="SE",
                         strandMode=0) {
    if(!.hasSlot(counts, ".ASpliVersion")){
      counts@.ASpliVersion = "1" #Last version before 2.0.0 was 1.14.0. 
    }
    if(counts@.ASpliVersion == "1"){
      stop("Your version of ASpliCounts can not be used with this version of ASpli, please run gbCounts first. See vignette for details on the new pipeline.")
    }
    as <- AsDiscover( counts = counts, 
                      targets = NULL, 
                      features = features, 
                      bam = NULL, readLength = minReadLength, 
                      threshold = threshold,
                      cores = 1, minAnchor = 10,
                      libType = libType,
                      strandMode = strandMode )
    as@.ASpliVersion = "2" #Marks ASpliCounts object with the ASpli update 2.0.0    
    return(as)
  }
)

setGeneric (
  name= "AsDiscover",
  def = function( counts, 
                  targets,
                  features, 
                  bam, 
                  readLength, 
                  threshold = 5,
                  cores = 1, 
                  minAnchor = 10,
                  libType=libType,
                  strandMode=strandMode
                  ) standardGeneric("AsDiscover") )

setMethod(
  f = "AsDiscover",
  signature = "ASpliCounts",
  definition = function( counts, 
                         targets,
                         features, 
                         bam, 
                         readLength, 
                         threshold = 5,
                         cores = 1, 
                         minAnchor = 10,
                         libType=libType,
                         strandMode=strandMode) {
    
    if(!.hasSlot(counts, ".ASpliVersion")){
      counts@.ASpliVersion = "1" #Last version before 2.0.0 was 1.14.0. 
    }
    if(counts@.ASpliVersion == "1"){
      #Version conflict
      if(is.null(bam)){
        stop("Counts object is ASpli v1 but no bam was loaded. Please see vignette for new pipeline.")
      }
      .Deprecated("jCounts")
    }else{
      targets <- counts@targets
    }
    minReadLength <- readLength
    cores <- 1
    libType=libType
    strandMode=strandMode
    
    as  <- new(Class = "ASpliAS")
    as@.ASpliVersion = "1" #Last version before 2.0.0 was 1.14.0.    
    as@targets <- targets
    
    df0 <- countsj(counts)[ countsj(counts)$multipleHit == "-", ]
    df0 <- df0[ df0$gene != "noHit" , ]
    
    targets <- .condenseTargetsConditions( targets )
    
    jcounts <- .filterJunctionBySample( df0=df0,
                                        targets=targets,
                                        threshold=threshold )
    
    # Junctions PSI:
    junctionsPSI    <- .junctionsPSI_SUM( df0, targets )
    as@junctionsPJU <- junctionsPSI
    message("Junctions PJU completed")
    
    # Junctions PIR:
    if(is.null(bam)) {
      ntargets <- nrow(targets)
      for(target in 1:ntargets){
         
         if(ntargets > 1){
          #Load bam from current target
          #agrego el libType y StrandMode
          bam <- loadBAM(targets[target, ], cores = NULL, 
                         libType=libType, strandMode=strandMode)
          junctionsPIR <- .junctionsDiscover( df=jcounts, 
                                              minReadLength=minReadLength, 
                                              targets=targets[target, ], 
                                              features=features,
                                              minAnchor = minAnchor,
                                              bam=bam)  
          
         }
        
        
        
        if(ncol(as@junctionsPIR) == 0){
          as@junctionsPIR <- junctionsPIR
        }else{
          as@junctionsPIR <- cbind(as@junctionsPIR, junctionsPIR[, 3:6])
        }
      }
      junctions.order <- c(1, 2, 
                           seq(from=3, to=ncol(as@junctionsPIR), by=4),
                           seq(from=4, to=ncol(as@junctionsPIR), by=4),
                           seq(from=5, to=ncol(as@junctionsPIR), by=4))
      as@junctionsPIR <- as@junctionsPIR[, junctions.order]
      colnames(as@junctionsPIR)[c(-2, -1)] <- rep(rownames(targets), times=3)
      inicio_j1 <- 3
      inicio_j2 <- inicio_j1+nrow(targets)
      inicio_j3 <- inicio_j2+nrow(targets)
      j1 <- .sumByCond( as@junctionsPIR[, inicio_j1:(inicio_j1+nrow(targets)-1)],     targets )
      j2 <- .sumByCond( as@junctionsPIR[, inicio_j2:(inicio_j2+nrow(targets)-1)],     targets )
      j3 <- .sumByCond( as@junctionsPIR[, inicio_j3:(inicio_j3+nrow(targets)-1)],     targets )
      pirValues <- ( j1 + j2 ) / ( j1 + j2 + 2 * j3 )
      as@junctionsPIR <- cbind(as@junctionsPIR, pirValues)
    }else{
      
      junctionsPIR <- .junctionsDiscover( df=jcounts, 
                                          minReadLength=minReadLength, 
                                          targets=targets, 
                                          features=features,
                                          minAnchor = minAnchor,
                                          bam=bam) 
      as@junctionsPIR <- junctionsPIR
    }
    message("Junctions PIR completed")
    
    jranges <- .createGRangesExpJunctions( rownames( jcounts ) )
    
    # : refactor this code to other functions 
    # ---------------------------------------------------------------------- #
    # Get all bins that are intronic or are associated to a Intron retention 
    # event
    ic <- rbind( countsb(counts)[countsb(counts)$feature == "I",], 
                 countsb(counts)[countsb(counts)$feature == "Io",], 
                 countsb(counts)[countsb(counts)$event   == "IR*",],
                 countsb(counts)[countsb(counts)$event   == "IR",])
    # Get A GRanges object for intron bins, ordered by ic
    intranges <- featuresb(features)[ rownames(ic) ]
    
    # get exclusion junction counts, and make and index to ordered by ic
    dfe1e2 <- .e1e2JPIR( intranges, jcounts, targets )
    colnames(dfe1e2)[c(-1,-2)] <- rownames(targets)
    indexOrder <- match( dfe1e2$jbin, rownames( ic ) )
    
    # Get counts of inclusion junctions
    e1i <- .extractCountColumns( countse1i( counts ), targets )[ rownames(ic) ,]
    ie2 <- .extractCountColumns( countsie2( counts ), targets )[ rownames(ic) ,]
    
    j3 <- data.frame( matrix( NA, 
                              nrow =  nrow( e1i ), 
                              ncol =  length( targets$condition ) ), 
                      stringsAsFactors = FALSE )
    colnames( j3 ) <- colnames( e1i )  
    
    j3bin <- rep( NA , nrow( j3 ) )
    j3bin[ indexOrder ] <- rownames( dfe1e2 )
    j3[ indexOrder, ] <- .extractCountColumns( dfe1e2, targets )
    
    # Sum exclusion and inclusion counts by condition
    sumE1i <- .sumByCond( e1i, targets )
    sumIe2 <- .sumByCond( ie2, targets )
    sumJ3  <- .sumByCond( j3,  targets )
    
    # Calculates pir
    pirValues <- ( sumE1i + sumIe2 ) / ( sumE1i + sumIe2 + 2 * sumJ3 )  
    
    # Creates result object
    result <- cbind( 
      data.frame( event = ic$event ), 
      data.frame( J1 = paste( rownames( e1i ), "E1I", sep="_") ), 
      e1i, 
      data.frame( J2 = paste( rownames( ie2 ), "IE2", sep="_") ), 
      ie2,
      data.frame( J3 = j3bin ),
      j3, 
      pirValues ) 
    
    
    message("Junctions IR PIR completed")
    
    as@irPIR <- result
    # ---------------------------------------------------------------------- #
    
    # ---------------------------------------------------------------------- #
    # Get all exons, except those that are associated to a intron retention
    # event
    ec <- countsb(counts)[countsb(counts)$feature == "E",]
    ec <- ec[ec$event != "IR",]
    ec <- ec[ec$event != "IR*",]
    
    exranges <- featuresb( features )[ rownames( ec ) ]
    
    fillAndReorderBy <- function( df , orderNames ) {
      indexOrder <- match( rownames( df ) , orderNames )
      result <- data.frame( 
        matrix( 
          NA,
          nrow = length( orderNames ),
          ncol = ncol( df ) ) )
      result[ indexOrder, ] <- df
      colnames( result ) <- colnames( df )
      rownames( result ) <- orderNames
      return( result )
    }
    
    dfstart  <- .getJPSIByOverlap( jranges, exranges, jcounts, targets, 'start' )
    dfstart  <- fillAndReorderBy( dfstart , rownames( ec ) )
    dfend    <- .getJPSIByOverlap( jranges, exranges, jcounts, targets, 'end' )   
    dfend    <- fillAndReorderBy( dfend , rownames( ec ) )
    dfwithin <- .getJPSIByOverlap( jranges, exranges, jcounts, targets, 'within' )
    dfwithin <- fillAndReorderBy( dfwithin , rownames( ec ) )
    
    events   <- mcols( exranges ) $ event
    # ---------------------------------------------------------------------- #
    
    # ---------------------------------------------------------------------- #
    # Get the subset of previosly selected exons and gets only those associated
    # with an alternative splicing site usage event 
    getAlternativeSS <- function( df, events ) {
      rbind(
        df[ events == "Alt3ss", ],
        df[ events == "Alt5ss", ],
        df[ events == "Alt3ss*", ],
        df[ events == "Alt5ss*", ] )
    }
    
    altJ1 <- getAlternativeSS( dfstart , events )
    altJ2 <- getAlternativeSS( dfend , events )
    altJ3 <- getAlternativeSS( dfwithin , events )
    colnames(altJ1)[-ncol(altJ1)] <- rownames(targets)
    colnames(altJ2)[-ncol(altJ2)] <- rownames(targets)
    colnames(altJ3)[-ncol(altJ3)] <- rownames(targets)
      
    sumAltJ1 <- .sumByCond( .extractCountColumns( altJ1, targets ), targets )
    sumAltJ1[is.na(sumAltJ1)] <- 0 
    sumAltJ2 <- .sumByCond( .extractCountColumns( altJ2, targets ), targets )
    sumAltJ2[is.na(sumAltJ2)] <- 0
    sumAltJ3 <- .sumByCond( .extractCountColumns( altJ3, targets ), targets )
    sumAltJ3[is.na(sumAltJ3)] <- 0
    
    altPsiValues <- ( sumAltJ1 + sumAltJ2 ) / ( sumAltJ1 + sumAltJ2 + sumAltJ3 )
    
    result <- cbind( 
      data.frame( event = mcols( exranges[ rownames( altJ1) ] )$ event ), 
      data.frame( J1 = altJ1$overlappedSubjectNames ), 
      .extractCountColumns( altJ1, targets ),
      data.frame( J2 = altJ2$overlappedSubjectNames ), 
      .extractCountColumns( altJ2, targets ),
      data.frame( J3 = altJ3$overlappedSubjectNames ),
      .extractCountColumns( altJ3, targets ), 
      altPsiValues )
    
    message("Junctions AltSS PSI completed")
    altPSI( as ) <- result
    # ---------------------------------------------------------------------- #
    
    # ---------------------------------------------------------------------- #
    # Get the subset of previosly selected exons and gets only those associated
    # with an exon skipping event and those not assigned to any splice event.  
    getES <- function( df, events ) {
      rbind(
        df[ events == "ES", ],
        df[ events == "-", ],
        df[ events == "ES*", ] )
    }
    
    esJ1 <- getES( dfstart , events )
    esJ2 <- getES( dfend , events )
    esJ3 <- getES( dfwithin , events )
    colnames(esJ1)[-ncol(esJ1)] <- rownames(targets)
    colnames(esJ2)[-ncol(esJ2)] <- rownames(targets)
    colnames(esJ3)[-ncol(esJ3)] <- rownames(targets)
    
    sumEsJ1 <- .sumByCond( .extractCountColumns( esJ1, targets ), targets )
    sumEsJ1[is.na(sumEsJ1)] <- 0 
    sumEsJ2 <- .sumByCond( .extractCountColumns( esJ2, targets ), targets )
    sumEsJ2[is.na(sumEsJ2)] <- 0
    sumEsJ3 <- .sumByCond( .extractCountColumns( esJ3, targets ), targets )
    sumEsJ3[is.na(sumEsJ3)] <- 0
    
    esPsiValues <- ( sumEsJ1 + sumEsJ2 ) / ( sumEsJ1 + sumEsJ2 + 2 * sumEsJ3 )
    
    result <- cbind( 
      data.frame( event = mcols( exranges[ rownames( esJ1) ] )$ event ), 
      data.frame( J1 = esJ1$overlappedSubjectNames ), 
      .extractCountColumns( esJ1, targets ), 
      data.frame( J2 = esJ2$overlappedSubjectNames ), 
      .extractCountColumns( esJ2, targets ),
      data.frame( J3 = esJ3$overlappedSubjectNames ),
      .extractCountColumns( esJ3, targets ),
      esPsiValues )
    
    message("Junctions ES PSI completed")
    
    esPSI( as ) <- result
    # ---------------------------------------------------------------------- #
    
    # TODO: joint podria ser un getter, pero no es necesario mantener toda
    # esta data repetida
    joint( as ) <- rbind( altPSI( as ), esPSI( as ), irPIR( as ) )
    
    return( as )
    
  })

setMethod( 
  f = 'subset',
  signature = 'ASpliAS',
  def = function( x, targets, select) .subset.ASpliAS( x, targets, select ) )

# ---------------------------------------------------------------------------- #
# writeAS
setGeneric (
  name = "writeAS",
  def  = function(as, output.dir = "as" )
    standardGeneric( "writeAS" ) )

setMethod(
  f = "writeAS",
  signature = "ASpliAS",
  definition = function( as, output.dir = "as" ) {
    
    # Creates output folder structure
    exonsFilePSI     <- file.path( output.dir, "exons", "exon.altPSI.tab" )       
    exonsFileES      <- file.path( output.dir, "exons", "exon.altES.tab" )       
    intronsFile      <- file.path( output.dir, "introns", "intron.irPIR.tab" )
    junctionsFilePIR <- file.path( output.dir, "junctions", "junction.PIR.tab" )
    junctionsFilePSI <- file.path( output.dir, "junctions", "junction.PSI.tab" )
    asDiscoverFile   <- file.path( output.dir, "as_discovery.tab" )
    
    file.exists( output.dir ) || dir.create( output.dir )
    for ( folder in unique( lapply( c( exonsFilePSI, exonsFileES, intronsFile, 
                                       junctionsFilePIR, junctionsFilePSI ), dirname ) ) ) {
      dir.create( folder )
    }
    
    # Export exons
    write.table( altPSI(as), exonsFilePSI, sep="\t", quote=FALSE, col.names=NA)
    write.table( esPSI(as), exonsFileES, sep="\t", quote=FALSE, col.names=NA)
    
    # Export Introns
    write.table( irPIR(as), intronsFile, sep="\t", quote=FALSE, col.names=NA)
    
    # Export Junctions
    write.table(junctionsPIR(as), junctionsFilePIR, sep="\t", quote=FALSE, col.names=NA)
    write.table(junctionsPJU(as), junctionsFilePSI, sep="\t", quote=FALSE, col.names=NA)
    
    # Export AS discovery table
    write.table( joint(as), asDiscoverFile, sep="\t", quote=FALSE, col.names=NA)
  }
)

# TODO:  Es necesario agregar todos los parametros con valores por default en
# la firma del metodo ? 
setGeneric (
  name = "DUreport.norm",
  def = function( counts, 
                  minGenReads  = 10,
                  minBinReads  = 5,
                  minRds = 0.05,
                  contrast = NULL,
                  ignoreExternal = TRUE,
                  ignoreIo = TRUE, 
                  ignoreI = FALSE,
                  filterWithContrasted = TRUE,
                  verbose = FALSE,
                  threshold = 5
  ) standardGeneric("DUreport.norm") )

#setGeneric (
#  name = "DUreport_DEXSeq",
#  def = function ( counts, ... ) standardGeneric("DUreport_DEXSeq") )

setMethod(
  f = "DUreport.norm",
  signature = "ASpliCounts",
  definition = function( counts, 
                         minGenReads  = 10,
                         minBinReads  = 5,
                         minRds = 0.05,
                         contrast = NULL,
                         ignoreExternal = TRUE,
                         ignoreIo = TRUE, 
                         ignoreI = FALSE,
                         filterWithContrasted = TRUE,
                         verbose = FALSE,
                         threshold = 5
  ) { 
    offset = FALSE
    offsetAggregateMode = c( "geneMode", "binMode" )[1]
    offsetUseFitGeneX = TRUE    
    .DUreport( counts, counts@targets, minGenReads, minBinReads, minRds, offset, 
               offsetAggregateMode, offsetUseFitGeneX, contrast, 
               ignoreExternal, ignoreIo, ignoreI, filterWithContrasted, verbose, threshold  )
  }
)

setGeneric (
  name = "DUreport.offset",
  def = function( counts, 
                  minGenReads  = 10,
                  minBinReads  = 5,
                  minRds = 0.05,
                  offsetAggregateMode = c( "geneMode", "binMode" )[1],
                  offsetUseFitGeneX = TRUE,
                  contrast = NULL,
                  ignoreExternal = TRUE,
                  ignoreIo = TRUE, 
                  ignoreI = FALSE,
                  filterWithContrasted = TRUE,
                  verbose = FALSE
  ) standardGeneric("DUreport.offset") )

setMethod(
  f = "DUreport.offset",
  signature = "ASpliCounts",
  definition = function( counts, 
                         minGenReads  = 10,
                         minBinReads  = 5,
                         minRds = 0.05,
                         offsetAggregateMode = c( "geneMode", "binMode" )[1],
                         offsetUseFitGeneX = TRUE,
                         contrast = NULL,
                         ignoreExternal = TRUE,
                         ignoreIo = TRUE, 
                         ignoreI = FALSE,
                         filterWithContrasted = TRUE,
                         verbose = FALSE
  ) { 
    offset = TRUE
    .DUreport( counts, counts@targets, minGenReads, minBinReads, minRds, offset, 
               offsetAggregateMode, offsetUseFitGeneX, contrast,                ignoreExternal, ignoreIo, ignoreI, filterWithContrasted, verbose  )
  }
)


setGeneric( name = 'gbDUreport',
            def = function( counts, 
                            minGenReads  = 10, minBinReads = 5, 
                            minRds = 0.05, contrast = NULL, 
                            ignoreExternal = TRUE, ignoreIo = TRUE, 
                            ignoreI = FALSE, 
                            filterWithContrasted = TRUE, 
                            verbose = TRUE, 
                            formula = NULL, coef = NULL ) 
              standardGeneric( 'gbDUreport'))

setMethod( 
  f = 'gbDUreport',
  signature = 'ASpliCounts',
  definition = function( counts, 
                         minGenReads  = 10, 
                         minBinReads = 5,
                         minRds = 0.05, 
                         contrast = NULL, 
                         ignoreExternal = TRUE,
                         ignoreIo = TRUE, 
                         ignoreI = FALSE, 
                         filterWithContrasted = TRUE,
                         verbose = TRUE, 
                         formula = NULL,
                         coef = NULL) {
  if(!.hasSlot(counts, ".ASpliVersion")){
    counts@.ASpliVersion = "1" #Last version before 2.0.0 was 1.14.0. 
  }
  if(counts@.ASpliVersion == "1"){
    stop("Your version of ASpliCounts can not be used with this version of ASpli, please run gbCounts first. See vignette for details on the new pipeline.")
  }
  
  du <- .DUreportBinSplice( counts, targets = NULL, minGenReads, minBinReads, minRds, 
                            contrast, forceGLM = NULL, ignoreExternal, ignoreIo, ignoreI, 
                            filterWithContrasted, verbose, formula, coef ) #forceGLM was deprecated
  du@.ASpliVersion <- "2"
  return(du)  
})

setGeneric( name = 'DUreportBinSplice',
            def = function( counts, targets, minGenReads  = 10, minBinReads = 5, 
                            minRds = 0.05, contrast = NULL, forceGLM = FALSE,  
                            ignoreExternal = TRUE, ignoreIo = TRUE, ignoreI = FALSE, 
                            filterWithContrasted = FALSE, verbose = TRUE ) 
              standardGeneric( 'DUreportBinSplice'))

setMethod( 
  f = 'DUreportBinSplice',
  signature = 'ASpliCounts',
  definition = function( counts, 
                         targets, 
                         minGenReads  = 10, 
                         minBinReads = 5,
                         minRds = 0.05, 
                         contrast = NULL, 
                         forceGLM = FALSE, 
                         ignoreExternal = TRUE,
                         ignoreIo = TRUE, 
                         ignoreI = FALSE, 
                         filterWithContrasted = TRUE,
                         verbose = TRUE ) {
    .Deprecated("gbDUreport")
    du <- .DUreportBinSplice( counts, targets, minGenReads, minBinReads, minRds, 
                              contrast, forceGLM, ignoreExternal, ignoreIo, ignoreI, 
                              filterWithContrasted, verbose = TRUE ) 
    du@.ASpliVersion <- "1"
    return(du)
  })


setGeneric( name = "junctionDUreport",
            def = function (  counts, 
                              targets, 
                              appendTo = NULL, 
                              minGenReads = 10,
                              minRds = 0.05,
                              threshold = 5,
                              offset   = FALSE,
                              offsetUseFitGeneX = TRUE,
                              contrast = NULL,
                              forceGLM = FALSE 
            ) standardGeneric("junctionDUreport") )


setMethod(
  f = "junctionDUreport",
  signature = "ASpliCounts",
  definition = function ( 
    counts, 
    targets, 
    appendTo = NULL,
    minGenReads = 10,
    minRds = 0.05,
    threshold = 5,
    offset = FALSE,
    offsetUseFitGeneX = TRUE,
    contrast = NULL,
    forceGLM = FALSE 
    # -------------------------------------------------------------------- #
    # Comment to disable priorcounts usage in bin normalization 
    # , priorCounts = 0 
    # -------------------------------------------------------------------- #
  ) {
    .Deprecated("jDUreport")
    .junctionDUreport( counts, targets, appendTo,  minGenReads,  minRds, 
                       threshold, offset, offsetUseFitGeneX, contrast, 
                       forceGLM ) 
  }
)

setGeneric( name = "jDUreport",
            def = function (asd, 
                            minAvgCounts                       = 5, 
                            contrast                           = NULL,
                            filterWithContrasted               = TRUE,
                            runUniformityTest                  = FALSE,
                            mergedBams                         = NULL,
                            maxPValForUniformityCheck          = 0.2,
                            strongFilter                       = TRUE,
                            maxConditionsForDispersionEstimate = 24,
                            formula                            = NULL,
                            coef                               = NULL,
                            maxFDRForParticipation             = 0.05,
                            useSubset                          = FALSE
            ) standardGeneric("jDUreport") )


setMethod(
  f = "jDUreport",
  signature = "ASpliAS",
  definition = function (
    asd,
    minAvgCounts                       = 5, 
    contrast                           = NULL,
    filterWithContrasted               = TRUE,
    runUniformityTest                  = FALSE,
    mergedBams                         = NULL,
    maxPValForUniformityCheck          = 0.2,
    strongFilter                       = TRUE,
    maxConditionsForDispersionEstimate = 24,
    formula                            = NULL,
    coef                               = NULL,
    maxFDRForParticipation             = 0.2,
    useSubset                          = FALSE
  ) {

    if(!.hasSlot(asd, ".ASpliVersion")){
      asd@.ASpliVersion = "1" #Last version before 2.0.0 was 1.14.0. 
    }
    if(asd@.ASpliVersion == "1"){
      stop("Your version of ASpliAS can not be used with this version of ASpli, please run jCounts first. See vignette for details on the new pipeline.")
    }
    jdu <- .junctionDUreportExt( asd, 
                                 minAvgCounts, 
                                 contrast, 
                                 filterWithContrasted, 
                                 runUniformityTest, 
                                 mergedBams, 
                                 maxPValForUniformityCheck, 
                                 strongFilter,
                                 maxConditionsForDispersionEstimate, 
                                 formula, 
                                 coef, 
                                 maxFDRForParticipation, 
                                 useSubset) 
    jdu@.ASpliVersion <- "2"
    return(jdu)
  }
)


setGeneric( name = "splicingReport",
            def = function (bdu, 
                            jdu, 
                            counts
            ) standardGeneric("splicingReport") )


setMethod(
  f = "splicingReport",
  signature = "ASpliDU",
  definition = function (
    bdu,
    jdu,
    counts
  ) {
    if(!.hasSlot(bdu, ".ASpliVersion")){
      bdu@.ASpliVersion = "1" #Last version before 2.0.0 was 1.14.0. 
    }
    if(bdu@.ASpliVersion == "1"){
      stop("Your version of ASpliDU can not be used with this version of ASpli, please run gbDUreport first. See vignette for details on the new pipeline.")
    }
    if(!.hasSlot(jdu, ".ASpliVersion")){
      jdu@.ASpliVersion = "1" #Last version before 2.0.0 was 1.14.0. 
    }
    if(jdu@.ASpliVersion == "1"){
      stop("Your version of ASpliJDU can not be used with this version of ASpli, please run jDUreport first. See vignette for details on the new pipeline.")
    }
    if(!.hasSlot(counts, ".ASpliVersion")){
      counts@.ASpliVersion = "1" #Last version before 2.0.0 was 1.14.0. 
    }
    if(counts@.ASpliVersion == "1"){
      stop("Your version of ASpliCounts can not be used with this version of ASpli, please run gbCounts first. See vignette for details on the new pipeline.")
    }    
    sr <- .splicingReport( bdu, jdu, counts ) 
    sr@.ASpliVersion = "2"
    return(sr)
  }
)

setGeneric( name = "integrateSignals",
            def = function (sr = NULL,
                            asd = NULL, 
                            bin.FC = 3,
                            bin.fdr = 0.05,
                            nonunif = 1,
                            usenonunif = FALSE,
                            bin.inclussion = 0.2,
                            bjs.inclussion = 10.3,
                            bjs.fdr = 0.01,
                            a.inclussion = 0.3,
                            a.fdr = 0.01,
                            l.inclussion = 0.3,
                            l.fdr = 0.01,
                            otherSources = NULL,
                            overlapType = "any"
            ) standardGeneric("integrateSignals") )


setMethod(
  f = "integrateSignals",
  signature = "ASpliSplicingReport",
  definition = function (
    sr = NULL,
    asd = NULL, 
    bin.FC = 3,
    bin.fdr = 0.05,
    nonunif = 1,
    usenonunif = FALSE,
    bin.inclussion = 0.2,
    bjs.inclussion = 10.3,
    bjs.fdr = 0.01,
    a.inclussion = 0.3,
    a.fdr = 0.01,
    l.inclussion = 0.3,
    l.fdr = 0.01,
    otherSources = NULL,
    overlapType = "any"
  ) {
    if(!.hasSlot(sr, ".ASpliVersion")){
      sr@.ASpliVersion = "1" #Last version before 2.0.0 was 1.14.0. 
    }
    if(sr@.ASpliVersion == "1"){
      stop("Your version of ASpliSplicingReport can not be used with this version of ASpli, please run splicingReport first. See vignette for details on the new pipeline.")
    }
    if(!.hasSlot(asd, ".ASpliVersion")){
      asd@.ASpliVersion = "1" #Last version before 2.0.0 was 1.14.0. 
    }
    if(asd@.ASpliVersion == "1"){
      stop("Your version of ASpliAS can not be used with this version of ASpli, please run jCounts first. See vignette for details on the new pipeline.")
    }    
    is <- .integrateSignals(sr, asd, bin.FC, bin.fdr, nonunif, usenonunif, bin.inclussion, bjs.inclussion, bjs.fdr, a.inclussion, a.fdr,
                      l.inclussion, l.fdr, otherSources, overlapType) 
    is@.ASpliVersion = "2"
    return(is)
  }
)

setMethod( f = 'subset',
           signature = 'ASpliCounts',
           def = function( x, targets, select ) { .subset.ASpliCounts( x, targets, select ) }  )

setGeneric( 
  name = 'filterDU',
  def = function( du, what = c( 'genes','bins','junctions'), fdr = 1, 
                  logFC = 0, absLogFC = TRUE, logFCgreater = TRUE ) standardGeneric('filterDU') )

setMethod(
  f = 'filterDU',
  signature = "ASpliDU",
  definition = function( du, what = c( 'genes','bins','junctions'), fdr = 1, 
                         logFC = 0, absLogFC = TRUE, logFCgreater = TRUE ) {
    .Deprecated("", msg = "filterDU is deprecated and is no longer needed. See ASpli vignette for new pipeline.")
    .filter.ASpliDU( du, what, fdr, logFC, absLogFC, logFCgreater ) } )


# ---------------------------------------------------------------------------- #
# writeDU

setGeneric( name = 'containsGenesAndBins', 
            def = function ( du ) standardGeneric("containsGenesAndBins") )

setMethod( f = 'containsGenesAndBins', 
           signature = "ASpliDU",
           definition = function ( du ) {
             nrow( genesDE( du ) ) > 0  & nrow( binsDU( du) ) > 0 
           } )

setGeneric( name = 'containsJunctions', 
            def = function ( du ) standardGeneric("containsJunctions") )

setMethod( f = 'containsJunctions', 
           signature = "ASpliDU",
           definition = function ( du ) {
             nrow( junctionsDU( du ) ) > 0
           } )

setGeneric( name = "writeDU", 
            def = function ( du, output.dir="du"  ) standardGeneric( "writeDU" ) )

setMethod(
  f = "writeDU",
  signature = "ASpliDU",
  definition = function( du, output.dir="du" ) {

    file.exists( output.dir ) || dir.create( output.dir , recursive = TRUE)
    output.dir <- paste(output.dir, paste(names(du@contrast)[du@contrast != 0], collapse="-"), sep="/")
    file.exists( output.dir ) || dir.create( output.dir , recursive = TRUE )

    if ( containsGenesAndBins( du ) ) {
      # Export Genes  
      write.table( genesDE( du ), paste(normalizePath(output.dir), "gene.de.tab", sep="/"), sep = "\t", quote = FALSE, 
                   col.names = NA )
      
      # Export Exons
      exonBins <- binsDU(du)[binsDU(du)$feature == "E",]
      exonBins <- exonBins[exonBins$event !="IR",]
      write.table( exonBins, paste(normalizePath(output.dir), "exon.du.tab", sep="/"), sep="\t", quote=FALSE, col.names=NA)
      
      # Export Introns 
      intronBins <- rbind( 
        binsDU(du)[binsDU(du)$feature == "I" ,], 
        binsDU(du)[binsDU(du)$feature == "Io",],
        binsDU(du)[binsDU(du)$event   == "IR",])
      write.table( intronBins, paste(normalizePath(output.dir), "intron.du.tab", sep="/"), sep = "\t", quote = FALSE, 
                   col.names = NA )
    }
    # Export Junctions
    if ( containsJunctions( du ) ) {
      write.table( junctionsDU( du ), paste(normalizePath(output.dir), "junction.du.tab", sep="/"), sep = "\t", quote = FALSE, 
                   col.names=NA )
    }
  }
)

setGeneric( name = "writeJDU", 
            def = function ( jdu, output.dir="jdu"  ) standardGeneric( "writeJDU" ) )

setMethod(
  f = "writeJDU",
  signature = "ASpliJDU",
  definition = function( jdu, output.dir="jdu" ) {
    
    file.exists( output.dir ) || dir.create( output.dir , recursive = TRUE )
    output.dir <- paste(output.dir, paste(names(jdu@contrast)[jdu@contrast != 0], collapse="-"), sep="/")
    file.exists( output.dir ) || dir.create( output.dir , recursive = TRUE )
    
    for(slotName in slotNames(jdu)){
      # Export Genes  
      write.table( slot(jdu, slotName), paste(output.dir, paste0(slotName, ".txt"), sep="/"), sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE )
    }
    
  }
)

setGeneric( name = "writeSplicingReport", 
            def = function ( sr, output.dir="sr"  ) standardGeneric( "writeSplicingReport" ) )

setMethod(
  f = "writeSplicingReport",
  signature = "ASpliSplicingReport",
  definition = function( sr, output.dir="sr" ) {
    
    file.exists( output.dir ) || dir.create( output.dir , recursive = TRUE)
    
    for(slotName in slotNames(sr)){
      # Export Genes  
      if(class(slot(sr, slotName)) == "data.frame"){
        write.table( slot(sr, slotName), paste(output.dir, paste0(slotName, ".txt"), sep="/"), sep = "\t", 
                     quote = FALSE, col.names = TRUE, row.names = FALSE )
      }
    }
    
  }
)


setGeneric( name = 'mergeBinDUAS',
            def = function( du, as, targets, contrast = NULL  ) 
              standardGeneric( 'mergeBinDUAS' ))

setMethod( f = 'mergeBinDUAS',
           signature = c( 'ASpliDU', 'ASpliAS' ),
           definition = function( du, as, targets, contrast = NULL  ) {
             .Deprecated("", msg = "mergeBinDUAS is deprecated and is no longer needed. See ASpli vignette for new pipeline.")
             .mergeBinDUAS( du, as, targets, contrast ) } )


setGeneric( name = "exportSplicingReports", 
            def = function ( sr, output.dir="sr" , openInBrowser = FALSE, maxBinFDR = 0.2, maxJunctionFDR = 0.2 ) standardGeneric( "exportSplicingReports" ) )

setMethod(
  f = "exportSplicingReports",
  signature = "ASpliSplicingReport",
  definition = function( sr, output.dir="sr" , openInBrowser = FALSE, maxBinFDR = 0.2, maxJunctionFDR = 0.2 ) {
    output.dir <- paste0(output.dir, "/", paste0(names(sr@contrast)[sr@contrast != 0], collapse="-"))
    file.exists( output.dir ) || dir.create( output.dir , recursive = TRUE)
    
    
    for(s in slotNames(sr)){
      if(class(slot(sr, s)) == "data.frame"){
        b <- slot(sr, s)
        if(s == "binbased"){
            b <- b[union(which(b$bin.fdr < maxBinFDR), which(b$junction.fdr < maxJunctionFDR)), ]
            b$feature <- as.factor(b$feature)
            b$bin.event <- as.factor(b$bin.event)
            #b[, c(10:12, 15:ncol(b))] <- apply(b[, c(10:12, 15:ncol(b))], 2, function(s){return(signif(as.numeric(s), digits = 4))})
        }else{
            b <- b[union(which(b$bin.fdr < maxBinFDR), which(b$cluster.fdr < maxJunctionFDR)), ] 
        }
        if(nrow(b)==0){
          cat("Not a single",s,"event passed the filtering step. You might consider relaxing thresholds.\n")
          next
        }
        columnas_numericas <- which(sapply(b, class) == "numeric")
	      b[, columnas_numericas] <- apply(b[, columnas_numericas], 2, function(s){return(signif(as.numeric(s), digits = 4))})
        titulo <- paste0('ASpli: ', s, ". Contrasts: ", paste(names(sr@contrast)[sr@contrast != 0], collapse = " - "))
        if(s != "localebased"){
          y <- datatable(b,
                         escape = TRUE,
                         filter ="top",
                         extensions = c('Buttons', 'KeyTable'), 
                         options = list(dom = 'lfrtBip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print', I('colvis')),
                                        columnDefs = list(
                                          #  list(visible = FALSE, targets = c(0, 2, 3)),
                                          list(orderable = FALSE, className =
                                                 'details-control', targets = 1)
                                        ),
                                        keys = TRUE
                         ),   caption = htmltools::tags$caption(
                           style = 'caption-side: top; text-align: left;',
                           htmltools::h1(titulo)
                         ))    
        }else{
          clusters  <- unique(b[, c("junction.cluster", "cluster.locus", "cluster.size", "cluster.LR", "cluster.pvalue", "cluster.fdr", "cluster.range", "cluster.participation")])
          
          b <- slot(sr, s)
          columnas_numericas <- which(sapply(b, class) == "numeric")
          b[, columnas_numericas] <- apply(b[, columnas_numericas], 2, function(s){return(signif(as.numeric(s), digits = 4))})
          
          junctions <- unique(b[, c("junction.cluster", colnames(b)[!colnames(b) %in% colnames(clusters)])])
          colnames(clusters)[1] <- "cluster"
          subtables <- "var subtables = [];"
          
          for (i in clusters$cluster) {
            cluster_junctions <- junctions[junctions$junction.cluster == i, ]
            datos_bin         <- data.frame(bin = as.character(aggregate(bin ~ junction, cluster_junctions, FUN=function(s){paste(s, collapse=";")}, na.action=na.pass)[, 2]),
                                            bin.pvalue = as.character(aggregate(bin.pvalue ~ junction, cluster_junctions, FUN=function(s){paste(s, collapse=";")}, na.action=na.pass)[, 2]),
                                            bin.fdr = as.character(aggregate(bin.fdr ~ junction, cluster_junctions, FUN=function(s){paste(s, collapse=";")}, na.action=na.pass)[, 2]),
                                            stringsAsFactors = FALSE)
            cluster_junctions <- unique(cluster_junctions[, !colnames(cluster_junctions) %in% c("bin", "bin.pvalue", "bin.fdr")])
            cluster_junctions <- data.frame(cluster_junctions, datos_bin)
            subtables <- paste0(subtables, "subtables[", i, "] = '<table><tr>")
            subtables <- paste0(subtables, paste0(paste0("<th bgcolor=\"#808080\">", colnames(cluster_junctions), "</th>"), collapse=""))
            color <- "#FFFFFF"
            subtables <- paste0(subtables, "</tr>")         
            for(j in 1:nrow(cluster_junctions)){
              subtables <- paste0(subtables, '<tr>')
              aux <- gsub("rcolor", color, paste0(paste0("<td bgcolor=\"rcolor\">", cluster_junctions[j, ], "</td>"), collapse=""))
              subtables <- paste0(subtables, aux)
              subtables <- paste0(subtables, "</tr>")
              color <- ifelse(color == "#FFFFFF", "#CDCDCD", "#FFFFFF")
            }
            subtables <- paste0(subtables, "</table>';")
          }
                      
          text <- ";var format = function(d) { text = '<div>' + subtables[d[1]] + '</div>'; return text;};"
          y <- datatable(cbind(' ' = '&oplus;', clusters),
                         rownames = FALSE,
                         escape = -1,
                         filter ="top",
                         fillContainer = FALSE,
                         extensions = c('Buttons', 'KeyTable'), 
                         options = list(dom = 'lfrtBip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print', I('colvis')),
                                        columnDefs = list(
                                          #  list(visible = FALSE, targets = c(0, 2, 3)),
                                          list(orderable = FALSE, className =
                                                 'details-control', targets = 0)
                                        ),
                                        keys = TRUE
                         ),   caption = htmltools::tags$caption(
                           style = 'caption-side: top; text-align: left;',
                           htmltools::h1(titulo)
                         ), 
                         callback = JS(paste0("
                                      table.order([5, 'asc']).draw();
                                      table.column(0).nodes().to$().css({cursor: 'pointer'});",
                                      subtables,
                                      text,
                                      " 
                                        table.on('click', 'td.details-control', function() {
                                         var td = $(this), row = table.row(td.closest('tr'));
                                         if (row.child.isShown()) {
                                            row.child.hide();
                                            td.html('&oplus;');
                                        } else {
                                            row.child(format(row.data())).show();
                                            td.html('&CircleMinus;');
                                        }
                                     });"))
                  )       
        }
        ffile <- paste0(normalizePath(output.dir), "/", s, "Report.html")
        suppressWarnings(saveWidget(y, file = ffile, selfcontained = FALSE, title = paste0(names(sr@contrast)[sr@contrast != 0], collapse="-")))
        if(openInBrowser == TRUE) browseURL(ffile)
      }
    }    
  }
)

setGeneric( name = "exportIntegratedSignals", 
            def = function ( is, output.dir="is", 
                             sr, counts, features, asd,
                             mergedBams, 
                             jCompletelyIncluded = FALSE, 
                             zoomRegion = 1.5, 
                             useLog = FALSE, 
                             tcex = 1, 
                             ntop = NULL, 
                             openInBrowser = FALSE, 
                             makeGraphs = TRUE,
                             bforce=FALSE
                             ) standardGeneric( "exportIntegratedSignals" ) )

setMethod(
  f = "exportIntegratedSignals",
  signature = "ASpliIntegratedSignals",
  definition = function( is, output.dir="is", sr, counts, features, asd, mergedBams, jCompletelyIncluded = FALSE, zoomRegion = 1.5, useLog = FALSE, tcex = 1, ntop = NULL, 
                         openInBrowser = FALSE, makeGraphs = TRUE,bforce=FALSE) {

    if(class(is) != "ASpliIntegratedSignals"){
      stop("is must be an ASpliIntegratedSignals object")
    }
        
    if(class(sr) != "ASpliSplicingReport"){
      stop("sr must be an ASpliSplicingReport object")
    }

    if(class(counts) != "ASpliCounts"){
      stop("counts must be an ASpliCounts object")
    }

    if(class(features) != "ASpliFeatures"){
      stop("features must be an ASpliFeatures object")
    }
    
    if(class(asd) != "ASpliAS"){
      stop("asd must be an ASpliAS object")
    }
    
    if(nrow(mergedBams) > 4 & makeGraphs == TRUE){
      continue <- ""
      while(!continue %in% c("y", "n")){
        continue <- readline(prompt=paste("Warning, we are about to generate", nrow(mergedBams), "images for every event. This could take too long, do you want to continue? (y/n)"))
      }
      if(continue == "n"){
        return()
      }
    }
    filters <- is@filters
    is <- is@signals
    
    saux <- names(sr@contrast)[sr@contrast != 0]
    ina <- which(is.na(match(saux,mergedBams$condition)))
    if(length(ina)>0){
      stop("Check merged-bams data.frame.\nContrast condition(s)", paste(saux[ina],collapse=", ")," not found in any merged-bam data.frame condition category.")  
    }
    
    mergedBams <- mergedBams[mergedBams$condition %in% saux, ]
    if(nrow(mergedBams) == 0 ){
      stop("Merged bams dont match with contrasts")  
    }
    
    output.dir <- paste0(output.dir, "/", paste0(names(sr@contrast)[sr@contrast != 0], collapse="-"))
    output.dir <- substr(output.dir, 1, 255) #maximum dir length is 255
    file.exists( output.dir ) || dir.create( output.dir, recursive = TRUE)
    file.exists( paste0(output.dir, "/img") ) || dir.create( paste0(output.dir, "/img") )
    
    is$feature[is.na(is$feature)] <- "-"
    is$b <- as.factor(is$b)
    is$bjs <- as.factor(is$bjs)
    is$ja <- as.factor(is$ja)
    is$jl <- as.factor(is$jl)
    is$feature <- as.factor(is$feature)
    is$bin.event <- as.factor(is$bin.event)
    
  
    is$b.fdr <- signif(as.numeric(is$b.fdr), 4)
    is$b.logfc <- signif(as.numeric(is$b.logfc), 4)    
    is$bjs.lr <- signif(as.numeric(is$bjs.lr), 4)
    is$bjs.fdr <- signif(as.numeric(is$bjs.fdr), 4)    
    is$bjs.nonuniformity <- signif(as.numeric(is$bjs.nonuniformity), 4)
    is$bjs.inclussion <- signif(as.numeric(is$bjs.inclussion), 4)    
    is$a.lr <- signif(as.numeric(is$a.lr), 4)
    is$a.fdr <- signif(as.numeric(is$a.fdr), 4)        
    is$a.nonuniformity <- signif(as.numeric(is$a.nonuniformity), 4)
    is$l.lr <- signif(as.numeric(is$l.lr), 4)    
    is$l.fdr <- signif(as.numeric(is$l.fdr), 4)
    is$l.participation <- signif(as.numeric(is$l.participation), 4)    
    is$l.dparticipation <- signif(as.numeric(is$l.dparticipation), 4)
    is$a.dpir <- signif(as.numeric(is$a.dpir), 4)        
     
    is$bjs.inclussion_sign <- as.factor(.my_replace_na(sign(as.numeric(is$bjs.inclussion)), 0))     
    is$a.dpir_sign <- as.factor(.my_replace_na(sign(as.numeric(is$a.dpir)), 0))    
    is$a.dpir      <- abs(is$a.dpir)
    

    message("Generating graphs...")
    if(!is.numeric(ntop)){
        ntop <- length(is$region)
    }else{
      if(ntop < 1){
        ntop <- length(is$region)
      }
    }
    if(ntop > 1){
      pb <- txtProgressBar(min=1,max=ntop,style=3)
      if(makeGraphs){
        for(i in 1:ntop){
          r <- is$region[i]
          #if(i %% 10 == 0){
          # message(paste0(signif(i/ntop, 2)*100, "% completed"))
          #}
          setTxtProgressBar(pb,i)
          tryCatch({
            if(bforce | !file.exists(paste0(normalizePath(output.dir), "/img/", r, "_gene.png"))){
              png(width = 1400, height=700, filename = paste0(normalizePath(output.dir), "/img/", gsub("-", "_", gsub(":", "_", r)), "_gene.png"))
              .plotSplicingPattern(r, is, counts, features, mergedBams, sr, asd, genePlot = TRUE, jCompletelyIncluded, zoomRegion, useLog, tcex)
              dev.off()
            }
          }, warning = function(warning_condition) {
              #message(warning_condition)   
              #dev.off()
          }, error = function(error_condition) {
              #message(error_condition)
              #dev.off()
          }, finally={
            
          })
        }
      }
    }
    #close(pb)
    titulo <- paste(paste('ASpli: integrated signals. Contrasts:', 
                     paste(names(sr@contrast)[sr@contrast != 0], collapse = " - ")))
    sketch = htmltools::withTags(table(
      class = 'display',
      tags$thead(
        tags$tr(
          tags$th(rowspan = 2, 'View'),
          tags$th(rowspan = 2, 'Region'),
          tags$th(rowspan = 2, 'Event'),
          tags$th(rowspan = 2, 'Locus'),
          tags$th(rowspan = 2, 'Locus overlap'),
          tags$th(rowspan = 2, 'Bin Evidence'),
          tags$th(rowspan = 2, 'Bin SJ Evidence'),
          tags$th(rowspan = 2, 'Anchor Evidence'),
          tags$th(rowspan = 2, 'Locale Evidence'),
          tags$th(rowspan = 2, 'Bin'),
          tags$th(rowspan = 2, 'Feature'),
          tags$th(colspan = 2, 'Bins', bgcolor="#DCDCDC"),
          tags$th(colspan = 5, 'Bin Supporting Junctions', bgcolor="#C0C0C0"),
          tags$th(colspan = 5, 'Anchor Junctions', bgcolor="#DCDCDC"),
          tags$th(colspan = 5, 'Locale Junctions', bgcolor="#C0C0C0")
        ),
        tags$tr(
          lapply(c("logFC", "FDR", "LR", "FDR", "Non Uniformity", "dInclussion", "Sign Inclussion", "LR", "FDR", "Non uniformity", "dInclussion", "Sign Inclussion", "LR", "FDR", "Participation", "dParticipation"), tags$th)
        )
      )
    ))
    y <- datatable(cbind('&oplus;', is[1:ntop,c("region", "bin.event", "locus", "locus_overlap", "b", "bjs", "ja", "jl", "bin", "feature", "b.logfc", "b.fdr", "bjs.lr", "bjs.fdr", "bjs.nonuniformity", "bjs.inclussion", "bjs.inclussion_sign", "a.lr", "a.fdr", "a.nonuniformity", "a.dpir", "a.dpir_sign", "l.lr", "l.fdr", "l.participation", "l.dparticipation")]),
              rownames = FALSE,
              escape = -1,
              filter ="top",
              fillContainer = FALSE,
              extensions = c('Buttons', 'KeyTable'), 
              options = list(dom = 'lfrtBip',
                             buttons = c('copy', 'csv', 'excel', 'pdf', 'print', I('colvis')),
                             columnDefs = list(
                               #  list(visible = FALSE, targets = c(0, 2, 3)),
                               list(orderable = FALSE, className =
                                      'details-control', targets = 0)
                             ),
                             keys = TRUE
              ),   caption = htmltools::tags$caption(
                style = 'caption-side: top; text-align: left;',
                htmltools::h1(titulo), htmltools::h4(paste("Filters:", paste(paste(names(filters), filters, sep="="), collapse="; ")))
              ), container = sketch,
              callback = JS("
            table.order([10, 'asc']).draw();
            table.column(0).nodes().to$().css({cursor: 'pointer'});
            var format = function(d) {
             return '<div><h4>Gene view</h4></br><img src=\"img/' + d[1].replace(\":\", \"_\").replace(\"-\", \"_\") + '_gene.png\" height=\"700\"></img></div>';
            };
            table.on('click', 'td.details-control', function() {
               var td = $(this), row = table.row(td.closest('tr'));
               if (row.child.isShown()) {
                  row.child.hide();
                  td.html('&oplus;');
              } else {
                  row.child(format(row.data())).show();
                  td.html('&CircleMinus;');
              }
           });")
    )    
    
    ffile <- paste0(normalizePath(output.dir), "/integratedSignals.html")
    suppressWarnings(saveWidget(y, file = ffile, title = paste(names(sr@contrast)[sr@contrast != 0], collapse = " - "), selfcontained = F))
    if(openInBrowser == T) browseURL(ffile)
    
  }
)

setGeneric( name = 'filterSignals',
            def = function( sr,
                            bin.FC = 3,
                            bin.fdr = 0.05,
                            nonunif=1,
                            bin.inclussion = 0.1,
                            bjs.inclussion = 0.2,
                            bjs.fdr = 0.1,
                            a.inclussion = 0.3,
                            a.fdr = 0.05,
                            l.inclussion = 0.3,
                            l.fdr = 0.05, 
                            bDetectionSummary=FALSE  ) 
              standardGeneric( 'filterSignals' ))

setMethod( f = 'filterSignals',
           signature = c( 'ASpliSplicingReport' ),
           definition = function( sr,
                                  bin.FC = 3,
                                  bin.fdr = 0.05,
                                  nonunif=1,
                                  bin.inclussion = 0.1,
                                  bjs.inclussion = 0.2,
                                  bjs.fdr = 0.1,
                                  a.inclussion = 0.3,
                                  a.fdr = 0.05,
                                  l.inclussion = 0.3,
                                  l.fdr = 0.05, 
                                  bDetectionSummary=FALSE ) {
             .filterSignals( sr, 
                             bin.FC, 
                             bin.fdr, 
                             nonunif, 
                             bin.inclussion, 
                             bjs.inclussion,
                             bjs.fdr,
                             a.inclussion,
                             a.fdr,
                             l.inclussion,
                             l.fdr, 
                             bDetectionSummary ) } )

# ---------------------------------------------------------------------------- #
# plotBins
setGeneric( name = "plotBins",
            def = function(
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
              outfolder       = NULL,
              outfileType     = c( 'png', 'bmp', 'jpeg', 'tiff', 'pdf')[1],
              deviceOpt       = NULL ) standardGeneric( 'plotBins' ) )

setMethod( 
  f = "plotBins",
  signature = 'ASpliCounts',
  definition = function( 
    counts, 
    as,
    bin, 
    factorsAndValues, 
    targets,
    main            = NULL,
    colors          = c( '#2F7955', '#79552F', '#465579', '#A04935', 
                         '#752020', '#A07C35') ,
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
    outfolder       = NULL,
    outfileType     = c( 'png', 'bmp', 'jpeg', 'tiff', 'pdf')[1],
    deviceOpt       = NULL ) {
    
    .plotBins( counts, as, bin, factorsAndValues, targets, main, colors, 
               panelTitleColors, panelTitleCex, innerMargins, outerMargins, 
               useBarplots, barWidth, barSpacer, las.x, useHCColors, legendAtSide,
               outfolder, outfileType, deviceOpt )
  } 
)
# ---------------------------------------------------------------------------- # 
# ---------------------------------------------------------------------------- #
# plotGenomicRegions
setGeneric( 
  name = "plotGenomicRegions", 
  def = function( 
    #counts,
    features,
    x, 
    genomeTxDb, 
    targets, 
    xIsBin = TRUE, 
    layout = 'auto', 
    colors = 'auto', 
    plotTitles = 'auto', 
    sashimi = FALSE, 
    zoomOnBins= FALSE, 
    deviceOpt = NULL, 
    highLightBin = TRUE, 
    outfolder = NULL, 
    outfileType = 'png', 
    mainFontSize = 24, 
    annotationHeight = 0.2,
    annotationCol = 'black', 
    annotationFill = 'gray', 
    annotationColTitle = 'black',
    preMergedBAMs = NULL,
    useTransparency = FALSE,
    tempFolder = 'tmp',
    avoidReMergeBams = FALSE,
    verbose = TRUE ) standardGeneric( "plotGenomicRegions" ) )

setMethod(
  f = "plotGenomicRegions",
  signature = "ASpliFeatures",
  definition = function ( 
    #        counts,
    features, 
    x, 
    genomeTxDb, 
    targets, 
    xIsBin = TRUE, 
    layout = 'auto', 
    colors = 'auto', 
    plotTitles = 'auto', 
    sashimi = FALSE, 
    zoomOnBins= FALSE, 
    deviceOpt = NULL, 
    highLightBin = TRUE, 
    outfolder = NULL, 
    outfileType = 'png', 
    mainFontSize = 24, 
    annotationHeight = 0.2, 
    annotationCol = 'black', 
    annotationFill = 'gray', 
    annotationColTitle = 'black',
    preMergedBAMs = NULL,
    useTransparency = FALSE,
    tempFolder = 'tmp',
    avoidReMergeBams = FALSE,
    verbose = TRUE ) {
    
    .Deprecated("", msg = "plotGenomicRegions is deprecated and is no longer needed. See ASpli vignette for new pipeline.")
    
    .plotGenomicRegions(
      x, 
      genomeTxDb, 
      #              counts,
      features,
      targets, 
      xIsBin, 
      layout, 
      colors, 
      plotTitles, 
      sashimi, 
      zoomOnBins, 
      deviceOpt, 
      highLightBin, 
      outfolder, 
      outfileType,
      mainFontSize, 
      annotationHeight, 
      annotationCol, 
      annotationFill, 
      annotationColTitle,
      preMergedBAMs,
      useTransparency,
      tempFolder,
      avoidReMergeBams ,
      verbose )
  }
)
# ---------------------------------------------------------------------------- #
#agrego para no romper el old pipeline
# TODO:  Es necesario agregar todos los parametros con valores por default en
# la firma del metodo ? 
setGeneric (
  name = "DUreport",
  def = function( counts, 
                  targets, 
                  minGenReads  = 10,
                  minBinReads  = 5,
                  minRds = 0.05,
                  offset = FALSE,
                  offsetAggregateMode = c( "geneMode", "binMode" )[1],
                  offsetUseFitGeneX = TRUE,
                  contrast = NULL,
                  forceGLM = FALSE,
                  ignoreExternal = TRUE,
                  ignoreIo = TRUE, 
                  ignoreI = FALSE,
                  filterWithContrasted = FALSE,
                  verbose = FALSE
  ) standardGeneric("DUreport") )

#setGeneric (
#  name = "DUreport_DEXSeq",
#  def = function ( counts, ... ) standardGeneric("DUreport_DEXSeq") )

setMethod(
  f = "DUreport",
  signature = "ASpliCounts",
  definition = function( counts, 
                         targets, 
                         minGenReads  = 10,
                         minBinReads  = 5,
                         minRds = 0.05,
                         offset = FALSE,
                         offsetAggregateMode = c( "geneMode", "binMode" )[1],
                         offsetUseFitGeneX = TRUE,
                         contrast = NULL,
                         forceGLM = FALSE,
                         ignoreExternal = TRUE,
                         ignoreIo = TRUE, 
                         ignoreI = FALSE,
                         filterWithContrasted = FALSE,
                         verbose = FALSE
  ) { 
    .Deprecated(c("DUreport.offset", "DUreport.norm"))
    .DUreport( counts, targets, minGenReads, minBinReads, minRds, offset, 
               offsetAggregateMode, offsetUseFitGeneX, contrast, forceGLM,
               ignoreExternal, ignoreIo, ignoreI, filterWithContrasted, verbose  )
  }
)
