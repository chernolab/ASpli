.counterGenes <- function( reads, feature ) {
  cores = 1
  hits <- mclapply( reads, mc.cores = cores, function( x ) { 
        co <- countOverlaps( feature, x, ignore.strand = TRUE ) 
        gc()
        return(co)
      } )  
  
  # Create result dataframe
  effectiveGeneLength <- sum( width( feature ) )
  geneStarts          <- sapply( start( feature ), min )
  geneEnds            <- sapply( end( feature ), max)
  geneWidths          <- geneEnds - geneStarts + 1
  strand              <- as.character( unlist(runValue( strand( feature ) ) ))
  
  result <- data.frame(
      symbol = feature@elementMetadata$symbol , 
      locus_overlap = feature@elementMetadata$locus_overlap , 
      gene_coordinates = feature@elementMetadata$gene_coordinates , 
      start = geneStarts ,
      end = geneEnds , 
      length = geneWidths , 
      effective_length = effectiveGeneLength)
  #Le ponemos los hits despues porque si los nombres empiezan con un numero se rompen los nombres de las columnas
  for(i in 1:length(hits)){
    result[, names(hits)[i]] <- hits[[i]]
  }
  
  #ACH 20190520: reassure rownames are not lost
  rownames(result) <- names(hits[[1]])
  
  return(result)
}

.counterBin <- function( reads, feature, genes  ) { 
  cores = 1
  hits <- mclapply( reads, mc.cores = cores, function(x) {
        co <- countOverlaps( feature, x, ignore.strand = TRUE)
        gc()
        return(co)
      } )
  
  # Create result dataframe
  binToGeneIndex  <- match( feature@elementMetadata$locus, rownames(genes) )
  geneCoordinates <- genes$gene_coordinates[ binToGeneIndex ]
  
  result <- data.frame( feature@elementMetadata[c ('feature', 'event', 'locus',
              'locus_overlap', 'symbol')], 
      geneCoordinates, 
      as.data.frame( feature@ranges )[ ,c("start", "end", "width")])
  
  #Le ponemos los hits despues porque si los nombres empiezan con un numero se rompen los nombres de las columnas
  for(i in 1:length(hits)){
    result[, names(hits)[i]] <- hits[[i]]
  }
  
  colnames(result)[1:9] <- c( "feature","event","locus","locus_overlap","symbol",
      "gene_coordinates","start","end","length" )
  
  #ACH 20190520: reassure rownames are not lost
  rownames(result) <- names(hits[[1]])
  
  return(result)
}

# ---------------------------------------------------------------------------- #
# .counterJbin Cuenta reads que atraviesan dos bins. 
# TODO: el argumento l es la longitud de una read.  Que pasa con datos que 
#       tienen un tamano de read variable ?
.counterJbin <- function(reads, feature, genes, l) {
  cores=1
  ungapped <- mclapply( reads, mc.cores = cores, function(x) { x[ njunc( x ) == 0 , ] } )
  
  hits <- mclapply( ungapped, mc.cores = cores, function(x) { 
        co <- countOverlaps( feature, x, ignore.strand = TRUE,  minoverlap = l) 
        gc()
        return(co)
      })  
  
  jbinToGeneIndex <- match( feature@elementMetadata$locus, rownames(genes) )
  
  gene_coordinates <- genes$gene_coordinates[jbinToGeneIndex]
  
  result <- data.frame(
      feature@elementMetadata[ c( 'event', 'locus', 'locus_overlap', 'symbol') ], 
      gene_coordinates, 
      as.data.frame(feature@ranges))
  
  #Le ponemos los hits despues porque si los nombres empiezan con un numero se rompen los nombres de las columnas
  for(i in 1:length(hits)){
    result[, names(hits)[i]] <- hits[[i]]
  }
  
  result$names <- NULL
  
  colnames(result)[1:8] <- c( "event", "locus", "locus_overlap", "symbol",
      "gene_coordinates", "start", "end", "length" )
  
  #ACH 20190520: reassure rownames are not lost
  rownames(result) <- names(hits[[1]])
  
  return(result)
  
}

.getGeneGRanges <- function( aspliFeatures ) {
  
  feature <- featuresg( aspliFeatures )
  
  aggregate_first <- function (data, by){
    d = b = NULL # due to NSE notes in R CMD check
    data <- data.table(d=data, b=by) 
    ans  <- data[,list(A = first(d)), by = b]
    return(ans$A)
  }
  
  # Search junctions within genes
  unlistedFeatures <- unlist(feature)
  
  geneAndChr <- paste( names(unlistedFeatures) , as.character( seqnames(unlistedFeatures)), sep="_" )
  geneChr <- aggregate_first( as.data.frame(seqnames(unlistedFeatures))[, 1], by = geneAndChr)
  
  geneStarts <- aggregate( as.data.frame(start(unlistedFeatures)), by = list(geneAndChr), FUN = min)[,2]
  geneEnds <- aggregate( as.data.frame(end(unlistedFeatures)), by = list(geneAndChr), FUN = max)[,2]
  strand <- aggregate_first( as.data.frame(strand(unlistedFeatures))[, 1], by = geneAndChr)
  
  geneCoordinates <- rep( feature@elementMetadata$gene_coordinates ,table(names(unlistedFeatures)))

  geneCoordinates <- aggregate_first( geneCoordinates, by = geneAndChr)
  
  symbols <- rep( feature@elementMetadata$symbol ,table(names(unlistedFeatures)))
  symbols <- aggregate_first( symbols, by = geneAndChr)
  
  geneNames <- aggregate_first( as.data.frame(names(unlistedFeatures),stringsAsFactors=FALSE)[, 1] , by = geneAndChr)
  
  genes <- GRanges(
      seqnames = geneChr ,
      strand = strand,
      ranges = IRanges(geneStarts,geneEnds), 
      gene_coordinates = geneCoordinates,
      symbol=symbols )
  
  mcols(genes)$names <- geneNames
  
  return ( genes )
}

.getJunctionOverlapGeneData <- function(jranges, genes) {
  
  hitGen <- rep("-", length(jranges))
  hitGenStrand <- rep("*", length(jranges))
  gene_coordinates <- rep("-", length(jranges))
  ambiguos <- rep("-", length(jranges))
  symbol <- rep("-", length(jranges))    
  
  overGene <- findOverlaps(jranges, genes, type="within")
  overGeneDF <- as.data.frame(overGene)
  posJrange <- overGeneDF$queryHits 
  posGene <- overGeneDF$subjectHits
  overGeneDF$queryHits <- names(jranges)[as.numeric(overGeneDF$queryHits)]
  overGeneDF$subjectHits <- mcols(genes)$names[as.numeric(overGeneDF$subjectHits)]
  table <- table(overGeneDF$queryHits)

# BUG FIX: aggregate fails with 0-rows dfCountsStart. 
  if ( nrow( overGeneDF  ) > 0 ) {
    ttG <- data.frame(aggregate(subjectHits ~ queryHits, data = overGeneDF, paste, collapse=";"))
  } else {
    ttG <- data.frame( names = character(0) ) 
    for ( i in 1:ncol( overGeneDF ) ) {
      ttG[, i+1] <- integer(0)
    }
    colnames( ttG )[2:ncol(ttG)] <- colnames( overGeneDF )
  }
  
  dd0 <- match(ttG$queryHits, names(jranges))
  hitGen[dd0] <- ttG$subjectHits
  dd <- match(ttG$queryHits, names(table) )
  ttG$undef <- table[dd]
  ttG$tag <- rep("-", nrow(ttG))
  ttG$tag[ttG$undef>1] <- "yes"
  
  ambiguos[dd0] <- ttG$tag
  
  hitGen[posJrange] <-  mcols(genes)$names[posGene]
  hitGen[-posJrange] <- "noHit"
  hitGenStrand[posJrange] <- as.character(strand(genes)[posGene])
  gene_coordinates[posJrange] <- mcols(genes)$gene_coordinates[posGene] 
  symbol[posJrange] <- as.character(genes@elementMetadata$symbol[posGene])
  
  return( list( ambiguos, hitGen, hitGenStrand, gene_coordinates, symbol ) )
}

.getJunctionMatchingBins <- function( features, jranges ) {
  hitBin <- rep("-", length(jranges))
  annJunctions <- featuresj(features)
  overJ <- findOverlaps(jranges, annJunctions, type="equal") #identify annotated junctions
  overJDF <- as.data.frame(overJ) #get a df
  namesJ <- as.numeric(overJDF[,1])  #get index of jrangs that hit against annJunctions
  namesAnnJ <- as.numeric(overJDF[,2]) #get index of annJunctions thta hit against jranges
  #jname[namesJ] <- names(jranges[namesJ])
  hitBin[namesJ] <- names(annJunctions[namesAnnJ]) #ok, metadata vector
  hitBin[-namesJ] <- "noHit" #ok,  metadata vector. 
  return(hitBin)
}

.getJunctionSpanningBins <- function( jranges, exonsBins) {
  over <- findOverlaps(jranges, exonsBins)
  overDF <- as.data.frame(over)
  namesJ <- as.numeric(overDF[,1])
  overDF[,1] <- names(jranges[namesJ])
  namesBins <- as.numeric(overDF[,2])
  overDF[,2] <- names(exonsBins[namesBins])
  
# BUG FIX: aggregate fails with 0-rows dfCountsStart. 
  if ( nrow( overDF  ) > 0 ) {
    tt <- data.frame(aggregate(subjectHits ~ queryHits, data = overDF, paste, collapse=";")) 
  } else {
    tt <- data.frame( names = character(0) ) 
    for ( i in 1:ncol( overDF ) ) {
      tt[, i+1] <- integer(0)
    }
    colnames( tt )[2:ncol(tt)] <- colnames( overDF )
  }
  
  span <- rep("-", length(jranges))
  te <- match(names(jranges), tt$queryHits)
  span[ ! is.na(te) ] <- tt$subjectHits[te][ ! is.na(te) ]
  return(span)
}

.getJunctionWithinBins <- function(jranges, exonsBins) {
  j_within_bin <- rep("-", length(jranges))
  overJunctionWithinBins <- findOverlaps(jranges, exonsBins, type="within")
  if (length (overJunctionWithinBins ) > 0) {
    overJunctionWithinBinsDF <- as.data.frame(overJunctionWithinBins)
    namesJ <- as.numeric(overJunctionWithinBinsDF[,1])
    namesB <- as.numeric(overJunctionWithinBinsDF[,2])
    overJunctionWithinBinsDF[,1] <- names(jranges[namesJ])
    overJunctionWithinBinsDF[,2] <- names(exonsBins[namesB])
    agtt <- data.frame( aggregate( subjectHits ~ queryHits,
            data = overJunctionWithinBinsDF, paste, collapse=";")) 
    tw <- match(names(jranges), agtt$queryHits) 
    j_within_bin[ ! is.na(tw) ] <- agtt$subjectHits[tw][ ! is.na(tw) ]
  }
  return( j_within_bin )
}

.ovBinJunction <- function( features, jranges ) {
  
# Creates GRanges for genes, including repeating genes in different
# chromosomes
  genes <- .getGeneGRanges( features )
  
  # Looks for data of junction overlapping genes
  ovGene <- .getJunctionOverlapGeneData( jranges, genes )
  ambiguos <- ovGene[[1]]
  hitGen <- ovGene[[2]]
  hitGenStrand <- ovGene[[3]]
  gene_coordinates <- ovGene[[4]]
  symbol <- ovGene[[5]]
  
  # Search junction equal to annotated junctions
  hitBin <- .getJunctionMatchingBins(features, jranges)
  
  # Search junctions spanning bins
  exonsBins <- featuresb(features)[featuresb(features)@elementMetadata$feature=="E",]
  span <- .getJunctionSpanningBins(jranges, exonsBins )
  
  # Search junctions within exons!
  j_within_bin <- .getJunctionWithinBins(jranges, exonsBins)

  mcols(jranges) <- append( mcols(jranges), DataFrame(
          hitBin=hitBin, 
          hitGen=hitGen, 
          hitGenStrand=hitGenStrand,
          gene_coordinates=gene_coordinates,
          undef=ambiguos,
          bin_spanned=span,
          j_within_bin=j_within_bin,
          symbol=symbol))
  return(jranges)
}  

.counterJunctions <- function(features, bam, maxISize) {
  cores=1
    ujunctions <- mclapply (bam, mc.cores = cores, function(x)    {  
          junctions <- unlist(junctions(x) )
          strand(junctions) <- "*"
          start(junctions) <- start(junctions)-1
          end(junctions) <- end(junctions)+1
          ujunctions <- unique(junctions)
          gc()
          return(ujunctions)   
        } )
  jranges <- unique( unlist( GRangesList( unlist( ujunctions ) ) ) )
  
  # Filter junctions by Intron size
  maxWidth <-  maxISize+2
  jranges <- jranges[width(jranges)<= maxISize]
  
  #Here I summarize hits agains the element
  fcoord <- paste(seqnames(jranges), 
      start(jranges),
      end(jranges) , sep="." )
  jranges@ranges@NAMES <- fcoord

  jcounts <- mclapply(bam, mc.cores = cores, function(x) {
        junctions <- unlist(junctions(x) )
        strand(junctions)<- "*"
        start(junctions) <- start(junctions)-1
        end(junctions) <- end(junctions)+1
        count <- countMatches(jranges, junctions)    
        jc <- data.frame(row.names=names(jranges), count)
        gc()
        return(jc)
      })
  
  df <- do.call("cbind", jcounts)
  colnames(df) <- names(jcounts)

  jranges <- .ovBinJunction(features, jranges)
  
  jrdf <- data.frame(
      as.data.frame(jranges@elementMetadata$hitBin), 
      as.data.frame(jranges@elementMetadata$hitGen), 
      as.data.frame(jranges@elementMetadata$hitGenStrand),
      as.data.frame(jranges@elementMetadata$undef),
      as.data.frame(jranges@elementMetadata$symbol), 
      as.data.frame(jranges@elementMetadata$gene_coordinates),
      as.data.frame(jranges@elementMetadata$bin_spanned),
      as.data.frame(jranges@elementMetadata$j_within_bin),
      row.names=names(jranges)  )
  
  colnames(jrdf) <- c("junction", 
      "gene", 
      "strand",
      "multipleHit",
      "symbol",
      "gene_coordinates",
      "bin_spanned",
      "j_within_bin")
  
  rownames(jrdf) <- names(jranges)
  aa <- merge(jrdf, df, by.x="row.names", by.y="row.names", sort=FALSE)
  rnames <- paste(start(jranges)-1,end(jranges)+1, sep="-" )
  rownames(aa) <- fcoord
  aa$Row.names <- NULL
  return(aa)
}
