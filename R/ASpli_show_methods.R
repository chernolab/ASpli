# ---------------------------------------------------------------------------- #
# Getters and setters for ASpliFeatures
setGeneric( name = "featuresg", def = function( x ) standardGeneric("featuresg") )
setMethod( f = "featuresg", signature = "ASpliFeatures",
            definition = function(x){  x@genes } )

setGeneric( name = "featuresg<-", def = function( x, value ) standardGeneric("featuresg<-") )
setReplaceMethod( f = "featuresg", signature = c( "ASpliFeatures", 'GRangesList' ),
    definition = function( x , value){  x@genes <- value; return( x ) } )

setGeneric( name = "featuresj", def = function( x ){ standardGeneric("featuresj") } )
setMethod( f = "featuresj", signature = "ASpliFeatures", 
    definition = function(x){ x@junctions })

setGeneric( name = "featuresj<-", def = function( x, value ){ standardGeneric("featuresj<-") } )
setReplaceMethod( f = "featuresj", signature = c("ASpliFeatures", "GRanges"), 
    definition = function( x, value ){ x@junctions <- value; return ( x )})

setGeneric( name = "featuresb", def = function( x ) standardGeneric("featuresb"))
setMethod( f = "featuresb",signature = "ASpliFeatures",
           definition = function( x ){ x@bins })

setGeneric( name = "featuresb<-", def = function( x, value ) standardGeneric("featuresb<-"))
setReplaceMethod( f = "featuresb",signature = c( "ASpliFeatures","GRanges" ),
   definition = function( x, value ){ x@bins <- value; return( x ) } )

setGeneric( name = "transcriptExons", def = function( x ) standardGeneric("transcriptExons"))
setMethod( f = "transcriptExons",signature = "ASpliFeatures",
           definition = function( x ){ x@transcriptExons })

setGeneric( name = "transcriptExons<-", def = function( x, value ) standardGeneric("transcriptExons<-"))
setReplaceMethod( f = "transcriptExons",signature = c( "ASpliFeatures","GRangesList" ),
                  definition = function( x, value ){ x@transcriptExons <- value; return( x ) } )
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Getters and Setters for ASpliCounts
setGeneric( name= "countsg", def=function( x ) standardGeneric("countsg") )
setMethod( f = "countsg", signature = "ASpliCounts",
    definition = function( x ){ x@gene.counts} )
       
setGeneric( name= "countsg<-", def=function( x, value ) standardGeneric("countsg<-") )
setReplaceMethod( f = "countsg", signature = c( "ASpliCounts", "data.frame") ,
   definition = function( x, value ){ x@gene.counts <- value; return ( x )} )
       
setGeneric ( name = "countsj", def = function( x ) standardGeneric("countsj") )
setMethod (f = "countsj", signature = "ASpliCounts", 
    definition = function( x ){ x@junction.counts })

setGeneric ( name = "countsj<-", def = function( x, value ) standardGeneric("countsj<-") )
setReplaceMethod (f = "countsj", signature = c( "ASpliCounts","data.frame"), 
    definition = function( x, value ){ x@junction.counts <- value ; return (x) })
       
setGeneric( name = "countsb", def = function( x ) standardGeneric("countsb") )
setMethod( f = "countsb", signature = "ASpliCounts", 
    definition = function( x ){ x@exon.intron.counts } )

setGeneric( name = "countsb<-", def = function( x, value ) standardGeneric("countsb<-") )
setReplaceMethod( f = "countsb", signature = c("ASpliCounts","data.frame"), 
    definition = function( x, value ){ x@exon.intron.counts <- value; return( x ) } )

setGeneric( name = "countse1i", def = function( x ) standardGeneric("countse1i") )
setMethod (f = "countse1i", signature = "ASpliCounts", 
    definition = function( x ){ x@e1i.counts })

setGeneric( name = "countse1i<-", def = function( x, value ) standardGeneric("countse1i<-") )
setReplaceMethod (f = "countse1i", signature = c("ASpliCounts", "data.frame"), 
    definition = function( x, value ){ x@e1i.counts <- value; return( x ) } )

setGeneric( name = "countsie2", def = function( x ) standardGeneric("countsie2"))
setMethod( f = "countsie2", signature = "ASpliCounts", 
    definition = function( x ){ x@ie2.counts })

setGeneric( name = "countsie2<-", def = function( x, value ) standardGeneric("countsie2<-"))
setReplaceMethod( f = "countsie2", signature = c("ASpliCounts","data.frame"), 
    definition = function( x, value ){ x@ie2.counts <- value; return( x ) })

setGeneric( name = "rdsg", def = function( x ) standardGeneric("rdsg"))
setMethod( f = "rdsg", signature = "ASpliCounts",
    definition = function( x ){ x@gene.rd })

setGeneric( name = "rdsg<-", def = function( x, value ) standardGeneric("rdsg<-"))
setReplaceMethod( f = "rdsg", signature = c("ASpliCounts","data.frame"),
    definition = function( x, value ){ x@gene.rd <- value; return( x )})

setGeneric( name = "rdsb", def = function( x ) standardGeneric("rdsb") )
setMethod( f = "rdsb", signature = "ASpliCounts", 
    definition = function( x ){ x@bin.rd })

setGeneric( name = "rdsb<-", def = function( x, value ) standardGeneric("rdsb<-") )
setReplaceMethod( f = "rdsb", signature = c("ASpliCounts","data.frame"), 
    definition = function( x, value ){ x@bin.rd <- value; return ( x )})

setGeneric( name = "condition.order", def = function( x ) standardGeneric("condition.order") )
setMethod( f = "condition.order", signature = "ASpliCounts", 
           definition = function( x ){ x@condition.order })

setGeneric( name = "targets", def = function( x ) standardGeneric("targets") )
setMethod( f = "targets", signature = "ASpliCounts", 
           definition = function( x ){ x@targets })
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Getters and Setters for ASpliAS
setGeneric( name =  "irPIR", def = function( x ) standardGeneric("irPIR") )
setMethod( f = "irPIR", signature =  "ASpliAS", 
    definition = function( x ){ x@irPIR } )

setGeneric( name =  "irPIR<-", def = function( x, value ) standardGeneric("irPIR<-") )
setReplaceMethod( f = "irPIR", signature = c("ASpliAS",'data.frame'), 
    definition = function( x, value ){ x@irPIR <- value; return( x ) } )       

setGeneric( name =  "altPSI", def = function( x ) standardGeneric("altPSI"))
setMethod(f = "altPSI", signature = "ASpliAS", 
    definition = function( x ) { x@altPSI })

setGeneric( name =  "altPSI<-", def = function( x, value ) standardGeneric("altPSI<-"))
setReplaceMethod( f = "altPSI", signature = c("ASpliAS","data.frame"), 
    definition = function( x, value ) { x@altPSI <- value; return( x ) })

setGeneric( name = "esPSI", def = function( x ) standardGeneric("esPSI"))
setMethod( f = "esPSI", signature = "ASpliAS",
    definition = function( x ) { x@esPSI })

setGeneric( name = "esPSI<-", def = function( x, value ) standardGeneric("esPSI<-"))
setReplaceMethod( f = "esPSI", signature = c("ASpliAS","data.frame"),
    definition = function( x, value ) { x@esPSI <- value; return( x ) })

setGeneric( name = "junctionsPIR", def = function( x ) standardGeneric("junctionsPIR"))
setMethod( f ="junctionsPIR", signature = "ASpliAS", 
    definition = function( x ){ x@junctionsPIR })

setGeneric( name = "junctionsPIR<-", def = function( x, value ) standardGeneric("junctionsPIR<-"))
setReplaceMethod( f ="junctionsPIR", signature = c("ASpliAS","data.frame"), 
    definition = function( x, value){ x@junctionsPIR <- value; return( x ) })

setGeneric( name = "junctionsPJU", def = function( x ) standardGeneric("junctionsPJU"))
setMethod( f = "junctionsPJU", signature = "ASpliAS",
    definition = function( x ){ x@junctionsPJU } )

setGeneric( name = "junctionsPJU<-", def = function( x, value ) standardGeneric("junctionsPJU<-"))
setReplaceMethod( f = "junctionsPJU", signature = c("ASpliAS","data.frame"),
    definition = function( x, value ){ x@junctionsPJU <- value; return( x )} )

setGeneric( name = "joint", def = function( x ) standardGeneric("joint") )
setMethod( f = "joint", signature = "ASpliAS", 
    definition = function( x ){ x@join })

setGeneric( name = "joint<-", def = function( x, value ) standardGeneric("joint<-") )
setReplaceMethod( f = "joint", signature = c("ASpliAS","data.frame"), 
    definition = function( x, value ){ x@join <- value; return( x ) })
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Getters and Setters for ASpliDU
setGeneric( name = "genesDE", def = function( x ) standardGeneric("genesDE"))
setMethod( f = "genesDE", signature = "ASpliDU", 
    definition = function( x ){ x@genes })

setGeneric( name = "genesDE<-", def = function( x, value ) standardGeneric("genesDE<-"))
setReplaceMethod( f = "genesDE", signature = c("ASpliDU","data.frame"), 
    definition = function( x,value ){ x@genes <- value; return( x )})

setGeneric( name = "binsDU", def = function( x ) standardGeneric("binsDU"))
setMethod( f = "binsDU", signature = "ASpliDU", 
    definition = function( x ){ x@bins })

setGeneric( name = "binsDU<-", def = function( x, value ) standardGeneric("binsDU<-"))
setReplaceMethod( f = "binsDU", signature = "ASpliDU", 
    definition = function( x, value ){ x@bins <- value; return( x ) })

setGeneric(  name = "junctionsDU", def = function( x ) standardGeneric("junctionsDU") )
setMethod( f = "junctionsDU", signature = "ASpliDU", 
    definition = function( x ){ x@junctions })

setGeneric(  name = "junctionsDU<-", def = function( x, value ) standardGeneric("junctionsDU<-") )
setReplaceMethod( f = "junctionsDU", signature = c( "ASpliDU", 'data.frame'), 
    definition = function( x, value ){ x@junctions <- value; return( x ) } )
#-----------------------------------------------------------------------

setGeneric( name = "filters", def = function( x ) standardGeneric("filters"))
setMethod( f = "filters", signature = "ASpliIntegratedSignals",
           definition = function( x ){ x@filters } )

setGeneric( name = "filters<-", def = function( x, value ) standardGeneric("filters<-"))
setReplaceMethod( f = "filters", signature = c("ASpliIntegratedSignals","data.frame"),
                  definition = function( x, value ){ x@filters <- value; return( x )} )

setGeneric( name = "signals", def = function( x ) standardGeneric("signals") )
setMethod( f = "signals", signature = "ASpliIntegratedSignals", 
           definition = function( x ){ x@signals })

setGeneric( name = "signals<-", def = function( x, value ) standardGeneric("signals<-") )
setReplaceMethod( f = "signals", signature = c("ASpliIntegratedSignals","data.frame"), 
                  definition = function( x, value ){ x@signals <- value; return( x ) })
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# show methods
setMethod( 'show', 'ASpliFeatures', function( object ) {
      cat("Object of class", class(object),"\n")
      cat("Genes: GRangesList of length ", 
          length(object@genes),
          "Access using featuresg(object)", "\n")
      cat("Bins: GRanges of length", 
          length(object@bins),
          "Access using featuresb(object)", "\n")
      cat("Junctions: GRanges of length", 
          length(object@junctions),
          "Access using featuresj(object)","\n")
      cat("Transcript Exons: GRanges of length", 
          length(object@transcriptExons),
          "Access using transcriptExons(object)","\n")
    })

setMethod( 'show', 'ASpliCounts', function( object ) {
      cat("Object of class", class(object),"\n")
      cat("Gene counts:", 
          dim(object@gene.counts)[1], "genes analysed.",
          "Access using countsg(object)", "\n")
      cat("Gene RD:", 
          dim(object@gene.rd)[1], "genes analysed.",
          "Access using rdsg(object)", "\n")
      cat("Bin counts:", 
          dim(object@exon.intron.counts)[1], "bins analysed.",
          "Access using countsb(object)", "\n")
      cat("Bin RD:", 
          dim(object@bin.rd)[1],"bins analysed.",
          "Access using rdsb(object)", "\n")
      cat("Junction counts:", 
          dim(object@junction.counts)[1], "junctions analysed.",
          "Access using countsj(object)", "\n")
    })
        
setMethod( 'show', 'ASpliAS', function( object )  {
      cat("Object of class", class(object),"\n")
      cat("IR PIR: ", 
          dim(object@irPIR)[1], "intron bins analysed.",
          sum(!is.na(object@irPIR$J3)), "intron bins passed the filters.",
          "Access using irPIR(object)", "\n")
      cat("ES PSI:", 
          dim(object@esPSI)[1], "exon bins analysed.",
          sum(!is.na(object@esPSI$J3)), "exon bins passed the filters.",
          " Access using esPSI(object)", "\n")
      cat("AltSS PSI:", 
          dim(object@altPSI)[1], "exon bins analysed.",
          sum(!is.na(object@altPSI$J3)), "exon bins passed the filters.",
          " Access using altPSI(object)", "\n")
      cat("Junctions PIR:", 
          dim(object@junctionsPIR)[1], "junctions analysed.",
          "Access using junctionsPIR(object)", "\n")
      cat("Junctions PJU:", 
          dim(object@junctionsPJU)[1], "junctions analysed.",
          "Access using junctionsPJU(object)", "\n")
    })

setMethod('show', 'ASpliDU', function( object ) {
      cat("Object of class", class(object),"\n")
      if ( containsGenesAndBins( object ) ) {
        cat("Gene DE:", 
            dim(object@genes)[1], "genes analysed.",
            "Access using genesDE(object)", "\n")
        
        cat("Bins DU:", 
            dim(object@bins)[1], "bins analysed.",
            "Access using binsDU(object)", "\n")
      }
      if ( containsJunctions( object ) ) {
        cat("Junctions DU:", 
            dim(object@junctions)[1],"junctions analysed.",
            "Access using junctionsDU(object)", "\n")
      }
    })

setMethod('show', 'ASpliJDU', function( object ) {
  cat("Object of class", class(object),"\n")
  
  cat("Locale Junctions:", 
      nrow(localej(object)), "junctions.",
      "Access using localej(object)", "\n")

  cat("Locale Clusters:", 
      nrow(localec(object)), "clusters",
      "Access using localec(object)", "\n")
  
  cat("Anchor Junctions:", 
      nrow(anchorj(object)), "junctions.",
      "Access using anchorj(object)", "\n")
  
  cat("Anchor Clusters:", 
      nrow(anchorc(object)), "clusters",
      "Access using anchorc(object)", "\n")
  
  cat("Intron Retention Junctions:",
      nrow(jir(object)), "junctions.",
      "Access using jir(object)", "\n")
  
  cat("Exon Skipping Junctions:", 
      nrow(jes(object)), "junctions.",
      "Access using jes(object)", "\n")
  
  cat("Alternative Splicing Junctions:", 
      nrow(jalt(object)), "junctions.",
      "Access using jalt(object)", "\n")
})

setMethod( 'show', 'ASpliSplicingReport', function( object ) {
  
  # my_table    <- function(data, jcol="junction"){
  #   binNA       <- is.na(data$bin.fdr)
  #   junctionNA  <- is.na(data$junction.fdr)
  #   
  #   binNNA      <- data[!binNA, ]
  #   binNA       <- data[binNA, ]
  #   junctionNNA <- data[!junctionNA, ]
  #   junctionNA  <- data[junctionNA, ]
  #   
  #   binT      <- binNNA$bin.fdr < 0.05
  #   binF      <- unique(binNNA[!binT, colnames(data) == jcol])
  #   binT      <- unique(binNNA[binT, colnames(data) == jcol])
  #   
  #   if("cluster.fdr" %in% colnames(data)){
  #     junctionT <- junctionNNA$junction.fdr < 0.05 | junctionNNA$cluster.fdr < 0.05
  #   }else{
  #     junctionT <- junctionNNA$junction.fdr < 0.05
  #   }
  #   junctionF <- unique(junctionNNA[!junctionT, colnames(data) == jcol])
  #   junctionT <- unique(junctionNNA[junctionT, colnames(data) == jcol])
  #   
  #   binNA     <- unique(binNA[, colnames(data) == jcol])
  #   binNNA    <- unique(binNNA[, colnames(data) == jcol])
  #   
  #   junctionNA  <- unique(junctionNA[, colnames(data) == jcol])
  #   junctionNNA <- unique(junctionNNA[, colnames(data) == jcol])
  #   
  #   m         <- matrix(c(length(intersect(binF, junctionF)), 
  #                         length(intersect(binF, junctionT)), 
  #                         length(intersect(binF, junctionNA)), 
  #                         length(intersect(binT, junctionF)), 
  #                         length(intersect(binT, junctionT)),
  #                         length(intersect(binT, junctionNA)),
  #                         length(intersect(binNA, junctionF)),
  #                         length(intersect(binNA, junctionT)),
  #                         length(intersect(binNA, junctionNA))),
  #                       nrow=3, ncol=3, byrow=T)
  #   rownames(m) <- c("Bin F", "Bin T", "Bin <NA>")                           
  #   colnames(m) <- c("Junction F", "Junction T", "Junction <NA>")
  #   return(m)
  # }
  
  if(FALSE){
    integratedSignals <- integrateSignals(object, asd = new(Class = "ASpliAS"))
    cat("Object of class", class(object),"\n")
    
    cat("Locale based: \n") 
    cat("De novo junctions: ", sum(is.na(integratedSignals$feature[integratedSignals$jl == 1])), "\n")
    cat("Annotated features:\n")
    print(table(integratedSignals[integratedSignals$jl == 1,c("feature","bin.event")]))
    cat("Access using localebased(object)\n\n")
    #print(table(bin=localebased(object)$bin.fdr < 0.05, junction=(localebased(object)$junction.fdr < 0.05 | localebased(object)$cluster.fdr < 0.05), useNA = "always"))
      
    cat("Anchor based: \n")
    cat("De novo junctions: ", sum(is.na(integratedSignals$feature[integratedSignals$ja == 1])), "\n")
    cat("Annotated features:\n")
    print(table(integratedSignals[integratedSignals$ja == 1,c("feature","bin.event")]))
    cat("Access using anchorbased(object)\n\n")
    #print(table(bin=anchorbased(object)$bin.fdr < 0.05, junction=(anchorbased(object)$junction.fdr < 0.05 | anchorbased(object)$cluster.fdr < 0.05), useNA = "always"))  
  
    cat("Bin based: \n") 
    cat("Annotated features:\n")  
    print(table(integratedSignals[integratedSignals$b == 1,c("feature","bin.event")]))
    #print(table(bin=binbased(object)$bin.fdr < 0.05, junction=binbased(object)$junction.fdr < 0.05, useNA = "always"))
    cat("Access using binbased(object)\n\n")  
  }else{
    cat("Object of class", class(object),"\n")
    cat("Locale based: \n") 
    cat("Access using localebased(object)\n\n")

    cat("Anchor based: \n")
    cat("Access using anchorbased(object)\n\n")

    cat("Bin based: \n") 
    cat("Access using binbased(object)\n\n")  
  }
})
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# Write methods
setGeneric ( name =  "writeCounts",
  def = function( counts,  output.dir = "counts" ) standardGeneric( "writeCounts" ) )

setMethod(
  f = "writeCounts",
  signature = "ASpliCounts",
  definition = function( counts, output.dir = "counts" ) {
    genesFile <- file.path( output.dir, "genes", "gene.counts.tab" )       
    exonsFile <- file.path( output.dir, "exons", "exon.counts.tab")       
    intronsFile    <- file.path( output.dir, "introns", "intron.counts.tab")
    intronsIE1File <- file.path( output.dir, "introns", "e1i.counts.tab")
    intronsEI2File <- file.path( output.dir, "introns", "ie2.counts.tab")
    junctionsFile <- file.path( output.dir, "junctions", "junction.counts.tab")              	     	     	     	     	    
    
    file.exists( output.dir ) || dir.create ( output.dir )
    for (folder in unique( lapply( c(genesFile, exonsFile, intronsIE1File, 
                intronsEI2File, junctionsFile ), 
            dirname ))) {
      dir.create( folder )
    }

    # Export genes 
    write.table(countsg(counts), genesFile, sep="\t", quote=FALSE, col.names=NA)

    # Export exons
    ec <- countsb(counts)[ countsb(counts)$feature=="E",]
    ec <- ec[ ec$event != "IR", ]
    ec <- ec[ ec$event != "IR*", ]
    write.table( ec, exonsFile, sep="\t", quote=FALSE, col.names=NA)

    # Export introns
    ic <- rbind( countsb( counts )[ countsb( counts )$feature == "I",], 
                 countsb( counts )[ countsb( counts )$feature == "Io",], 
                 countsb( counts )[ countsb( counts )$event == "IR",],
                 countsb( counts )[ countsb( counts )$event == "IR*",])
    write.table(ic, intronsFile, sep="\t", quote=FALSE, col.names=NA)
    write.table(countse1i(counts), intronsIE1File, sep="\t", quote=FALSE, 
        col.names=NA)
    write.table(countsie2(counts), intronsEI2File, sep="\t", quote=FALSE, 
        col.names=NA)
    
    # Export junctions 
    write.table(countsj(counts), junctionsFile, sep="\t", quote=FALSE, 
        col.names=NA)
  }
)

setGeneric ( name = "writeRds", 
    def = function( counts, output.dir="rds" ) standardGeneric( "writeRds" ) )

setMethod(
  f = "writeRds",
  signature = "ASpliCounts",
  definition = function( counts, output.dir="rds") {
    
    # Create output folder structure
    file.exists( output.dir ) || dir.create( output.dir )
    
    genesFile   <- file.path( output.dir, "genes", "gene.rd.tab" )       
    exonsFile   <- file.path( output.dir, "exons", "exon.rd.tab" )       
    intronsFile <- file.path( output.dir, "introns", "intron.rd.tab" )

    # Export gene read densities
    write.table( rdsg(counts), genesFile, sep="\t", quote = FALSE, 
        col.names = NA )

    # Export exon read densities
    erd <- rdsb(counts)[rdsb(counts)$feature == "E",]
    erd <- erd[erd$event != "IR",]
    erd <- erd[erd$event != "IR*",]
    write.table( erd, exonsFile, sep="\t", quote = FALSE, col.names = NA )

    # Export Intron read densities
    ird <- rbind( rdsb(counts)[rdsb(counts)$feature == "I",], 
                  rdsb(counts)[rdsb(counts)$feature == "Io",], 
                  rdsb(counts)[rdsb(counts)$eventJ  == "IR",])
    write.table( ird, intronsFile, sep="\t", quote = FALSE, col.names = NA )
  }
)

setGeneric ( name = "writeAll", 
    def = function( counts, du, as, output.dir="output" ) 
      standardGeneric( "writeAll" ) )

setMethod(
  f = "writeAll",
  signature = 'ASpliCounts',
  definition = function( counts, du, as, output.dir="output" ) {
    
    # Exports Counts, read densities, AS y DU
    suppressWarnings( writeCounts( counts, output.dir ) )
    suppressWarnings( writeRds( counts, output.dir ) )
    suppressWarnings( writeAS( as, output.dir ) )
    suppressWarnings( writeDU( du, output.dir ) )
    
    # Export as+du table
    colnames( irPIR( as ) ) <- colnames( altPSI( as ) )
    conP <- rbind( altPSI(as), esPSI(as), irPIR(as) )
    ii   <- match( rownames( binsDU( du ) ), row.names(conP) )
    bins.join <- data.frame( binsDU(du), conP[ii,])
    bins.join$feature <- NULL
    bins.join$event.1 <- NULL
    file <- file.path( output.dir, "bins_du_psi_pir.tab" )
    write.table( bins.join, file, sep="\t", quote=FALSE, col.names=NA )
    
    # Export Summary
    group <- colnames(countsg(counts))[ 8 : ncol( countsg(counts) ) ]
    summary <- bins.join[ , c(1:11) ]
    summary <- cbind( summary, bins.join[ , colnames( bins.join ) == levels( group ) ] )
    file <- file.path( output.dir, "summary.tab")
    write.table( summary, file, sep="\t", quote=FALSE, col.names = NA)  
    
  })
# End of write methods
# ---------------------------------------------------------------------------- # 
#
# ---------------------------------------------------------------------------- #
# Getters and Setters for ASpliJDU
setGeneric( name = "localej", def = function( x ) standardGeneric("localej"))
setMethod( f = "localej", signature = "ASpliJDU", 
           definition = function( x ){ x@localej })

setGeneric( name = "localej<-", def = function( x, value ) standardGeneric("localej<-"))
setReplaceMethod( f = "localej", signature = c("ASpliJDU","data.frame"), 
                  definition = function( x,value ){ x@localej <- value; return( x )})

setGeneric( name = "localec", def = function( x ) standardGeneric("localec"))
setMethod( f = "localec", signature = "ASpliJDU", 
           definition = function( x ){ x@localec })

setGeneric( name = "localec<-", def = function( x, value ) standardGeneric("localec<-"))
setReplaceMethod( f = "localec", signature = c("ASpliJDU","data.frame"), 
                  definition = function( x,value ){ x@localec <- value; return( x )})

setGeneric( name = "anchorj", def = function( x ) standardGeneric("anchorj"))
setMethod( f = "anchorj", signature = "ASpliJDU", 
           definition = function( x ){ x@anchorj })

setGeneric( name = "anchorj<-", def = function( x, value ) standardGeneric("anchorj<-"))
setReplaceMethod( f = "anchorj", signature = c("ASpliJDU","data.frame"), 
                  definition = function( x,value ){ x@anchorj <- value; return( x )})

setGeneric( name = "anchorc", def = function( x ) standardGeneric("anchorc"))
setMethod( f = "anchorc", signature = "ASpliJDU", 
           definition = function( x ){ x@anchorc })

setGeneric( name = "anchorc<-", def = function( x, value ) standardGeneric("anchorc<-"))
setReplaceMethod( f = "anchorc", signature = c("ASpliJDU","data.frame"), 
                  definition = function( x,value ){ x@anchorc <- value; return( x )})

setGeneric( name = "jir", def = function( x ) standardGeneric("jir"))
setMethod( f = "jir", signature = "ASpliJDU", 
           definition = function( x ){ x@jir })

setGeneric( name = "jir<-", def = function( x, value ) standardGeneric("jir<-"))
setReplaceMethod( f = "jir", signature = c("ASpliJDU","data.frame"), 
                  definition = function( x,value ){ x@jir <- value; return( x )})

setGeneric( name = "jes", def = function( x ) standardGeneric("jes"))
setMethod( f = "jes", signature = "ASpliJDU", 
           definition = function( x ){ x@jes })

setGeneric( name = "jes<-", def = function( x, value ) standardGeneric("jes<-"))
setReplaceMethod( f = "jes", signature = c("ASpliJDU","data.frame"), 
                  definition = function( x,value ){ x@jes <- value; return( x )})

setGeneric( name = "jalt", def = function( x ) standardGeneric("jalt"))
setMethod( f = "jalt", signature = "ASpliJDU", 
           definition = function( x ){ x@jalt })

setGeneric( name = "jalt<-", def = function( x, value ) standardGeneric("jalt<-"))
setReplaceMethod( f = "jalt", signature = c("ASpliJDU","data.frame"), 
                  definition = function( x,value ){ x@jalt <- value; return( x )})

# ---------------------------------------------------------------------------- #

#
# ---------------------------------------------------------------------------- #
# Getters and Setters for ASpliSplicingReport
setGeneric( name = "binbased", def = function( x ) standardGeneric("binbased"))
setMethod( f = "binbased", signature = "ASpliSplicingReport", 
           definition = function( x ){ x@binbased })

setGeneric( name = "binbased<-", def = function( x, value ) standardGeneric("binbased<-"))
setReplaceMethod( f = "binbased", signature = c("ASpliSplicingReport","data.frame"), 
                  definition = function( x,value ){ x@binbased <- value; return( x )})

setGeneric( name = "localebased", def = function( x ) standardGeneric("localebased"))
setMethod( f = "localebased", signature = "ASpliSplicingReport", 
           definition = function( x ){ x@localebased })

setGeneric( name = "localebased<-", def = function( x, value ) standardGeneric("localebased<-"))
setReplaceMethod( f = "localebased", signature = c("ASpliSplicingReport","data.frame"), 
                  definition = function( x,value ){ x@localebased <- value; return( x )})

setGeneric( name = "anchorbased", def = function( x ) standardGeneric("anchorbased"))
setMethod( f = "anchorbased", signature = "ASpliSplicingReport", 
           definition = function( x ){ x@anchorbased })

setGeneric( name = "anchorbased<-", def = function( x, value ) standardGeneric("anchorbased<-"))
setReplaceMethod( f = "anchorbased", signature = c("ASpliSplicingReport","data.frame"), 
                  definition = function( x,value ){ x@anchorbased <- value; return( x )})

# ---------------------------------------------------------------------------- #
