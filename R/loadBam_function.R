loadBAM <- function( targets, cores = 1 ,libType, strandMode) {

  if(!is.null(cores)){ #Uses cores as flag for ASpli version. NULL means v2
    .Deprecated("", msg = "loadBAM is deprecated and is no longer needed. See ?gbCounts.")
    }
  datac <- lapply( as.character(targets$bam), function(x){

#option to load as SE o PE, default PE:
  if(libType=="SE"){ 
    r <- readGAlignments(x) 
  }
    else {
      r <- readGAlignmentPairs(x, strandMode=strandMode) 
    }
  
if(length(grep("[.]", seqlevels(r)) > 0)){
        seqlevels(r) <- gsub("[.]", "_", seqlevels(r))
        warning("Some seqnames had a '.' present in their names. ASpli had to normalize them using '_'.")
      }
      gc()
      return(r)
    })
  names( datac ) <- rownames(targets)
  
  return(datac)
  
}
