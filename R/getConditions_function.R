# This function return the names of the unique conditions corresponing to the 
# targets in the order that they are assigned to define the constrat.
getConditions <- function( targets ) {
  return( unique ( .condenseTargetsConditions( targets )$condition ))
}