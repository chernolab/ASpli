\name{DU accessors}
\alias{genesDE}
\alias{genesDE<-}
\alias{binsDU}
\alias{binsDU<-}
\alias{junctionsDU}
\alias{junctionsDU<-}

\title{Accessors for ASpliDU object}

\description{ Accessors for ASpliDU object }
\usage{ 
  genesDE( x )
  binsDU( x )
  junctionsDU( x )
}

\arguments{
  \item{ x }{ An ASpliDU object }
}

\value{ Returns dataframes with genomic metadata and logFC and pvalue }


\author{ Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz  }

\examples{
  
  # Get data tables from an ASpliDU object
  
  #du <- aspliDUexample1()
  
  #gde <- genesDE( du )
  #bdu <- binsDU( du )
  #jdu <- junctionsDU( du )
  
  # Set data tables to an ASpliDU object 
  
  #genesDE( du )     <- gde
  #binsDU( du )      <- bdu
  #junctionsDU( du ) <- jdu

}


