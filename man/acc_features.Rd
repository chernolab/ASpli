\name{Features accesors}
\alias{featuresb}
\alias{featuresg}
\alias{featuresj}
\alias{transcriptExons}
\alias{featuresb<-}
\alias{featuresg<-}
\alias{featuresj<-}
\alias{transcriptExons<-}

\title{ Accessors for ASpliFeatures object }

\description{ Accessors for ASpliFeatures object }

\usage{
  featuresg( x )
  featuresb( x )
  featuresj( x )
  transcriptExons( x )
}

\arguments{
  \item{x}{ An ASpliFeatures object }
}

\value{
  Returns a GenomicRanges object. Function featuresg returns a GRangesList 
  object containing exon ranges for each gene. Functions featuresb and 
  featuresj, returns GRanges object for all bins and junctions. 
}

\author{ Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz}

\examples{
  # Get data from an ASpliFeatures object
  
  features <- aspliFeaturesExample()
  
  fg <- featuresg( features )
  fb <- featuresb( features )
  fj <- featuresj( features )

  
  # Set data to an ASpliFeatures object
  
  featuresg( features ) <- fg 
  featuresb( features ) <- fb 
  featuresj( features ) <- fj 

}


