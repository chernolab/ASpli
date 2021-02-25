\name{ASpliDU-class}

\Rdversion{1.1}

\docType{class}

\alias{ASpliDU}
\alias{ASpliDU-class}
\alias{binsDU,ASpliDU-method}
\alias{genesDE,ASpliDU-method}
\alias{junctionsDU,ASpliDU-method}
\alias{binsDU<-,ASpliDU-method}
\alias{genesDE<-,ASpliDU,data.frame-method}
\alias{junctionsDU<-,ASpliDU,data.frame-method}
\alias{writeDU,ASpliDU-method}
\alias{filterDU,ASpliDU-method}
\alias{mergeBinDUAS,ASpliDU,ASpliAS-method}
\alias{containsGenesAndBins,ASpliDU-method}
\alias{containsJunctions,ASpliDU-method}
\alias{splicingReport,ASpliDU-method}

\title{Class \code{"ASpliDU"}}

\description{ Contains results of differential expression at gene level and 
  differential usage at bin and junction level estimation using DEreport method.
}

\section{Slots}{
  \describe{
    \item{\code{genes}}{}
    \item{\code{bins}}{}
    \item{\code{junctions}}{}
    \item{\code{contrast}}{}
    \item{\code{.ASpliVersion}:}{
      ASpli version when this object was created. It should not be modified by the user.
    } 
  }
}


\author{ Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz  }
