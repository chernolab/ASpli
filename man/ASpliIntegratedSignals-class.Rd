\name{ASpliIntegratedSignals-class}

\Rdversion{1.1}

\docType{class}

\alias{ASpliIntegratedSignals}
\alias{ASpliIntegratedSignals-class}
\alias{signals,ASpliIntegratedSignals-method}
\alias{filters,ASpliIntegratedSignals-method}
\alias{signals<-,ASpliIntegratedSignals,data.frame-method}
\alias{filters<-,ASpliIntegratedSignals,data.frame-method}
\alias{splicingReport,ASpliIntegratedSignals-method}
\alias{show,ASpliIntegratedSignals-method}
\alias{writeIntegratedSignals,ASpliIntegratedSignals-method}
\alias{exportIntegratedSignals,ASpliIntegratedSignals-method}


\title{Class \code{"ASpliIntegratedSignals"}}

\description{ Contains results of differential expression at junction level.
}

\section{Slots}{
  \describe{
    \item{\code{signals}}{}
    \item{\code{filters}}{}
    \item{\code{.ASpliVersion}:}{
      ASpli version when this object was created. It should not be modified by the user.
    }     
  }
}

\author{ Andres Rabinovich, Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz  }