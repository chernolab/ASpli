\name{ASpliSplicingReport-class}

\Rdversion{1.1}

\docType{class}

\alias{ASpliSplicingReport}
\alias{ASpliSplicingReport-class}
\alias{anchorbased,ASpliSplicingReport-method}
\alias{binbased,ASpliSplicingReport-method}
\alias{localebased,ASpliSplicingReport-method}
\alias{anchorbased<-,ASpliSplicingReport,data.frame-method}
\alias{localebased<-,ASpliSplicingReport,data.frame-method}
\alias{binbased<-,ASpliSplicingReport,data.frame-method}
\alias{splicingReport,ASpliSplicingReport-method}
\alias{show,ASpliSplicingReport-method}
\alias{writeSplicingReport,ASpliSplicingReport-method}
\alias{integrateSignals,ASpliSplicingReport-method}
\alias{exportSplicingReports,ASpliSplicingReport-method}
\alias{filterSignals,ASpliSplicingReport-method}


\title{Class \code{"ASpliSplicingReport"}}

\description{ Contains results of differential expression at junction level.
}

\section{Slots}{
  \describe{
    \item{\code{binbased}}{}
    \item{\code{localebased}}{}
    \item{\code{anchorbased}}{}
    \item{\code{contrast}}{}
    \item{\code{.ASpliVersion}:}{
      ASpli version when this object was created. It should not be modified by the user.
    }     
  }
}

\author{ Andres Rabinovich, Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz  }