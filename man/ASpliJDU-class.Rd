\name{ASpliJDU-class}

\Rdversion{1.1}

\docType{class}

\alias{ASpliJDU}
\alias{ASpliJDU-class}
\alias{anchorc,ASpliJDU-method}
\alias{anchorj,ASpliJDU-method}
\alias{localec,ASpliJDU-method}
\alias{localej,ASpliJDU-method}
\alias{jir,ASpliJDU-method}
\alias{jes,ASpliJDU-method}
\alias{jalt,ASpliJDU-method}
\alias{anchorc<-,ASpliJDU,data.frame-method}
\alias{anchorj<-,ASpliJDU,data.frame-method}
\alias{localec<-,ASpliJDU,data.frame-method}
\alias{localej<-,ASpliJDU,data.frame-method}
\alias{jir<-,ASpliJDU,data.frame-method}
\alias{jes<-,ASpliJDU,data.frame-method}
\alias{jalt<-,ASpliJDU,data.frame-method}
\alias{jDUreport,ASpliJDU-method}
\alias{show,ASpliJDU-method}
\alias{writeJDU,ASpliJDU-method}

\title{Class \code{"ASpliJDU"}}

\description{ Contains results of differential expression at junction level.
}

\section{Slots}{
  \describe{
    \item{\code{localec}}{}
    \item{\code{localej}}{}
    \item{\code{anchorc}}{}
    \item{\code{anchorj}}{}
    \item{\code{jir}}{}
    \item{\code{jes}}{}
    \item{\code{jalt}}{}
    \item{\code{contrast}}{}   
    \item{\code{.ASpliVersion}:}{
      ASpli version when this object was created. It should not be modified by the user.
    }     
  }
}

\author{ Andres Rabinovich, Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz  }