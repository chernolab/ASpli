\name{plotBins}

\alias{plotBins}

\title{
  Draw plots of gene counts, bin counts, PSI/PIR value, inclusion and exclusion
  junctions for selected bins. 
}

\description{ Creates a plot with gene counts, bin counts, PSI/PIR value, 
  inclusion and exclusion junctions for selected bins and conditions. }

\usage{
  plotBins( counts, 
            as,
            bin, 
            factorsAndValues, 
            targets,
            main            = NULL,
            colors          = c( '#2F7955', '#79552F', '#465579', 
                                 '#A04935', '#752020', '#A07C35') ,
            panelTitleColors = '#000000',
            panelTitleCex   = 1,
            innerMargins    = c( 2.1, 3.1, 1.1, 1.1 ),
            outerMargins    = c( 0, 0, 2.4, 0 ), 
            useBarplots     = NULL,
            barWidth        = 0.9,
            barSpacer       = 0.4,
            las.x           = 2,
            useHCColors     = FALSE,
            legendAtSide    = TRUE,
            outfolder       = NULL,
            outfileType     = c( 'png', 'bmp', 'jpeg', 'tiff', 'pdf')[1],
            deviceOpt       = NULL ) 
}

\arguments{
  
  \item{counts}{ An object of class \code{ASpliCounts} }
  
  \item{as}{ An object of class \code{ASpliAS} }
  
  \item{bin}{ A character vector with the names of the bins to be plotted. }
  
  \item{factorsAndValues}{ A list containing the factor and the values for each
    factor to be plotted. The order of the factors will modify how the 
    conditions are grouped in the plot. \code{factorsAndValues} must be a named
    list, where the name of each element is a factor and the list element itself
    is a character vector of the values of this factor in the order to be 
    plotted. See examples for more details.
  }
    
  \item{targets}{ A data frame containing sample, bam files and experimental 
    factor columns 
  }
  
  \item{main}{ Main title of the plot. If \code{NULL} the bin name is used as
    title.
  }
  
  \item{colors}{ A vector of character colors for lines and bar plots. }
  
  \item{panelTitleColors}{ A vector of character colors for the titles of each
    plot panel.
  }
  
  \item{panelTitleCex}{ Character size expansion for panel titles.}
  
  \item{innerMargins}{ A numerical vector of the form c(bottom, left, top, 
    right) which gives the size of each plot panel margins. Defaults to c( 2.1, 
    3.1, 1.1, 1.1 ) }
  
  \item{outerMargins}{ A numerical vector of the form c(bottom, left, top, 
    right) which gives the size of margins. Defaults to c( 0, 0, 2.4, 0 ) }

  \item{useBarplots}{ A logical value that indicates the type of plot to be 
    used. If \code{TRUE} bar plots are used, if \code{FALSE} lines are used. If
    \code{NULL} the type is bar plot if there just two conditions and lines if 
    there are more than two conditions.}
  
  \item{barWidth}{ The width of the bars in bar plots. \code{barWidth} must be 
    in (0,1] range. Default value is 0.9.}
  
  \item{barSpacer}{ Fraction of \code{barwidth} used as spacer between bar plot 
    groups. Defaule value is 0.4.}
  
  \item{las.x}{ Text orientation of x-axis labels. }
  
  \item{useHCColors}{ A logical value. If \code{TRUE} \code{panelTitleColors}} 
    are not used, instead panel title are automatically chosen to have high 
    contrast against \code{colors}.

  \item{legendAtSide}{ A logical value that forces panel title to be shown on 
    the y-axis, instead of over the plot. }

  \item{outfolder}{ Path to output folder to write plot images. Is \code{NULL},
    plot are rendered on the default device }
    
  \item{outfileType}{ File format of the output files used if \code{outfolder} 
    is not \code{NULL}. Accepted values are 'png', 'jpeg', 'tiff', 'pdf'. Each
    value selects the graphic device of the same name. The name of the image
    file is the name of bin with the corresponding extension given by the chosen
    type }
  
  \item{deviceOpt}{A list of named options to be passed to the graphic device 
    selected in \code{outfileType}}
    
}

\value{
Returns a png for each selected bin }

\author{ Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz  }

\seealso{
  \code{\link{plotGenomicRegions}}, 
}

\examples{

  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  #library(GenomicFeatures)
  #genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
  #                               package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  #features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  #bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
  #                   "A_D_0.bam", "A_D_1.bam", "A_D_2.bam",
  #                   "B_C_0.bam", "B_C_1.bam", "B_C_2.bam", 
  #                   "B_D_0.bam", "B_D_1.bam", "B_D_2.bam" )
                     
  #targets <- data.frame( 
  #             row.names = paste0('Sample_',c(1:12)),
  #             bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
  #             factor1 = c( 'A','A','A','A','A','A','B','B','B','B','B','B'),
  #             factor2 = c( 'C','C','C','D','D','D','C','C','C','D','D','D') )
  
  # Load reads from bam files
  #bams <- loadBAM( targets )
  
  # Read counts from bam files
  #counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
  #                        maxISize = 50000 )
  
  # Calculate differential usage of genes, bins and junctions 
  #du       <- DUreport.norm( counts, targets , contrast = c(1,-1,-1,1))

  # Calculate PSI / PIR for bins and junction.
  #as       <- AsDiscover( counts, targets, features, bams, readLength = 100, 
  #                        threshold = 5, cores = 1 )
  
  # Plot bin data. Factor2 is the main factor for graphic representation in
  # this example as it is the first in factorsAndValues argument.
  # This makes a bar plot comparing four conditions, grouped by factor1.
  #plotBins( counts, as, 'GENE03:E002', 
  #  factorsAndValues = list( 
  #    factor2 = c('C','D'), 
  #    factor1 = c('A','B') ),
  #  las.x = 1,
  #  legendAtSide = TRUE,
  #  useHCColors = TRUE,   
  #  targets = targets,
  #  barWidth = 0.95,
  #  innerMargins = c( 2.1, 4.1, 1.1, 1.1 ) )
    
    
  # Redefine targets  
  #targets <- data.frame( 
  #             row.names = paste0('Sample_',c(1:12)),
  #             bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
  #             factor1 = c( 'A','A','B','B','C','C','D','D','E','E','F','F') )
  
  #as       <- AsDiscover( counts, targets, features, bams, readLength = 100, 
  #                        threshold = 5, cores = 1 )
  
  # This makes a line plot for six conditions, grouped by factor1.                       
  #plotBins( counts, as, 'GENE03:E002', 
  #  factorsAndValues = list( 
  #    factor1 = c('A','B','C','D','E','F') ),
  #  las.x = 1,
  #  legendAtSide = FALSE,
  #  targets = targets,
  #  innerMargins = c( 2.1, 4.1, 1.1, 1.1 ) )                        
  
}
