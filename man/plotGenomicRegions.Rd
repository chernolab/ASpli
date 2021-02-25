\name{plotGenomicRegions}

\alias{plotGenomicRegions}

\title{ Create genomic regions coverage plots }

\description{
  Graphic representation of coverage and junctions is useful to complement 
  the results of differential usage of bins and junction and differential 
  expression analysis. 
  
  Function plotGenomicRegions allow to create plots for multiple conditions for
  bins and genes. Each individual plot can only correspond to a single gene or 
  bin, but can contain many panels for different experimental conditions. 
}

\usage{
  plotGenomicRegions( features, x, genomeTxDb, targets, xIsBin = TRUE, 
    layout = 'auto', colors = 'auto', plotTitles = 'auto', sashimi = FALSE, 
    zoomOnBins = FALSE, deviceOpt = NULL, highLightBin = TRUE, outfolder = NULL, 
    outfileType = 'png', mainFontSize = 24, annotationHeight = 0.2,
    annotationCol = 'black', annotationFill = 'gray', annotationColTitle = 'black',
    preMergedBAMs = NULL, useTransparency = FALSE, tempFolder = 'tmp', 
    avoidReMergeBams= FALSE, verbose = TRUE )
}

\arguments{
  \item{features}{
    An ASpliFeatures object, generated with \code{binGenoms} function.
  }
  \item{x}{ 
    A character vector with the names of bins or genes to plot. To plot into a
    window is recommended that the length of x be one.
  }
  \item{genomeTxDb}{
    A TxDb object with the annotation of reference genome.
  }
  \item{targets}{A data frame containing sample, bam files and experimental 
    factor columns }
  \item{xIsBin}{
    A logical value that indicates if values in x corresponds to gene names or
    bin names.
  }
  \item{layout}{ 
    A character with value 'auto' or a character matrix with condition names 
    arranged with the desired layout of the panels in plots. The dimensions of 
    layout matrix, colors matrix and plotTitles matrices must be the same. 
    Matrix can have NA values, however, the height of the panels corresponding
    to that column are modified to occupy the total height.  
    The default value is 'auto'.
  }
  \item{colors}{ 
    A character containing value 'auto' or containing colors strings, or a 
    character matrix with color names arranged with the layout specified in 
    layout argument. The dimensions of layout matrix, colors matrix and 
    plotTitles matrices must be the same. The default value is 'auto'. 
  }
  \item{plotTitles}{ 
    A character containing value 'auto', or a character matrix with titles for
    each panel arranged with the layout specified in 
    layout argument. The dimensions of layout matrix, colors matrix and 
    plotTitles matrices must be the same. The default value is 'auto'. 
  }
  \item{sashimi}{ 
    A logical value that specifies that a sashimi plot for junctions must be 
    included into each panel. The default value is FALSE.
  }
  \item{zoomOnBins}{ 
    A FALSE logical value or a double value between 0 and 1. If value is FALSE
    then the genomic range to be plotted correspond to the complete gene, 
    otherwise the genomic range is that the size of the bin corresponds to a 
    fraction equals to zoomOnBins value of the total and is centered in the bin. 
    Is used only when xIsBin is TRUE. The default value is FALSE.}
  
  \item{deviceOpt}{ 
    A named list of arguments to be passed to the graphic device used to plot.
    This allow to further customization of the plot. The default value is an
    empty list.}
  
  \item{highLightBin}{
    A logical value that indicates if the bin should be highlighted. 
    The default value is TRUE.
  }
  \item{outfolder}{
    NULL or a character vector representing a folder path that will be used to 
    save the plot images. If the folder doesn't exists it is created. 
    If NULL, the plot is made in a window. The default value is NULL.
    For each bin, a single image file is generated. The name of the file is the
    name of the bin, added with a '.gr.' string, and the file extension at the 
    end. If the name of the bin contains invalid character for a file name, 
    those will be replaced by an underscore character.
  }
  \item{outfileType}{
    A character value the specifies the file format of the plot to be created.
    Is used only when outfolder is not NULL. Accepted values are 'png', 'jpeg',
    'bmp', 'tiff', 'pdf'. Each value is the graphic device used to create the 
    image. The default value is 'png'. 
  }
  \item{mainFontSize}{
    A numeric value specifying the size of the main title. The default value is 
    24.
  }
  \item{annotationHeight}{
    A double value specifying the proportion of the total height used to 
    represent the gene model. The default value is 0.2. 
  }
  \item{annotationCol}{
    A character value that specifies the color of the borders of bars in gene 
    model representation. The default value is 'black'.
  }
  \item{annotationFill}{
    A character value that specifies the color of the filling of bars in gene 
    model representation. The default value is 'gray'.
  }
  \item{annotationColTitle}{
    A character value that specifies the color of text in gene model 
    representation. The default value is 'black'.
  }
  \item{preMergedBAMs}{
    A one column data frame that associates a condition, specified in the row 
    name, with a character value representing the path of bam file with the 
    reads for that condition. This bam file is typically generated by merging
    the bam files of all replicates for that condition. The default value is 
    NULL, this specifies that not merged bam files are used, instead on-the-fly
    read extraction and merging is done from the bam files specified in the
    targets argument.   
  }
  \item{useTransparency}{
    A logical value that specifies if transparency will be used to generate the
    plots, this leads to better looking plot, however not all graphic devices
    support transparency. The default value is FALSE.
  }
  \item{tempFolder}{
    A character value specifying the path to store intermediate files produced
    while extracting and merging reads from bam files. It is only used when 
    preMergedBAMs arguments is NULL. The files created are not automatically 
    removed after plotting because can be reused to create a new plot of the
    same genes or bins with different graphic options. 
    The default value is 'tmp', that means a new 'tmp' folder in the current
    working folder.
  }
  \item{avoidReMergeBams}{
    A logical value specifying that extraction and merging of bam files will be
    avoided. This is only meaningful when the extraction and merging of the same
    set of genes and bins was done in the previous execution of 
    plotGenomicRegions function. The default value is FALSE.  
  }
  \item{verbose}{
    A logical value specifying that detailed information about the execution
    will be informed to the user. The default value is TRUE. 
  }
}

\value{ Returns a png for each selected bin }

\author{ Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz  }


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
  #                   "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  #targets <- data.frame( 
  #             row.names = paste0('Sample_',c(1:6)),
  #             bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
  #             factor1 = c( 'C','C','C','D','D','D'),
  #             stringsAsFactors = FALSE )
  
  # Plot a single bin to a window
  #plotGenomicRegions( 
  #  features, 
  #  'GENE01:E002', 
  #  genomeTxDb, 
  #  targets, 
  #  sashimi = FALSE,
  #  colors = '#AA4444', 
  #  annotationHeight = 0.1, 
  #  tempFolder = 'tmp', 
  #  verbose = TRUE , 
  #  avoidReMergeBams = FALSE, 
  #  useTransparency = FALSE ) 
  #  
  # plot two bins to pdf files.
  #plotGenomicRegions( 
  #  features, c( 'GENE01:E002', 'GENE02:E002' ), 
  #  genomeTxDb, 
  #  targets, 
  #  layout = matrix( c( 'C', 'D'), ncol = 1),
  #  colors = matrix( c( '#663243', '#363273'), ncol = 1),
  #  plotTitles = matrix( c( 'C condition', 'D condition'), ncol = 1),
  #  sashimi = FALSE,
  #  mainFontSize = 12,
  #  annotationHeight = 0.1, 
  #  tempFolder = 'tmp', 
  #  verbose = TRUE , 
  #  avoidReMergeBams = FALSE, 
  #  useTransparency = TRUE,
  #  outfolder = '.',
  #  outfileType = 'pdf',
  #  deviceOpt = list( height = 6, width = 5, paper = 'a4r' ) ) 
}

\seealso{
  \code{\link{Devices}},
  \code{\link{pdf}}, 
  \code{\link{png}},
  \code{\link{bmp}},
  \code{\link{jpeg}},
  \code{\link{tiff}}
}
