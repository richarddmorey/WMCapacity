     \name{wommbatGUI}
     \Rdversion{1.1}
     \alias{wommbatGUI}
     \title{Start GUI for working memory capacity estimation.}
     \description{
       \code{wommbatGUI} starts the graphical user interface for building and
     estimating parameters of working memory models. This is the main
     interface for the \code{WMCapacity} package.
     }
     \usage{
     wommbatGUI(project = NULL, projectFile= NULL, CSVfile = NULL, dataFrame = NULL, devel=FALSE)
     }
     \arguments{
       \item{dataFrame}{a data frame containing the trial-by-trial data to be
     analyzed.}
       \item{CSVfile}{the location of a CSV file containing the
     trial-by-trial data to be analyzed.}
       \item{projectFile}{a saved WMCapacity analysis (.Rdata file), containing
	 all the necessary information necessary to start an analysis.}
       \item{project}{a environment of a WMCapacity analysis, containing
	 all the necessary information necessary to start an analysis.}
	   \item{devel}{Turn on (undocumented) testing features.}
	 
     }
     \details{This function starts the GUI for estimating working memory
     capacity from change detection data, using the hierarchical
     Bayesian models described in Morey (2011).

     There are a number of ways of passing data to the function. If you
     have the trial-by-trial data already loaded in R, you can pass the
     data via the \code{dataFrame} argument. If it is in a CSV file, you can load
     it by passing the path via the \code{CSVfile} argument. If you saved the analysis
     and wish to reload it, you
     can use the pass the saved file name to the \code{projectFile} argument.
     argument. Alternatively, you may load data via the GUI.

     For further details, click on the Help button in the GUI, or see
     the user's manual at \url{http://wmcapacity.r-forge.r-project.org/}.
   }

     \value{This function returns nothing. The defined models and results may be extracted my means of the \code{\link{womExtractModel}} function.}
     \seealso{
       \code{\link{wommbatNoGUI}}, for the non-GUI interface
     (useful for simulations).
     }
     \examples{
     \dontrun{
     ## load Visual Array data set (Rouder et al., 2008)
     data(VisualArray)

     wommbatGUI(dataFrame=VisualArray)
     }
     }
     \keyword{nonlinear}
