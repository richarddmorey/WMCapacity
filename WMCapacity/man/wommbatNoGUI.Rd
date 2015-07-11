     \name{wommbatNoGUI}
     \Rdversion{1.1}
     \alias{wommbatNoGUI}
     \title{Analyze working memory models with no GUI.}
     \description{
       The \code{wommbatNoGUI} function analyzes the models defined in a saved wommbat analysis file, without using the GUI. This is useful for simulations.
     }
     \usage{
     wommbatNoGUI(project=NULL, projectFile= NULL, settings)
     }
     \arguments{
       \item{settings}{a list containing the MCMC settings for the analysis. See \code{\link{womExtractModel}}.}
       \item{projectFile}{a saved WMCapacity analysis (.Rdata file), containing
	 all the necessary information necessary to start an analysis.}
       \item{project}{a environment of a WMCapacity analysis, containing
	 all the necessary information necessary to start an analysis.}
     }
     \details{This function analyzes the models specified in the saved analysis file, using no GUI.

     For further details, click on the Help button in the GUI, or see
     the user's manual at \url{http://wmcapacity.r-forge.r-project.org/}.
   }

     \value{This function returns nothing. The defined models and results may be extracted my means of the \code{\link{womExtractModel}} function.}
     \seealso{
       \code{\link{wommbatGUI}}, for the GUI interface
     (useful for simulations).
     }
     \examples{
     }
     \keyword{nonlinear}
