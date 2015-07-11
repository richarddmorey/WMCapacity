     \name{womExtractModel}
     \Rdversion{1.1}
     \alias{womExtractModel}
     \title{Extract WMCapacity model specification and results.}
     \description{
       The \code{womExtractModel} function allows the user to extract model specification and results that were defined in the \code{wommbatGUI}
	   graphical user interface. This allows the user to extract 
	   results for the purposes of making plots, et cetera, in R.} 
     
     \usage{
     womExtractModel(name=1)
     }
     
     \arguments{
       \item{name}{the name (or number) of the model defined in the GUI.}
       }
     \details{
     The \code{womExtractModel} function allows the user access to the model specifications and analysis results from the GUI, once
	 they've been defined. 

     For further details, see the user's manual at \url{http://wmcapacity.r-forge.r-project.org/}.
   }

     \value{
	 A list containing the following elements (if applicable):
	 \item{modelName}{The name of the model.}
	 \item{model}{A list containing the model specification.}
	 \item{priors}{A list containing the prior specification.}
	 \item{settings}{A list containing the MCMC settings, if an analysis has been performed.}
	 \item{results}{A list containing the results, if an analysis has been performed.}
	 
	 }
     \seealso{
       \code{\link{wommbatGUI}}, for fitting the working memory models.
     }
     \examples{
     \dontrun{
     ## load Visual Array data set (Rouder et al., 2008)
     data(VisualArray)
     
	 ## Define the model in the GUI
	 wommbatGUI(dataFrame=VisualArray)
	 
     # extract the first model. Replace 'Model' with the model name (in quotes)
     myModel = womExtractModel(name='Model')
     
	 # examine the posterior means (if an analysis has been performed)
	 myModel$results$pointEst
	 }
     
     }
     \keyword{nonlinear}
