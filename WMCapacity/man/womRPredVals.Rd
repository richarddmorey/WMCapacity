     \name{womRPredVals}
     \Rdversion{1.1}
     \alias{womRPredVals}
     \title{Compute predicted values for working memory models in the package WMCapacity}
     \description{
       The \code{RPredVals} function computes predicted response probabilities for a
     specified working memory model and design, given values of the effect
     parameters. This function is useful for model checking.} 
     
     \usage{
     womRPredVals(x,setup)
     }
     

     \arguments{
       \item{x}{the vector containing values at which to evaluate the function.}
       \item{setup}{a list object, created by means of the \code{\link{womExtractModel}} function, containing the model specification.}
       }
     \details{
     The \code{womRPredVals} function computes predicted probabilities for the
       model and design specified by the \code{setup} argument. The
       probability of responding ``change'' in a change detection task
       is predicted for each trial.

     For further details, see the user's manual at \url{http://wmcapacity.r-forge.r-project.org/}.
   }

     \value{A vector of predicted probabilities; each element
     corresponds to a trial.}
     \seealso{
       \code{\link{wommbatGUI}}, for fitting the working memory models, and \code{\link{womExtractModel}} for extracting a model specification.
     }
     \examples{
     \dontrun{
     ## load Visual Array data set (Rouder et al., 2008)
     data(VisualArray)

     wommbatGUI(dataFrame = VisualArray)
     
	 ## Once the model is defined via the GUI, 
     ## compute log-likelihood at posterior mean
	 ## replace 'Model' with model name (in quotes)
	 myModel = womExtractModel('Model')
	 posteriorMean = myModel$results$pointEst[,5]
	 
	 womRPredVals(posteriorMean,myModel)
     }
     
     }
     \keyword{nonlinear}
