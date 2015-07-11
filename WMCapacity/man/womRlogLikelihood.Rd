     \name{RlogLikelihood}
     \Rdversion{1.1}
     \alias{womRlogLikelihood}
     \alias{womRlogPosterior}
     \alias{womRgradLogPosterior}
     \alias{womRlogPosteriorWithCov}
     \alias{womRgradLogPosteriorWithCov}
     \title{Compute various quantities relating to the working memory models in  the package WMCapacity}
     \description{
       These functions will compute the log-likelihood, log-posterior, and the gradient of the log-posterior for the multinomial tree models of working memory capacity in the package \code{WMCapacity}. These functions are not meant to be called by the user.} 
     
     \usage{
     womRlogLikelihood(x,setup)
     
     womRlogPosterior(x,setup)
    
     womRgradLogPosterior(x,setup)
     
     womRlogPosteriorWithCov(x,setup,precList,means)
         
     womRgradLogPosteriorWithCov(x,setup,precList,means)
     }
     

     \arguments{
       \item{x}{the vector containing values at which to evaluate the function.}
       \item{setup}{a list object, created by means of the \code{\link{womExtractModel}} function, containing the model specification.}
       \item{precList}{a list of length equal to the number of covariance matrices modelled. Each element of the list must be a square matrix, representing a precision matrix of appropriate size. The sizes are determined by the model specification.}
       \item{means}{A matrix with rows corresponding to the precision matrices, and columns correpsonding to the parameters involved in the precision matrices. The number of rows will be the number of precision matrices modelled, and the number of columns will be equal to the number of columns in the largest precision matrix.}
       }
     \details{
     The \code{womRlogLikelihood} function computes the log-likelihood of the model at the given values. This is used internally for computing the deviance information criterion (DIC) value.

     The \code{womRlogPosterior} function computes the (negative) log-posterior of the model at the given values. This is used internally for obtaining starting values for MCMC.

     The \code{womRgradLogPosterior} function computes the (negative) gradient of the log-posterior of the model at the given values. This is used internally for obtaining starting values for MCMC.

     The \code{womRlogPosteriorWithCov} and \code{RgradLogPosteriorWithCov} compute the (negative) log-posterior and (negative) gradient, respectively. The difference between these two functions and the two above (without covariances) is that these two functions are designed for models with covariance matrices defined. The log-posterior and gradient will be conditioned on the given covariance matrices and means.

     For further details, see the user's manual at \url{http://wmcapacity.r-forge.r-project.org/}.
   }

     \value{For the likelihood and posterior functions, a single value. For the gradient functions, a vector of the same length as the parameter vector. Each element repesents the partial derivative at the given value.}
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
	 womRlogLikelihood(posteriorMean,myModel)
     }
     
     }
\keyword{nonlinear}
\keyword{internal}
