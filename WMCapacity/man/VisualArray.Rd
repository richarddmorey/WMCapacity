\name{VisualArray}
\alias{VisualArray}
\docType{data}
\title{Working memory visual array change detection responses}
\description{
This dataset contains the (cleaned) data which was analyzed in Rouder,
Morey, Cowan, Zwilling, Morey, and Pratte (2008).   
}
\usage{data(VisualArray)}
\format{
  A data frame with 12345 observations on the following 36 variables.
  \describe{
    \item{\code{sub}}{the subject number}
    \item{\code{blk}}{the block number}
    \item{\code{trl}}{the trial number within a block}
    \item{\code{prch}}{the log-odds of a change trial within that block}
    \item{\code{N}}{the visual array set size}
    \item{\code{ischange}}{whether the array changed or not (1=changed)}
    \item{\code{resp}}{whether the participant responded ``change'' (1=changed)}
    \item{\code{RT}}{the response time}
    \item{\code{probesq}}{which square number was probed}
    \item{\code{newcol}}{the color of the probed square}
    \item{\code{c1}}{the color of square 1}
    \item{\code{x1}}{the x location of square 1}
    \item{\code{y1}}{the y location of square 1}
    \item{\code{c2}}{the color of square 2}
    \item{\code{x2}}{the x location of square 2}
    \item{\code{y2}}{the y location of square 2}
    \item{\code{c3}}{the color of square 3}
    \item{\code{x3}}{the x location of square 3}
    \item{\code{y3}}{the y location of square 3}
    \item{\code{c4}}{the color of square 4}
    \item{\code{x4}}{the x location of square 4}
    \item{\code{y4}}{the y location of square 4}
    \item{\code{c5}}{the color of square 5}
    \item{\code{x5}}{the x location of square 5}
    \item{\code{y5}}{the y location of square 5}
    \item{\code{c6}}{the color of square 6}
    \item{\code{x6}}{the x location of square 6}
    \item{\code{y6}}{the y location of square 6}
    \item{\code{c7}}{the color of square 7}
    \item{\code{x7}}{the x location of square 7}
    \item{\code{y7}}{the y location of square 7}
    \item{\code{c8}}{the color of square 8}
    \item{\code{x8}}{the x location of square 8}
    \item{\code{y8}}{the y location of square 8}
    \item{\code{cor}}{whether the participant was correct or not}
    \item{\code{oldcol}}{the color of the probed square in the study array}
  }
}
\details{

  Each trial consisted of a presentation of N colored squares in an
  array, followed by a mask, and then a single square. The square was
  either the same color or a different color than the one in the same
  location in the first array. The participant responded ``change'' if
  they believed the square was different.
}
\source{
Rouder, J. N., Morey, R. D., Cowan, N., Zwilling, C. E., Morey, C. C. &
Pratte, M. S. (2008). An Assessment of Fixed-Capacity Models of Visual
Working Memory. Proceedings of the National Academy of Sciences, 105, 5976-5979.
}
\examples{
data(VisualArray)
## Aggregated hit and false alarm rates
tapply(VisualArray$resp,list(VisualArray$ischange,
       VisualArray$N,VisualArray$prch),mean)

}
\keyword{datasets}
