\name{extract.lmeDesign}
\alias{extract.lmerDesign}
\alias{extract.lmeDesign}
\title{Extract the Design of a linear mixed model}
\description{
 These functions extract various elements of the design of a fitted \code{lme}- or \code{lmer}-Object.
 They are called by \code{exactRLRT} and \code{exactLRT}.
}
\usage{
extract.lmeDesign(m)
extract.lmerDesign(m)
}
\arguments{
  \item{m}{fitted \code{lme}- or \code{lmer}-Object }
}
\value{a list with components
  \item{Vlambda}{estimated covariance of the response divided by the estimated variance of the residuals}
  \item{V}{estimated covariance of the error terms}
  \item{Vr}{estimated covariance of the random effects divided by the estimated variance of the residuals}
  \item{X}{design of the fixed effects}
  \item{Z}{design of the random effects}
  \item{sigmasq}{variance of the residuals}
  \item{lambda}{variance ratios of the random effects and the residuals}
  \item{y}{response variable}
}
\author{Fabian Scheipl}
\note{\code{extract.lmeDesign} is an adaptation of \code{mgcv::extract.lme.cov} by Simon Wood}
\examples{
library(nlme)
d <- extract.lmeDesign(lme(distance ~ age + Sex, data = Orthodont, random = ~ 1))
}
\keyword{utilities}