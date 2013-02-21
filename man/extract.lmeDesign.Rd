\name{extract.lmeDesign}
\alias{extract.lmerModDesign}
\alias{extract.lmerDesign}
\alias{extract.lmeDesign}
\title{Extract the Design of a linear mixed model}
\description{
 These functions extract various elements of the design of a fitted \code{lme}-, \code{mer} or \code{lmerMod}-Object.
 They are called by \code{exactRLRT} and \code{exactLRT}.
}
\usage{
extract.lmeDesign(m)
extract.lmerDesign(m)
extract.lmerModDesign(m)
}
\arguments{
  \item{m}{fitted \code{lme}- or \code{mer}-Object }
}
\value{a list with components
  \item{Vr}{estimated covariance of the random effects divided by the estimated variance of the residuals}
  \item{X}{design of the fixed effects}
  \item{Z}{design of the random effects}
  \item{sigmasq}{variance of the residuals}
  \item{lambda}{ratios of the variances of the random effects and the variance of the residuals}
  \item{y}{response variable}
}
\author{Fabian Scheipl, \code{extract.lmerModDesign} supplied by Ben Bolker. Many thanks to Andrzej Galecki and Tomasz Burzykowski for bug fixes.}
\note{\code{extract.lmeDesign} is based on \code{mgcv::extract.lme.cov} by Simon Wood}
\examples{
library(nlme)
d <- extract.lmeDesign(lme(distance ~ age + Sex, data = Orthodont, random = ~ 1))
}
\keyword{utilities}