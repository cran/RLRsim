\name{exactRLRT}
\alias{exactRLRT}
\title{Restricted Likelihood Ratio Tests for linear mixed models}
\description{
This function provides an (exact) restricted likelihood ratio test based on
simulated values from the finite sample distribution for testing whether the
variance of a random effect is 0 in a linear mixed model with known
correlation structure of the tested random effect and i.i.d. errors.
}
\usage{
exactRLRT(m, mA = NULL, m0 = NULL, seed = NA, nsim = 10000,
         log.grid.hi = 8, log.grid.lo = -10, gridlength = 200, print.p = TRUE,
         return.sample = FALSE, ...)
}
\arguments{
  \item{m}{The fitted model under the alternative or, for testing in models with
           multiple variance components, the reduced model containing only the
           random effect to be tested (see Details), an \code{lme}, \code{lmer} or  \code{spm} object}
  \item{mA}{The full model under the alternative for testing in models with multiple variance components}
  \item{m0}{The model under the null for testing in models with multiple variance components}
  \item{seed}{input for \code{set.seed}}
  \item{nsim}{Number of values to simulate}
  \item{log.grid.hi}{Lower value of the grid on the log scale. See \code{\link{exactRLRT}}. }
  \item{log.grid.lo}{Lower value of the grid on the log scale. See \code{\link{exactRLRT}}. }
  \item{gridlength}{Length of the grid. See \code{\link{exactLRT}}.  }
  \item{print.p}{print table with observed variance ratio, observed test statistic \eqn{-2\log RLR} and
                 p-value?}
  \item{return.sample}{return simulated sample?}
  \item{...}{Further parameters to be passed to \code{\link{RLRTSim}}.}
}
\details{
Tests in models with only a single variance component require only the first argument \code{m}.
For testing in models with multiple variance components, \code{m} contains only the random effect
whose variance is to be tested, while \code{mA} and \code{m0} are the models under the alternative
and the null, respectively. For models with a single variance component,
the simulated distribution is exact if the number
of parameters (fixed and random) is smaller than the number of observations.
For models with multiple components, theoretical results are still outstanding. Fairly extensive
simulation studies suggest that the application of the test is safe and the simulated distribution is
correct as long as the number
of parameters (fixed and random) is smaller than the number of observations and the
nuisance variance components are not superfluous or very small.
The exact finite sample distribution of the restricted likelihood ratio test statistic
that is simulated from was derived by Crainiceanu & Ruppert (2004).
}
\value{
If \code{return.sample = FALSE}, the p-value.
Else a list with components
  \item{p}{p-value for the observed test statistic}
  \item{sample}{the simulated sample of test statistics under the null}
}

\references{Crainiceanu, C. and Ruppert, D. (2004)
             Likelihood ratio tests in linear mixed models with one variance component,
             \emph{Journal of the Royal Statistical Society: Series B},\bold{66},165--185.}
\author{ Fabian Scheipl }

\seealso{\code{\link{RLRTSim}} for the underlying simulation algorithm;
    \code{\link{exactLRT}} for likelihood based tests}

}
\examples{
library(lme4)
data(sleepstudy)
mA <- lmer(Reaction ~ I(Days-4.5) + (1|Subject) + (0 + I(Days-4.5)|Subject), sleepstudy)
m0 <- update(mA, . ~ . - (0 + I(Days-4.5)|Subject))
m.slope  <- update(mA, . ~ . - (1|Subject))
#test for subject specific slopes:
exactRLRT(m.slope, mA, m0)

detach(package:lme4) #avoid conflicts
library(mgcv)
data(trees)
#test quadratic trend vs. smooth alternative
m.q<-gamm(I(log(Volume)) ~ Height + s(Girth, m = 3), data = trees, method = "REML")$lme
exactRLRT(m.q)
#test linear trend vs. smooth alternative
m.l<-gamm(I(log(Volume)) ~ Height + s(Girth, m = 2), data = trees, method = "REML")$lme
exactRLRT(m.l)
}
\keyword{htest}
