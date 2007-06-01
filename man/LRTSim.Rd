\name{LRTSim}
\alias{LRTSim}
\alias{RLRTSim}
\title{Simulation of the (Restricted) Likelihood Ratio Statistic}
\description{
These functions simulate values from the (exact) finite sample distribution of the
(restricted) likelihood
ratio statistic for testing the presence of the variance component
(and restrictions of the fixed effects) in a simple linear mixed model with known
correlation structure of the random effect and i.i.d. errors. They are usually called by
\code{exactLRT} or \code{exactRLRT}.}
\usage{
LRTSim(X, Z, q, sqrt.Sigma, seed = NA, nsim = 5000, log.grid.hi = 8,
       log.grid.lo=-10, gridlength=200)
RLRTSim(X, Z, sqrt.Sigma, lambda0 = NA, seed = NA, nsim = 5000, use.approx=0,
       log.grid.hi=8, log.grid.lo=-10, gridlength=200)
}
\arguments{
  \item{X}{ The fixed effects design matrix of the model under the alternative }
  \item{Z}{ The random effects design matrix of the model under the alternative }
  \item{q}{ The number of parameters restrictions on the fixed effects (see Details)}
  \item{sqrt.Sigma}{  The upper triangular cholesky factor of  the correlation matrix of the random effect}
  \item{lambda0}{The value of the ratio of the variances of the random effect and the errors under the null}
  \item{seed}{ Specify a seed for \code{set.seed}}
  \item{nsim}{ Number of values to simulate }
  \item{use.approx}{If 0, the exact distribution is simulated. If between 0 and 1, only the largest
  eigenvalues whose sum represents at least \code{use.approx}*(sum of all eigenvalues) are used.}
  \item{log.grid.hi}{ Lower value of the grid on the log scale. See \bold{Details}}
  \item{log.grid.lo}{ Lower value of the grid on the log scale. See \bold{Details}}
  \item{gridlength}{ Length of the grid.  }
}
\details{
The model under the alternative must be a linear mixed model
\eqn{y=X\beta+Zb+\varepsilon}{y=X*beta+Z*b+epsilon} with
a single random effect \eqn{b} with known correlation
structure \eqn{Sigma} and error terms that are i.i.d.
The simulated distribution of the likelihood ratio statistic
was derived by Crainiceanu & Ruppert (2004).
The simulation algorithm uses a gridsearch over a log-regular grid
of values of \eqn{\lambda=\frac{Var(b)}{Var(\varepsilon)}}{lambda=Var(b)/Var(epsilon)} to maximize the
likelihood under the alternative for \code{nsim} realizations of \eqn{y}
drawn under the null hypothesis.
\code{log.grid.hi} and \code{log.grid.lo} are  the lower and upper
limits of this grid on the log scale.
\code{gridlength} is the number of points on the grid.
}
\value{A dataframe with \code{nsim} rows and 2 columns:
  \item{lambda}{estimates of \eqn{\lambda} under the alternative}
  \item{lrt}{simulated realizations of the test statistic}
}
\references{ Crainiceanu, C. and Ruppert, D. (2004)
             Likelihood ratio tests in linear mixed models with one variance component,
             \emph{Journal of the Royal Statistical Society: Series B},\bold{66},165--185.
}
\author{Fabian Scheipl}
\seealso{ \code{\link{exactLRT}}, \code{\link{exactRLRT}} for tests}
\examples{
library(lme4)
g <- rep(1:10, e = 10)
x <- rnorm(100)
y <- 0.1 * x + rnorm(100)
m <- lmer(y ~ x + (1|g), method="ML")
m0 <- lm(y ~ 1)
X <- m@X
Z <- t(as.matrix(m@Zt))
lrt <- LRTSim(X, Z, 1, diag(10))
}
\keyword{datagen}
\keyword{distribution}