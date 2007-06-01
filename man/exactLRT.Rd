\name{exactLRT}
\alias{exactLRT}
\title{Likelihood Ratio Tests for simple linear mixed models}
\description{
This function provides an exact likelihood ratio test based on
simulated values from the finite sample distribution for simultaneous testing of the presence of the variance component
and some restrictions of the fixed effects in a simple linear mixed model with known
correlation structure of the random effect and i.i.d. errors.
}
\usage{
exactLRT(m, m0, seed = NA, nsim = 5000, log.grid.hi = 8,
    log.grid.lo = -10, gridlength = 200, print.p = TRUE, return.sample = FALSE)
}
\arguments{
  \item{m}{The fitted model under the alternative; of class \code{lme}, \code{lmer} or  \code{spm}}
  \item{m0}{The fitted model under the null hypothesis; of class \code{lm}}
  \item{seed}{ Specify a seed for \code{set.seed} }
  \item{nsim}{ Number of values to simulate }
  \item{log.grid.hi}{ Lower value of the grid on the log scale. See \code{\link{exactLRT}}. }
  \item{log.grid.lo}{ Lower value of the grid on the log scale. See \code{\link{exactLRT}}. }
  \item{gridlength}{ Length of the grid. See \code{\link{exactLRT}}.  }
  \item{print.p}{print table with observed variance ratio, observed test statistic \eqn{-2\log LR} and
                 p-value?}
  \item{return.sample}{return simulated sample?}
}
\details{
The model under the alternative must be a linear mixed model
\eqn{y=X\beta+Zb+\varepsilon}{y=X*beta+Z*b+epsilon} with
a \emph{single} random effect \eqn{b} with known correlation structure
and error terms that are i.i.d.
The hypothesis to be tested must be of the form
\deqn{H_0: \beta_{p+1-q}=\beta^0_{p+1-q},\dots,\beta_{p}=\beta^0_{p};\quad
Var(b)=0}{H0: beta_1=beta0_1,..,beta_q=beta0_q, Var(b)=0}
versus
\deqn{H_A:\; \beta_{p+1-q}\neq \beta^0_{p+1-q}\;\mbox{or}\dots
\mbox{or}\;\beta_{p}\neq \beta^0_{p}\;\;\mbox{or}\;Var(b)>0}{H0: beta_1 \neq beta0_1,..or..,beta_q \neq beta0_q ot Var(b)>0}
The exact finite sample  distribution of the likelihood ratio test statistic
that is simulated from was derived by Crainiceanu & Ruppert (2004).
}
\value{  If \code{return.sample = FALSE}, the p-value.
Else a list with components
  \item{p}{p-value for the observed test statistic}
  \item{sample}{the simulated sample of test statistics under the null}
}
\references{Crainiceanu, C. and Ruppert, D. (2004)
             Likelihood ratio tests in linear mixed models with one variance component,
             \emph{Journal of the Royal Statistical Society: Series B},\bold{66},165--185.}
\author{ Fabian Scheipl }

\seealso{\code{\link{LRTSim}} for the underlying simulation algorithm;
   \code{\link{RLRTSim}} and \code{\link{exactRLRT}} for restricted likelihood based tests}
\examples{
library(nlme);
data(Orthodont);

##test for Sex:Age interaction and Subject-Intercept
mA<-lme(distance ~ Sex * I(age - 11), random = ~ 1| Subject,
    data = Orthodont, method = "ML")
m0<-lm(distance ~ Sex + I(age - 11), data = Orthodont)
summary(mA)
summary(m0)
exactLRT(m = mA, m0 = m0)


library(SemiPar);
data(janka);
attach(janka);
mA <- spm(I(log(hardness)) ~ f(dens, basis = "trunc.poly", degree = 2),
    spar.method = "ML")
m0 <- lm(I(log(hardness))~ dens)

####test for linear trend vs. smooth alternative
\dontrun{exactLRT(m = mA, m0 = m0)}

}
\keyword{ htest }