`exactRLRT` <-
function(m,mA=NULL,m0=NULL, seed= NA, nsim=10000,
                   log.grid.hi=8, log.grid.lo=-10, gridlength=200)
{
if(class(m)=="spm") {m<-m$fit; class(m)<-"lme"}
if (!(c.m<-(class(m))) %in% c("mer","lme"))  stop("Invalid m specified. \n")
if("REML"!=switch(c.m,lme=m$method,mer=switch(m@dims['REML']+1,"ML","REML")))
cat("Using restricted likelihood evaluated at ML estimators.\n Refit with method=\"REML\" for exact results.\n")


d<-switch(c.m,lme=extract.lmeDesign(m),mer=extract.lmerDesign(m))
X<-d$X; Z<-d$Z; y<-d$y
Vr<-d$Vr

K<-ncol(Z)     # no. of random effects
n<-nrow(X)     # no. of obs
p<-ncol(X)     # no. of fixed effects

if(is.null(mA) && is.null(m0)) #test for simple model
{
  #2*restricted ProfileLogLik under H0: lambda=0
  res<-qr.resid(qr(X),y)
  R<-qr.R(qr(X)); detXtX<-det(t(R)%*%R)

  reml.H0<- -((n-p)*log(2*pi) + (n-p)*log(sum(res^2)) +
              log(detXtX) + (n-p) - (n-p)*log(n-p))


  #observed value of the test-statistic
  reml.obs<- 2*logLik(m,REML=TRUE)[1]   #2*observed restricted ML
  rlrt.obs<-max(0,reml.obs - reml.H0)
  lambda<-d$lambda
}
else     #test for model with multiple var.comp.
{
  if(c.m == 'lme')
  {
    if(any(mA$fixDF$terms!=m0$fixDF$terms)) 		stop("Fixed effects structures of mA and m0 not identical.\n REML-based inference not appropriate.")
  } else {
  		if(any(mA@X!=m0@X)) stop("Fixed effects structures of mA and m0 not identical.\n REML-based inference not appropriate.") 
  }

  
		if(diff(anova(mA,m0)$Df)>1)
  {
    stop("Random effects not independent - covariance(s) set to 0 under the null hypothesis.\n Approximation not appropriate.\n")
  }

  rlrt.obs<-max(0,2*(logLik(mA,REML=TRUE)[1]-logLik(m0,REML=TRUE)[1]))
}

if(rlrt.obs!=0)
{
  sample<-RLRTSim(X,Z,sqrt.Sigma=chol(cov2cor(Vr)),lambda0=0, seed=seed, nsim=nsim,
                   log.grid.hi=log.grid.hi, log.grid.lo=log.grid.lo,
                   gridlength=gridlength)
  p<-mean(rlrt.obs<sample)
}
else p=1

RVAL <- list(statistic = c(RLRT=rlrt.obs), p.value = p, 
             method=paste('simulated finite sample distribution of RLRT.\n (p-value based on',nsim,'simulated values)'))
    class(RVAL) <- "htest"
    return(RVAL)
}

