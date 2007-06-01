`exactRLRT` <-
function(m,mA=NULL,m0=NULL, seed= NA, nsim=10000,
                   log.grid.hi=8, log.grid.lo=-10, gridlength=200,
                   print.p=TRUE,return.sample=FALSE,...)
{
if(class(m)=="spm") {m<-m$fit; class(m)<-"lme"}
if (!(c.m<-(class(m))) %in% c("lmer","lme"))  stop("Invalid m specified. \n")
if("REML"!=switch(c.m,lme=m$method,lmer=switch(m@status[2]+1,"ML","REML")))
cat("Using restricted likelihood evaluated at ML estimators.\n Refit with method=\"REML\" for exact results.\n")


d<-switch(c.m,lme=extract.lmeDesign(m),lmer=extract.lmerDesign(m))
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
  dA<-switch(class(mA),lme=extract.lmeDesign(mA),lmer=extract.lmerDesign(mA))
  d0<-switch(class(m0),lme=extract.lmeDesign(m0),lmer=extract.lmerDesign(m0))
  check4cor<-switch(class(mA),lme=any(as.matrix(mA$modelStruct$reStruct[[1]])!=diag(diag(as.matrix(mA$modelStruct$reStruct[[1]]))))
                           ,lmer=any(dim(as.matrix(tail(mA@Omega, 1)[[1]]))!=c(1,1)))   #is cov(ranef) a diagonal matrix?
  if(check4cor)
  {
    cat("Random effects not independent - covariance(s) set to 0 under the null hypothesis.\n Approximation not appropriate.\n")
    return(invisible())
  }
  if(any(dim(dA$X)!=dim(d0$X)))
  {
    stop("Fixed effects structures of mA and m0 not identical.\n REML-based inference not appropriate.")
  }

  rlrt.obs<-max(0,2*(logLik(mA,REML=TRUE)[1]-logLik(m0,REML=TRUE)[1]))
  lambda<-tail(dA$lambda,1)
}

if((rlrt.obs!=0)||return.sample)
{
  sample<-RLRTSim(X,Z,sqrt.Sigma=chol(cov2cor(Vr)),lambda0=0, seed=seed, nsim=nsim,
                   log.grid.hi=log.grid.hi, log.grid.lo=log.grid.lo,
                   gridlength=gridlength,...)[,2]
  p<-mean(rlrt.obs<sample)
}
else p=1

if(print.p)
{
  ans<-matrix(c(lambda,rlrt.obs,p),1,3)
  colnames(ans)<-c("variance ratio", "RLRT", "p-value")
  rownames(ans)<-c(" ")
  print(ans,digits=4)
  cat(paste("      p-value based on",nsim, "simulated values. \n"))
}
if(return.sample) return(list(p=p,sample=sample))
invisible(p)
}

