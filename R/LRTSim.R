`LRTSim` <-
function(X,Z,q,sqrt.Sigma, seed=NA, nsim=5000,
                       log.grid.hi=8, log.grid.lo=-10, gridlength=200)
{
K<-NCOL(Z)     # no. of random effects
n<-NROW(X)     # no. of obs
p<-NCOL(X)     # no of fixed effects

#compute eigenvalues
mu<- (svd(sqrt.Sigma%*%t(qr.resid(qr(X), Z)), nu = 0, nv = 0)$d)^2
ksi<-(svd(sqrt.Sigma%*%t(Z),nu=0,nv=0)$d)^2

#norm eigenvalues
mu<-mu/max(mu)
ksi<-ksi/max(ksi)


#generate random samples
if(!is.na(seed)) set.seed(seed)
w.k.sq.mat<-matrix(rchisq(nsim*K,1),nrow=K)  #K*nsim matrix of ChiSq(1)
w.sum2<-rchisq(nsim,n-p-K)                   #1*nsim vector of ChiSq(n-p-K)
w.sum<-colSums(w.k.sq.mat)+w.sum2            #1*nsim of ChiSq(n-p)
u<-ifelse(rep(q,nsim)==0,rep(0,nsim),rchisq(nsim,q))   #1*nsim of ChiSq(q)

f.lambdafix.vec<-function(lambda,mu.w,w.k.sq.mat,w.sum2,mu,ksi,n,p)
#vectorized version: lambda as scalar
{
     num<-  colSums((lambda*mu.w)/(1+lambda*mu) )
     den<-  colSums( w.k.sq.mat / (1+lambda*mu) ) + t(w.sum2)
     return( n*log( 1+ (num/den) ) - sum(log(1+lambda*ksi)))
}

    #matrix (K*nsim) (mu_i*w_ij) i=1:K, j=1:nsim
    mu.w<-mu*w.k.sq.mat
    
    lambda.grid    <-c(0,exp(seq(log.grid.lo,log.grid.hi,length=gridlength-1)))
    lrt.array     <-sapply(lambda.grid,f.lambdafix.vec,mu.w,w.k.sq.mat,w.sum2,mu,ksi,n,p) + n*log(1+(u/w.sum)) #grid.length * nsim
    max.index.lrt.array  <-apply(lrt.array,1,which.max) #column-indices of rowmaxima
    lrt.max       <-lrt.array[cbind(1:nsim,max.index.lrt.array)]
    lambda.max     <-lambda.grid[max.index.lrt.array]
    res            <-as.data.frame(cbind(lambda.max,lrt.max))

    colnames(res)<-c("lambda","lrt")

#calculate p-value
return(res)
}#end LRT1Sim

