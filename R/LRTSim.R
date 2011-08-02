`LRTSim` <-
        function(X,Z,q, sqrt.Sigma, seed=NA, nsim=10000,
                log.grid.hi=8, log.grid.lo=-10, gridlength=200)
{
    K<-NCOL(Z)     # no. of random effects
    n<-NROW(X)     # no. of obs
    p<-NCOL(X)     # no of fixed effects
    
#compute eigenvalues
    mu<- (svd(sqrt.Sigma%*%t(qr.resid(qr(X), Z)), nu = 0, nv = 0)$d)^2
    xi<-(svd(sqrt.Sigma%*%t(Z),nu=0,nv=0)$d)^2
    
#norm eigenvalues
    mu<-mu/max(mu,xi)
    xi<-xi/max(mu,xi)
    
    lambda.grid    <-c(0,exp(seq(log.grid.lo,log.grid.hi,length=gridlength-1)))
    
    if (!is.na(seed)) 
        set.seed(seed)
    
    res <- .C("RLRsim", 
            p=as.integer(p),
            k=as.integer(K),
            n=as.integer(n),
            s=as.integer(nsim),
            g=as.integer(gridlength),
            q=as.integer(q),
            mu =as.double(mu),
            lambda =as.double(lambda.grid),
            lambda0 =as.double(0),
            xi = as.double(xi),
            REML =as.logical(FALSE),
            res = double(nsim),
            lambdaind=as.integer(rep(1,nsim)))
    lambda <- lambda.grid[res$lambdaind] 

    ret <- res$res
    attr(ret, "lambda") <- lambda.grid[res$lambdaind] 
    return(ret)
}

