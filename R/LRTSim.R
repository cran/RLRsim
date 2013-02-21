`LRTSim` <- function(X,Z,q, sqrt.Sigma, seed=NA, nsim=10000,
           log.grid.hi=8, log.grid.lo=-10, gridlength=200,
           parallel = c("no", "multicore", "snow"), 
           ncpus = 1L, cl = NULL){
    
    parallel <- match.arg(parallel)
    have_mc <- have_snow <- FALSE
    if (parallel != "no" && ncpus > 1L) {
      if (parallel == "multicore") 
        have_mc <- .Platform$OS.type != "windows"
      else if (parallel == "snow") 
        have_snow <- TRUE
      if (!have_mc && !have_snow) 
        ncpus <- 1L
    }      
    
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
    
    res <- if (ncpus > 1L && (have_mc || have_snow)) {
      nsim. <- as.integer(ceiling(nsim/ncpus))
      if (have_mc) {
        tmp <- parallel::mclapply(seq_len(ncpus), function(i){
          .C("RLRsim", p = as.integer(p), k = as.integer(K), 
             n = as.integer(n), s = as.integer(nsim.), g = as.integer(gridlength), 
             q = as.integer(0), mu = as.double(mu), lambda = as.double(lambda.grid), 
             lambda0 = as.double(0), xi = as.double(mu), REML = as.logical(FALSE), 
             res = double(nsim.), lambdaind=as.integer(rep(1,nsim.)))
        }, mc.cores = ncpus)
        do.call(mapply, c(tmp, FUN=c))
      } else { 
        if (have_snow) {
          if (is.null(cl)) {
            cl <- parallel::makePSOCKcluster(rep("localhost", 
                                                 ncpus))
            if (RNGkind()[1L] == "L'Ecuyer-CMRG") {
              parallel::clusterSetRNGStream(cl)
            }
            tmp <- parallel::parLapply(cl, seq_len(ncpus), function(i){
              .C("RLRsim", p = as.integer(p), k = as.integer(K), 
                 n = as.integer(n), s = as.integer(nsim.), g = as.integer(gridlength), 
                 q = as.integer(0), mu = as.double(mu), lambda = as.double(lambda.grid), 
                 lambda0 = as.double(0), xi = as.double(mu), REML = as.logical(FALSE), 
                 res = double(nsim.), lambdaind=as.integer(rep(1,nsim.)))
            })
            parallel::stopCluster(cl)
            do.call(mapply, c(tmp, FUN=c))
          } else {
            tmp <- parallel::parLapply(cl, seq_len(ncpus), function(i){
              .C("RLRsim", p = as.integer(p), k = as.integer(K), 
                 n = as.integer(n), s = as.integer(nsim.), g = as.integer(gridlength), 
                 q = as.integer(0), mu = as.double(mu), lambda = as.double(lambda.grid), 
                 lambda0 = as.double(0), xi = as.double(mu), REML = as.logical(FALSE), 
                 res = double(nsim.), lambdaind=as.integer(rep(1,nsim.)))
            })
            do.call(mapply, c(tmp, FUN=c))
          }  
        }
      }
    } else {
      .C("RLRsim", p = as.integer(p), k = as.integer(K), 
         n = as.integer(n), s = as.integer(nsim), g = as.integer(gridlength), 
         q = as.integer(0), mu = as.double(mu), lambda = as.double(lambda.grid), 
         lambda0 = as.double(0), xi = as.double(mu), REML = as.logical(FALSE), 
         res = double(nsim), lambdaind=as.integer(rep(1,nsim)))
    }
    
    
    
    lambda <- lambda.grid[res$lambdaind] 
    
    ret <- res$res
    attr(ret, "lambda") <- lambda.grid[res$lambdaind] 
    return(ret)
  }

