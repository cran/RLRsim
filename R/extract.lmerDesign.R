`extract.lmerDesign` <-
  function(m)
  {
    ## code by Martin Maechler 
    ## dl'ed @ https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
    ## 21 Mar 2012
    packageLoaded <- function(name)
    {
      (paste("package:", name, sep="") %in% search()) ||
        (name %in% loadedNamespaces())
    }   
    
    if(packageLoaded("lme4")) {
      varcorr <- lme4::VarCorr
    } else {
      if(packageLoaded("lme4.0")) {
         stop("RLRsim no longer supports models fit with lme4.0.")
      }    
    } 
    
    
    X<-m@X
    Z<-t(as.matrix(m@Zt))
    #Cov(b)/ Var(Error)
    Sigma.l<- lapply(.Call(lme4:::mer_ST_chol, m),crossprod)
    k <- length(Sigma.l) #how many grouping factors
    q <- lapply(Sigma.l,NROW) #how many variance components in each grouping factor
    nlevel<-sapply(m@flist, function(x) length(levels(x))) #how many inner blocks in Sigma_i
    Vr <- matrix(0,NCOL(Z),NCOL(Z)) #Cov(RanEf)/Var(Error)
    from <- 1        
    for(i in 1:k)
    {
      ii<-nlevel[i]
      inner.block<-as.matrix(Sigma.l[[i]])
      to<-from-1+ii*NCOL(inner.block)
      Vr[from:to,from:to]<- inner.block %x% diag(ii)  
      from<-to+1
    }
    return(list(
      Vr=Vr, #Cov(RanEf)/Var(Error)
      X=X,
      Z=Z,
      sigmasq=attributes(varcorr(m))$sc^2,
      lambda=unique(diag(Vr)),
      y=as.numeric(m@y),
      k=k
    )
    )
  }

