`extract.lmerDesign` <-
function(m)
{
  X<-m@X
  Z<-as.matrix(t(m@Zt))
  Sigma.l<- lapply(.Call(lme4:::mer_ST_chol, m),crossprod) #Cov(b)/ Var(Error)
  k <- length(Sigma.l) #how many grouping factors
  q <- lapply(Sigma.l,NROW) #how many variance components in each grouping factor
  
  nlevel<-18#numeric(ngroup) #how many levels per group / inner blocks in Sigma_i
  Vr<-matrix(0,NCOL(Z),NCOL(Z)) #Cov(RanEf)/Var(Error)
  from<-1
  for(i in 1:k)
      {
          ii<-nlevel[i]#<-length(attr(m@flist[[i]],"l"))
          inner.block<-as.matrix(Sigma.l[[i]])
          to<-from-1+ii*NCOL(inner.block)
          Vr[from:to,from:to]<- inner.block %x% diag(ii)  
          from<-to+1
      }
   return(list(
             Vr=Vr, #Cov(RanEf)/Var(Error)
             X=X,
             Z=Z,
             sigmasq=attributes(VarCorr(m))$sc^2,
             lambda=unique(diag(Vr)),
             y=as.numeric(m@y)
           )
      )
}

