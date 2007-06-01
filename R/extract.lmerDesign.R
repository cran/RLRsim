`extract.lmerDesign` <-
function(m)
{
  X<-m@X
  Z<-as.matrix(t(m@Zt))
  Om<-m@Omega
  ngroup<-length(m@flist) #how many groups/outer blocks in Sigma
  nlevel<-numeric(ngroup) #how many levels per group / inner blocks in Sigma_i
  Vr<-matrix(0,NCOL(Z),NCOL(Z)) #Cov(RanEf)/Var(Error)
  from<-1
  for(i in 1:ngroup)
      {
          ii<-nlevel[i]<-length(attr(m@flist[[i]],"l"))
          inner.block<-solve(as.matrix(Om[[i]]))
          to<-from-1+ii*NCOL(inner.block)
          Vr[from:to,from:to]<-diag(ii) %x% inner.block
          from<-to+1
      }
   Vlambda<-diag(nrow(X))+ Z %*% Vr %*% t(Z)
   return(list(
             Vlambda=Vlambda, #Cov(y)/Var(Error)
             Vr=Vr, #Cov(RanEf)/Var(Error)
             X=X,
             Z=Z,
             sigmasq=attributes(VarCorr(m))$sc^2,
             lambda=unique(diag(Vr)),
             y=as.numeric(m@y)
           )
      )
}

