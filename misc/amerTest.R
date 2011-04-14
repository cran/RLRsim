# generate data with group-specific smooth trajectories
# similar to data in M. Durban et al 2005

install.packages("~/RLRsim_source/RLRsim_2.0-3.tar.gz")

library(lme4)
library(RLRsim)
nsubjects <- 2000
nmeasures <- 5
n <- nsubjects  * nmeasures
subject <- gl(n=nsubjects, k=n/nsubjects, labels=1:nsubjects)
x <- jitter(rep(seq(-1, 1, l=n/nsubjects), times=nsubjects))
i <- 1
b_subj <- rep(as.numeric(levels(subject))/10, e=nmeasures)
y <- b_subj + x + rnorm(n)
m <- lmer(y ~ x + (1|subject))
exactRLRT(m)

X <- matrix(rnorm(20000), 1000, 20)
replicate(100, 	system.time((svd(X, nu = 0, nv = 0)$d)^2)[3])
replicate(100, 	system.time((eigen(tcrossprod(X), symmetric= TRUE, only.values = TRUE)))[3])

if(require(amer)){
	nsubjects <- 60
	nmeasures <- 15
	n <- nsubjects  * nmeasures
	ntreatments <- 4
	snr <- 15
	
	x <- jitter(rep(seq(-1, 1, l=n/nsubjects), times=nsubjects))
	subject <- gl(n=nsubjects, k=n/nsubjects, labels=1:nsubjects)
	treatment <- rep(1:ntreatments, each = n/ntreatments)
	f_treatment <- rep(NA, n)
	f_treatment[treatment == 1] <- 0
	f_treatment[treatment == 2] <- 2*x[treatment == 2]+2
	f_treatment[treatment == 3] <- (x[treatment == 3]+2)^2
	f_treatment[treatment == 4] <- exp(1.5*x[treatment == 4])
 	treatment <- factor(treatment)
	
	i <- 1
	f_subj <- unlist(tapply(x, subject, function(x){
						ff <- rnorm(1) + rnorm(1)*x + 2*sin(2*pi*x + i/10)
						i <<- i+1
						return(ff)
					}))  
	
	
	f_subj <- scale(f_subj)/2
	f_treatment <- 2*scale(f_treatment)
	
	xyplot(f_subj~x, groups=subject, type="l")
	xyplot(f_treatment~x, groups=treatment)
	
	
	f <- f_subj + f_treatment
	
	y <-  f + sqrt(var(f)/snr)*rnorm(n)
	d <- data.frame(y, x=scale(x), subject, treatment)
		
	(m1 <- amer(y ~ tp(x, by = treatment, k = 10) + tp(x, by = subject, k = 10, allPen = T), data = d))
	plotF(m1, legendPos="none")
	
	(m0 <- amer(y ~ tp(x, by = treatment, k = 10), data = d))
	(m <- amer(y ~ tp(x, by = subject, k = 10, allPen = T), data = d))
	
	require(RLRsim)
	debug(exactRLRT)
	
}

