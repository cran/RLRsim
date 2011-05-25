'exactRLRT' <- function(m, mA = NULL, m0 = NULL, seed = NA, 
		nsim = 10000, log.grid.hi = 8, log.grid.lo = -10, gridlength = 200) {
	if (class(m) == "spm") {
		m <- m$fit
		class(m) <- "lme"
	}
	if (class(m) == "amer") class(m) <- "mer"
	if (!(c.m <- (class(m))) %in% c("mer", "lme")) 
		stop("Invalid m specified. \n")
	if(c.m == "mer"){
		if(length(m@muEta))
			stop("exactRLRT can only be used for mixed models for Gaussian responses.")
	}
	if ("REML" != switch(c.m, 
			lme = m$method, 
			mer = switch(m@dims["REML"] + 1, "ML", "REML"))){
		cat("Using restricted likelihood evaluated at ML estimators.\n Refit with method=\"REML\" for exact results.\n")	
	}
	
	d <- switch(c.m, lme = extract.lmeDesign(m), mer = extract.lmerDesign(m))
	X <- d$X
	qrX <- qr(X)
	Z <- d$Z
	y <- d$y
	Vr <- d$Vr
	K <- ncol(Z)
	n <- nrow(X)
	p <- ncol(X)
	if (is.null(mA) && is.null(m0)) {
		if(length(d$lambda) != 1 || d$k != 1) 
			stop("multiple random effects in model - exactRLRT needs 'm' with only a single random effect.")
		#2*restricted ProfileLogLik under H0: lambda=0
		res <- qr.resid(qrX, y)
		R <- qr.R(qrX)
		detXtX <- det(t(R) %*% R)
		reml.H0 <- -((n - p) * log(2 * pi) + (n - p) * log(sum(res^2)) + 
					log(detXtX) + (n - p) - (n - p) * log(n - p))
		#observed value of the test-statistic
		reml.obs <- 2 * logLik(m, REML = TRUE)[1]
		rlrt.obs <- max(0, reml.obs - reml.H0)
		lambda <- d$lambda
	}
	else {
		if (c.m == "lme") {
			if (any(mA$fixDF$terms != m0$fixDF$terms)) 
				stop("Fixed effects structures of mA and m0 not identical.\n REML-based inference not appropriate.")
		}
		else {
			if (any(mA@X != m0@X)) 
				stop("Fixed effects structures of mA and m0 not identical.\n REML-based inference not appropriate.")
		}
		## bug fix submitted by Andrzej Galecki 3/10/2009
		DFx <- switch(c.m, lme = anova(mA,m0)$df, mer = anova(mA,m0)$Df) 
		if (abs(diff(DFx)) > 1) {
			stop("Random effects not independent - covariance(s) set to 0 under the null hypothesis.\n exactRLRT can only test a single variance.\n")
		}
		rlrt.obs <- max(0, 2 * (logLik(mA, REML = TRUE)[1] - 
							logLik(m0, REML = TRUE)[1]))
	}
	if (rlrt.obs != 0) {
		sample <- RLRTSim(X, Z, qrX=qrX, sqrt.Sigma = chol(cov2cor(Vr)), 
				lambda0 = 0, seed = seed, nsim = nsim, log.grid.hi = log.grid.hi, 
				log.grid.lo = log.grid.lo, gridlength = gridlength)
		p <- mean(rlrt.obs < sample)
	}
	else p = 1
	RVAL <- list(statistic = c(RLRT = rlrt.obs), p.value = p, 
			method = paste("simulated finite sample distribution of RLRT.\n (p-value based on", 
					nsim, "simulated values)"))
	class(RVAL) <- "htest"
	return(RVAL)
} 