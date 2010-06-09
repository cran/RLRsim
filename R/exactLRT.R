`exactLRT` <-
function(m, m0, seed = NA, nsim = 10000, 
				log.grid.hi = 8, log.grid.lo = -10, gridlength = 200) 
{
	if (class(m0) != "lm") 
		stop("m0 not an lm-object. \n")
	if (class(m) == "spm") {
		m <- m$fit
		class(m) <- "lme"
	}
	if (class(m) == "amer") 
		class(m) <- "mer"
	if (!((c.m <- class(m)) %in% c("mer", "lme"))) 
		stop("Invalid m specified. \n")
	d <- switch(c.m, lme = extract.lmeDesign(m), mer = extract.lmerDesign(m))
	X <- d$X
	Z <- d$Z
	y <- d$y
	Vr <- d$Vr
	K <- NCOL(Z)
	n <- NROW(X)
	p <- NCOL(X)
	q <- p - length(coefficients(m0)[!is.na(coefficients(m0))])
	if (n != length(m0$fitted)) 
		stop("different data under the null and alternative. \n")
	if (q < 0) 
		stop("m0 not nested in m. \n")
	if (n - p - K < 1) 
		stop("No. of effects greater than n. Reduce model complexity.\n")
	if (q == 0) 
		cat("No restrictions on fixed effects. REML-based inference preferable. \n")
	if ("ML" != switch(c.m, lme = m$method, lmer = switch(m@status[2] + 
							1, "ML", "REML"))) {
		cat("Using likelihood evaluated at REML estimators.\nPlease refit model with method=\"ML\" for exact results.\n")
	}
	#observed value of the LRT
	lrt.obs <- max(0, 2 * logLik(m, REML = FALSE)[1] - 2 * logLik(m0, 
					REML = FALSE)[1])
	sample <- LRTSim(X, Z, q, sqrt.Sigma = chol(cov2cor(Vr)), 
			seed = seed, nsim = nsim, log.grid.hi = log.grid.hi, 
			log.grid.lo = log.grid.lo, gridlength = gridlength)
	if (quantile(sample, 0.9) == 0) 
		(cat("Warning: Null distribution has", mean(sample == 
											0), "mass at zero.\n"))
	p <- mean(lrt.obs < sample)
	RVAL <- list(statistic = c(LRT = lrt.obs), p.value = p, method = paste("simulated finite sample distribution of LRT. (p-value based on", 
					nsim, "simulated values)"))
	class(RVAL) <- "htest"
	return(RVAL)
} 

