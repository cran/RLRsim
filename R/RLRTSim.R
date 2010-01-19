'RLRTSim' <- function(X, Z, qrX, sqrt.Sigma, lambda0 = NA, 
		seed = NA, nsim = 10000, use.approx = 0, log.grid.hi = 8, 
		log.grid.lo = -10, gridlength = 200) {
	if (is.na(lambda0)) 
		lambda0 <- 0
	#checking args:
	if (!is.numeric(lambda0) | (lambda0 < 0) | length(lambda0) != 
			1) 
		stop("Invalid lambda0 specified. \n")
	if (lambda0 > exp(log.grid.hi)) {
		log.grid.hi <- log(10 * lambda0)
		warning("lambda0 smaller than upper end of grid: \n Setting log.grid.hi to ln(10*lambda0).\n", 
				immediate. = T)
	}
	if ((lambda0 != 0) && (lambda0 < exp(log.grid.lo))) {
		log.grid.lo <- log(-10 * lambda0)
		warning("lambda0 > 0 and larger than lower end of grid: \n Setting log.grid.lo to ln(-10*lambda0).\n", 
				immediate. = T)
	}
	n <- NROW(X)
	p <- NCOL(X)
	K <- min(n, NCOL(Z))
	if (any(is.na(sqrt.Sigma))) 
		sqrt.Sigma <- diag(NCOL(Z))
	mu <- (svd(sqrt.Sigma %*% t(qr.resid(qrX, Z)), nu = 0, 
						nv = 0)$d)^2
	#normalize
	mu <- mu/max(mu)
	if (!is.na(seed)) 
		set.seed(seed)
	if (use.approx) {
		#eigenvalue pattern of balanced ANOVA: mu_s=const for s=1,..,K-1, mu_K approx. 0
		if ((length(unique(round(mu, 6))) == 2) & (1000 * mu[K] < 
					mu[1])) {
			cat("using simplified distribution for balanced ANOVA \n")
			approx.constantmu <- function(nsim, n, p, K, mu) #simplified distribution for balanced ANOVA:
			#mu_s=const for s=1,..,K-1 and mu_K=0
			{
				w.K <- rchisq(nsim, (K - 1))
				w.n <- rchisq(nsim, (n - p - K + 1))
				lambda <- pmax(rep(0, nsim), ((((n - p - K + 
												1)/(K - 1)) * w.K/w.n - 1)/mu[1]))
				rlrt <- rep(0, nsim)
				rlrt[lambda != 0] = ((n - p) * log((w.K + w.n)/(n - 
										p)) - (n - p - K + 1) * log(w.n/(n - p - K + 
										1)) - (K - 1) * log(w.K/(K - 1)))[lambda != 
								0]
				return(cbind(lambda, rlrt))
			}
			res <- approx.constantmu(nsim, n, p, K, mu)
			return(res)
		}
		#eigenvalue pattern for P-splines: exponential decrease
		if (mu[1]/sum(mu) > use.approx) {
			cat("using simplified distribution for 1 single dominating eigenvalue \n")
			approx.scalarmu <- function(nsim, n, p, K, mu) #simplified distribution for B-splines:
			#mu_1 >>> mu_s for s=2,..,K
			{
				mu <- mu[1]
				w.1 <- rchisq(nsim, 1)
				w.n <- rchisq(nsim, (n - p - 1))
				lambda <- pmax(rep(0, nsim), ((((n - p - 1) * 
											w.1)/w.n) - 1)/mu)
				rlrt <- rep(0, nsim)
				rlrt[lambda != 0] = log(((w.1 + w.n)/(n - p))^(n - 
									p)/(w.1 * (w.n/(n - p - 1))^(n - p - 1)))[lambda != 
								0]
				return(cbind(lambda, rlrt))
			}
			res <- approx.scalarmu(nsim, n, p, K, mu)
			return(res)
		}
		#use only first k elements of mu, adapt K<-k accordingly
		#how many eigenvalues are needed to represent at least approx.ratio of the sum of all eigenvalues  (at least 1, of course)
		new.K <- max(sum((cumsum(mu)/sum(mu)) < use.approx), 
				1)
		if (new.K < K) 
			cat(paste("Approximation used:", new.K, "biggest eigenvalues instead of", 
							K, "\n"))
		mu <- mu[1:new.K]
		K <- new.K
	}
	make.lambdagrid <- function(lambda0, gridlength, log.grid.lo, 
			log.grid.hi) #generate symmetric grid around lambda0 that is log-equidistant to the right,
	{
		if (lambda0 == 0) 
			return(c(0, exp(seq(log.grid.lo, log.grid.hi, length = gridlength - 
													1))))
		else {
			leftratio <- min(max((log(lambda0)/((log.grid.hi) - 
											(log.grid.lo))), 0.2), 0.8)
			leftlength <- max(round(leftratio * gridlength) - 
							1, 2)
			leftdistance <- lambda0 - exp(log.grid.lo)
			#make sure leftlength doesn't split the left side into too small parts:
			if (leftdistance < (leftlength * 10 * .Machine$double.eps)) {
				leftlength <- max(round(leftdistance/(10 * .Machine$double.eps)), 
						2)
			}
			#leftdistance approx. 1 ==> make a regular grid, since (1 +- epsilon)^((1:n)/n) makes a too concentrated grid
			if (abs(leftdistance - 1) < 0.3) {
				leftgrid <- seq(exp(log.grid.lo), lambda0, length = leftlength + 
								1)[-(leftlength + 1)]
			}
			else {
				leftdiffs <- ifelse(rep(leftdistance > 1, leftlength - 
										1), leftdistance^((2:leftlength)/leftlength) - 
								leftdistance^(1/leftlength), leftdistance^((leftlength - 
										1):1) - leftdistance^(leftlength))
				leftgrid <- lambda0 - rev(leftdiffs)
			}
			rightlength <- gridlength - leftlength
			rightdistance <- exp(log.grid.hi) - lambda0
			rightdiffs <- rightdistance^((2:rightlength)/rightlength) - 
					rightdistance^(1/rightlength)
			rightgrid <- lambda0 + rightdiffs
			return(c(0, leftgrid, lambda0, rightgrid))
		}
	}
	lambda.grid <- make.lambdagrid(lambda0, gridlength, log.grid.lo = log.grid.lo, 
			log.grid.hi = log.grid.hi)
	res <- .C("RLRsim", p = as.integer(p), k = as.integer(K), 
			n = as.integer(n), s = as.integer(nsim), g = as.integer(gridlength), 
			q = as.integer(0), mu = as.double(mu), lambda = as.double(lambda.grid), 
			lambda0 = as.double(lambda0), xi = as.double(mu), REML = as.logical(TRUE), 
			res = double(nsim))
	return(res$res)
} 