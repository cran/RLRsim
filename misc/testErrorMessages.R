# TODO: Add comment
# 
# Author: fabians
###############################################################################

n <- 100
x <- runif(n)
g1 <- sample(gl(10, n/10))
g2 <- sample(gl(10, n/10))
y1 <- x+ rnorm(n)
y2 <- rbinom(n, 1, p = plogis(x))



require(lme4)
require(RLRsim)
## multiple random effects
m1 <- lmer(y1 ~ x + (1|g1) + (1|g2))
exactRLRT(m1)
m2 <- lmer(y1 ~ x + (x|g1))
exactRLRT(m2)

## non_Gaussian response

mB <- lmer(y2 ~ x + (1|g1), family="binomial")
exactRLRT(mB)


#################
detach(package::lme4)
require(nlme)
data(Orthodont)


m<-lme(distance ~ Sex, random = ~ age | Subject,
		data = Orthodont, method = "ML")
exactRLRT(m)

g4 <- gl(10, n/10)
g5 <- gl(20, n/20)

m <-lme(y1 ~ x, random = ~ 1| g5/g4, data= data.frame(y1, x, g4, g5))
exactRLRT(m)

m2 <- nlme(y1 ~ x, random = ~1|g1)
exactRLRT(m2)