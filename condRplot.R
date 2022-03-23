########################################################################################################
# Plotting function for conditional repeatability
########################################################################################################
#
# *Usage* 
# condRplot(condR, Rmar=NULL, minVix=NULL, xmin=NULL, ...) {
#
# *Arguments*
# Parameters of the phenotypic equation (fixed or estimated)
# condR:  A data.frame with one column x for covariate values and one column condR for conditional repeatabilities
# Vi:     Optional value for the marginalized between group variance (to be drawn as a horizontal line)
# xmin:	  Optional value for the value of covariate x with minimal between-group variance
# minVix: Optional value for the minimum 
# ...:    Range of covariate values x for which the conditional repeatability is to be calculated
#
# Reference:
# Schielzeth, H. & Nakagawa, S. (2022). Conditional repeatability and the variance explained by reaction 
# norm variation in random slope models. Methods in Ecology and Evolution (in press)
#
# Example:
# Dummy data generation
# ngr    <- 50
# nrep   <- 5
# dat    <- expand.grid(groupID=1:ngr, Rep=1:nrep)
# dat$x  <- runif(ngr*nobs, -2.2, 2.3)
# Sigma  <- matrix(c(1.5, 0.2, 0.2, 0.3), ncol=2)
# rand   <- MASS::mvrnorm(ngr, c(Intercept=0,Slope=0), Sigma)
# dat$y  <- 5 + (1+rep(rand[,"Slope"], nrep)) * dat$x + rep(rand[,"Intercept"], nrep) + rnorm(ngr*nrep, 0, sqrt(2))
# Model fit using lme4::lmer
# mod    <- lme4::lmer(y ~ x + (x|groupID), data=dat)
# Extraction of components
# randest<- data.frame(lme4::VarCorr(mod))
# betaX  <- mod@beta[2]
# meanX  <- mean(dat$x)
# Vx     <- var(dat$x)
# Vr     <- attr(lme4::VarCorr(mod), "sc")
# Sigma  <- lme4::VarCorr(mod)$groupID
# res    <- condR(betaX, meanX, Vx, Vr, Sigma=Sigma, xrange=range(dat$x))
# Plotting
# condRplot(res$condR, res$Vi, res$minVix, res$xmin, lwd=2)

condRplot <- function(condR, Vi=NULL, minVix=NULL, xmin=NULL, ...) {
	Covariate <- condR$x
	conditional_Variance <- condR$condR
	plot(Covariate, conditional_Variance, type="l", ...)
	if(!is.null(Vi)) abline(h=Vi)
	if(!is.null(minVix) & !is.null(xmin)) {
		segments(x0=xmin, x1=xmin, y0=0, y1=minVix)
		segments(x0=-100, x1=xmin, y0=minVix, y1=minVix)
	}
}
