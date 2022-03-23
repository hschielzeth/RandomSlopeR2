########################################################################################################
# Function for calculations of conditional repeatability, random-slope R2 and conditional repeatability
########################################################################################################
#
# *Usage* 
# condR(betaX, meanX, Vx, Vr, Vu, Vv, Cuv, Sigma=NULL, Vo=NULL, Vp=NULL, xrange=c(-2,2))
#
# *Arguments*
# Parameters of the phenotypic equation (fixed or estimated)
# betaX:  Slope of y on x
# meanX:  Mean of covariate x
# Vx:     Variance in covariate x
# Vr:     Residual variance
# Vu:     Random-intercept variance
# Vv:     Random-slope variance
# Cuv:    Random intercept-slope covariance
# Sigma:  A positive definite matrix of Vu (in cell [1,1]), Vv (in cell [2,2]) and Cuv (in cells [2,1] and [1,2])
#         will be prioritized over values given for Vu, Vv and Cuv
# Vo:     Optional value or vector of other variance components 
# Vp:     Total variance in response y (if NULL, Vp will be estimated from the components provided
# xrange: Range of covariate values x for which the conditional repeatability is to be plotted (calculations are
#         done based on Vx
##
# *Value* 
# A list of variance components and ratios
# Vf:     Variance explained by covariate x
# Vi:     Average between-group variance across x
# Vs:     Variance explaind by random slopes
# R2s:    R2 for random slopes (in proportion of Vp)
# Rmar:   Marginalized repeatability (between-group variance across x in propoportion of Vp)
# minVix: Minimum value of the between-group variance (at xmin)
# xmin:	  Value of covariate x where the between-group variance is at its minimum
# CondR:  A data.frame of covariate values and conditional repeatabilities
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
# res


condR <- function(betaX, meanX, Vx, Vr, Vu, Vv, Cuv, Sigma=NULL, Vp=NULL, Vo=NULL, xrange=c(-2,2)) {
	if(!is.null(Sigma)) {
		Vu  <- Sigma[1,1]
		Cuv <- Sigma[1,2]
		Vv  <- Sigma[2,2]
	}
	if(is.null(Vo)) Vo <- 0
	Vo     <- sum(Vo)
	Vf     <- betaX^2 * Vx
	Vi     <- Vu + Vv*Vx + meanX^2*Vv + 2*meanX*Cuv
	Vs     <- Vv*Vx + meanX^2 * Vv
	if(is.null(Vp))
		Vp <- Vf + Vi + Vr + Vo
	R2s    <- Vs / Vp
	Rmar   <- Vi / Vp
	xmin   <- -Cuv/Vv
	minVix <- Vu - Cuv^2/Vv
	x 	   <- seq(xrange[1], xrange[2], length.out=501)
	condR  <- data.frame(x=x, condR=Vu + 2*x*Cuv + x^2*Vv)
	cat("Variance explained by fixed effects\nVf = ", Vf, "\n",
	    "Average between-group variance\nVi = ", Vi, "\n",
	    "Variance explaind by random slopes\nVs = ", Vs, "\n",
	    "R2 for random slopes\nR2s = ", R2s, "\n",
	    "Marginalized repeatability\nRmar = ", Rmar, "\n",
	    "Minimum value of the between-group variance\nminVix = ", minVix, "\n",
	    "Covariate value at minimum between-group variance\nxmin = ", xmin, "\n", sep = "")
	return(list(Vf=Vf, Vi=Vi, Vs=Vs, R2s=R2s, Rmar=Rmar, minVix=minVix, xmin=xmin, condR=condR))
}



