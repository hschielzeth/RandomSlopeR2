#####################################################################################
# varDecomp
# A function for calculating variance components from fixed parameters or model estimates
#####################################################################################
# Arguments: Parameters of the phenotypic equation (fixed or estimated)
# Value:     Vector of variance components
varDecomp = function(alpha, betaX, betaX2, meanX, meanX2, Vx, Vx2, Vr, Vu, Vv, Cuv, Xmat=NULL, vect=NULL, msteps=1) {
	if(!is.null(vect)) {
		alpha  = as.numeric(vect["alpha"])
		betaX  = as.numeric(vect["betaX"])
		betaX2 = as.numeric(vect["betaX2"])
		meanX  = as.numeric(vect["meanX"])
		meanX2 = as.numeric(vect["meanX2"])
		Vx     = as.numeric(vect["Vx"])
		Vx2    = as.numeric(vect["Vx2"])
		Vr     = as.numeric(vect["Vr"])
		Vu     = as.numeric(vect["Vu"])
		Vv     = as.numeric(vect["Vv"])
		Cuv    = as.numeric(vect["Cuv"])
	}
	Sigma  = matrix(c(Vu, Cuv, Cuv, Vv), ncol=2)
	Vf     = betaX^2 * Vx
	Vi     = Vu + Vv*Vx + meanX^2 * Vv
	Vr     = Vr
	Vo     = betaX2^2 * Vx2
	Vp     = Vf + Vi + Vr + Vo
	Vs     = Vv*Vx + meanX^2 * Vv
	R2s    = Vs / Vp
	Rmar   = Vi / Vp
	xmin   = -Cuv/Vv
	minVix = Vu - Cuv^2/Vv
	if(!is.null(Xmat)) { # Johnson 2014, in two steps, because of memory issues with large datasets
		if(msteps==1) meanVi = mean(diag(Xmat%*%Sigma%*%t(Xmat))) 
		if(msteps==2) {		
			temp   = Xmat%*%Sigma
			meanVi = mean(diag(temp%*%t(Xmat))) 
		}
	}
	if( is.null(Xmat)) meanVi = NULL
	res = c(Vf=Vf, Vi=Vi, Vo=Vo, Vp=Vp, Vs=Vs, R2s=R2s, Rmar=Rmar, xmin=xmin, minVix=minVix, meanVi=meanVi)
	return(res)
}

#####################################################################################
# conditionalR
# A function for calculating conditional between-individual variation
#####################################################################################
# Arguments: Parameters of the phenotypic equation (fixed or estimated) 
#            x: values of x to be evaluated 
# Value:     Vector of conditional between-individual variation for all x values
conditionalR = function(x, meanX, Vu, Vv, Cuv, vect=NULL) {
	if(!is.null(vect)) {
		Vu  = as.numeric(vect["Vu"])
		Vv  = as.numeric(vect["Vv"])
		Cuv = as.numeric(vect["Cuv"])
	}
	Vix = Vu + 2*x*Cuv + x^2*Vv 
	return(Vix)
}

#####################################################################################
# estExtr
# A function for extracting parameter estimates from a fitted model
#####################################################################################
# Arguments: A fitted model and the original dataframe
#            randSlopeOmitted: specifies if random-slope have been omitted in model fitting
# Value:     Vector of extracted parameter estimates
estExtr = function(mod, md, randSlopeOmitted) {
	ests = summary(mod)$coefficients
	alpha  = ests["(Intercept)", "Estimate"]
	betaX  = ests["x", "Estimate"]
	betaX2 = ests["x2", "Estimate"]
	meanX  = mean(md$x)
	meanX2 = mean(md$x2)
	Vx     = var(md$x)
	Vx2    = var(md$x2)
	if(!randSlopeOmitted) {
		varcors = as.data.frame(summary(mod)$varcor)
		Vu  = varcors$vcov[varcors$grp=="IndID" & varcors$var1=="(Intercept)" & is.na(varcors$var2)]
		Vv  = varcors$vcov[varcors$grp=="IndID" & varcors$var1=="x"]
		Cuv = varcors$vcov[varcors$grp=="IndID" & varcors$var1=="(Intercept)" & !is.na(varcors$var2) & varcors$var2=="x"]
		Vr  = varcors$vcov[varcors$grp=="Residual"]
	}
	if(randSlopeOmitted) {
		varcors = as.data.frame(summary(mod)$varcor)
		Vu  = varcors$vcov[varcors$grp=="IndID" & varcors$var1=="(Intercept)" & is.na(varcors$var2)]
		Vv  = 0
		Cuv = 0
		Vr  = varcors$vcov[varcors$grp=="Residual"]
	}
	res = c(alpha=alpha, betaX=betaX, betaX2=betaX2, meanX=meanX, meanX2=meanX2, Vx=Vx, Vx2=Vx2, Vr=Vr, Vu=Vu, Vv=Vv, Cuv=Cuv)
	return(res)
}

#####################################################################################
# simData
# A function for simulating data from a phenotypic equation
#####################################################################################
# Arguments: Parameters of the phenotypic equation and target sample sizes
#            errX: gives the proportional error in covariate x
#            balanced: states if sampling of observations across individuals is balanced
#            xdist: gives the target distribution of fixed effects (three currently implemented)
# Details:   (a) x values are first generated from the target distribution specified by xdist
#                in a next step they are transformed to target meanX and Vx
#            (b) Covariate x is an observation-level predictor
#            (c) Covariate x2 is an individual-level predictor
#            (d) Error in x is introduced such that the target outcome var(x) is Vx, but the signal
#                in x with respect to the response is reduced to Vx*(1-errX)
#                the error is thus proportional to the variance in x
#                Response values are generated from the x without error (variance Vx*(1-errX)) and
#                error to x is added in a separate step
# Value:     A data frame with simulated values (for the response, covariates and random effect levels)
#            xOrg in the return object refers to the error-free covariate used for generating the response
#            sloInd in the return object refers to the simulated random slope for each individual
simData = function(nind, nrep, alpha, betaX, betaX2, meanX, meanX2, Vx, Vx2, Vr, Vu, Vv, Cuv, errX=0, balanced=FALSE, xdist=c("normal", "uniform", "lognorm"), vect=NULL) {
	if(!is.null(vect)) {
		nind     = as.numeric(vect["nind"])
		nrep     = as.numeric(vect["nrep"])
		alpha    = as.numeric(vect["alpha"])
		betaX    = as.numeric(vect["betaX"])
		betaX2   = as.numeric(vect["betaX2"])
		meanX    = as.numeric(vect["meanX"])
		meanX2   = as.numeric(vect["meanX2"])
		Vx       = as.numeric(vect["Vx"])
		Vx2      = as.numeric(vect["Vx2"])
		Vr       = as.numeric(vect["Vr"])
		Vu       = as.numeric(vect["Vu"])
		Vv       = as.numeric(vect["Vv"])
		Cuv      = as.numeric(vect["Cuv"])
		errX     = as.numeric(vect["errX"])
		xdist    = as.character(vect["xdist"])
		balanced = as.logical(vect["balanced"])
	}
	if(balanced) {
		md = expand.grid(IndID=1:nind,Obs=1:nrep)
		md$ObsID = paste(md$IndID, md$Obs, sep="_")
	}
	if(!balanced) {
		IndID = c(1:nind, sample(1:nind, nind*nrep-nind, replace=TRUE))
		ObsID = 1:length(IndID)
		md = data.frame(IndID=IndID, ObsID=ObsID)
	}
	Sigma    = matrix(c(Vu, Cuv, Cuv, Vv), ncol=2)
	randcomp = MASS::mvrnorm(nind, c(0,0), Sigma)
	randint  = randcomp[,1]
	randslo  = randcomp[,2]
	md$x2 = rnorm(nind, meanX2, sqrt(Vx2))[md$IndID] # between-individual covariate
	if(xdist=="normal") {
		md$xOrg = rnorm(nind*nrep, meanX, sqrt(Vx*(1-errX)))
	}
	if(xdist=="uniform") {
		xasifnormal = rnorm(nind*nrep, meanX, sqrt(Vx*(1-errX)))
		md$xOrg = runif(nind*nrep, 0, 1)
		md$xOrg = scale(md$xOrg) * sd(xasifnormal) + mean(xasifnormal)

	}
	if(xdist=="lognorm") {
		xasifnormal = rnorm(nind*nrep, meanX, sqrt(Vx*(1-errX)))
		md$xOrg = rlnorm(nind*nrep, 0, 1)
		md$xOrg = scale(md$xOrg) * sd(xasifnormal) + mean(xasifnormal)
	}
	md$x            = md$xOrg + rnorm(nind*nrep, 0, sqrt(Vx*errX))
	md$sloInd       = randslo[md$IndID]
	md$yMeanIndComp = alpha + randint[md$IndID]
	md$ySlopeCompX  = md$xOrg * betaX
	md$ySlopeCompX2 = md$x2 * betaX2
	md$yRandSloComp = md$xOrg * md$sloInd
	md$yResComp     = rnorm(nind*nrep, 0, sqrt(Vr))
	md$yPred        = md$yMeanIndComp + md$ySlopeCompX + md$ySlopeCompX2 + md$yRandSloComp
	md$y            = md$yPred + md$yResComp
	return(md)
}

#####################################################################################
# modFit
# Fist a mixed effects model to the data
#####################################################################################
# Arguments: A data frame and a binary indicatior whether random slope should be fitted
# Value:     A fitted mixed effects model
modFit = function(md, randSlopeOmitted=FALSE) {
	if( randSlopeOmitted) mod = lme4::lmer(y ~ x + x2 + (1|IndID), data=md)
	if(!randSlopeOmitted) mod = lme4::lmer(y ~ x + x2 + (x|IndID), data=md)
	return(mod)
}
