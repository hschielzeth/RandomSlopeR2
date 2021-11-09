####################################################################################################################
# Simulation
####################################################################################################################

# R version 4.1.2 (2021-11-01)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Matrix products: default
# attached base packages:
# stats     graphics  grDevices utils     datasets  methods   base     
# lme4_1.1-27.1    Matrix_1.3-4     matrixcalc_1.0-5 psych_2.1.9     
# loaded via a namespace (and not attached):
# Rcpp_1.0.7      lattice_0.20-45 MASS_7.3-54     grid_4.1.2      nlme_3.1-153    minqa_1.2.4    
# nloptr_1.2.2.2  boot_1.3-28     splines_4.1.2   tools_4.1.2     parallel_4.1.2  compiler_4.1.2 
# mnormt_2.0.2    tmvnsim_1.0-2  

rm(list=ls())
require(psych)
require(matrixcalc)
require(lme4)
sessionInfo()

path = "C:\\Users\\Holger\\Documents\\A_PaperDrafts\\B_Submitted\\HS_RandomSlopeR2\\Submitted\\Submitted to MethEcolEvol\\"

source(paste0(path, "RandSlopeR2_SimulationFunctions.R"))
source(paste0(path, "RandSlopeR2_PlottingFunctions.R"))

scenario = 0
#scenario = 1 # Random slopes not fitted
#scenario = 2 # Small sample size
#scenario = 3 # Balanced sampling of observations within individuals
#scenario = 4 # 30% error in X
#scenario = 5 # normal x distribution
#scenario = 6 # lognormal x distribution

baseScenario = function() {
	# Sampling
	assign("nind",    60, .GlobalEnv)
	assign("nrep",    10, .GlobalEnv)
	# Fixed effects
	assign("alpha",    3, .GlobalEnv)
	assign("betaX", -0.5, .GlobalEnv)
	assign("betaX2", 0.5, .GlobalEnv)
	# Variances
	assign("Vr",       1, .GlobalEnv)
	assign("Vu",       1, .GlobalEnv)
	assign("Vv",     0.5, .GlobalEnv)
	assign("Cuv",    0.3, .GlobalEnv) 
	# Means and variance in covarate
	assign("meanX",  0.5, .GlobalEnv) 
	assign("meanX2",-0.5, .GlobalEnv) 
	assign("Vx",     1.2, .GlobalEnv) 
	assign("Vx2",      1, .GlobalEnv) 
	# Other settings
	assign("errX",     0, .GlobalEnv) 
	assign("xdist", "uniform", .GlobalEnv) 
	assign("balanced", FALSE, .GlobalEnv) 
	assign("randSlopeOmitted", FALSE, .GlobalEnv) 
}

# Set parameters for simulation
baseScenario()
if(scenario == 1) { # Random slopes not fitted
	randSlopeOmitted = TRUE }
if(scenario == 2) { # Small sample size
	nind = 30; nrep = 3}
if(scenario == 3) { # Balanced
	balanced = TRUE }
if(scenario == 4) { # error in X
	errX=0.3 }
if(scenario == 5) { # normal x distribution
	xdist="normal" }
if(scenario == 6) { # lognormal x distribution
	xdist="lognorm" }

# Derived quantities
Sigma = matrix(c(Vu, Cuv, Cuv, Vv), ncol=2)
if(!matrixcalc::is.positive.definite(Sigma)) cat("Sigma is not positive definite!\n")
nobs     = nind*nrep
settings = c(nind=nind, nrep=nrep, alpha=alpha, betaX=betaX, betaX2=betaX2, meanX=meanX, meanX2=meanX2, Vx=Vx, Vx2=Vx2, Vr=Vr, Vu=Vu, Vv=Vv, Cuv=Cuv, 
			 randSlopeOmitted=randSlopeOmitted, errX=errX, balanced=balanced, xdist=xdist)

# Sim calculations
condRsim    = varDecomp(vect=settings)
allsim      = c(settings, condRsim, covxmin=NA, covxmax=NA)

# Simulation, model fittings and extraction
md        = simData(vect=settings)
mod       = modFit(md, randSlopeOmitted)
ests      = estExtr(mod, md, randSlopeOmitted)
condRests = varDecomp(vect=ests)
allests   = c(ests, condRests, covxmin=min(md$x), covxmax=max(md$x))

# Visual check for whether condition R is a good fit to 
# between-individual variance estimated for bins of x
windows()
plotCheck(md, settings, ests, binfact=3, ptscale=5, las=1)

# Preparation for looping
nsim=500
resSim = allsim
resEst = allests
for(i in 2:nsim) resEst = rbind(resEst, rep(NA, length(allests))) 
rownames(resEst) = NULL

# Setup of pregressive plot for simulation
finalOrd = rev(c("alpha", "betaX", "meanX", "Vx", "betaX2", "meanX2", "Vx2", 
	"Vu", "Vv", "Cuv", "Vp", "Vr", "Vo", "Vf", "Vi", "Vs", "Rmar", "R2s", "xmin", "minVix"))
ylabs = rev(c(expression(alpha), expression(beta[X]), expression(mu[X]), expression("V"[X]), 
		expression(beta[X2]), expression(mu[X2]), expression("V"[X2]), 
		expression("V"[u]), expression("V"[v]), expression("C"[uv]), 
		expression("V"[P]), expression("V"[R]), expression("V"[phi]), expression("V"[F]), expression("V"[I]), expression("V"[S]), 
		expression("R"[mar]), expression(paste("R"^2)[S]), expression("x"[min]), expression(paste("min(V"["I,x"],") "))))
windows()
par(mfrow=c(1,1), mar=c(5,5,1,1))
plot((allests[finalOrd]-as.numeric(allsim[finalOrd]))/abs(as.numeric(allsim[finalOrd])), jitter(1:length(finalOrd)), ylab="", yaxt="n", ylim=c(1,length(finalOrd)), xlab="", xlim=c(-0.4, 0.4))
title(xlab="Proportional difference to simulated values\n(Est-Sim)/Sim", line=3.5)
axis(2, at=1:length(finalOrd), labels=ylabs, las=1)
abline(v=0)
abline(h=1:I(length(finalOrd)+1)-0.5)
abline(h=c(1,3,5,10,11,14,17,20,21)-0.5, lwd=2)
for(i in 2:nsim) {
	print(i)
	md = simData(vect=settings)
	mod = modFit(md, randSlopeOmitted)
	ests = estExtr(mod, md, randSlopeOmitted)
	condRests = varDecomp(vect=ests)
	allests = c(ests, condRests, covxmin=min(md$x), covxmax=max(md$x))
	resEst[i,] = allests
	points((allests[finalOrd]-as.numeric(allsim[finalOrd]))/abs(as.numeric(allsim[finalOrd])), jitter(1:length(finalOrd)))
}

# Figure 3
windows()
par(mfrow=c(1,1), mar=c(5,5,1.4,4))
plot((allests[finalOrd]-as.numeric(allsim[finalOrd]))/abs(as.numeric(allsim[finalOrd])), jitter(1:length(finalOrd)), ylab="", type="n", yaxt="n", ylim=c(1,length(finalOrd)), xlab="", xlim=c(-1.1, 1.4))
title(xlab="Proportional difference to simulated values\n(Est-Sim)/Sim", line=3.5)
axis(2, at=1:length(finalOrd), labels=ylabs, las=1)
axis(4, at=1:length(finalOrd), labels=paste0("= ", round(as.numeric(allsim[finalOrd]),2)), las=1, tick=FALSE, line=-0.3)
mtext("Sim", 3, at=1.385, line=-0.1, adj=2)
abline(v=0, lwd=3)
abline(h=1:I(length(finalOrd)+1)-0.5)
abline(h=c(1,3,5,10,11,14,17,20,21)-0.5, lwd=2)
for(i in 1:length(finalOrd)) {
	devs = (resEst[,finalOrd[i]]-as.numeric(allsim[finalOrd[i]])) / as.numeric(allsim[finalOrd[i]])
	rect(xleft=min(devs), xright=max(devs), ybottom=i-0.5, ytop=i+0.5, col=rgbCol("blue",15), lwd=c(1,1))
	rect(xleft=quantile(devs, p=c(0.025)), xright=quantile(devs, p=c(0.975)), ybottom=i-0.5, ytop=i+0.5, col=rgbCol("blue",30), lwd=c(1,1))
	rect(xleft=quantile(devs, p=c(0.25)), xright=quantile(devs, p=c(0.75)), ybottom=i-0.5, ytop=i+0.5, col=rgbCol("blue",60), lwd=c(1,1))
	points(devs, jitter(rep(i, length(devs)),amount=0.2))
	segments(x0=mean(devs), x1=mean(devs), y0=i-0.5, y1=i+0.5, col="red", lwd=3)
}

# Figure 2 and 4 (depending on scenarios)
windows(7,7)
par(mfrow=c(1,1))
source(paste0(path, "RandSlopeR2_PlottingFunctions.R"))
plotFig2(resEst, allsim, main=paste0("Scenario ", scenario))
