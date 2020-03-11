#####################################################################################
# plotCheck
# Check for whether condition R is a good fit to between-individual variance estimated for bins of x
# points are estimates from raw data with size indicating relative sample size
#####################################################################################
# Arguments: Original data frame, settings for simulation settings, ests for model estimates
#            binfact: a factor for efficient binning
#            ptscale: a factor for point size scaling
# Details:   Black curve shows simulated conditional between-individual variance
#            Blue curve shows estimated conditional between-individual variance
#            Black line shows simulated Vi
#            Blue line shows estiamted Vi
#            Orange line shows estiamted Vi using Johnson's method
# Value:     No specific return object
plotCheck = function(md, settings, ests, binfact=3, ptscale=5, ...) {
	nobs     = nrow(md)
	Xmat     = cbind(rep(1,nobs), md$x)
	vars     = tapply(md$yPred, round(md$x*binfact,0)/binfact, var)
	weight   = tapply(md$yPred, round(md$x*binfact,0)/binfact, length) 
	xmin     = min(as.numeric(names(vars)))
	xmax     = max(as.numeric(names(vars)))
	x        = seq(xmin, xmax, length.out=1001)
	condRsim = conditionalR(x, vect=settings)
	condRest = conditionalR(x, vect=ests)
	simVi    = varDecomp(vect=settings)["Vi"]
	estVi    = varDecomp(vect=ests)["Vi"]
	meanVi   = varDecomp(vect=ests, Xmat=Xmat)["meanVi"]
	xvals    = as.numeric(names(vars))
	plot(xvals, vars, type="n", ylab="", xlab="", ylim=c(0,max(vars, na.rm=TRUE)), xlim=c(xmin, xmax), yaxt="n", ...) #  
	title(ylab="Between-individual variance", xlab="Covariate", cex.lab=1.2)
	points(xvals, vars, cex=weight/max(weight)*ptscale)  # between-individual variances estimated from raw data
	points(x, condRsim, type="l")                        # between-individual variances from simulation settings
	points(x, condRest, type="l", col="blue", lwd=2)   	 # between-individual variances estimated from model fit
	abline(h = simVi)                                    # Vi from simulation settings
	abline(h = estVi, col="blue")                        # Vi estimated from model fit
	abline(h = meanVi, col="orange")                     # Vi estimated from using Johnson's method
}


#####################################################################################
# plotFig2a
# A function for plotting
#####################################################################################
# Arguments: ....
# Value:     No specific return object
plotFig2a = function(resEst, allsim, finalOrd, xlim=c(-0.4, 0.4)) {
	par(mar=c(5,5,1,1))
	plot((resEst[1,finalOrd]-allsim[finalOrd])/abs(allsim[finalOrd]), jitter(1:length(finalOrd)), ylab="", yaxt="n", ylim=c(1,length(finalOrd)), xlab="", xlim=xlim)
	title(xlab="Proportional difference to simulated values\n(Est-Sim)/Sim", line=3.5)
	axis(2, at=1:length(finalOrd), labels=finalOrd, las=1)
	abline(v=0)
	abline(h=1:I(length(finalOrd)+1)-0.5)
	abline(h=c(1,4,7,8,12,16,19)-0.5, lwd=2)
	for(i in 2:nsim) points((resEst[i,finalOrd]-allsim[finalOrd])/abs(allsim[finalOrd]), jitter(1:length(finalOrd)))
}

#####################################################################################
# plotFig2b
# A function for plotting 
#####################################################################################
# Arguments: ....
# Value:     No specific return object
plotFig2b = function(resEst, allsim, finalOrd, xlim=c(-2, 3.8)) {
	par(mar=c(5,5,1,1))
	plot(allsim[finalOrd], jitter(1:length(finalOrd)), ylab="", yaxt="n", ylim=c(1,length(finalOrd)), xlab="", xlim=xlim)
	title(xlab="Absolute estimates", line=3.5)
	axis(2, at=1:length(finalOrd), labels=finalOrd, las=1)
	abline(v=0)
	abline(h=1:I(length(finalOrd)+1)-0.5)
	abline(h=c(1,4,7,8,12,16,19)-0.5, lwd=2)
	for(i in 1:nrow(resEst)) points(resEst[i,finalOrd], jitter(1:length(finalOrd)))
	segments(x0=allsim[finalOrd], x1=allsim[finalOrd], y0=1:length(finalOrd)-0.5, y1=1:length(finalOrd)+0.5, lwd=2, col="red")
}

#####################################################################################
# plotFig2c
# A function for plotting 
#####################################################################################
# Arguments: ....
# Value:     No specific return object
plotFig2c = function(resEst, allsim, main="") {
	if(main=="") par(mar=c(5,5,1,1))
	if(main!="") par(mar=c(4.5,5,2,1))
	xmin = allsim["meanX2"] -2*sqrt(allsim["Vx"])
	xmax = allsim["meanX2"] +2*sqrt(allsim["Vx"])
	x = seq(xmin,xmax,length.out=1001)
	condR = allsim["Vu"] + 2*x*allsim["Cuv"] + x^2*allsim["Vv"]
	plot(x, condR, type="l", ylab="Between-individual variance", xlab="Covariate", main=main, ylim=c(0,6), cex.lab=1.4, yaxt="n", xaxt="n", xaxs="i") #  
	abline(h=0)
	# Minimum value
	col = rgbCol("darkolivegreen2")
	for(i in 1:nrow(resEst)) {
		segments(x0=-3, x1=resEst[i,"xmin"], y0=resEst[i,"minVix"], y1=resEst[i,"minVix"], col=col)
		segments(x0=resEst[i,"xmin"], x1=resEst[i,"xmin"], y0=-3, y1=resEst[i,"minVix"], col=col)
	}
	# random intercepts
	col = rgbCol("lightblue")
	for(i in 1:nrow(resEst)) {
		segments(x0=-3, x1=0, y0=resEst[i,"Vu"], y1=resEst[i,"Vu"], col=col)
		segments(x0=0, x1=0, y0=-3, y1=resEst[i,"Vu"], col=col)
	}
	# Individual component
	col = rgbCol("tan1")
	for(i in 1:nrow(resEst)) {
		abline(h=resEst[i,"Vi"], lty=1, col=col)
	}
	# Conditional repeatability
	col = rgbCol("grey50")
	for(i in 1:nrow(resEst)) {
		condRest = resEst[i, "Vu"] + 2*x*resEst[i,"Cuv"] + x^2*resEst[i,"Vv"]
		points(x, condRest, type="l", col=col) 
	}
	# Average values
	segments(x0=-3, x1=mean(resEst[,"xmin"]), y0=mean(resEst[,"minVix"]), y1=mean(resEst[,"minVix"]), col="green", lwd=4)
	segments(x0=mean(resEst[,"xmin"]), x1=mean(resEst[,"xmin"]), y0=-3, y1=mean(resEst[,"minVix"]), col="green", lwd=4)
	segments(x0=-3, x1=0, y0=mean(resEst[,"Vu"]), y1=mean(resEst[,"Vu"]), col="blue", lwd=4)
	segments(x0=0, x1=0, y0=-3, y1=mean(resEst[,"Vu"]), col="blue", lwd=4)
	abline(h=mean(resEst[,"Vu"]+resEst[,"Vv"]*resEst[,"Vx"]), lty=1, col="orangered", lwd=4)
	# Simulated values
	Vu = settings["Vu"]
	Vv = settings["Vv"]
	Cuv = settings["Cuv"]
	Vx = settings["Vx"]
	# Conditional repeatability
	points(x, condR, type="l", lwd=2) 
	# random intercepts
	segments(x0=-3, x1=0, y0=allsim["Vu"], y1=allsim["Vu"], lwd=2)
	segments(x0=0, x1=0, y0=-3, y1=allsim["Vu"], lwd=2)
	axis(2, at=Vu, labels=expression("V"[u]), las=1, cex.axis=1.2)
	axis(1, at=0, labels=0, las=1, cex.axis=1.2)
	# Minimum value
	segments(x0=-3, x1=allsim["xmin"], y0=allsim["minVix"], y1=allsim["minVix"], lwd=2)
	segments(x0=allsim["xmin"], x1=allsim["xmin"], y0=-3, y1=allsim["minVix"], lwd=2)
	axis(1, at=allsim["xmin"], labels=expression("x"[min]), las=1, cex.axis=1.2)
	axis(2, at=allsim["minVix"], labels=expression(paste("min(V"["I,x"],") ")), las=1, cex.axis=1.2)
	# Indivudal component
	abline(h=allsim["Vi"], lty=1, lwd=2)
	axis(2, at=allsim["Vi"], labels=expression("V"[I]), las=2, cex.axis=1.2)
	box(lwd=2)
}

#####################################################################################
# rgbCol
# A helper function that converts character value colors to rbg scores
#####################################################################################
# Arguments: Color given as a character vector and alpha for transparency
# Value:     An rbg valued color with transparency
rgbCol = function(color, alpha=100) {
	rgb(col2rgb(color)["red",1], col2rgb(color)["green",1], col2rgb(color)["blue",1], alpha=alpha, maxColorValue = 255)
}
