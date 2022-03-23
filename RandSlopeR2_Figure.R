####################################################################################################################
# Figure 1
####################################################################################################################

rm(list=ls())
require(lme4)
require(psych)

# Settings
# Variances
Vu = 1
Vv = 1
Cuv = 0.8
# Variance in covarate
Vx = 1
covBalancedYN = TRUE

# Derived quantities
Sigma = matrix(c(Vu, Cuv, Cuv, Vv), ncol=2)
matrixcalc::is.positive.definite(Sigma)

# Plotting
xmin = -2
xmax = 2
x = seq(xmin,xmax,length.out=1001)
condR = Vu+2*x*Cuv+x^2*Vv

windows(7, 14)
par(mfrow=c(1, 2), mar=c(5,5,1,1))

# plot a
plot(x, condR, type="l", lwd=2, las=1, ylab="Between-individual variance", xlab="Covariate\n(environment, context, time or age)", ylim=c(0,6), xaxt="n", yaxt="n", cex.lab=1.2, xaxs="i")
#abline(h=0)
segments(x0=-3, x1=0, y0=Vu, y1=Vu)
segments(x0=0, x1=0, y0=-3, y1=Vu)
axis(2, at=Vu, labels=expression("V"[u]), las=1)
axis(1, at=0, labels=0, las=1)
segments(x0=-3, x1=-Cuv/Vv, y0=Vu-2*Cuv^2/Vv+Cuv^2, y1=Vu-2*Cuv^2/Vv+Cuv^2)
segments(x0=-Cuv/Vv, x1=-Cuv/Vv, y0=-3, y1=Vu-2*Cuv^2/Vv+Cuv^2)
axis(1, at=-Cuv/Vv, labels=expression("x"[min]), las=1)
axis(2, at=Vu-Cuv^2/Vv, labels=expression(paste("min(V"["I,x"],") ")), las=1)
abline(h=Vu+Vv*Vx, lty=2)
axis(2, at=Vu+Vv*Vx, labels=expression("V"[I]), las=2)
#text(min(x)+0.2, 6*0.98, "(a)", cex=1.4)
mtext("(a)", 2, line=2, at=6*1.05, cex=1.6, las=1)
box(lwd=2)


# plot b
mu = 10
beta = 0.5
nind=500
randid = MASS::mvrnorm(nind, c(mu, beta), Sigma)
ymin = randid[,1] + min(x)*randid[,2]
ymax = randid[,1] + max(x)*randid[,2]
plot(x, seq(min(ymin), max(ymax), length.out=length(x)), type="n", xlab="Covariate", ylab="", xaxt="n", yaxt="n", ylim=c(5,19), cex.lab=1.2, xaxs="i")
title(ylab="Trait value", line=2, cex.lab=1.2)
axis(1, at=0, labels="0")
for(i in 1:nind)
	abline(a=randid[i,1], b=randid[i,2])
#text(min(x)+0.2, 5+(19-5)*0.98, "(b)", cex=1.4)
mtext("(b)", 2, line=2, at=5+(19-5)*1.05, cex=1.6, las=1)
box(lwd=2)