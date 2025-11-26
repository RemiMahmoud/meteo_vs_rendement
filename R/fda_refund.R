library(refund)

data(DTI)
DTI1 <- DTI[DTI$visit==1 & complete.cases(DTI),]
par(mfrow=c(1,2))

# Fit model with linear functional term for CCA
fit.lf <- pfr(pasat ~ lf(cca, k=30, bs="ps"), data=DTI1)
plot(fit.lf, ylab=expression(paste(beta(t))), xlab="t")
## Not run: 
# Alternative way to plot
bhat.lf <- coef(fit.lf, n=101)
bhat.lf$upper <- bhat.lf$value + 1.96*bhat.lf$se
bhat.lf$lower <- bhat.lf$value - 1.96*bhat.lf$se
matplot(bhat.lf$cca.argvals, bhat.lf[,c("value", "upper", "lower")],
        type="l", lty=c(1,2,2), col=1,
        ylab=expression(paste(beta(t))), xlab="t")

# Fit model with additive functional term for CCA, using tensor product basis
fit.af <- pfr(pasat ~ af(cca, Qtransform=TRUE, k=c(7,7)), data=DTI1)
plot(fit.af, scheme=2, xlab="t", ylab="cca(t)", main="Tensor Product")
plot(fit.af, scheme=2, Qtransform=TRUE,
     xlab="t", ylab="cca(t)", main="Tensor Product")

# Change basistype to thin-plate regression splines
fit.af.s <- pfr(pasat ~ af(cca, basistype="s", Qtransform=TRUE, k=50),
                data=DTI1)
plot(fit.af.s, scheme=2, xlab="t", ylab="cca(t)", main="TPRS", rug=FALSE)
plot(fit.af.s, scheme=2, Qtransform=TRUE,
     xlab="t", ylab="cca(t)", main="TPRS", rug=FALSE)

# Visualize bivariate function at various values of x
par(mfrow=c(3,2))
vis.pfr(fit.af, xval=.2)
vis.pfr(fit.af, xval=.4)
vis.pfr(fit.af, xval=.45)
vis.pfr(fit.af, xval=.5)
vis.pfr(fit.af, xval=.6)
vis.pfr(fit.af, xval=.8)

# Include random intercept for subject
DTI.re <- DTI[complete.cases(DTI$cca),]
DTI.re$ID <- factor(DTI.re$ID)
fit.re <- pfr(pasat ~ lf(cca, k=30) + re(ID), data=DTI.re)
coef.re <- coef(fit.re)
par(mfrow=c(1,2))
plot(fit.re, select = 1)
summary(fit.re)
plot(coef.re, select = "cca")

# FPCR_R Model
fit.fpc <- pfr(pasat ~ fpc(cca), data=DTI.re)
plot(fit.fpc)

# PEER Model with second order difference penalty
DTI.use <- DTI[DTI$case==1,]
DTI.use <- DTI.use[complete.cases(DTI.use$cca),]
fit.peer <- pfr(pasat ~ peer(cca, argvals=seq(0,1,length=93),
                             integration="riemann", pentype="D"), data=DTI.use)
plot(fit.peer)

## End(Not run)