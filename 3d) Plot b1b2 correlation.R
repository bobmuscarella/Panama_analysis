# Plot species-level correlations between beta1 and beta2 parameters...

setwd("/Users/Bob/Projects/Postdoc/Panama/RESULTS/")

d <- read.csv('Species_Correlation_beta1_beta2.csv')
d1 <- d[d$metric=='growth',]
d2 <- d[d$metric=='survival',]

col <- brewer.pal(11, 'BrBG')[c(3,9,11)]

##### TRADE-OFFS (TABLE 2) #####
pdf(file="/Users/Bob/Projects/Postdoc/Panama/Figures/Tradeoffs_table2.pdf")

par(mfrow=c(2,1), mar=c(0,4,2,0), oma=c(3,5,4,10))
xs <- c(1,1.25, 2, 2.25, 3, 3.25)
colreps <- rep(rep(col, each=2), times=2)
bg <- ifelse(sign(d1$q2.5)==sign(d1$q97.5), colreps, 'white')
ylim <- c(min(d1$q2.5, na.rm=T) - (0.5*sd(d1$q50)), max(d1$q97.5, na.rm=T) + (0.5*sd(d1$q50)))
plot(xs, d1$q50, ylim=ylim, xlim=c(0.75,3.5), col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d1$q2.5, xs, d1$q97.5, len=0.075, angle=90, code=3, col=colreps, lwd=1.5)
points(xs, d1$q50, pch=21, bg=bg, col=colreps, cex=1.5)
axis(2); box()
mtext('Pearsons R', 2, 2, font=2)
mtext('(A) Growth: average vs. sensitivity to NCI', 3, 0.25, font=2, at=0.6, adj=0)
xs2 <- seq(1.125, 3.125, by=1)
mtext(rep(c('Cocoli','BCI','Sherman')), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
axis(1, labels=rep(c('',''), 3), at=xs, cex.axis=1.5, las=2, tck = 0.03)
bg <- ifelse(sign(d2$q2.5)==sign(d2$q97.5), colreps, 'white')
xs <- rep(1:(nrow(d2)/2), each=2) + rep(c(0, 0.25), times=(nrow(d2)/2))
ylim <- c(min(d2$q2.5, na.rm=T) - (0.5*sd(d2$q50)), max(d2$q97.5, na.rm=T) + (0.5*sd(d2$q50)))
plot(xs, d2$q50, ylim=ylim, xlim=c(0.75,3.5), col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d2$q2.5, xs, d2$q97.5, len=0.075, angle=90, code=3, col=colreps, lwd=1.5)
points(xs, d2$q50, pch=21, bg=bg, col=colreps, cex=1.5)
axis(2); box()
mtext('Pearsons R', 2, 2, font=2)
mtext('(B) Survival: average vs. sensitivity to NCI', 3, 0.25, font=2, at=0.6, adj=0)
axis(1, labels=rep(c('',''), 3), at=xs, cex.axis=1.5, las=2, tck = 0.03)
mtext(rep(c('Small','Large'), 6), side=1, line=0.3, at=xs, cex=0.8, las=2)
mtext(c('Cocoli','BCI','Sherman'), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)

dev.off()