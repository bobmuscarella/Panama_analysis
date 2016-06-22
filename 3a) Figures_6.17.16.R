setwd("/Users/Bob/Projects/Postdoc/Panama/RESULTS/_6.14.16")
library(jagsUI)
library(rjags)
library(RColorBrewer)

col <- brewer.pal(11, 'BrBG')[c(3,9,11)]

res <- data.frame()
for(j in 1:2){
m <- c('growth','survival')[j]
for(i in 1:6){
file <- list.files(m)[i]
r <- readRDS(paste(m,file,sep='/'))
samp <- do.call(rbind, r$samples)
samp <- samp[,!grepl('sigma', colnames(samp))]
samp <- samp[,!grepl('deviance', colnames(samp))]
q <- round(apply(samp, 2, quantile, c(0.025, 0.05, 0.5, 0.95, 0.975)), 3)
s <- substring(file, 5, 6)
p <- substring(file, 1, 3)
param <- colnames(q)
rownames(q) <- paste('q',c(2.5, 5, 50, 95, 97.5),sep='')
res <- rbind(res, data.frame(m, s, p, param, t(q)))
}}

rownames(res) <- NULL
res$p <- ifelse(res$p=='coc',1,ifelse(res$p=='bci',2,3))

# SD(GROWTH) OF ORIGINAL DATA SCALED WITHIN SIZE CLASS...
load("/Users/Bob/Projects/Postdoc/Panama/DATA/Panama_AnalysisData_6.14.16.RDA")
tdata <- droplevels(tdata[!is.na(tdata$growth) & !is.na(tdata$All.NCI) & tdata$Not.Edge==1,])
x <- split(tdata, paste(tdata$plot, ifelse(tdata$dbh>=100,'large','small')))
xsd5 <- lapply(x, function(x) sd(x$growth)*5)
for(i in 1:length(x)){
  x[[i]] <- x[[i]][abs(x[[i]]$growth) < xsd5[i] & x[[i]]$Growth.Include.3,]
}
f <- unlist(lapply(x, function(x) sd(x$growth)))


########################################################
########## TRAIT EFFECT ON SENSITIVITY TO NCI ##########
########################################################

pdf(file="/Users/Bob/Projects/Postdoc/Panama/Figures/beta.t2.pdf")
##### GET GROWTH RESULTS
g <- res[res$m=='growth',]
b <- c('beta.wd[2]', 'beta.lma[2]')
d <- g[as.character(g$param) %in% b,]
d <- d[order(rev(d$param), d$p, rev(d$s)),]

par(mfrow=c(2,1), mar=c(0,4,0.25,0), oma=c(5,3,4,1))
colreps <- rep(rep(col, each=2), times=2)
bg <- ifelse(sign(d$q2.5)==sign(d$q97.5), colreps, 'white')
xs <- rep(1:(nrow(d)/2), each=2) + rep(c(0, 0.25), times=(nrow(d)/2))
xlim <- c(min(xs)-0.1875, max(xs)+0.1875)
ylim <- c(min(d$q2.5, na.rm=T) - (3*sd(d$q50)), max(d$q97.5, na.rm=T) + (3*sd(d$q50)))
plot(xs, d$q50, ylim=ylim, xlim=xlim, col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d$q2.5, xs, d$q97.5, len=0.1, angle=90, code=3, col=colreps, lwd=1.5)
segments(xs, d$q5, xs, d$q95, col=colreps, lwd=3)
points(xs, d$q50, pch=21, bg=bg, col=colreps, cex=1.5)
abline(v=mean(xs))
axis(2); box()
mtext('Standard Effect Size', 2, 2, font=2)
mtext('(A) Growth: WSG', 3, -1.25, font=2, at=0.6875, adj=0)
mtext('(B) Growth: LMA', 3, -1.25, font=2, at=3.75, adj=0)
#legend(2.5, -0.08, legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
axis(1, labels=rep(c('',''), 6), at=xs, cex.axis=1.5, las=2, tck = 0.02)
xs2 <- seq(1.125, 6.125, by=1)
mtext(rep(c('Cocoli','BCI','Sherman'), 2), side=1, line=-1.25, at=xs2, cex.axis=1.5, col=col, font=2)
abline(v=c(1.625,2.625,4.625,5.625), lty=3)


##### GET SURVIVAL RESULTS
g <- res[res$m=='survival',]
b <- c('beta.wd[2]', 'beta.lma[2]')
d <- g[as.character(g$param) %in% b,]
d <- d[order(rev(d$param), d$p, rev(d$s)),]

colreps <- rep(rep(col, each=2), times=2)
bg <- ifelse(sign(d$q2.5)==sign(d$q97.5), colreps, 'white')
xs <- rep(1:(nrow(d)/2), each=2) + rep(c(0, 0.25), times=(nrow(d)/2))
ylim <- c(min(d$q2.5, na.rm=T) - (3*sd(d$q50)), max(d$q97.5, na.rm=T) + (3*sd(d$q50)))
xlim <- c(min(xs)-0.1875, max(xs)+0.1875)
plot(xs, d$q50, ylim=ylim, xlim=xlim, col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d$q2.5, xs, d$q97.5, len=0.1, angle=90, code=3, col=colreps, lwd=1.5)
segments(xs, d$q5, xs, d$q95, col=colreps, lwd=3)
points(xs, d$q50, pch=21, bg=bg, col=colreps, cex=1.5)
abline(v=mean(xs))
axis(2); box()
mtext('Standard Effect Size', 2, 2, font=2)
mtext('(C) Survival: WSG', 3, -1.25, font=2, at=0.6875, adj=0)
mtext('(D) Survival: LMA', 3, -1.25, font=2, at=3.75, adj=0)
#legend(5.25, -0.1, legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
#legend(2.5, -0.1, legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
axis(1, labels=rep(c('',''), 6), at=xs, cex.axis=1.5, las=2, tck = 0.02)
#mtext(rep(c('Sm','Lg'), 6), side=1, line=0.3, at=xs, cex.axis=1.5)
mtext(rep(c('Small','Large'), 6), side=1, line=0.3, at=xs, cex=0.8, las=2)
abline(v=c(1.625,2.625,4.625,5.625), lty=3)
xs2 <- seq(1.125, 6.125, by=1)
mtext(rep(c('Cocoli','BCI','Sherman'), 2), side=1, line=-1.25, at=xs2, cex.axis=1.5, col=col, font=2)


dev.off()



####################################################
########## TRAIT EFFECT ON AVERAGE RATES  ##########
####################################################

##### GET GROWTH RESULTS
g <- res[res$m=='growth',]
b <- c('beta.wd[1]', 'beta.lma[1]')
d <- g[as.character(g$param) %in% b,]
d <- d[order(rev(d$param), d$p, rev(d$s)),]

par(mfrow=c(2,1), mar=c(0,4,0.25,0), oma=c(5,3,4,1))
colreps <- rep(rep(col, each=2), times=2)
bg <- ifelse(sign(d$q2.5)==sign(d$q97.5), colreps, 'white')
xs <- rep(1:(nrow(d)/2), each=2) + rep(c(0, 0.25), times=(nrow(d)/2))
ylim <- c(min(d$q2.5, na.rm=T) - (2*sd(d$q50)), max(d$q97.5, na.rm=T) + (2*sd(d$q50)))
plot(xs, d$q50, ylim=ylim, xlim=range(xs), col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d$q2.5, xs, d$q97.5, len=0.1, angle=90, code=3, col=colreps, lwd=2)
segments(xs, d$q5, xs, d$q95, col=colreps, lwd=5)
points(xs, d$q50, pch=21, bg=bg, col=colreps, cex=1.5)
abline(v=mean(xs))
axis(2); box()
mtext('Standard Effect Size', 2, 3, font=2)
mtext('(A) Growth, WSG', 3, -1.25, font=2, at=1.75)
mtext('(B) Growth, LMA', 3, -1.25, font=2, at=4.5)
legend(0.8, -0.225, legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
legend(3.8, -0.225, legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
axis(1, labels=rep(c('',''), 6), at=xs, cex.axis=1.5, las=2, tck = 0.03)

##### GET SURVIVAL RESULTS
g <- res[res$m=='survival',]
b <- c('beta.wd[1]', 'beta.lma[1]')
d <- g[as.character(g$param) %in% b,]
d <- d[order(rev(d$param), d$p, rev(d$s)),]

colreps <- rep(rep(col, each=2), times=2)
bg <- ifelse(sign(d$q2.5)==sign(d$q97.5), colreps, 'white')
xs <- rep(1:(nrow(d)/2), each=2) + rep(c(0, 0.25), times=(nrow(d)/2))
ylim <- c(min(d$q2.5, na.rm=T) - (2*sd(d$q50)), max(d$q97.5, na.rm=T) + (2*sd(d$q50)))
plot(xs, d$q50, ylim=ylim, xlim=range(xs), col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d$q2.5, xs, d$q97.5, len=0.1, angle=90, code=3, col=colreps, lwd=2)
segments(xs, d$q5, xs, d$q95, col=colreps, lwd=5)
points(xs, d$q50, pch=21, bg=bg, col=colreps, cex=1.5)
abline(v=mean(xs))
axis(2); box()
mtext('Standard Effect Size', 2, 3, font=2)
mtext('(A) Survival, LMA', 3, -1.25, font=2, at=1.75)
mtext('(B) Survival, WSG', 3, -1.25, font=2, at=4.5)
legend(2.5, -0.4, legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
legend(5.25, -0.4, legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
axis(1, labels=rep(c('',''), 6), at=xs, cex.axis=1.5, las=2, tck = 0.03)
mtext(rep(c('S','L'), 6), side=1, line=0.3, at=xs, cex.axis=1.5)




#########################################
########## AVERAGE NCI EFFECT  ##########
#########################################

##### GET GROWTH RESULTS
g <- res[res$m=='growth',]
b <- c('mu.beta[2]')
d <- g[as.character(g$param) %in% b,]
d <- d[order(rev(d$param), d$p, rev(d$s)),]

par(mfrow=c(2,1), mar=c(0,4,0.25,0), oma=c(3,5,4,10))
xs <- c(1,1.25, 2, 2.25, 3, 3.25)
colreps <- rep(rep(col, each=2), times=2)
bg <- ifelse(sign(d$q2.5)==sign(d$q97.5), colreps, 'white')
ylim <- c(min(d$q2.5, na.rm=T) - (2*sd(d$q50)), max(d$q97.5, na.rm=T) + (2*sd(d$q50)))
plot(xs, d$q50, ylim=ylim, xlim=c(0.5,3.75), col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d$q2.5, xs, d$q97.5, len=0.1, angle=90, code=3, col=colreps, lwd=2)
segments(xs, d$q5, xs, d$q95, col=colreps, lwd=5)
points(xs, d$q50, pch=21, bg=bg, col=colreps, cex=1.5)
axis(2); box()
mtext('Standard Effect Size', 2, 3, font=2)
mtext('(A) Growth', 3, -1.5, font=2, at=1)
legend('topright', legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
axis(1, labels=rep(c('',''), 3), at=xs, cex.axis=1.5, las=2, tck = 0.03)

##### GET SURVIVAL RESULTS
g <- res[res$m=='survival',]
b <- c('mu.beta[2]')
d <- g[as.character(g$param) %in% b,]
d <- d[order(rev(d$param), d$p, rev(d$s)),]

xs <- c(1,1.25, 2, 2.25, 3, 3.25)
colreps <- rep(rep(col, each=2), times=2)
bg <- ifelse(sign(d$q2.5)==sign(d$q97.5), colreps, 'white')
ylim <- c(min(d$q2.5, na.rm=T) - (2*sd(d$q50)), max(d$q97.5, na.rm=T) + (2*sd(d$q50)))
plot(xs, d$q50, ylim=ylim, xlim=c(0.5,3.75), col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d$q2.5, xs, d$q97.5, len=0.1, angle=90, code=3, col=colreps, lwd=2)
segments(xs, d$q5, xs, d$q95, col=colreps, lwd=5)
points(xs, d$q50, pch=21, bg=bg, col=colreps, cex=1.5)
axis(2); box()
mtext('Standard Effect Size', 2, 3, font=2)
mtext('(B) Survival', 3, -1.5, font=2, at=1)
legend('bottomright', legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
axis(1, labels=rep(c('',''), 3), at=xs, cex.axis=1.5, las=2, tck = 0.03)
mtext(rep(c('S','L'), 3), side=1, line=0.3, at=xs, cex.axis=1.5)


#########################################
########## AVERAGE SIZE EFFECT  ##########
#########################################

##### GET GROWTH RESULTS
g <- res[res$m=='growth',]
b <- c('mu.beta[3]')
d <- g[as.character(g$param) %in% b,]
d <- d[order(rev(d$param), d$p, rev(d$s)),]

par(mfrow=c(2,1), mar=c(0,4,0.25,0), oma=c(3,5,4,10))
xs <- c(1,1.25, 2, 2.25, 3, 3.25)
colreps <- rep(rep(col, each=2), times=2)
bg <- ifelse(sign(d$q2.5)==sign(d$q97.5), colreps, 'white')
ylim <- c(min(d$q2.5, na.rm=T) - (2*sd(d$q50)), max(d$q97.5, na.rm=T) + (2*sd(d$q50)))
plot(xs, d$q50, ylim=ylim, xlim=c(0.5,3.75), col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d$q2.5, xs, d$q97.5, len=0.1, angle=90, code=3, col=colreps, lwd=2)
segments(xs, d$q5, xs, d$q95, col=colreps, lwd=5)
points(xs, d$q50, pch=21, bg=bg, col=colreps, cex=1.5)
axis(2); box()
mtext('Standard Effect Size', 2, 3, font=2)
mtext('(A) Growth', 3, -1.5, font=2, at=1)
legend('topright', legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
axis(1, labels=rep(c('',''), 3), at=xs, cex.axis=1.5, las=2, tck = 0.03)

##### GET SURVIVAL RESULTS
g <- res[res$m=='survival',]
b <- c('mu.beta[3]')
d <- g[as.character(g$param) %in% b,]
d <- d[order(rev(d$param), d$p, rev(d$s)),]

xs <- c(1,1.25, 2, 2.25, 3, 3.25)
colreps <- rep(rep(col, each=2), times=2)
bg <- ifelse(sign(d$q2.5)==sign(d$q97.5), colreps, 'white')
ylim <- c(min(d$q2.5, na.rm=T) - (2*sd(d$q50)), max(d$q97.5, na.rm=T) + (2*sd(d$q50)))
plot(xs, d$q50, ylim=ylim, xlim=c(0.5,3.75), col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d$q2.5, xs, d$q97.5, len=0.1, angle=90, code=3, col=colreps, lwd=2)
segments(xs, d$q5, xs, d$q95, col=colreps, lwd=5)
points(xs, d$q50, pch=21, bg=bg, col=colreps, cex=1.5)
axis(2); box()
mtext('Standard Effect Size', 2, 3, font=2)
mtext('(B) Survival', 3, -1.5, font=2, at=1)
legend('bottomright', legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
axis(1, labels=rep(c('',''), 3), at=xs, cex.axis=1.5, las=2, tck = 0.03)
mtext(rep(c('S','L'), 3), side=1, line=0.3, at=xs, cex.axis=1.5)



#####################################
########## AVERAGE RATES   ##########
#####################################

##### GET GROWTH RESULTS
g <- res[res$m=='growth',]
b <- c('mu.beta[1]')
d <- g[as.character(g$param) %in% b,]
d <- d[order(rev(d$param), d$p, rev(d$s)),]

# RESCALE TO MM/YR BASED ON SD(GROWTH) OF ORIGINAL DATA SCALED WITHIN SIZE CLASS...
d[,5:9] <- d[,5:9] * f[c(2,1,4,3,6,5)]

par(mfrow=c(2,1), mar=c(0,4,0.25,0), oma=c(3,5,4,10))
xs <- c(1,1.25, 2, 2.25, 3, 3.25)
colreps <- rep(rep(col, each=2), times=2)
bg <- ifelse(sign(d$q2.5)==sign(d$q97.5), colreps, 'white')
ylim <- c(min(d$q2.5, na.rm=T) - (2*sd(d$q50)), max(d$q97.5, na.rm=T) + (2*sd(d$q50)))
plot(xs, d$q50, ylim=ylim, xlim=c(0.5,3.75), col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d$q2.5, xs, d$q97.5, len=0.1, angle=90, code=3, col=colreps, lwd=2)
segments(xs, d$q5, xs, d$q95, col=colreps, lwd=5)
points(xs, d$q50, pch=21, bg=bg, col=colreps, cex=1.5)
axis(2); box()
mtext('Predicted growth 
(mm / yr)', 2, 3, font=2)
mtext('(A) Growth', 3, -1.5, font=2, at=1)
legend('topright', legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
axis(1, labels=rep(c('',''), 3), at=xs, cex.axis=1.5, las=2, tck = 0.03)

##### GET SURVIVAL RESULTS
g <- res[res$m=='survival',]
b <- c('mu.beta[1]')
d <- g[as.character(g$param) %in% b,]
d <- d[order(rev(d$param), d$p, rev(d$s)),]

d[,5:9] <- d[,5:9]/(1+d[,5:9])

xs <- c(1,1.25, 2, 2.25, 3, 3.25)
colreps <- rep(rep(col, each=2), times=2)
bg <- ifelse(sign(d$q2.5)==sign(d$q97.5), colreps, 'white')
ylim <- c(min(d$q2.5, na.rm=T) - (2*sd(d$q50)), max(d$q97.5, na.rm=T) + (2*sd(d$q50)))
plot(xs, d$q50, ylim=ylim, xlim=c(0.5,3.75), col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d$q2.5, xs, d$q97.5, len=0.1, angle=90, code=3, col=colreps, lwd=2)
segments(xs, d$q5, xs, d$q95, col=colreps, lwd=5)
points(xs, d$q50, pch=21, bg=bg, col=colreps, cex=1.5)
axis(2); box()
mtext('Predicted annual survival
(%)', 2, 3, font=2, cex=1)
mtext('(B) Survival', 3, -1.5, font=2, at=1)
legend('topright', legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
axis(1, labels=rep(c('',''), 3), at=xs, cex.axis=1.5, las=2, tck = 0.03)
mtext(rep(c('S','L'), 3), side=1, line=0.3, at=xs, cex.axis=1.5)

