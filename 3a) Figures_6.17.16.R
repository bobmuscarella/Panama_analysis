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
q50 <- round(unlist(r$q50), 3)
q50 <- q50[!grepl('sigma',names(q50))]
q50 <- q50[!grepl('deviance',names(q50))]

q2.5 <- round(unlist(r$q2.5), 3)
q2.5 <- q2.5[!grepl('sigma',names(q2.5))]
q2.5 <- q2.5[!grepl('deviance',names(q2.5))]

q97.5 <- round(unlist(r$q97.5), 3)
q97.5 <- q97.5[!grepl('sigma',names(q97.5))]
q97.5 <- q97.5[!grepl('deviance',names(q97.5))]

s <- substring(file, 5, 6)
p <- substring(file, 1, 3)

param <- names(q50)
res <- rbind(res, data.frame(m, s, p, param, q50, q2.5, q97.5))
}

rownames(res) <- NULL
}


res$p <- ifelse(res$p=='coc',1,ifelse(res$p=='bci',2,3))


########################################################
########## TRAIT EFFECT ON SENSITIVITY TO NCI ##########
########################################################

##### GET GROWTH RESULTS
g <- res[res$m=='growth',]
b <- c('beta.wd2', 'beta.lma2')
d <- g[as.character(g$param) %in% b,]
d <- d[order(d$param, d$p, d$s),]

par(mfrow=c(2,1), mar=c(0,4,0.25,0), oma=c(5,3,4,1))
colreps <- rep(rep(col, each=2), times=2)
bg <- ifelse(sign(d$q2.5)==sign(d$q97.5), colreps, 'white')
xs <- rep(1:(nrow(d)/2), each=2) + rep(c(0, 0.25), times=(nrow(d)/2))
ylim <- c(min(d$q2.5, na.rm=T) - (2*sd(d$q50)), max(d$q97.5, na.rm=T) + (2*sd(d$q50)))
plot(xs, d$q50, ylim=ylim, xlim=range(xs), col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d$q2.5, xs, d$q97.5, len=0.1, angle=90, code=3, col=colreps, lwd=2)
points(xs, d$q50, pch=21, bg=bg, col=colreps, cex=1.5)
abline(v=mean(xs))
axis(2); box()
mtext('Standard Effect Size', 2, 3, font=2)
mtext('(A) Growth, LMA', 3, -1.25, font=2, at=1.75)
mtext('(B) Growth, WSG', 3, -1.25, font=2, at=4.5)
legend(2.5, -0.08, legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
legend(5.25, -0.08, legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
axis(1, labels=rep(c('',''), 6), at=xs, cex.axis=1.5, las=2, tck = 0.03)

##### GET SURVIVAL RESULTS
g <- res[res$m=='survival',]
b <- c('beta.wd2', 'beta.lma2')
d <- g[as.character(g$param) %in% b,]
d <- d[order(d$param, d$p, d$s),]

#par(mfrow=c(2,1), mar=c(0,4,0.25,0), oma=c(5,3,4,1))
colreps <- rep(rep(col, each=2), times=2)
bg <- ifelse(sign(d$q2.5)==sign(d$q97.5), colreps, 'white')
xs <- rep(1:(nrow(d)/2), each=2) + rep(c(0, 0.25), times=(nrow(d)/2))
ylim <- c(min(d$q2.5, na.rm=T) - (2*sd(d$q50)), max(d$q97.5, na.rm=T) + (2*sd(d$q50)))
plot(xs, d$q50, ylim=ylim, xlim=range(xs), col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d$q2.5, xs, d$q97.5, len=0.1, angle=90, code=3, col=colreps, lwd=2)
points(xs, d$q50, pch=21, bg=bg, col=colreps, cex=1.5)
abline(v=mean(xs))
axis(2); box()
mtext('Standard Effect Size', 2, 3, font=2)
mtext('(A) Survival, LMA', 3, -1.25, font=2, at=1.75)
mtext('(B) Survival, WSG', 3, -1.25, font=2, at=4.5)
legend(5.25, -0.1, legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
legend(2.5, -0.1, legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
axis(1, labels=rep(c('',''), 6), at=xs, cex.axis=1.5, las=2, tck = 0.03)




####################################################
########## TRAIT EFFECT ON AVERAGE RATES  ##########
####################################################

##### GET GROWTH RESULTS
g <- res[res$m=='growth',]
b <- c('beta.wd1', 'beta.lma1')
d <- g[as.character(g$param) %in% b,]
d <- d[order(d$param, d$p, d$s),]

par(mfrow=c(2,1), mar=c(0,4,0.25,0), oma=c(5,3,4,1))
colreps <- rep(rep(col, each=2), times=2)
bg <- ifelse(sign(d$q2.5)==sign(d$q97.5), colreps, 'white')
xs <- rep(1:(nrow(d)/2), each=2) + rep(c(0, 0.25), times=(nrow(d)/2))
ylim <- c(min(d$q2.5, na.rm=T) - (2*sd(d$q50)), max(d$q97.5, na.rm=T) + (2*sd(d$q50)))
plot(xs, d$q50, ylim=ylim, xlim=range(xs), col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d$q2.5, xs, d$q97.5, len=0.1, angle=90, code=3, col=colreps, lwd=2)
points(xs, d$q50, pch=21, bg=bg, col=colreps, cex=1.5)
abline(v=mean(xs))
axis(2); box()
mtext('Standard Effect Size', 2, 3, font=2)
mtext('(A) Growth, LMA', 3, -1.25, font=2, at=1.75)
mtext('(B) Growth, WSG', 3, -1.25, font=2, at=4.5)
legend(0.8, -0.225, legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
legend(3.8, -0.225, legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
axis(1, labels=rep(c('',''), 6), at=xs, cex.axis=1.5, las=2, tck = 0.03)

##### GET SURVIVAL RESULTS
g <- res[res$m=='survival',]
b <- c('beta.wd1', 'beta.lma1')
d <- g[as.character(g$param) %in% b,]
d <- d[order(d$param, d$p, d$s),]

colreps <- rep(rep(col, each=2), times=2)
bg <- ifelse(sign(d$q2.5)==sign(d$q97.5), colreps, 'white')
xs <- rep(1:(nrow(d)/2), each=2) + rep(c(0, 0.25), times=(nrow(d)/2))
ylim <- c(min(d$q2.5, na.rm=T) - (2*sd(d$q50)), max(d$q97.5, na.rm=T) + (2*sd(d$q50)))
plot(xs, d$q50, ylim=ylim, xlim=range(xs), col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d$q2.5, xs, d$q97.5, len=0.1, angle=90, code=3, col=colreps, lwd=2)
points(xs, d$q50, pch=21, bg=bg, col=colreps, cex=1.5)
abline(v=mean(xs))
axis(2); box()
mtext('Standard Effect Size', 2, 3, font=2)
mtext('(A) Survival, LMA', 3, -1.25, font=2, at=1.75)
mtext('(B) Survival, WSG', 3, -1.25, font=2, at=4.5)
legend(2.5, -0.4, legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
legend(5.25, -0.4, legend=c('Cocoli','BCI','Sherman'), bty='n', pch=16, col=col, pt.cex=2)
axis(1, labels=rep(c('',''), 6), at=xs, cex.axis=1.5, las=2, tck = 0.03)

