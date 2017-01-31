setwd("/Users/Bob/Projects/Postdoc/Panama/RESULTS/_6.14.16")
library(jagsUI)
library(rjags)
library(RColorBrewer)

res <- data.frame()
for(j in 1:2){
  m <- c('growth','survival')[j]
  for(i in 1:6){
    filenames <- list.files(m)[!grepl('notrait',list.files(m))]
    file <- filenames[i]
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


col <- brewer.pal(11, 'BrBG')[c(3,9,11)]

plot_average_rates <- function(d1, d2, mtxt=NULL){
  
  par(mfrow=c(2,1), mar=c(0,4,2,0), oma=c(3,5,4,10))
  xs <- c(1,1.25, 2, 2.25, 3, 3.25)
  colreps <- rep(rep(col, each=2), times=2)
  bg <- ifelse(sign(d1$q2.5)==sign(d1$q97.5), colreps, 'white')
  ylim <- c(min(d1$q2.5, na.rm=T) - (0.5*sd(d1$q50)), max(d1$q97.5, na.rm=T) + (0.5*sd(d1$q50)))
  plot(xs, d1$q50, ylim=ylim, xlim=c(0.75,3.5), col=NA, axes=F, xlab='', ylab='')
  #abline(h=0, lty=2)
  arrows(xs, d1$q2.5, xs, d1$q97.5, len=0.075, angle=90, code=3, col=colreps, lwd=1.5)
  #segments(xs, d1$q5, xs, d1$q95, col=colreps, lwd=3)
  #points(xs, d1$q50, pch=21, bg=bg, col=colreps, cex=1.5)
  points(xs, d1$q50, pch=c(21,23), bg=bg, col=colreps, cex=1.5)
  axis(2); box()
  mtext('Predicted Avg. Growth', 2, 3.5, font=2)
  mtext(mtxt, 2, 4)
  mtext(bquote('( mm yr'^-1~')'), 2, 2)
  mtext('(A) Growth', 3, 0.25, font=2, at=0.75, adj=0)
  xs2 <- seq(1.125, 3.125, by=1)
  mtext(rep(c('Dry','Intermediate','Wet')), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
  #  mtext(rep(c('Cocoli','BCI','Sherman')), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
  axis(1, labels=rep(c('',''), 3), at=xs, cex.axis=1.5, las=2, tck = 0.03)
  #abline(v=c(1.625,2.625), lty=3)
  
  bg <- ifelse(sign(d2$q2.5)==sign(d2$q97.5), colreps, 'white')
  xs <- rep(1:(nrow(d2)/2), each=2) + rep(c(0, 0.25), times=(nrow(d2)/2))
  d2 <- d2*100
  ylim <- c(min(d2$q2.5, na.rm=T) - (0.5*sd(d2$q50)), max(d2$q97.5, na.rm=T) + (0.5*sd(d2$q50)))
  plot(xs, d2$q50, ylim=ylim, xlim=c(0.75,3.5), col=NA, axes=F, xlab='', ylab='')
  abline(h=0, lty=2)
  arrows(xs, d2$q2.5, xs, d2$q97.5, len=0.075, angle=90, code=3, col=colreps, lwd=1.5)
  #segments(xs, d2$q5, xs, d2$q95, col=colreps, lwd=3)
  points(xs, d2$q50, pch=c(21,23), bg=bg, col=colreps, cex=1.5)
  axis(2); box()
  mtext('Predicted Avg. Survival', 2, 3.5, font=2)
  mtext(bquote('( % yr'^-1~')'), 2, 2, font=2)
  mtext(mtxt, 2, 4)
  mtext('(B) Survival', 3, 0.25, font=2, at=0.75, adj=0)
  axis(1, labels=rep(c('',''), 3), at=xs, cex.axis=1.5, las=2, tck = 0.03)
  mtext(rep(c('Small','Large'), 6), side=1, line=0.3, at=xs, cex=0.8, las=2)
  #abline(v=c(1.625,2.625), lty=3)
  xs2 <- seq(1.125, 3.125, by=1)
  mtext(rep(c('Dry','Intermediate','Wet')), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
  # mtext(c('Cocoli','BCI','Sherman'), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
}

plot_nontrait_effects <- function(d1, d2, mtxt=NULL){
  par(mfrow=c(2,1), mar=c(0,4,2,0), oma=c(3,5,4,10))
  xs <- c(1,1.25, 2, 2.25, 3, 3.25)
  colreps <- rep(rep(col, each=2), times=2)
  bg <- ifelse(sign(d1$q2.5)==sign(d1$q97.5), colreps, 'white')
  ylim <- c(min(d1$q2.5, na.rm=T) - (1*sd(d1$q50)), max(d1$q97.5, na.rm=T) + (1*sd(d1$q50)))
  plot(xs, d1$q50, ylim=ylim, xlim=c(0.75,3.5), col=NA, axes=F, xlab='', ylab='')
  abline(h=0, lty=2)
  arrows(xs, d1$q2.5, xs, d1$q97.5, len=0.075, angle=90, code=3, col=colreps, lwd=1.5)
  #segments(xs, d1$q5, xs, d1$q95, col=colreps, lwd=3)
  points(xs, d1$q50, pch=c(21,23), bg=bg, col=colreps, cex=1.5)
  axis(2); box()
  mtext(mtxt, 2, 4)
  mtext('Standard Effect Size', 2, 2, font=2)
  mtext('(A) Growth', 3, 0.25, font=2, at=0.75, adj=0)
  xs2 <- seq(1.125, 3.125, by=1)
  mtext(rep(c('Dry','Intermediate','Wet')), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
  #   mtext(rep(c('Cocoli','BCI','Sherman')), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
  axis(1, labels=rep(c('',''), 3), at=xs, cex.axis=1.5, las=2, tck = 0.03)
  #abline(v=c(1.625,2.625), lty=3)
  
  bg <- ifelse(sign(d2$q2.5)==sign(d2$q97.5), colreps, 'white')
  xs <- rep(1:(nrow(d2)/2), each=2) + rep(c(0, 0.25), times=(nrow(d2)/2))
  ylim <- c(min(d2$q2.5, na.rm=T) - (1*sd(d2$q50)), max(d2$q97.5, na.rm=T) + (2*sd(d2$q50)))
  plot(xs, d2$q50, ylim=ylim, xlim=c(0.75,3.5), col=NA, axes=F, xlab='', ylab='')
  abline(h=0, lty=2)
  arrows(xs, d2$q2.5, xs, d2$q97.5, len=0.075, angle=90, code=3, col=colreps, lwd=1.5)
  #segments(xs, d2$q5, xs, d2$q95, col=colreps, lwd=3)
  points(xs, d2$q50, pch=c(21,23), bg=bg, col=colreps, cex=1.5)
  axis(2); box()
  mtext(mtxt, 2, 4)
  mtext('Standard Effect Size', 2, 2, font=2)
  mtext('(B) Survival', 3, 0.25, font=2, at=0.75, adj=0)
  axis(1, labels=rep(c('',''), 3), at=xs, cex.axis=1.5, las=2, tck = 0.03)
  mtext(rep(c('Small','Large'), 6), side=1, line=0.3, at=xs, cex=0.8, las=2)
  xs2 <- seq(1.125, 3.125, by=1)
  mtext(rep(c('Dry','Intermediate','Wet')), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
#   mtext(c('Cocoli','BCI','Sherman'), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
  #abline(v=c(1.625,2.625), lty=3)
}

plot_trait_effects <- function(d1, d2, mtxt=NULL){
  
  par(mfrow=c(2,1), mar=c(0,4,2,0), oma=c(5,3,4,1))
  colreps <- rep(rep(col, each=2), times=2)
  bg <- ifelse(sign(d1$q2.5)==sign(d1$q97.5), colreps, 'white')
  xs <- rep(1:(nrow(d1)/2), each=2) + rep(c(0, 0.25), times=(nrow(d1)/2))
  xlim <- c(min(xs)-0.1875, max(xs)+0.1875)
  ylim <- c(min(d1$q2.5, na.rm=T) - (0.5 * sd(d1$q50)), max(d1$q97.5, na.rm=T) + (0.75 * sd(d1$q50)))
  plot(xs, d1$q50, ylim=ylim, xlim=xlim, col=NA, axes=F, xlab='', ylab='')
  abline(h=0, lty=2)
  arrows(xs, d1$q2.5, xs, d1$q97.5, len=0.075, angle=90, code=3, col=colreps, lwd=1.5)
  #segments(xs, d1$q5, xs, d1$q95, col=colreps, lwd=3)
  points(xs, d1$q50, pch=c(21,23), bg=bg, col=colreps, cex=1.5)
  abline(v=mean(xs))
  axis(2); box()
  mtext('Standard Effect Size', 2, 2, font=2)
  mtext('(A) Growth: Wood density', 3, 0.25, font=2, at=0.6875, adj=0)
  mtext('(B) Growth: LMA', 3, 0.25, font=2, at=3.75, adj=0)
  axis(1, labels=rep(c('',''), 6), at=xs, cex.axis=1.5, las=2, tck = 0.03)
  xs2 <- seq(1.125, 6.125, by=1)
  mtext(rep(c('Dry','Intermediate','Wet'), 2), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
#   mtext(rep(c('Cocoli','BCI','Sherman'), 2), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
#  abline(v=c(1.625,2.625,4.625,5.625), lty=3)
  mtext(mtxt, 2, 4)

  bg <- ifelse(sign(d2$q2.5)==sign(d2$q97.5), colreps, 'white')
  xs <- rep(1:(nrow(d2)/2), each=2) + rep(c(0, 0.25), times=(nrow(d2)/2))
  ylim <- c(min(d2$q2.5, na.rm=T) - (1 * sd(d2$q50)), max(d2$q97.5, na.rm=T) + (1 * sd(d2$q50)))
  xlim <- c(min(xs)-0.1875, max(xs)+0.1875)
  plot(xs, d2$q50, ylim=ylim, xlim=xlim, col=NA, axes=F, xlab='', ylab='')
  abline(h=0, lty=2)
  arrows(xs, d2$q2.5, xs, d2$q97.5, len=0.075, angle=90, code=3, col=colreps, lwd=1.5)
  #segments(xs, d2$q5, xs, d2$q95, col=colreps, lwd=3)
  points(xs, d2$q50, pch=c(21,23), bg=bg, col=colreps, cex=1.5)
  abline(v=mean(xs))
  axis(2); box()
  mtext('Standard Effect Size', 2, 2, font=2)
  mtext('(C) Survival: Wood density', 3, 0.25, font=2, at=0.6875, adj=0)
  mtext('(D) Survival: LMA', 3, 0.25, font=2, at=3.75, adj=0)
  axis(1, labels=rep(c('',''), 6), at=xs, cex.axis=1.5, las=2, tck = 0.03)
  mtext(rep(c('Small','Large'), 6), side=1, line=0.3, at=xs, cex=0.8, las=2)
  xs2 <- seq(1.125, 6.125, by=1)
  mtext(rep(c('Dry','Intermediate','Wet'), 2), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
  mtext(mtxt, 2, 4)
#   mtext(rep(c('Cocoli','BCI','Sherman'), 2), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
#  abline(v=c(1.625,2.625,4.625,5.625), lty=3)  
}

plot_conceptual <- function(d1, d2){
  par(mfrow=c(2,1), mar=c(0,4,2,0), oma=c(3,5,4,10))
  xs <- c(1,1.25, 2, 2.25, 3, 3.25)
  colreps <- rep(rep(col, each=2), times=2)
  bg <- colreps
  ylim <- c(min(d1$q50, na.rm=T) - (1*sd(d1$q50)), max(d1$q50, na.rm=T) + (1*sd(d1$q50)))
  plot(xs, d1$q50, ylim=ylim, xlim=c(0.75,3.5), col=NA, axes=F, xlab='', ylab='')
  abline(h=0, lty=2)
  points(xs, d1$q50, pch=c(21,23), bg=bg, col=colreps, cex=1.5)
  axis(2, labels=c('-','+'), at=range(d1$q50), cex.axis=2); box()
  mtext('Standard Effect Size', 2, 2, font=2)
  mtext('(A) Growth', 3, 0.25, font=2, at=0.75, adj=0)
  xs2 <- seq(1.125, 3.125, by=1)
  mtext(rep(c('Dry','Intermediate','Wet')), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
  axis(1, labels=rep(c('',''), 3), at=xs, cex.axis=1.5, las=2, tck = 0.03)
  
  xs <- rep(1:(nrow(d2)/2), each=2) + rep(c(0, 0.25), times=(nrow(d2)/2))
  ylim <- c(min(d2$q50, na.rm=T) - (1*sd(d2$q50)), max(d2$q50, na.rm=T) + (2*sd(d2$q50)))
  plot(xs, d2$q50, ylim=ylim, xlim=c(0.75,3.5), col=NA, axes=F, xlab='', ylab='')
  abline(h=0, lty=2)
  points(xs, d2$q50, pch=c(21,23), bg=bg, col=colreps, cex=1.5)
  axis(2, labels=c('-','+'), at=range(d1$q50), cex.axis=2); box()
  mtext('Standard Effect Size', 2, 2, font=2)
  mtext('(B) Survival', 3, 0.25, font=2, at=0.75, adj=0)
  axis(1, labels=rep(c('',''), 3), at=xs, cex.axis=1.5, las=2, tck = 0.03)
  mtext(rep(c('Small','Large'), 6), side=1, line=0.3, at=xs, cex=0.8, las=2)
  xs2 <- seq(1.125, 3.125, by=1)
  mtext(rep(c('Dry','Intermediate','Wet')), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
}



##### AVERAGE RATES #####
#pdf(file="/Users/Bob/Projects/Postdoc/Panama/Figures/mu.beta1.pdf")

pdf(file="/Users/Bob/Projects/Postdoc/Panama/Figures/FigS2_average_rates.pdf")
g <- res[res$m=='growth',]
b <- c('mu.beta[1]')
d1 <- g[as.character(g$param) %in% b,]
d1 <- d1[order(rev(d1$param), d1$p, rev(d1$s)),]
d1[,5:9] <- d1[,5:9] * f[c(2,1,4,3,6,5)]
g <- res[res$m=='survival',]
b <- c('mu.beta[1]')
d2 <- g[as.character(g$param) %in% b,]
d2 <- d2[order(rev(d2$param), d2$p, rev(d2$s)),]
d2 <- (exp(d2[,5:9])/(1 + exp(d2[,5:9])))
plot_average_rates(d1,d2)
dev.off()


##### TRAIT EFFECT ON AVERAGE RATES #####
#pdf(file="/Users/Bob/Projects/Postdoc/Panama/Figures/beta.t1.pdf")
pdf(file="/Users/Bob/Projects/Postdoc/Panama/Figures/FigS1_trait_effects_on_avg_rates.pdf")
g <- res[res$m=='growth',]
b <- c('beta.wd[1]', 'beta.lma[1]')
d1 <- g[as.character(g$param) %in% b,]
d1 <- d1[order(rev(d1$param), d1$p, rev(d1$s)),]
g <- res[res$m=='survival',]
b <- c('beta.wd[1]', 'beta.lma[1]')
d2 <- g[as.character(g$param) %in% b,]
d2 <- d2[order(rev(d2$param), d2$p, rev(d2$s)),]
plot_trait_effects(d1, d2, mtxt='Trait effect 
on Average Rates')
dev.off()


##### TRAIT EFFECT ON NCI EFFECT #####
#pdf(file="/Users/Bob/Projects/Postdoc/Panama/Figures/beta.t2.pdf")
pdf(file="/Users/Bob/Projects/Postdoc/Panama/Figures/Fig4_trait_effect_on_NCIsensitivity.pdf")
g <- res[res$m=='growth',]
b <- c('beta.wd[2]', 'beta.lma[2]')
d1 <- g[as.character(g$param) %in% b,]
d1 <- d1[order(rev(d1$param), d1$p, rev(d1$s)),]
g <- res[res$m=='survival',]
b <- c('beta.wd[2]', 'beta.lma[2]')
d2 <- g[as.character(g$param) %in% b,]
d2 <- d2[order(rev(d2$param), d2$p, rev(d2$s)),]

plot_trait_effects(d1, d2, mtxt='Trait effect
on sensitivity to NCI')


dev.off()


##### AVERAGE NCI EFFECT #####
#pdf(file="/Users/Bob/Projects/Postdoc/Panama/Figures/mu.beta2.pdf")
pdf(file="/Users/Bob/Projects/Postdoc/Panama/Figures/Fig3_avg_effect_of_NCI.pdf")
g <- res[res$m=='growth',]
b <- c('mu.beta[2]')
d1 <- g[as.character(g$param) %in% b,]
d1 <- d1[order(rev(d1$param), d1$p, rev(d1$s)),]
g <- res[res$m=='survival',]
b <- c('mu.beta[2]')
d2 <- g[as.character(g$param) %in% b,]
d2 <- d2[order(rev(d2$param), d2$p, rev(d2$s)),]
plot_nontrait_effects(d1,d2,mtxt='Average effect of NCI')
dev.off()



##### AVERAGE SIZE EFFECT #####
#pdf(file="/Users/Bob/Projects/Postdoc/Panama/Figures/mu.beta3.pdf")
pdf(file="/Users/Bob/Projects/Postdoc/Panama/Figures/FigS3_avg_size_effect.pdf")
g <- res[res$m=='growth',]
b <- c('mu.beta[3]')
d1 <- g[as.character(g$param) %in% b,]
d1 <- d1[order(rev(d1$param), d1$p, rev(d1$s)),]
g <- res[res$m=='survival',]
b <- c('mu.beta[3]')
d2 <- g[as.character(g$param) %in% b,]
d2 <- d2[order(rev(d2$param), d2$p, rev(d2$s)),]
plot_nontrait_effects(d1,d2,mtxt='Average size effect')
dev.off()



##### Conceptual Figure #####
# pdf(file="/Users/Bob/Projects/Postdoc/Panama/Figures/Conceptual_Fig1.pdf")
d1 <- data.frame(m='growth', 
                 p=rep(1:3, each=2), 
                 q50=c(1,0.6,0,0,-1,-0.6))
d2 <- data.frame(m='growth', 
                 p=rep(1:3, each=2), 
                 q50=c(1,0.6,0,0,-1,-0.6))
plot_conceptual(d1,d2)

# dev.off()





##### AVERAGE NCI EFFECT (with new y-axis labels) #####
pdf(file="/Users/Bob/Projects/Postdoc/Panama/Figures/Fig3_avg_effect_of_NCI_V2.pdf")

g <- res[res$m=='growth',]
b <- c('mu.beta[2]')
d1 <- g[as.character(g$param) %in% b,]
d1 <- d1[order(rev(d1$param), d1$p, rev(d1$s)),]
g <- res[res$m=='survival',]
b <- c('mu.beta[2]')
d2 <- g[as.character(g$param) %in% b,]
d2 <- d2[order(rev(d2$param), d2$p, rev(d2$s)),]

par(mfrow=c(2,1), mar=c(0,4,2,0), oma=c(3,5,4,10))
  xs <- c(1,1.25, 2, 2.25, 3, 3.25)
  colreps <- rep(rep(col, each=2), times=2)
  bg <- ifelse(sign(d1$q2.5)==sign(d1$q97.5), colreps, 'white')
  ylim <- c(min(d1$q2.5, na.rm=T) - (1*sd(d1$q50)), max(d1$q97.5, na.rm=T) + (1*sd(d1$q50)))
  plot(xs, d1$q50, ylim=ylim, xlim=c(0.75,3.5), col=NA, axes=F, xlab='', ylab='')
  abline(h=0, lty=2)
  arrows(xs, d1$q2.5, xs, d1$q97.5, len=0.075, angle=90, code=3, col=colreps, lwd=1.5)
  #segments(xs, d1$q5, xs, d1$q95, col=colreps, lwd=4)
  points(xs, d1$q50, pch=c(21,23), bg=bg, col=colreps, cex=1.5)
  axis(2); box()
  mtext('Standard Effect Size', 2, 2, font=2)
  mtext('(A) Growth', 3, 0.25, font=2, at=0.75, adj=0)
  xs2 <- seq(1.125, 3.125, by=1)
  mtext(rep(c('Dry','Intermediate','Wet')), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
  axis(1, labels=rep(c('',''), 3), at=xs, cex.axis=1.5, las=2, tck = 0.03)
  mtext('(+) Crowding 
        increases 
        performance', 2, las=2, 4, cex=0.75, at=ylim[2]/1.3)
  mtext('(-) Crowding 
        decreases 
        performance', 2, las=2, 4, cex=0.75, at=ylim[1]/1.3)
  bg <- ifelse(sign(d2$q2.5)==sign(d2$q97.5), colreps, 'white')
  xs <- rep(1:(nrow(d2)/2), each=2) + rep(c(0, 0.25), times=(nrow(d2)/2))
  ylim <- c(min(d2$q2.5, na.rm=T) - (1*sd(d2$q50)), max(d2$q97.5, na.rm=T) + (2*sd(d2$q50)))
  plot(xs, d2$q50, ylim=ylim, xlim=c(0.75,3.5), col=NA, axes=F, xlab='', ylab='')
  abline(h=0, lty=2)
  arrows(xs, d2$q2.5, xs, d2$q97.5, len=0.075, angle=90, code=3, col=colreps, lwd=1.5)
  #segments(xs, d2$q5, xs, d2$q95, col=colreps, lwd=4)
  points(xs, d2$q50, pch=c(21,23), bg=bg, col=colreps, cex=1.5)
  axis(2); box()
  mtext('Standard Effect Size', 2, 2, font=2)
  mtext('(B) Survival', 3, 0.25, font=2, at=0.75, adj=0)
  axis(1, labels=rep(c('',''), 3), at=xs, cex.axis=1.5, las=2, tck = 0.03)
  mtext(rep(c('Small','Large'), 6), side=1, line=0.3, at=xs, cex=0.8, las=2)
  xs2 <- seq(1.125, 3.125, by=1)
  mtext(rep(c('Dry','Intermediate','Wet')), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
  mtext('(+) Crowding 
        increases 
        performance', 2, las=2, 4, cex=0.75, at=ylim[2]/1.3)
  mtext('(-) Crowding 
        decreases 
        performance', 2, las=2, 4, cex=0.75, at=ylim[1]/1.3)

dev.off()



##### TRAIT EFFECT ON NCI EFFECT (with new y-axis labels) #####
pdf(file="/Users/Bob/Projects/Postdoc/Panama/Figures/Fig4_trait_effect_on_NCIsensitivity_V2.pdf")

g <- res[res$m=='growth',]
b <- c('beta.wd[2]', 'beta.lma[2]')
d1 <- g[as.character(g$param) %in% b,]
d1 <- d1[order(rev(d1$param), d1$p, rev(d1$s)),]
g <- res[res$m=='survival',]
b <- c('beta.wd[2]', 'beta.lma[2]')
d2 <- g[as.character(g$param) %in% b,]
d2 <- d2[order(rev(d2$param), d2$p, rev(d2$s)),]

par(mfrow=c(2,1), mar=c(0,5,2,0), oma=c(5,3,4,1))
colreps <- rep(rep(col, each=2), times=2)
bg <- ifelse(sign(d1$q2.5)==sign(d1$q97.5), colreps, 'white')
xs <- rep(1:(nrow(d1)/2), each=2) + rep(c(0, 0.25), times=(nrow(d1)/2))
xlim <- c(min(xs)-0.1875, max(xs)+0.1875)
ylim <- c(min(d1$q2.5, na.rm=T) - (0.5 * sd(d1$q50)), max(d1$q97.5, na.rm=T) + (0.75 * sd(d1$q50)))
plot(xs, d1$q50, ylim=ylim, xlim=xlim, col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d1$q2.5, xs, d1$q97.5, len=0.075, angle=90, code=3, col=colreps, lwd=1.5)
#segments(xs, d1$q5, xs, d1$q95, col=colreps, lwd=4)
points(xs, d1$q50, pch=c(21,23), bg=bg, col=colreps, cex=1.5)
abline(v=mean(xs))
axis(2); box()
mtext('Standard Effect Size', 2, 2, font=2)
mtext('(A) Growth: Wood density', 3, 0.25, font=2, at=0.6875, adj=0)
mtext('(B) Growth: LMA', 3, 0.25, font=2, at=3.75, adj=0)
axis(1, labels=rep(c('',''), 6), at=xs, cex.axis=1.5, las=2, tck = 0.03)
xs2 <- seq(1.125, 6.125, by=1)
mtext(rep(c('Dry','Intermediate','Wet'), 2), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
mtext('(+) Acquisitive
      species more
      sensitive to
      crowding', 2, las=2, 3.25, cex=0.75, at=ylim[2]/1.3)
mtext('(-) Conservative
      species more
      sensitive to
      crowding', 2, las=2, 3.25, cex=0.75, at=ylim[1]/1.3)

bg <- ifelse(sign(d2$q2.5)==sign(d2$q97.5), colreps, 'white')
xs <- rep(1:(nrow(d2)/2), each=2) + rep(c(0, 0.25), times=(nrow(d2)/2))
ylim <- c(min(d2$q2.5, na.rm=T) - (3 * sd(d2$q50)), max(d2$q97.5, na.rm=T) + (3 * sd(d2$q50)))
xlim <- c(min(xs)-0.1875, max(xs)+0.1875)
plot(xs, d2$q50, ylim=ylim, xlim=xlim, col=NA, axes=F, xlab='', ylab='')
abline(h=0, lty=2)
arrows(xs, d2$q2.5, xs, d2$q97.5, len=0.075, angle=90, code=3, col=colreps, lwd=1.5)
#segments(xs, d2$q5, xs, d2$q95, col=colreps, lwd=4)
points(xs, d2$q50, pch=c(21,23), bg=bg, col=colreps, cex=1.5)
abline(v=mean(xs))
axis(2); box()
mtext('Standard Effect Size', 2, 2, font=2)
mtext('(C) Survival: Wood density', 3, 0.25, font=2, at=0.6875, adj=0)
mtext('(D) Survival: LMA', 3, 0.25, font=2, at=3.75, adj=0)
axis(1, labels=rep(c('',''), 6), at=xs, cex.axis=1.5, las=2, tck = 0.03)
mtext(rep(c('Small','Large'), 6), side=1, line=0.3, at=xs, cex=0.8, las=2)
xs2 <- seq(1.125, 6.125, by=1)
mtext(rep(c('Dry','Intermediate','Wet'), 2), side=3, line=-1, at=xs2, cex.axis=1.5, col=col, font=2)
mtext('(+) Acquisitive
      species more
      sensitive to
      crowding', 2, las=2, 3.25, cex=0.75, at=ylim[2]/1.3)
mtext('(-) Conservative
species more
sensitive to
crowding', 2, las=2, 3.25, cex=0.75, at=ylim[1]/1.3)

dev.off()




