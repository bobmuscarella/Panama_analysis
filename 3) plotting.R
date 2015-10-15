setwd("/Users/Bob/Projects/Postdoc/Panama/Results")
library(rjags)
load("coda_APR2_notrait_allplots.RDA")
plot(coda.results, ask=T)
load("summary_APR2_notrait_allplots.RDA")
res






library(RColorBrewer)
col <- brewer.pal(11, 'BrBG')[c(3,9,11)]
plot((1:3), cocoli.res[[2]][,3], ylim=c(-1.25,1.25), pch=21, bg=col[1], xlim=c(0.9,3.3), cex=2, axes=F, xlab='Parameter', ylab='Standardized Effect Size')
axis(1, at=c(1.1,2.1,3.1), labels=c('Avg. growth', 'DBH effect', 'NCI effect'), las=1, cex.axis=1)
axis(2)
segments((1:3), cocoli.res[[2]][,1], (1:3), cocoli.res[[2]][,5], col=col[1], lwd=4)
segments((1:3)+.1, bci.res[[2]][,1], (1:3)+.1, bci.res[[2]][,5], col=col[2], lwd=4)
points(1:3, cocoli.res[[2]][,3], pch=21, bg=col[1], cex=2)
points((1:3)+.1, bci.res[[2]][,3], pch=21, bg=col[2], cex=2)
segments((1:3)+.2, sherman.res[[2]][,1], (1:3)+.2, sherman.res[[2]][,5], col=col[3], lwd=4)
points((1:3)+.2, sherman.res[[2]][,3], pch=21, bg=col[3], cex=2)
abline(h=0, lty=3)
box()
legend('topleft', legend=c('cocoli (dry: 1950 mm/yr)', 'bci (moist: 2500 mm/yr)', 'sherman (moist: 2700-3000 mm/yr)'), pch=21, pt.bg=col, bty='n', pt.cex=2)









library(RColorBrewer)
col <- brewer.pal(11, 'BrBG')[c(1,3,9,11)]
col <- terrain.colors(8)[c(6,3,2,1)]


plot((1:3), res[[2]][c('mu.beta.1[2]','mu.beta.2[2]','mu.beta.3[2]'),3], ylim=c(-2,1.5), pch=21, bg=col[1], xlim=c(0.9,3.3), cex=1.5, axes=F, xlab='Parameter', ylab='Standardized Effect Size')
segments((1:3), res[[2]][c('mu.beta.1[2]','mu.beta.2[2]','mu.beta.3[2]'),1], (1:3), res[[2]][c('mu.beta.1[2]','mu.beta.2[2]','mu.beta.3[2]'),5], col=col[1], lwd=4)
points((1:3), res[[2]][c('mu.beta.1[2]','mu.beta.2[2]','mu.beta.3[2]'),3], pch=21, bg=col[1], cex=1.5)

segments((1:3)+.1, res[[2]][c('mu.beta.1[1]','mu.beta.2[1]','mu.beta.3[1]'),1], (1:3)+.1, res[[2]][c('mu.beta.1[1]','mu.beta.2[1]','mu.beta.3[1]'),5], col=col[2], lwd=4)
points((1:3)+.1, res[[2]][c('mu.beta.1[1]','mu.beta.2[1]','mu.beta.3[1]'),3], pch=21, bg=col[2], cex=1.5)

segments((1:3)+.2, res[[2]][c('mu.beta.1[4]','mu.beta.2[4]','mu.beta.3[4]'),1], (1:3)+.2, res[[2]][c('mu.beta.1[4]','mu.beta.2[4]','mu.beta.3[4]'),5], col=col[3], lwd=4)
points((1:3)+.2, res[[2]][c('mu.beta.1[4]','mu.beta.2[4]','mu.beta.3[4]'),3], pch=21, bg=col[3], cex=1.5)

segments((1:3)+.3, res[[2]][c('mu.beta.1[3]','mu.beta.2[3]','mu.beta.3[3]'),1], (1:3)+.3, res[[2]][c('mu.beta.1[3]','mu.beta.2[3]','mu.beta.3[3]'),5], col=col[4], lwd=4)
points((1:3)+.3, res[[2]][c('mu.beta.1[3]','mu.beta.2[3]','mu.beta.3[3]'),3], pch=21, bg=col[4], cex=1.5)

abline(h=0, lty=3)
box()

axis(1, at=1:3+.15, labels=c('Avg. growth', 'NCI effect', 'DBH effect'), las=1, cex.axis=1)
axis(2)

legend('topleft', legend=c('cocoli (dry: 1950 mm/yr)', 'bci (moist: 2500 mm/yr)', 'sherman (moist: 2850 mm/yr)','luquillo (wet: 3500 mm/yr)'), pch=21, pt.bg=col, bty='n', pt.cex=2)








setwd("/Users/Bob/Projects/Postdoc/Panama/Results")
library(rjags)
load("coda_APR4_trait_allplots.RDA")
#plot(coda.results, ask=T)
load("summary_APR4_trait_allplots.RDA")
res


library(RColorBrewer)
col <- brewer.pal(11, 'BrBG')[c(1,3,9,11)]
col <- terrain.colors(8)[c(6,3,2,1)]


plot((1:2), res[[2]][c('beta.t.1[2]','beta.t.2[2]'),3], ylim=c(-4,4), pch=21, bg=col[1], xlim=c(0.9,2.3), cex=1.5, axes=F, xlab='Parameter', ylab='Standardized Effect Size')
segments((1:2), res[[2]][c('beta.t.1[2]','beta.t.2[2]'),1], (1:2), res[[2]][c('beta.t.1[2]','beta.t.2[2]'),5], col=col[1], lwd=4)
points((1:2), res[[2]][c('beta.t.1[2]','beta.t.2[2]'),3], pch=21, bg=col[1], cex=1.5)

segments((1:2)+.1, res[[2]][c('beta.t.1[1]','beta.t.2[1]'),1], (1:2)+.1, res[[2]][c('beta.t.1[1]','beta.t.2[1]'),5], col=col[2], lwd=4)
points((1:2)+.1, res[[2]][c('beta.t.1[1]','beta.t.2[1]'),3], pch=21, bg=col[2], cex=1.5)

segments((1:2)+.2, res[[2]][c('beta.t.1[4]','beta.t.2[4]'),1], (1:2)+.2, res[[2]][c('beta.t.1[4]','beta.t.2[4]'),5], col=col[3], lwd=4)
points((1:2)+.2, res[[2]][c('beta.t.1[4]','beta.t.2[4]'),3], pch=21, bg=col[3], cex=1.5)

segments((1:2)+.3, res[[2]][c('beta.t.1[3]','beta.t.2[3]'),1], (1:2)+.3, res[[2]][c('beta.t.1[3]','beta.t.2[3]'),5], col=col[4], lwd=4)
points((1:2)+.3, res[[2]][c('beta.t.1[4]','beta.t.2[4]'),3], pch=21, bg=col[4], cex=1.5)

abline(h=0, lty=3)
box()

axis(1, at=1:2+.15, labels=c('Trait effect on avg. growth', 'Trait effect on NCI effect'), las=1, cex.axis=1)
axis(2)

legend('topleft', legend=c('cocoli (dry: 1950 mm/yr)', 'bci (moist: 2500 mm/yr)', 'sherman (moist: 2850 mm/yr)','luquillo (wet: 3500 mm/yr)'), pch=21, pt.bg=col, bty='n', pt.cex=2)




library(RColorBrewer)
col <- brewer.pal(11, 'BrBG')[c(1,3,9,11)]
col <- terrain.colors(8)[c(6,3,2,1)]


plot((1:2), res[[2]][c('beta.t.1[2]','beta.t.2[2]'),3], ylim=c(-4,4), pch=21, bg=col[1], xlim=c(0.9,3.3), cex=1.5, axes=F, xlab='Parameter', ylab='Standardized Effect Size')
segments((1:2), res[[2]][c('beta.t.1[2]','beta.t.2[2]'),1], (1:2), res[[2]][c('beta.t.1[2]','beta.t.2[2]'),5], col=col[1], lwd=4)
points((1:2), res[[2]][c('beta.t.1[2]','beta.t.2[2]'),3], pch=21, bg=col[1], cex=1.5)

segments((1:2)+.1, res[[2]][c('beta.t.1[1]','beta.t.2[1]'),1], (1:2)+.1, res[[2]][c('beta.t.1[1]','beta.t.2[1]'),5], col=col[2], lwd=4)
points((1:2)+.1, res[[2]][c('beta.t.1[1]','beta.t.2[1]'),3], pch=21, bg=col[2], cex=1.5)

segments((1:2)+.2, res[[2]][c('beta.t.1[4]','beta.t.2[4]'),1], (1:2)+.2, res[[2]][c('beta.t.1[4]','beta.t.2[4]'),5], col=col[3], lwd=4)
points((1:2)+.2, res[[2]][c('beta.t.1[4]','beta.t.2[4]'),3], pch=21, bg=col[3], cex=1.5)

segments((1:2)+.3, res[[2]][c('beta.t.1[3]','beta.t.2[3]'),1], (1:2)+.3, res[[2]][c('beta.t.1[3]','beta.t.2[3]'),5], col=col[4], lwd=4)
points((1:2)+.3, res[[2]][c('beta.t.1[4]','beta.t.2[4]'),3], pch=21, bg=col[4], cex=1.5)

segments(3+c(.1,.2,.3), res[[2]][c('mu.beta[1]','mu.beta[2]','mu.beta[3]'),1], 3+c(.1,.2,.3), res[[2]][c('mu.beta[1]','mu.beta[2]','mu.beta[3]'),5], col=1, lwd=4)
points(3+c(.1,.2,.3), res[[2]][c('mu.beta[1]','mu.beta[2]','mu.beta[3]'),3], pch=21, bg='white', cex=1.5)

abline(h=0, lty=3)
box()
axis(1, at=1:2+.15, labels=c('Trait effect on avg. growth', 'Trait effect on NCI effect'), las=1, cex.axis=1)
axis(2)

text(3+c(.1,.2,.3), -1, labels=c('intercept for avg. growth', 'intercept for NCI effect', 'intercept for DBh effect'), srt=90, adj=1)

legend('topleft', legend=c('cocoli (dry: 1950 mm/yr)', 'bci (moist: 2500 mm/yr)', 'sherman (moist: 2850 mm/yr)','luquillo (wet: 3500 mm/yr)'), pch=21, pt.bg=col, bty='n', pt.cex=2)

