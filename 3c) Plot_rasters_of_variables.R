# Plot rasters of variables

#########################################
####  START HERE WITH PROCESSED DATA ####
#########################################
library(sp)
library(colorRamps)
library(fields)
library(raster)

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")
load("Panama_AnalysisData_6.14.16.RDA")

head(tdata)

g <- 15

tdata$x <- ifelse(tdata$x==g, tdata$x+0.01, tdata$x)
tdata$y <- ifelse(tdata$y==g, tdata$y+0.01, tdata$y)

dat <- tdata[tdata$census %in% c(1,4) & tdata$Not.Edge==1,]

dat <- dat[order(dat$x, dat$y),]

xquad <- as.numeric(cut(dat$x, seq(0, 1000, g), include.lowest=T))
yquad <- as.numeric(cut(dat$y, seq(0, 1000, g), include.lowest=T))


for(i in 1:4){

  variable <- list(dat$wd, dat$lma, dat$log.lma, log(dat$All.NCI))[[i]]
  
  val <- tapply(variable, paste(xquad, yquad, dat$plot), mean, na.rm=T)
  xq <- tapply(xquad, paste(xquad, yquad, dat$plot), mean)
  yq <- tapply(yquad, paste(xquad, yquad, dat$plot), mean)
  p <- tapply(dat$plot, paste(xquad, yquad, dat$plot), unique)
  df <- data.frame(val=val, x=xq, y=yq, p=p)
  
fileend <- c("WD_grid_15.pdf","LMA_grid_15.pdf","log_LMA_grid_15.pdf","NCI_grid_15.pdf")[i]
file <- paste("/Users/Bob/Projects/Postdoc/Panama/Figures/", fileend, sep='')
pdf(file=file)

### PLOTTING
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

b <- seq(min(df$val, na.rm=T)-0.01, max(df$val, na.rm=T)+0.02, length.out=10)
#col <- rev(blue2green(length(b)-1))
#col <- rev(blue2yellow(length(b)-1))
col <- brewer.pal(9, 'YlGnBu')

par(mar=c(5,4,2,8))
x <- df[df$p==2,]
pts <- SpatialPointsDataFrame(coords=x[,c(2,3)]*g, data=as.data.frame(x$val))
r <- rasterFromXYZ(as.data.frame(pts)[, c(2,3,1)])
image(r, breaks=b, col=col, axes=F, ylim=c(0,500), xlim=c(0,1000), xlab='', ylab='')
axis(1, labels=seq(0,1000,50), at=seq(0,1000,50))
axis(2, labels=seq(0,1000,50), at=seq(0,1000,50))
mtext('BCI (Intermediate)', 3, 0, at=500, outer=F, font=2, cex=1.5)
mtext('X (m)', 1, 2.5, outer=F, font=2)
mtext('Y (m)', 2, 2.5, outer=F, font=2)
if(i==1) { mtext(bquote('Wood Density ( g cm'^-3~')'), 2, -32.5, font=2) }
if(i==2) { mtext(bquote('LMA ( g m'^-2~')'), 2, -32.5, font=2) }
if(i==3) { mtext(bquote('LMA log( g m'^-2~')'), 2, -32.5, font=2) }
if(i==4) { mtext('log(NCI)', 2, -32.5, font=2) }         
image.plot(r, legend.only=T, breaks=b, col=col, add=T, legend.shrink = 0.9, legend.mar=4.1)

plot(0,0,col=NA, axes=F, ylab='', xlab='')

par(mar=c(4,4,1,1))
x <- df[df$p==1,]
pts <- SpatialPointsDataFrame(coords=x[,c(2,3)]*g, data=as.data.frame(x$val))
r <- rasterFromXYZ(as.data.frame(pts)[, c(2,3,1)])
image(r, breaks=b, col=col, axes=F, ylim=c(0,300), xlim=c(0,300), xlab='', ylab='')
axis(1, labels=seq(0,200,50), at=seq(0,200,50))
axis(2, labels=seq(0,500,50), at=seq(0,500,50))
mtext('Cocoli (Dry)', 3, 0, at=100, outer=F, font=2, cex=1.5)
mtext('X (m)', 1, 2.5, at=100, outer=F, font=2)
mtext('Y (m)', 2, 2.5, outer=F, font=2)

x <- df[df$p==3,]
pts <- SpatialPointsDataFrame(coords=x[,c(2,3)]*g, data=as.data.frame(x$val))
r <- rasterFromXYZ(as.data.frame(pts)[, c(2,3,1)])
image(r, breaks=b, col=col, axes=F, ylim=c(0,350), xlim=c(0,350), xlab='', ylab='')
axis(1, labels=seq(0,250,50), at=seq(0,250,50))
axis(2, labels=seq(0,500,50), at=seq(0,500,50))
mtext('Sherman (Wet)', 3, 0, at=125, outer=F, font=2, cex=1.5)
mtext('X (m)', 1, 2.5, at=125, outer=F, font=2)
mtext('Y (m)', 2, 2.5, outer=F, font=2)

dev.off()
}








#======
#########################################
####  PLAYGROUND ####
#########################################
head(tdata)

tdata$x <- ifelse(tdata$x==20, tdata$x+0.01, tdata$x)
tdata$y <- ifelse(tdata$y==20, tdata$y+0.01, tdata$y)
dat <- tdata[tdata$census %in% c(1,4),]
dat <- dat[order(dat$x, dat$y),]
g <- 5
xquad <- as.numeric(cut(dat$x, seq(0, 1000, g), include.lowest=T))
yquad <- as.numeric(cut(dat$y, seq(0, 1000, g), include.lowest=T))

var1 <- log(dat$nci)
var2 <- dat$dbh
val1 <- tapply(var1, paste(xquad, yquad, dat$plot), mean, na.rm=T)
val2 <- tapply(var2, paste(xquad, yquad, dat$plot), length)
val3 <- tapply(var2, paste(xquad, yquad, dat$plot), mean, na.rm=T)

xq <- tapply(xquad, paste(xquad, yquad, dat$plot), mean)
yq <- tapply(yquad, paste(xquad, yquad, dat$plot), mean)
p <- tapply(dat$plot, paste(xquad, yquad, dat$plot), unique)

df <- data.frame(v1=val1, v2=val2, v3=val3, x=xq, y=yq, p=p)





par(mfcol=c(3,2))

xlim <- range(df$v1)
ylim <- range(df$v2)


plot(df$v1[df$p=='cocoli'], df$v2[df$p=='cocoli'], xlab='nci', ylab='stem density', pch=16, col=rgb(0,0,0,.4), xlim=xlim, ylim=ylim)
mtext(round(cor(df$v1[df$p=='cocoli'], df$v2[df$p=='cocoli']), 2), 3, -2, col=4)
abline(v=unlist(tapply(df$v1, df$p, range)), col=rep(1:3, each=2))
abline(h=unlist(tapply(df$v2, df$p, range)), col=rep(1:3, each=2))

plot(df$v1[df$p=='bci'], df$v2[df$p=='bci'], xlab='nci', ylab='stem density', pch=16, col=rgb(0,0,0,.4), xlim=xlim, ylim=ylim)
mtext(round(cor(df$v1[df$p=='bci'], df$v2[df$p=='bci']), 2), 3, -2, col=4)
abline(v=unlist(tapply(df$v1, df$p, range)), col=rep(1:3, each=2))
abline(h=unlist(tapply(df$v2, df$p, range)), col=rep(1:3, each=2))

plot(df$v1[df$p=='sherman'], df$v2[df$p=='sherman'], xlab='nci', ylab='stem density', pch=16, col=rgb(0,0,0,.4), xlim=xlim, ylim=ylim)
mtext(round(cor(df$v1[df$p=='sherman'], df$v2[df$p=='sherman']), 2), 3, -2, col=4)
abline(v=unlist(tapply(df$v1, df$p, range)), col=rep(1:3, each=2))
abline(h=unlist(tapply(df$v2, df$p, range)), col=rep(1:3, each=2))




ylim <- range(df$v3)

plot(df$v1[df$p=='cocoli'], df$v3[df$p=='cocoli'], xlab='nci', ylab='avg dbh', pch=16, col=rgb(0,0,0,.4), xlim=xlim, ylim=ylim)
mtext(round(cor(df$v1[df$p=='cocoli'], df$v3[df$p=='cocoli']), 2), 3, -2, col=4)
abline(v=unlist(tapply(df$v1, df$p, range)), col=rep(1:3, each=2))
abline(h=unlist(tapply(df$v3, df$p, range)), col=rep(1:3, each=2))

plot(df$v1[df$p=='bci'], df$v3[df$p=='bci'], xlab='nci', ylab='avg dbh', pch=16, col=rgb(0,0,0,.4), xlim=xlim, ylim=ylim)
mtext(round(cor(df$v1[df$p=='bci'], df$v3[df$p=='bci']), 2), 3, -2, col=4)
abline(v=unlist(tapply(df$v1, df$p, range)), col=rep(1:3, each=2))
abline(h=unlist(tapply(df$v3, df$p, range)), col=rep(1:3, each=2))

plot(df$v1[df$p=='sherman'], df$v3[df$p=='sherman'], xlab='nci', ylab='avg dbh', pch=16, col=rgb(0,0,0,.4), xlim=xlim, ylim=ylim)
mtext(round(cor(df$v1[df$p=='sherman'], df$v3[df$p=='sherman']), 2), 3, -2, col=4)
abline(v=unlist(tapply(df$v1, df$p, range)), col=rep(1:3, each=2))
abline(h=unlist(tapply(df$v3, df$p, range)), col=rep(1:3, each=2))




d <- tdata[tdata$plot=='sherman',]
wsg <- tapply(d$WSG, d$spcode, mean)
nci <- tapply(log(d$nci), d$spcode, mean)
d <- droplevels(d)
plot(d$WSG, log(d$nci))
points(wsg, nci, pch=21, bg=2)
cor(d$WSG, log(d$nci), use='p')
cor(wsg, nci, use='p')






d <- tdata[tdata$Growth.Include.3 == T,]
d <- d[!duplicated(d$id),]

dim(d)
table(is.na(d$WSG), d$plot)[2,] / colSums(table(is.na(d$WSG), d$plot))

table(is.na(d$LMALEAF_AVI), d$plot)[2,] / colSums(table(is.na(d$LMALEAF_AVI), d$plot))


d <- tdata[tdata$plot=='sherman',]

x <- table(d$spcode[is.na(d$WSG)])
sort(x[x>0])

x <- table(d$spcode[is.na(d$LMALEAF_AVI)])
sort(x[x>0])




unique(d$latin[d$spcode=='UNONPA'])
unique(d$latin[d$spcode=='TOVOLO'])
unique(d$latin[d$spcode=='TOVOST'])
unique(d$latin[d$spcode=='MANIBI'])



head(tdata[tdata$spcode=='MANIBI',])



sla <- read.csv('/Users/Bob/Desktop/sla_ordonez.csv')

head(sla)

sla$lma <- 1/(as.numeric(as.character(sla$SLA))) * 10000
sla$lma <- ifelse(sla$SLA=='NT', -9999, sla$lma)

sla <- sla[!is.na(sla$lma),]

unique(tdata$latin[is.na(tdata$LMALEAF_AVI)]) %in% as.character(sla$Species)

she <- tdata[tdata$plot=='sherman',]
#d <- sort(table(tdata$latin[is.na(tdata$LMALEAF_AVI)]))
d <- sort(table(she$latin[is.na(she$LMALEAF_AVI)]))
sum(d[names(d) %in% as.character(sla$Species[sla$SLA!=''])])

d[names(d) %in% as.character(sla$Species[sla$lma > 0])]


r <- (sla$Ref[as.character(sla$Species) %in% names(d)])

r <- sort(unique(trimws(unlist(strsplit(as.character(r), ';')))))

r

