##########################################
## GET SUMMARY STATS FOR PANAMA ANALYSIS
##########################################

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")

load("Panama_AnalysisData_6.14.16.RDA")
load("panama_ITV_traits_6.7.16.RDA")

cutoff <- 100
tdata$sc <- ifelse(tdata$dbh <= cutoff, 1, 2)

d <- tdata

gdat <- d[!is.na(d$growth) & d$Not.Edge==1,]
sd5 <- tapply(gdat$growth, paste(gdat$plot, gdat$sc), sd, na.rm=T)
sd5 <- sd5[match(paste(gdat$plot, gdat$sc), names(sd5))]
gdat <- gdat[gdat$Growth.Include.3==T & abs(gdat$growth) < sd5,]

# GROWTH
mean.growth.1 <- round(tapply(gdat$growth[gdat$sc==1], gdat$plot[gdat$sc==1], mean, na.rm=T), 2)
sd.growth.1 <- round(tapply(gdat$growth[gdat$sc==1], gdat$plot[gdat$sc==1], sd, na.rm=T), 2)
mean.growth.2 <- round(tapply(gdat$growth[gdat$sc==2], gdat$plot[gdat$sc==2], mean, na.rm=T), 2)
sd.growth.2 <- round(tapply(gdat$growth[gdat$sc==2], gdat$plot[gdat$sc==2], sd, na.rm=T), 2)
mean.growth.3 <- round(tapply(gdat$growth, gdat$plot, mean, na.rm=T), 2)
sd.growth.3 <- round(tapply(gdat$growth, gdat$plot, sd, na.rm=T), 2)

# SIZES
mean.dbh.1 <- round(tapply(d$dbh[d$sc==1], d$plot[d$sc==1], mean, na.rm=T), 2)/10
sd.dbh.1 <- round(tapply(d$dbh[d$sc==1], d$plot[d$sc==1], sd, na.rm=T), 2)/10
mean.dbh.2 <- round(tapply(d$dbh[d$sc==2], d$plot[d$sc==2], mean, na.rm=T), 2)/10
sd.dbh.2 <- round(tapply(d$dbh[d$sc==2], d$plot[d$sc==2], sd, na.rm=T), 2)/10
mean.dbh.3 <- round(tapply(d$dbh, d$plot, mean, na.rm=T), 2)/10
sd.dbh.3 <- round(tapply(d$dbh, d$plot, sd, na.rm=T), 2)/10

# NCI
mean.nci.1 <- round(tapply(log(d$All.NCI)[d$sc==1], d$plot[d$sc==1], mean, na.rm=T), 2)
sd.nci.1 <- round(tapply(log(d$All.NCI)[d$sc==1], d$plot[d$sc==1], sd, na.rm=T), 2)
mean.nci.2 <- round(tapply(log(d$All.NCI)[d$sc==2], d$plot[d$sc==2], mean, na.rm=T), 2)
sd.nci.2 <- round(tapply(log(d$All.NCI)[d$sc==2], d$plot[d$sc==2], sd, na.rm=T), 2)
mean.nci.3 <- round(tapply(log(d$All.NCI), d$plot, mean, na.rm=T), 2)
sd.nci.3 <- round(tapply(log(d$All.NCI), d$plot, sd, na.rm=T), 2)

# STEMS AND DEATHS
total.stems.1 <- round(tapply(d$survival[d$sc==1] > (-1), d$plot[d$sc==1], sum, na.rm=T), 0)
dead.1 <- round(tapply(!d$survival[d$sc==1], d$plot[d$sc==1], sum, na.rm=T), 2)
prop.dead.1 <- round(dead.1/total.stems.1, 4) * 100

total.stems.2 <- round(tapply(d$survival[d$sc==2] > (-1), d$plot[d$sc==2], sum, na.rm=T), 0)
dead.2 <- round(tapply(!d$survival[d$sc==2], d$plot[d$sc==2], sum, na.rm=T), 2)
prop.dead.2 <- round(dead.2/total.stems.2, 4) * 100

total.stems.3 <- round(tapply(d$survival > (-1), d$plot, sum, na.rm=T), 0)
dead.3 <- round(tapply(!d$survival, d$plot, sum, na.rm=T), 2)
prop.dead.3 <- round(dead.3/total.stems.3, 4) * 100

# MEAN TRAITS
head(d)
head(traits)

d$WSG <- traits$WD.mean[match(d$spplot, traits$sp)]
d$LMA <- traits$LMA.mean[match(d$spplot, traits$sp)]

wsg.1 <- round(tapply(d$WSG[d$sc==1], d$plot[d$sc==1], mean, na.rm=T), 2)
wsg.2 <- round(tapply(d$WSG[d$sc==2], d$plot[d$sc==2], mean, na.rm=T), 2)
wsg.3 <- round(tapply(d$WSG, d$plot, mean, na.rm=T), 2)

lma.1 <- round(tapply(d$LMA[d$sc==1], d$plot[d$sc==1], mean, na.rm=T), 2)
lma.2 <- round(tapply(d$LMA[d$sc==2], d$plot[d$sc==2], mean, na.rm=T), 2)
lma.3 <- round(tapply(d$LMA, d$plot, mean, na.rm=T), 2)

tab <- rbind(total.stems.1, total.stems.2, total.stems.3,
      prop.dead.1, prop.dead.2, prop.dead.3,
      mean.dbh.1, sd.dbh.1,
      mean.dbh.2, sd.dbh.2,
      mean.dbh.3, sd.dbh.3,
      mean.growth.1, sd.growth.1,
      mean.growth.2, sd.growth.2,
      mean.growth.3, sd.growth.3,
      mean.nci.1, sd.nci.1,
      mean.nci.2, sd.nci.2,
      mean.nci.3, sd.nci.3,
      wsg.1, wsg.2, wsg.3,
      lma.1, lma.2, lma.3)


write.csv(tab, '/Users/Bob/Projects/Postdoc/Panama/RESULTS/summary.table.csv')



