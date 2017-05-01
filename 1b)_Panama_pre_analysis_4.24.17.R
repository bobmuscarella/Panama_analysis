library(sp)
library(reshape)

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")

# Load wood density data for AGB growth metric (this data is derived below...)
load("panama_ITV_traits_6.7.16.RDA") # traits
get.agb <- function(wd, dbh){
  wd * (exp(-1.499 + (2.1481 * log(dbh)) + (0.207 * log(dbh^2)) - (0.0281 * log(dbh^3))))
}
spcodes <- read.csv('bci/spcodes.csv')
############################################
### CHUNK 1 : BUILD PANAMA DATA FROM RAW ####
############################################

### BCI ###
# bci4 <- read.csv("bci/BCI_census4.csv")
# bci5 <- read.csv("bci/BCI_census5.csv")
# save(bci4, file='bci/bci4.RDA')
# save(bci5, file='bci/bci5.RDA')
load('bci/bci4.RDA')
load('bci/bci5.RDA')

bci4$Date <- as.Date(bci4$Date, format='%m/%d/%y')
bci5$Date <- as.Date(bci5$Date, format='%m/%d/%y')

# subset to main stems only
bci4 <- bci4[bci4$Stem == "main",]
bci5 <- bci5[bci5$Stem == "main",]

# subset to alive stems in C4 with DBH measured and alive or dead stems in C5
bci4 <- bci4[ bci4$Status %in% c('alive'),]
#bci4 <- bci4[ ! is.na(bci4$DBH),]
bci4 <- droplevels(bci4)
bci5 <- bci5[ ! bci5$Status %in% c('missing'),]

# SET BORDER TO x METERS
r <- 15
bci4$Not.edge <- ifelse(bci4$PX > r & bci4$PX < (1000-r) & bci4$PY > r & bci4$PY < (500-r), 1, 0)
bci5$Not.edge <- ifelse(bci5$PX > r & bci5$PX < (1000-r) & bci5$PY > r & bci5$PY < (500-r), 1, 0)

bci4$plot <- 'bci'
bci5$plot <- 'bci'

bci4$census <- 4
bci5$census <- 5

bci4 <- bci4[,c('Tag','Latin','PX','PY','Not.edge','plot', 'census', 'DBH', 'HOM','Date','Status')]
bci5 <- bci5[,c('Tag','Latin','PX','PY','Not.edge','plot', 'census', 'DBH', 'HOM','Date','Status')]

bci4$match <- paste(bci4$Tag, bci4$Latin, bci4$HOM)
bci5$match <- paste(bci5$Tag, bci5$Latin, bci5$HOM)
bci4$matchdeads <- paste(bci4$Tag, bci4$Latin)
bci5$matchdeads <- paste(bci5$Tag, bci5$Latin)

bci4$DBH2 <- bci5$DBH[match(bci4$match, bci5$match)]
bci4$Status2 <- bci5$Status[match(bci4$match, bci5$match)]
bci4$Date2 <- bci5$Date[match(bci4$match, bci5$match)]

bci4$DBH3 <- bci5$DBH[match(bci4$matchdeads, bci5$matchdeads)]
bci4$Status3 <- bci5$Status[match(bci4$matchdeads, bci5$matchdeads)]
bci4$DateFinal <- bci5$Date[match(bci4$matchdeads, bci5$matchdeads)]

bci4$StatusStart <- bci4$Status
bci4$StatusFinal <- ifelse(!is.na(bci4$DBH2), bci4$Status2, bci4$Status3)
bci4$StatusFinal <- ifelse(bci4$StatusFinal == 1, 'alive', ifelse(bci4$StatusFinal == 2, 'dead', NA))

bci4$dbhStart <- bci4$DBH
bci4$dbhFinal <- bci4$DBH2

# Remove cases where ambiguous DBHs were assigned (emerging from multi HOM strangeness)
bci4 <- bci4[bci4$DBH2==bci4$DBH3 | is.na(bci4$DBH2) | is.na(bci4$DBH3),]

# Add wood density to calculate biomass growth
bci4$wd <- traits$WD.mean[match(paste0('2.',bci4$Latin), traits$sp)]
bci4$agbStart <- get.agb(bci4$wd, bci4$dbhStart)
bci4$agbFinal <- get.agb(bci4$wd, bci4$dbhFinal)

bci4$survival <- ifelse(bci4$StatusFinal == 'alive', 1, 0)
bci4$days <- as.numeric(bci4$DateFinal - bci4$Date)
bci4$growth <- ((bci4$dbhFinal - bci4$dbhStart)/bci4$days) * 365
bci4$log.growth <- ((log(bci4$dbhFinal) - log(bci4$dbhStart))/(bci4$days)) * 365
bci4$RGR <- log(bci4$dbhFinal/bci4$dbhStart)/(bci4$days/365)

bci4$AGB.growth <- ((bci4$agbFinal - bci4$agbStart)/bci4$days) * 365

bci <- bci4

# NOTE THAT I CHANGED 3 SPECIES CODES FROM THE ONLINE BCI SPCODE DATA 
# TO BE CONSISTENT WITH THE BCI TRAIT DATA...
spcodes <- read.csv("bci/spcodes.csv")
bci$spcode <- spcodes$spcode[match(bci$Latin, paste(spcodes$Genus, spcodes$SpeciesName))]

bci <- droplevels(bci)

bci <- bci[,c('Tag','Latin','spcode','PX','PY','Not.edge','plot', 'census', 
              'dbhStart', 'growth','log.growth','RGR','survival','days','AGB.growth')]
names(bci) <- c('tag','latin','spcode','x','y','Not.Edge','plot','census','dbh',
                'growth','log.growth','RGR','survival','days','AGB.growth')
bci$id <- as.numeric(rownames(bci))

bci <- bci[!is.na(bci$survival),]
bci <- bci[!is.na(bci$x),]
bci <- bci[!is.na(bci$y),]
bci$tag <- NULL

# save(bci, file='bci/bci_preNCI.RDA')
# save(bci, file='bci/bci_preNCI_10.23.15.RDA')
# save(bci, file='bci/bci_preNCI_12.9.15.RDA')
# save(bci, file='bci/bci_preNCI_6.1.16.RDA')

################################
### ADD OTHER PANAMA PLOTS ###
################################
# load(file='bci/bci_preNCI_6.1.16.RDA')
head(bci)

coc <- read.csv("cocoli/cocoli.csv")
cocsp <- read.table("cocoli/cocolisp.txt", header=T)

she <- read.csv("sherman/sherman.csv")
shesp <- read.table("sherman/shermansp.txt", header=T)

warning(paste('setting edge to', r, 'meters'))
# plot(she$x, she$y, xlim=c(0,450), ylim=c(0,450), pch=16, cex=.2)
she$Not.Edge <- point.in.polygon(she$x, she$y, c(r,(140-r),(140-r),(240-r),(240-r),(140+r),(140+r),r), 
                                 c(r,r,(40+r),(40+r),(340-r),(340-r),(140-r),(140-r)))
she$Not.Edge <- ifelse(she$Not.Edge > 0, 1, 0)
# points(she$x, she$y, col= she$Not.Edge +1, pch=16, cex=.2)
# sort(table(she$spcode[she$Not.Edge ==1]))

# new <- point.in.polygon(she$x, she$y, c(140,240,240,140), c(340,340,440,440))
# col <- ifelse(new==1, 3, she$Not.Edge+1)
# points(she$x, she$y, col= col, pch=16, cex=.2)

# plot(coc$x, coc$y, xlim=c(0,200), ylim=c(0,300), pch=16, cex=.5)
coc$Not.Edge <- point.in.polygon(coc$x, coc$y, c(r,(200-r),(200-r),(100-r),(100-r),r), 
                                 c(r,r,(100-r),(100-r),(300-r),(300-r)))
coc$Not.Edge <- ifelse(coc$Not.Edge > 0, 1, 0)
# points(coc$x, coc$y, col= coc$Not.Edge +1, pch=16, cex=.5)
# sum(table(coc$spcode[coc$Not.Edge==1]))

she$plot <- "sherman"
coc$plot <- "cocoli"

coc <- coc[!duplicated(coc$tag),] # remove a duplicate tag
she <- she[!duplicated(she$tag),] # remove a duplicate tag

tmp <- rbind(coc, she)
head(tmp)

Date.Convert <- function(x){
	x <- as.Date(x, format='%m/%d/%y')
	julian(x)
	}
tmp$date1 <- Date.Convert(tmp$date1)
tmp$date2 <- Date.Convert(tmp$date2)
tmp$date3 <- Date.Convert(tmp$date3)

head(tmp)

# GROWTH
tmp$dbh1 <- ifelse(tmp$dbh1 < 0, NA, tmp$dbh1)
tmp$dbh2 <- ifelse(tmp$dbh2 < 0, NA, tmp$dbh2)
tmp$dbh3 <- ifelse(tmp$dbh3 < 0, NA, tmp$dbh3)

# agb growth
tmp$spbinom <- paste(ifelse(tmp$plot=='cocoli',1,3), paste(spcodes$Genus, spcodes$SpeciesName)[match(tmp$spcode, spcodes$spcode)], sep='.')
tmp$wd <- traits$WD.mean[match(tmp$spbinom, traits$sp)]
tmp$agb1 <- ifelse(tmp$dbh1 < 0, NA, get.agb(tmp$wd, tmp$dbh1))
tmp$agb2 <- ifelse(tmp$dbh2 < 0, NA, get.agb(tmp$wd, tmp$dbh2))
tmp$agb3 <- ifelse(tmp$dbh3 < 0, NA, get.agb(tmp$wd, tmp$dbh3))

tmp$grow12 <- ((tmp$dbh2 - tmp$dbh1) / (tmp$date2 - tmp$date1))* 365
tmp$log.grow12 <- ((log(tmp$dbh2) - log(tmp$dbh1)) / (tmp$date2 - tmp$date1))* 365
tmp$RGR12 <- log(tmp$dbh2/tmp$dbh1) / ((tmp$date2 - tmp$date1) / 365)
tmp$agb.grow12 <- ((tmp$agb2 - tmp$agb1) / (tmp$date2 - tmp$date1))* 365

valid1 <- ifelse(tmp$recr1=='A' & as.character(tmp$recr1) == as.character(tmp$recr2) & tmp$pom1 == tmp$pom2, 1, 0)
tmp$grow12 <- ifelse(valid1==1, tmp$grow12, NA)
tmp$log.grow12 <- ifelse(valid1==1, tmp$log.grow12, NA)
tmp$RGR12 <- ifelse(valid1==1, tmp$RGR12, NA)
tmp$agb.grow12 <- ifelse(valid1==1, tmp$agb.grow12, NA)

tmp$grow23 <- ((tmp$dbh3 - tmp$dbh2) / (tmp$date3 - tmp$date2))* 365
tmp$log.grow23 <- ((log(tmp$dbh3) - log(tmp$dbh2) / (tmp$date3 - tmp$date2)))* 365
tmp$RGR23 <- log(tmp$dbh3/tmp$dbh2) / ((tmp$date3 - tmp$date2) / 365)
tmp$agb.grow23 <- ((tmp$agb3 - tmp$agb2) / (tmp$date3 - tmp$date2))* 365

valid2 <- ifelse(tmp$recr2=='A' & as.character(tmp$recr2) == as.character(tmp$recr3) & tmp$pom2 == tmp$pom3, 1, 0)
tmp$grow23 <- ifelse(valid2==1, tmp$grow23, NA)
tmp$log.grow23 <- ifelse(valid2==1, tmp$log.grow23, NA)
tmp$RGR23 <- ifelse(valid2==1, tmp$RGR23, NA)
tmp$agb.grow23 <- ifelse(valid2==1, tmp$agb.grow23, NA)


# SURVIVAL
tmp$survive12 <- ifelse(tmp$recr1=='A' & as.character(tmp$recr1) == as.character(tmp$recr2), 1, NA)
tmp$survive12 <- ifelse(tmp$recr1=='A' & tmp$recr2=='D' , 0, tmp$survive1)

tmp$survive23 <- ifelse(tmp$recr2=='A' & as.character(tmp$recr2) == as.character(tmp$recr3), 1, NA)
tmp$survive23 <- ifelse(tmp$recr2=='A' & tmp$recr3=='D' , 0, tmp$survive2)

tmp$int12 <- tmp$date2 - tmp$date1
tmp$int23 <- tmp$date3 - tmp$date2
head(tmp)


# Clean up
tmp <- tmp[,c("tag","spcode","x","y",'dbh1','dbh2',"Not.Edge","plot",
						"grow12","grow23","log.grow12","log.grow23","RGR12","RGR23",
            "survive12","survive23","int12","int23","agb.grow12","agb.grow23")]

tmp2 <- reshape(tmp, varying=c(list(5:6), list(9:10), list(11:12), list(13:14), list(15:16), list(17:18), list(19:20)), direction='long')
tmp2 <- tmp2[order(tmp2$tag),]
names(tmp2)[7:14] <- c('census','dbh','growth','log.growth','RGR','survival','days','AGB.growth')
#tmp2 <- tmp2[,c(1:15)]
tmp2$tag <- NULL

rownames(tmp2) <- NULL
tdata <- tmp2

tdata$latin[tdata$plot=='cocoli'] <- paste(cocsp$genus, cocsp$species)[match(tdata$spcode, cocsp$spcode)][tdata$plot=='cocoli']
tdata$latin[tdata$plot=='sherman'] <- paste(shesp$genus, shesp$species)[match(tdata$spcode, shesp$spcode)][tdata$plot=='sherman']

tdata <- tdata[,c("latin","spcode","x","y","Not.Edge","plot","census","dbh",
                  "growth","log.growth","RGR","survival","days","id","AGB.growth")]

tdata <- rbind(tdata, bci)

tdata <- tdata[ ! is.na(tdata$dbh),]

tdata <- tdata[tdata$dbh != 0 , ]

# save(tdata, file='panama_preNCI.RDA')
# save(tdata, file='panama_preNCI_10.23.15.RDA')
# save(tdata, file='panama_preNCI_12.9.15.RDA')
# save(tdata, file='panama_preNCI_6.1.16.RDA')

##################################
### CHUNK 2 : IDENTIFY PALMS ####
##################################
# load('panama_preNCI_6.1.16.RDA')

bcisp <- read.csv("bci/spcodes.csv")
cocsp <- read.table("cocoli/cocolisp.txt", header=T)
shesp <- read.table("sherman/shermansp.txt", header=T)
bcisp <- bcisp[,c('spcode','Genus','SpeciesName','Family')]
names(bcisp) <- names(cocsp)

bcisp$plot <- 'bci'
cocsp$plot <- 'cocoli'
shesp$plot <- 'sherman'

spcodes <- rbind(bcisp, cocsp, shesp)
spcodes$spcodeplot <- paste(spcodes$spcode, spcodes$plot)
spcodeplot <- paste(tdata$spcode, tdata$plot)

tdata$palm <-  spcodeplot %in% spcodes$spcodeplot[spcodes$family == "Arecaceae" ]

####################################
### CHUNK 3 : ATTACH TRAIT DATA ####
####################################

#######################################################
#### THIS WAS TO USE THE PUBLIC WD DATA FROM CHAVE 2006
# wd <- read.csv("traits/ChaveWoodDensity.csv")[,c(1:4)]
# wd$binom <- paste(wd$GENUS, wd$SPECIES, sep=' ')
# tdata$chave.wd <- wd$WSG[match(tdata$latin, wd$binom)]
#######################################################
# 
# bci <- read.csv("traits/BCITRAITS_20101220.csv")
# 
# # GET WSG 100C IF AVAILABLE, OTHERWISE WSG 60C AND THEN CHAVE WSG:
# bci$WSG <- ifelse(!is.na(bci$SG100C_AVG), bci$SG100C_AVG, 
#                   ifelse(!is.na(bci$SG60C_AVG), bci$SG60C_AVG, bci$WSG_CHAVE))
# 
# # SUBSET TO TRAITS OF PRIMARY INTEREST...
# bci <- bci[,c('SP.', 'WSG','LMALEAF_AVI')]
# 
# ### MATCH SP CODES WITH TRAITS..... 
# traits <- bci[match(tdata$spcode, bci$SP.),]
# tdata <- cbind(tdata, traits[,-1])
# 
# ### TOTAL NUMBER OF INDIVIDUALS WITH TRAIT DATA
# sum(!is.na(tdata$WSG))/nrow(tdata) # WD data for ~95% of observations here...
# sum(!is.na(tdata$LMALEAF_AVI))/nrow(tdata) # WD data for ~95% of observations here...
# 
# ### NUMBER OF INDIVIDUALS PER PLOT WITH TRAIT DATA
# table(tdata$plot[!is.na(tdata$WSG)]) / table(tdata$plot)
# table(tdata$plot[!is.na(tdata$LMALEAF_AVI)]) / table(tdata$plot)
# 
# ### NUMBER OF SPECIES PER PLOT WITH TRAIT DATA
# plotsp <- paste(tdata$plot, tdata$spcode)
# tmp <- tdata[!duplicated(plotsp),]
# table(tmp$plot[!is.na(tmp$WSG)])/table(tmp$plot)
# table(tmp$plot[!is.na(tmp$LMALEAF_AVI)])/table(tmp$plot)
# 
# ### ADD TRAIT DATA FROM PR
# setwd('/Users/Bob/Projects/Thesis/DATA') 
# pr <- read.csv("traits/Fixed_traits.csv", row.names=1)
# prsp <- read.csv("census/PR_SP_LIST.csv")
# setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")
# pr <- cbind(pr, prsp[match(rownames(pr), prsp$CODE),])
# pr$name <- sub('_', ' ', pr$BINOM)
# head(pr)
# 
# # ADD DATA FROM PR TRAIT DATA
# unique(tdata$latin[is.na(tdata$WSG)])[unique(tdata$latin[is.na(tdata$WSG)]) %in% pr$name]
# unique(tdata$latin[is.na(tdata$LMALEAF_AVI)])[unique(tdata$latin[is.na(tdata$LMALEAF_AVI)]) %in% pr$name]
# #tdata$WSG[is.na(tdata$WSG)] <- pr$WD[match(tdata$latin, pr$name)][is.na(tdata$WSG)]
# #tdata$LMALEAF_AVI[is.na(tdata$LMALEAF_AVI)] <- ((1/pr$SLA.wp)*10000)[match(tdata$latin, pr$name)][is.na(tdata$LMALEAF_AVI)]
# tdata$o1.WSG <- pr$WD[match(tdata$latin, pr$name)]
# tdata$o1.LMA <- ((1/pr$SLA.wp)*10000)[match(tdata$latin, pr$name)]
# 
# # ADD DATA FROM ORDONEZ TRAIT DATA
# ord <- read.csv('traits/sla_ordonez.csv')
# ord$LMA <- (1/as.numeric(as.character(ord$SLA)) * 10000)
# ord <- ord[!is.na(ord$LMA),]
# # tdata$LMALEAF_AVI[is.na(tdata$LMALEAF_AVI)] <- ord$LMA[match(tdata$latin, ord$Species)][is.na(tdata$LMALEAF_AVI)]
# tdata$o2.LMA <- ord$LMA[match(tdata$latin, ord$Species)]
# 
# # ADD LMA DATA FROM NEOTROPICAL TRAIT DATABASE
# nt <- read.csv('traits/Neotropics traits longform.csv')
# nt$LMA <- (1/as.numeric(as.character(nt$SLAcm2g)) * 10000)
# ntlma <- nt[!is.na(nt$LMA),]
# ntpan <- ntlma[ntlma$SITE=='PANAMA',]
# # tdata$LMALEAF_AVI[is.na(tdata$LMALEAF_AVI)] <- ntpan$LMA[match(tdata$latin, ntpan$Species)][is.na(tdata$LMALEAF_AVI)]
# # ntcr <- ntlma[ntlma$SITE=='COSTARICA',]
# # tdata$LMALEAF_AVI[is.na(tdata$LMALEAF_AVI)] <- ntcr$LMA[match(tdata$latin, ntcr$Species)][is.na(tdata$LMALEAF_AVI)]
# # ntbo <- ntlma[ntlma$SITE=='BOLIVIA',]
# # tdata$LMALEAF_AVI[is.na(tdata$LMALEAF_AVI)] <- ntbo$LMA[match(tdata$latin, ntbo$Species)][is.na(tdata$LMALEAF_AVI)]
# # ntmx <- ntlma[ntlma$SITE=='MEXICO',]
# # tdata$LMALEAF_AVI[is.na(tdata$LMALEAF_AVI)] <- ntmx$LMA[match(tdata$latin, ntmx$Species)][is.na(tdata$LMALEAF_AVI)]
# tdata$o3.LMA <- ntpan$LMA[match(tdata$latin, ntpan$Species)]
# ntcr <- ntlma[ntlma$SITE=='COSTARICA',]
# tdata$o3.LMA <- ntcr$LMA[match(tdata$latin, ntcr$Species)]
# ntbo <- ntlma[ntlma$SITE=='BOLIVIA',]
# tdata$o3.LMA <- ntbo$LMA[match(tdata$latin, ntbo$Species)]
# ntmx <- ntlma[ntlma$SITE=='MEXICO',]
# tdata$o3.LMA <- ntmx$LMA[match(tdata$latin, ntmx$Species)]
# 
# 
# # ADD WSG DATA FROM NEOTROPICAL TRAIT DATABASE
# nt <- read.csv('traits/Neotropics traits longform.csv')
# ntwd <- nt[!is.na(nt$WDkgm3),]
# ntwd$WSG <- ntwd$WDkgm3/1000
# for(i in 1:6){
#   site <- c('PANAMA','COSTARICA','BOLIVIA','MEXICO','ECUADOR','PERU')[i]  
#   tmp <- ntwd[ntwd$SITE %in% site,]
# #  tdata$WSG[is.na(tdata$WSG)] <- tmp$WSG[match(tdata$latin, tmp$Species)][is.na(tdata$WSG)]  
#   tdata$o2.WSG <- tmp$WSG[match(tdata$latin, tmp$Species)]
# }

#############
### MAKE DATA REQUEST FORM FOR TRY
# TITLE FOR TRY REQUEST
# Trait-mediated neighborhood interactions across a regional precipitation gradient

### DESCRIPTION FOR TRY REQUEST
# To better understand how the influence of local biotic interactions varies with respect to abiotic conditions, we are building trait-based neighborhood models of individual tree growth and survival across a regional gradient encompassing variation in precipitation, soil conditions, and pest pressure. We are quantifying (i) how average performance and the effect of crowding vary across the gradient and (ii) how two traits (wood density and leaf mass per area [LMA]) mediate these processes.  We are making this request to fill gaps in our existing trait dataset.

# sp <- read.table("traits/TRY_species_list.txt", sep='\t')
# names(sp) <- as.matrix(sp[1,])
# sp <- sp[-1,]
# 
# needlma <- unique(tdata$latin[is.na(tdata$LMALEAF_AVI)])
# lmaspnumb <- sp$AccSpeciesID[sp$AccSpeciesName %in% needlma]
# 
# tlist <- read.table("traits/TRY_trait_table.txt", sep='\t')
# names(tlist) <- as.matrix(tlist[1,])
# tlist <- tlist[-1,]
# 
# slanumb <- as.numeric(as.character(tlist$TraitID[which(tlist$Trait == 'Leaf area per leaf dry mass (specific leaf area, SLA)')]))
# lmaspnumb <- paste(lmaspnumb, collapse=', ')
# lmaspnumb
# slanumb
# 
# wsgnumb <- as.numeric(as.character(tlist$TraitID[which(tlist$Trait == 'Stem dry mass per stem fresh volume (stem specific density, SSD, wood density)')]))
# needwsg <- unique(tdata$latin[is.na(tdata$WSG)])
# wsgspnumb <- sp$AccSpeciesID[sp$AccSpeciesName %in% needwsg]
# wsgspnumb <- paste(wsgspnumb, collapse=', ')
# wsgspnumb
# wsgnumb
#############
# 
# 
# #### ADD DATA FROM TRY
# sla <- read.csv('traits/TRY_SLA_2069/SLA_2069_reduced.csv')
# sla$lma <- 1/(sla$StdValue*0.001)
# lma <- tapply(sla$lma, sla$AccSpeciesName, mean)
# #tdata$LMALEAF_AVI[is.na(tdata$LMALEAF_AVI)] <- lma[match(tdata$latin, names(lma))][is.na(tdata$LMALEAF_AVI)]
# tdata$o4.LMA <- lma[match(tdata$latin, names(lma))]
# 
# wd <- read.csv('traits/TRY_WSG_2071/WSG_2071_reduced.csv')
# wd <- tapply(wd$StdValue, wd$AccSpeciesName, mean)
# #tdata$WSG[is.na(tdata$WSG)] <- wd[match(tdata$latin, names(wd))][is.na(tdata$WSG)]
# tdata$o3.WSG <- wd[match(tdata$latin, names(wd))]
# 
# ### Log transform skewed traits *before* calculating tNCI...
# #tdata[,'log.LMA'] <- log(tdata[,'LMALEAF_AVI'])
# 
# #save(tdata, file='panama_traits_preNCI_12.9.15.RDA')
# #save(tdata, file='panama_traits_preNCI_5.20.15.RDA')
# 

##############################
### THIS IS THE NEW WAY TO ADD TRAIT DATA TO INCORPORATE INTRASPECIFIC VARIABILITY...
##############################
tdata$plot <- ifelse(tdata$plot=='bci', 2, ifelse(tdata$plot=='cocoli', 1, 3))
tdata$spplot <- paste(tdata$plot, tdata$latin, sep='.')
splist <- data.frame(sp=sort(unique(tdata$spplot)), stringsAsFactors=F)
splist$plot <- substring(splist$sp,1,1)
splist$latin <- substring(splist$sp,3,nchar(splist$sp))

### ADD INDIVIDUAL-LEVEL LMA DATA FROM BCI DATABASE
lma <- read.csv('traits/TRY_LVS_20101127.csv')
nom <- read.csv('traits/RAW_BCI_WSG/nomenclature_R_20120305_Rready.csv')
nom$spplot <- paste('2.', nom$genus, ' ', nom$species, sep='')
lma <- lma[lma$MICROHABIT.=='SHADE',]
lma$LMA <- lma$LMA_DRY*10000
lma$spplot <- nom$spplot[match(lma$SP., nom$sp6)]
lma.mean <- tapply(lma$LMA, lma$spplot, mean, na.rm=T)
lma.sd <- tapply(lma$LMA, lma$spplot, sd, na.rm=T)
lma.n <- tapply(lma$LMA, lma$spplot, length)
log.lma.sd <- tapply(log(lma$LMA), lma$spplot, sd, na.rm=T)
splist$bci.LMA.mean <- lma.mean[match(splist$sp, names(lma.mean))]
splist$bci.LMA.sd <- lma.sd[match(splist$sp, names(lma.sd))]
splist$bci.log.LMA.sd <- log.lma.sd[match(splist$sp, names(log.lma.sd))]
splist$bci.n.LMA <- lma.n[match(splist$sp, names(lma.n))]

### Give BCI values across the gradient (to use if no other values are available...)
gen.lma.mean <- rep(lma.mean, times=3)
names(gen.lma.mean) <- paste(rep(1:3, each=length(lma.mean)), substring(names(gen.lma.mean), 2, nchar(names(gen.lma.mean))), sep='')
gen.lma.sd <- rep(lma.sd, times=3)
names(gen.lma.sd) <- names(gen.lma.mean)
gen.log.lma.sd <- rep(log.lma.sd, times=3)
names(gen.log.lma.sd) <- names(gen.lma.mean)
splist$bcigen.LMA.mean <- gen.lma.mean[match(splist$sp, names(gen.lma.mean))]
splist$bcigen.LMA.sd <- gen.lma.sd[match(splist$sp, names(gen.lma.sd))]
splist$bcigen.log.LMA.sd <- gen.log.lma.sd[match(splist$sp, names(gen.log.lma.sd))]

### ADD LMA DATA FROM MESSIER
mess <- read.csv("traits/Messier/lft - sharing - BobMuscarella.csv")
mess$plot <- ifelse(mess$Site == 'MET', 1, ifelse(mess$Site == 'BCI', 2, 3))
mess$spplot <- paste(mess$plot, mess$Species, sep='.')
mess.mean <- tapply(mess$LMA[mess$Strata==0], mess$spplot[mess$Strata==0], mean, na.rm=T)
mess.sd <- tapply(mess$LMA[mess$Strata==0], mess$spplot[mess$Strata==0], sd, na.rm=T)
mess.n <- tapply(mess$LMA[mess$Strata==0], mess$spplot[mess$Strata==0], length)
mess.log.sd <- tapply(log(mess$LMA)[mess$Strata==0], mess$spplot[mess$Strata==0], sd, na.rm=T)
splist$mess.LMA.mean <- mess.mean[match(splist$sp, names(mess.mean))]
splist$mess.LMA.sd <- mess.sd[match(splist$sp, names(mess.sd))]
splist$mess.log.LMA.sd <- mess.log.sd[match(splist$sp, names(mess.log.sd))]
splist$mess.n.LMA <- mess.n[match(splist$sp, names(mess.n))]

### ADD LMA DATA FROM CANOPY CRANE
crane <- read.csv('traits/SLA crane species.csv')
crane$plot <- ifelse(crane$SITE. == 'FTS', 3, ifelse(crane$SITE. == 'PNM', 1, NA))
crane <- crane[crane$STRATA.=='UNDERSTORY',]
crane$LMA <- (1/crane$SLA_LEAF)*10000
crane$spplot <- paste(crane$plot, '.', crane$GENUS., ' ', crane$SPECIES., sep='')
splist$crane.LMA.mean <- crane$LMA[match(splist$sp, crane$spplot)]

### ADD LMA DATA FROM PR
setwd('/Users/Bob/Projects/Thesis/DATA') 
pr <- read.csv("traits/Fixed_traits.csv", row.names=1)
prsp <- read.csv("census/PR_SP_LIST.csv")
setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")
pr <- cbind(pr, prsp[match(rownames(pr), prsp$CODE),])
pr$name <- sub('_', ' ', pr$BINOM)
prfull <- rbind(pr,pr,pr)
prfull$plot <- rep(1:3, each=nrow(pr))
prfull$spplot <- paste(prfull$plot, '.', prfull$name, sep='')
splist$pr.LMA.mean <- ((1/prfull$SLA.wp)*10000)[match(splist$sp, prfull$spplot)]

### ADD DATA FROM ORDONEZ TRAIT DATA
ord <- read.csv('traits/sla_ordonez.csv')
ord$LMA <- (1/as.numeric(as.character(ord$SLA)) * 10000)
ord <- ord[!is.na(ord$LMA),]
ordfull <- rbind(ord,ord,ord)
ordfull$plot <- rep(1:3, each=nrow(ord))
ordfull$spplot <- paste(ordfull$plot, '.', ordfull$Species, sep='')
splist$ord.LMA.mean <- ordfull$LMA[match(splist$sp, ordfull$spplot)]

### ADD LMA DATA FROM NEOTROPICAL TRAIT DATABASE
nt <- read.csv('traits/Neotropics traits longform.csv')
nt$LMA <- (1/as.numeric(as.character(nt$SLAcm2g)) * 10000)
ntlma <- nt[!is.na(nt$LMA),]
splist$nt.LMA.mean <- NA
for(i in 1:4){
  site <- c('PANAMA','COSTARICA','BOLIVIA','MEXICO')[i]  
  tmp <- ntlma[ntlma$SITE %in% site,]
  ntlmafull <- rbind(tmp,tmp,tmp)
  ntlmafull$plot <- rep(1:3, each=nrow(tmp))
  ntlmafull$spplot <- paste(ntlmafull$plot, '.', ntlmafull$Species, sep='')
  splist$nt.LMA.mean[is.na(splist$nt.LMA)] <- ntlmafull$LMA[match(splist$sp, ntlmafull$spplot)][is.na(splist$nt.LMA)]
}

#### ADD LMA DATA FROM TRY
sla <- read.csv('traits/TRY_SLA_2069/SLA_2069_reduced.csv')
sla$lma <- 1/(sla$StdValue*0.001)
lma <- tapply(sla$lma, sla$AccSpeciesName, mean)
trylma <- c(lma,lma,lma)
names(trylma) <- paste(rep(1:3, each=length(lma)), names(lma), sep='.')
splist$try.LMA.mean <- trylma[match(splist$sp, names(trylma))]

### ADD INDIVIDUAL-LEVEL WOOD DENSITY DATA FROM BCI DATABASE
wdnom <- read.csv('traits/RAW_BCI_WSG/nomenclature_R_20120305_Rready.csv', 1)
wd <- read.csv('traits/RAW_BCI_WSG/WSG_ForIndividualPlants_tempstep9b.csv', 1)
wdnom$spplot <- paste('2.', wdnom$genus, ' ', wdnom$species, sep='')
wd$spplot <- wdnom$spplot[match(wd$SP, wdnom$sp6)]
wd <- wd[!is.na(wd$spplot),]
wdmean <- tapply(wd$SG100C, wd$spplot, mean, na.rm=T)
wdsd <- tapply(wd$SG100C, wd$spplot, sd, na.rm=T)
wdn <- tapply(wd$SG100C, wd$spplot, length)
splist$bci.WD.mean <- wdmean[match(splist$sp, names(wdmean))]
splist$bci.WD.sd <- wdsd[match(splist$sp, names(wdsd))]
wd60mean <- tapply(wd$SG60C, wd$spplot, mean, na.rm=T)
wd60sd <- tapply(wd$SG60C, wd$spplot, sd, na.rm=T)
splist$bci.WD.mean[is.na(splist$bci.WD.mean)] <- wd60mean[match(splist$sp, names(wd60mean))][is.na(splist$bci.WD.mean)]
splist$bci.WD.sd[is.na(splist$bci.WD.sd)] <- wd60sd[match(splist$sp, names(wd60sd))][is.na(splist$bci.WD.sd)]
splist$bci.WD.n <- wdn[match(splist$sp, names(wdn))]

### Give BCI values across the gradient (to use if no other values are available...)
gen.wd.mean <- rep(wdmean, times=3)
names(gen.wd.mean) <- paste(rep(1:3, each=length(wdmean)), substring(names(gen.wd.mean), 2, nchar(names(gen.wd.mean))), sep='')
gen.wd.sd <- rep(wdsd, times=3)
names(gen.wd.sd) <- paste(rep(1:3, each=length(wdsd)), substring(names(gen.wd.sd), 2, nchar(names(gen.wd.sd))), sep='')
splist$bcigen.WD.mean <- gen.wd.mean[match(splist$sp, names(gen.wd.mean))]
splist$bcigen.WD.sd <- gen.wd.sd[match(splist$sp, names(gen.wd.sd))]

### ADD WD COLUMNS FOR TRAITS FROM OTHER LOCATIONS...
bci <- read.csv("traits/BCITRAITS_20101220.csv")
# get mean WSG 100C IF AVAILABLE, OTHERWISE WSG 60C AND THEN CHAVE WSG:
bci$WSG <- ifelse(!is.na(bci$SG100C_AVG), bci$SG100C_AVG, 
                  ifelse(!is.na(bci$SG60C_AVG), bci$SG60C_AVG, bci$WSG_CHAVE))
bci <- rbind(bci, bci, bci)
bci$plot <- rep(1:3, each=nrow(bci)/3)
bci$spplot <- paste(bci$plot, '.', bci$GENUS., ' ', bci$SPECIES, sep='')
splist$chave.WD.mean <- bci$WSG[match(splist$sp, bci$spplot)]

### ADD TRAIT DATA FROM PR
splist$pr.WD.mean <- prfull$WD[match(splist$sp, prfull$spplot)]

### ADD WSG DATA FROM NEOTROPICAL TRAIT DATABASE
nt <- read.csv('traits/Neotropics traits longform.csv')
ntwd <- nt[!is.na(nt$WDkgm3),]
ntwd$WSG <- ntwd$WDkgm3/1000
splist$nt.WD.mean <- NA
for(i in 1:6){
  site <- c('PANAMA','COSTARICA','BOLIVIA','MEXICO','ECUADOR','PERU')[i]  
  tmp <- ntwd[ntwd$SITE %in% site,]
  ntwdfull <- rbind(tmp,tmp,tmp)
  ntwdfull$plot <- rep(1:3, each=nrow(tmp))
  ntwdfull$spplot <- paste(ntwdfull$plot, '.', ntwdfull$Species, sep='')
  splist$nt.WD.mean[is.na(splist$nt.WD.mean)] <- ntwdfull$WSG[match(splist$sp, ntwdfull$spplot)][is.na(splist$nt.WD.mean)]
}

#### ADD WD DATA FROM TRY
wd <- read.csv('traits/TRY_WSG_2071/WSG_2071_reduced.csv')
wd <- tapply(wd$StdValue, wd$AccSpeciesName, mean)
trywd <- c(wd,wd,wd)
names(trywd) <- paste(rep(1:3, each=length(wd)), names(wd), sep='.')
splist$try.WD.mean <- trywd[match(splist$sp, names(trywd))]

### CHOOSE WHICH LMA VALUE TO USE
splist$LMA.mean <- ifelse(!is.na(splist$mess.LMA.mean), splist$mess.LMA.mean, 
                          ifelse(!is.na(splist$bci.LMA.mean), splist$bci.LMA.mean, 
                                 ifelse(!is.na(splist$bcigen.LMA.mean), splist$bcigen.LMA.mean, 
                                        ifelse(!is.na(splist$crane.LMA.mean), splist$crane.LMA.mean, 
                                               ifelse(!is.na(splist$nt.LMA.mean), splist$nt.LMA.mean, 
                                                      ifelse(!is.na(splist$ord.LMA.mean), splist$ord.LMA.mean, 
                                                             ifelse(!is.na(splist$try.LMA.mean), splist$try.LMA.mean, 
                                                                    ifelse(!is.na(splist$pr.LMA.mean), splist$pr.LMA.mean, NA))))))))

splist$LMA.sd <- ifelse(!is.na(splist$mess.LMA.sd), splist$mess.LMA.sd, 
                        ifelse(!is.na(splist$bci.LMA.sd), splist$bci.LMA.sd, 
                               ifelse(!is.na(splist$bci.LMA.sd), splist$bcigen.LMA.sd, NA)))

splist$log.LMA.sd <- ifelse(!is.na(splist$mess.log.LMA.sd), splist$mess.log.LMA.sd, 
                        ifelse(!is.na(splist$bci.log.LMA.sd), splist$bci.log.LMA.sd, 
                               ifelse(!is.na(splist$bcigen.log.LMA.sd), splist$bcigen.log.LMA.sd, NA)))

splist$LMA.source <- ifelse(!is.na(splist$mess.LMA.mean), 'MESS', 
                            ifelse(!is.na(splist$bci.LMA.mean), 'BCI', 
                                   ifelse(!is.na(splist$bcigen.LMA.mean), 'BCIGEN', 
                                          ifelse(!is.na(splist$crane.LMA.mean), 'CRANE', 
                                                 ifelse(!is.na(splist$nt.LMA.mean), 'NT', 
                                                        ifelse(!is.na(splist$ord.LMA.mean), 'ORD', 
                                                               ifelse(!is.na(splist$try.LMA.mean), 'TRY', 
                                                                      ifelse(!is.na(splist$pr.LMA.mean), 'PR', NA))))))))

### CHOOSE WHICH WD VALUE TO USE
splist$WD.mean <- ifelse(!is.na(splist$bci.WD.mean), splist$bci.WD.mean, 
                         ifelse(!is.na(splist$bcigen.WD.mean), splist$bcigen.WD.mean, 
                                ifelse(!is.na(splist$nt.WD.mean), splist$nt.WD.mean, 
                                       ifelse(!is.na(splist$chave.WD.mean), splist$chave.WD.mean, 
                                              ifelse(!is.na(splist$try.WD.mean), splist$try.WD.mean, 
                                                     ifelse(!is.na(splist$pr.WD.mean), splist$pr.WD.mean, NA))))))

splist$WD.sd <- ifelse(!is.na(splist$bci.WD.sd), splist$bci.WD.sd, 
                       ifelse(!is.na(splist$bcigen.WD.sd), splist$bcigen.WD.sd, NA))

splist$WD.source <- ifelse(!is.na(splist$bci.WD.mean), 'BCI', 
                           ifelse(!is.na(splist$bcigen.WD.mean), 'BCIGEN', 
                                  ifelse(!is.na(splist$nt.WD.mean), 'NT', 
                                         ifelse(!is.na(splist$chave.WD.mean), 'CHAVE', 
                                                ifelse(!is.na(splist$chave.WD.mean), 'TRY', 
                                                       ifelse(!is.na(splist$pr.WD.mean), 'PR', NA))))))

tdata$wd <- splist$WD.mean[match(tdata$spplot, splist$sp)]
tdata$lma <- splist$LMA.mean[match(tdata$spplot, splist$sp)]
tdata$log.lma <- log(splist$LMA.mean)[match(tdata$spplot, splist$sp)]
tdata$wd.source <- splist$WD.source[match(tdata$spplot, splist$sp)]
tdata$lma.source <- splist$LMA.source[match(tdata$spplot, splist$sp)]

range(c(splist$mess.n.LMA, splist$bci.n.LMA), na.rm=T)
  
  
### HOW MANY INDIVIDUALS ARE COVERED IN TOTAL?
sum(!is.na(tdata$wd))/nrow(tdata)
sum(!is.na(tdata$lma))/nrow(tdata)

### HOW MANY INDIVIDUALS ARE COVERED WITH BCI DATA?
sum((tdata$latin) %in% splist$latin[!is.na(splist$bcigen.WD.mean)] ) / nrow(tdata)
sum((tdata$latin) %in% splist$latin[!is.na(splist$bcigen.LMA.mean) | !is.na(splist$mess.LMA.mean)] ) / nrow(tdata)

### HOW MANY INDIVIDUALS ARE COVERED?
table(is.na(tdata$wd), tdata$plot)[1,]/table(tdata$plot)
table(is.na(tdata$lma), tdata$plot)[1,]/table(tdata$plot)

### HOW MANY INDIVIDUALS ARE COVERED BY DIFFERENT DATA SOURCES?
table(tdata$plot, tdata$wd.source)
table(tdata$plot, tdata$lma.source)

### HOW MANY SPECIES IN TOTAL ARE COVERED IN THE WOOD DENSITY DATA? 
length(unique(tdata$latin[!is.na(tdata$wd)])) / length(unique(tdata$latin))


length(unique(tdata$latin[!is.na(tdata$lma)])) / length(unique(tdata$latin))

### WHAT PERCENTAGE OF SPECIES ARE INCLUDED IN BCI OR MESS DATA?
length(unique(splist$latin[!is.na(splist$bcigen.WD.mean)])) / length(unique(splist$latin))
length(unique(splist$latin[!is.na(splist$bcigen.LMA.mean) | !is.na(splist$mess.LMA.mean)])) / length(unique(splist$latin))

length(unique(splist$latin[!is.na(splist$LMA.mean)])) / length(unique(splist$latin))

length(unique(splist$latin[!is.na(splist$LMA.mean)]))
length(unique(splist$latin[!is.na(splist$bcigen.LMA.mean) | !is.na(splist$mess.LMA.mean)]))

length(unique(splist$latin[!is.na(splist$WD.mean)]))
length(unique(splist$latin[!is.na(splist$bcigen.WD.mean)]))




# save(tdata, file='panama_traits_preNCI_6.1.16.RDA')

#######################################
###  START HERE WITH PROCESSED DATA ###
#######################################
# load('panama_traits_preNCI_6.1.16.RDA') # tdata
# load('NCI_for_panama_traits_preNCI_6.14.16.RDA') # nci

### Current version from PC, post-NCI incorporation
# load('panama_NCI_Traits_6.14.16.RDA') # tdata

tdata.tmp <- tdata
load('panama_NCI_Traits_6.14.16.RDA') # tdata
tdata.tmp$All.NCI <- tdata$All.NCI
tdata.tmp$Con.NCI <- tdata$Con.NCI
tdata <- tdata.tmp

# load("alldata_NCI.RDA")
# load("panama_NCI_Traits_tNCI_uNCI_10.26.15.RDA")
# load("panama_traits_preNCI_12.9.15.RDA")
# load("panama_traits_preNCI_5.20.15.RDA")
# load("panama_traits_preNCI_6.1.16.RDA")

# library(data.table)
# tdata <- data.table(tdata)
# 
# Full.Neigh.Fun <- function(i, tdata){   
#   # If stem is not on the edge
#   if(tdata$Not.Edge[i] == 1){
# 
#     # Focal stem in row i
#     foc.tree <- as.data.frame(tdata[i,])
#     
#     # Gets neighboring stems within 20 meters and from the same census year
#     neighz <- as.data.frame(tdata[(1:nrow(tdata))!=i 
#                     & sqrt((tdata$x - foc.tree$x)^2 + (tdata$y - foc.tree$y)^2) <= 20
#                     & tdata$census == foc.tree$census 
#                     & tdata$plot == foc.tree$plot,])
#     
#     dist2 <- (neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2    
# 
#     # If some of the neighbors have same coords as focal.tree, add an offset
#     if(0 %in% dist2){
#       neighz[dist2==0 , c('x','y')] <- neighz[dist2==0, c('x','y')] + 0.25
#     }
#         
#     diam2 <- (neighz$dbh/10)^2
#     dist2 <- (neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2    
#     
#     cdist2 <- dist2[neighz$spcode %in% foc.tree$spcode]
#     cdiam2 <- diam2[neighz$spcode %in% foc.tree$spcode]
#     
#     bigneighz <- neighz[neighz$dbh >= foc.tree$dbh,]
#     bigdiam2 <- diam2[neighz$dbh >= foc.tree$dbh]
#     bigdist2 <- dist2[neighz$dbh >= foc.tree$dbh]
#     
#     # Calculate NCI using DBH^2 and dist^-2
#     All.NCI <- sum(diam2/dist2)
#     Con.NCI <- sum(cdiam2/cdist2)
#     Big.NCI <- sum(bigdiam2/bigdist2)
#     
#     wsg.diff <- foc.tree[,'WSG'] - neighz[,'WSG']
#     lma.diff <- foc.tree[,'log.LMALEAF_AVI'] - neighz[,'log.LMALEAF_AVI']
#     ldmc.diff <- foc.tree[,'log.LDMC_AVI'] - neighz[,'log.LDMC_AVI']
#     ss.diff <- foc.tree[,'log.SEED_DRY'] - neighz[,'log.SEED_DRY']
#     hmax.diff <- foc.tree[,'HEIGHT_AVG'] - neighz[,'HEIGHT_AVG']
#     
#     All.NCI.wsg.h <- sum(wsg.diff * (diam2/dist2), na.rm=T)
#     All.NCI.wsg.a <- sum(abs(wsg.diff) * (diam2/dist2), na.rm=T)
#     All.NCI.lma.h <- sum(lma.diff * (diam2/dist2), na.rm=T)
#     All.NCI.lma.a <- sum(abs(lma.diff) * (diam2/dist2), na.rm=T)
#     All.NCI.ldmc.h <- sum(ldmc.diff * (diam2/dist2), na.rm=T)
#     All.NCI.ldmc.a <- sum(abs(ldmc.diff) * (diam2/dist2), na.rm=T)
#     All.NCI.ss.h <- sum(ss.diff * (diam2/dist2), na.rm=T)
#     All.NCI.ss.a <- sum(abs(ss.diff) * (diam2/dist2), na.rm=T)
#     All.NCI.hmax.h <- sum(hmax.diff * (diam2/dist2), na.rm=T)
#     All.NCI.hmax.a <- sum(abs(hmax.diff) * (diam2/dist2), na.rm=T)
# 
#     wsg.diff.big <- foc.tree[,'WSG'] - bigneighz[,'WSG']
#     lma.diff.big <- foc.tree[,'log.LMALEAF_AVI'] - bigneighz[,'log.LMALEAF_AVI']
#     ldmc.diff.big <- foc.tree[,'log.LDMC_AVI'] - bigneighz[,'log.LDMC_AVI']
#     ss.diff.big <- foc.tree[,'log.SEED_DRY'] - bigneighz[,'log.SEED_DRY']
#     hmax.diff.big <- foc.tree[,'HEIGHT_AVG'] - bigneighz[,'HEIGHT_AVG']
#     
#     Big.NCI.wsg.h <- sum(wsg.diff.big * (bigdiam2/bigdist2), na.rm=T)
#     Big.NCI.wsg.a <- sum(abs(wsg.diff.big) * (bigdiam2/bigdist2), na.rm=T)
#     Big.NCI.lma.h <- sum(lma.diff.big * (bigdiam2/bigdist2), na.rm=T)
#     Big.NCI.lma.a <- sum(abs(lma.diff.big) * (bigdiam2/bigdist2), na.rm=T)
#     Big.NCI.ldmc.h <- sum(ldmc.diff.big * (bigdiam2/bigdist2), na.rm=T)
#     Big.NCI.ldmc.a <- sum(abs(ldmc.diff.big) * (bigdiam2/bigdist2), na.rm=T)
#     Big.NCI.ss.h <- sum(ss.diff.big * (bigdiam2/bigdist2), na.rm=T)
#     Big.NCI.ss.a <- sum(abs(ss.diff.big) * (bigdiam2/bigdist2), na.rm=T)
#     Big.NCI.hmax.h <- sum(hmax.diff.big * (bigdiam2/bigdist2), na.rm=T)
#     Big.NCI.hmax.a <- sum(abs(hmax.diff.big) * (bigdiam2/bigdist2), na.rm=T)
#     
#     All.NCI.wsg.a <- ifelse(is.na(foc.tree[,'WSG']), NA, All.NCI.wsg.a)
#     All.NCI.wsg.h <- ifelse(is.na(foc.tree[,'WSG']), NA, All.NCI.wsg.h)
#     Big.NCI.wsg.a <- ifelse(is.na(foc.tree[,'WSG']), NA, Big.NCI.wsg.a)
#     Big.NCI.wsg.h <- ifelse(is.na(foc.tree[,'WSG']), NA, Big.NCI.wsg.h)
# 
#     All.NCI.lma.a <- ifelse(is.na(foc.tree[,'log.LMALEAF_AVI']), NA, All.NCI.lma.a)
#     All.NCI.lma.h <- ifelse(is.na(foc.tree[,'log.LMALEAF_AVI']), NA, All.NCI.lma.h)
#     Big.NCI.lma.a <- ifelse(is.na(foc.tree[,'log.LMALEAF_AVI']), NA, Big.NCI.lma.a)
#     Big.NCI.lma.h <- ifelse(is.na(foc.tree[,'log.LMALEAF_AVI']), NA, Big.NCI.lma.h)
# 
#     All.NCI.ldmc.a <- ifelse(is.na(foc.tree[,'log.LDMC_AVI']), NA, All.NCI.ldmc.a)
#     All.NCI.ldmc.h <- ifelse(is.na(foc.tree[,'log.LDMC_AVI']), NA, All.NCI.ldmc.h)
#     Big.NCI.ldmc.a <- ifelse(is.na(foc.tree[,'log.LDMC_AVI']), NA, Big.NCI.ldmc.a)
#     Big.NCI.ldmc.h <- ifelse(is.na(foc.tree[,'log.LDMC_AVI']), NA, Big.NCI.ldmc.h)
# 
#     All.NCI.ss.a <- ifelse(is.na(foc.tree[,'log.SEED_DRY']), NA, All.NCI.ss.a)
#     All.NCI.ss.h <- ifelse(is.na(foc.tree[,'log.SEED_DRY']), NA, All.NCI.ss.h)
#     Big.NCI.ss.a <- ifelse(is.na(foc.tree[,'log.SEED_DRY']), NA, Big.NCI.ss.a)
#     Big.NCI.ss.h <- ifelse(is.na(foc.tree[,'log.SEED_DRY']), NA, Big.NCI.ss.h)
# 
#     All.NCI.hmax.a <- ifelse(is.na(foc.tree[,'HEIGHT_AVG']), NA, All.NCI.hmax.a)
#     All.NCI.hmax.h <- ifelse(is.na(foc.tree[,'HEIGHT_AVG']), NA, All.NCI.hmax.h)
#     Big.NCI.hmax.a <- ifelse(is.na(foc.tree[,'HEIGHT_AVG']), NA, Big.NCI.hmax.a)
#     Big.NCI.hmax.h <- ifelse(is.na(foc.tree[,'HEIGHT_AVG']), NA, Big.NCI.hmax.h)
# 
# return(
#       as.data.frame(cbind(All.NCI, Con.NCI, Big.NCI, 
#           All.NCI.wsg.h, All.NCI.wsg.a,
#           Big.NCI.wsg.h, Big.NCI.wsg.a,
#           All.NCI.lma.h, All.NCI.lma.a,
#           Big.NCI.lma.h, Big.NCI.lma.a,
#           All.NCI.ldmc.h, All.NCI.ldmc.a,
#           Big.NCI.ldmc.h, Big.NCI.ldmc.a,
#           All.NCI.ss.h, All.NCI.ss.a,
#           Big.NCI.ss.h, Big.NCI.ss.a,
#           All.NCI.hmax.h, All.NCI.hmax.a,
#           Big.NCI.hmax.h, Big.NCI.hmax.a
#           ))
#     )    
#     # For edge stems   
#   } else {
#      rep(NA, 22)
#   }  
# }
# 
# 
# Neigh.Fun <- function(i, tdata, r=15){   
#   # If stem is not on the edge
#   if(tdata$Not.Edge[i] == 1){
#     
#     # Focal stem in row i
#     foc.tree <- as.data.frame(tdata[i,])
#     
#     # Gets neighboring stems within 20 meters and from the same census year
#     neighz <- as.data.frame(tdata[(1:nrow(tdata))!=i 
#                     & sqrt((tdata$x - foc.tree$x)^2 + (tdata$y - foc.tree$y)^2) <= r
#                     & tdata$census == foc.tree$census 
#                     & tdata$plot == foc.tree$plot,])
#     
#     dist2 <- (neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2    
#     
#     # If some of the neighbors have same coords as focal.tree, add an offset
#     if(0 %in% dist2){
#       neighz[dist2==0 , c('x','y')] <- neighz[dist2==0 , c('x','y')] + 0.25
#     }
#   
#     diam2 <- (neighz$dbh/10)^2        
#     dist2 <- (neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2    
#     cdist2 <- dist2[neighz$spcode %in% foc.tree$spcode]
#     cdiam2 <- diam2[neighz$spcode %in% foc.tree$spcode]
#     
#     # Calculate NCI using DBH^2 and dist^-2
#     all.NCI <- sum(diam2/dist2)
#     con.NCI <- sum(cdiam2/cdist2)
#     
#     return(
#       as.data.frame(cbind(all.NCI, con.NCI))
#     )    
#     # For edge stems
#   } else {
#     rep(NA, 2)
#   }  
# }
# 
# 
# # TEST THE FUNCTION
# tdata <- as.data.table(tdata)
# samp <- sample(which(tdata$Not.Edge==T), 50)
# test <- do.call('rbind', lapply(samp, Neigh.Fun, tdata))
# 
# # do.call('rbind', lapply(samp, Full.Neigh.Fun, tdata))$All.NCI
# # unlist(lapply(samp[7], Neigh.Fun))
# 
# # DO THIS ON THE PC...
# library(parallel)
# library(snow)
# cl <- makeCluster(8, type='SOCK')
# nci <- parLapply(cl, 1:nrow(tdata), Full.Neigh.Fun, tdata)
# nci <- do.call('rbind', nci)
# stopCluster(cl)
# 
# tdata <- cbind(tdata, nci)

# save(tdata, file="panama_NCI_Traits_12.9.15.RDA")

# traits <- splist
# save(traits, file="panama_ITV_traits_6.7.16.RDA")



###############################
### CURRENT CONSIDERATIONS: ###
###############################
#load("panama_NCI_Traits_6.14.16.RDA") # tdata
load("panama_ITV_traits_6.7.16.RDA") # traits

head(tdata)
head(traits)

### Remove palms from growth analyses...
tdata$growth[tdata$palm==T] <- NA
tdata$RGR[tdata$palm==T] <- NA

### What to do about growth outliers? 
### One option is the remove stems that grew more than a fixed amount (e.g. > 5 sd per plot, per size class)
growth.include <- abs(tdata$growth) < sd(tdata$growth, na.rm=T) * 5
growth.include <- ifelse(is.na(growth.include), FALSE, growth.include)
tdata$Growth.Include <- growth.include

### Another is to follow Condit et al. 2004...
# 1. Negative growth must be smaller than  AND
# 2. Postive growth must be less than 75 mm / yr
tdata$Growth.Include.2 <- ((tdata$growth > 4 * (-(0.0062 * tdata$dbh + 0.904))) & (tdata$growth < 75))
tdata$Growth.Include.2 <- ifelse(is.na(tdata$Growth.Include.2), FALSE, tdata$Growth.Include.2)

### Another CONDIT METHOD TO EXCLUDE OUTLIERS
# 1. Negative growth must be smaller than 25% DBH AND
# 2. Postive growth must be less than 75 mm / yr
tdata$Growth.Include.3 <- (!is.na(tdata$growth) & (tdata$growth) < 75 & tdata$growth > (tdata$dbh * -0.25))

# sum(tdata$Growth.Include, na.rm=T)
# sum(tdata$Growth.Include.2, na.rm=T)
# sum(tdata$Growth.Include.3, na.rm=T)

# hist(tdata$growth[tdata$Growth.Include==T & tdata$Growth.Include.2==F], breaks=100)
# hist(tdata$growth[tdata$Growth.Include==F & tdata$Growth.Include.2==T], breaks=100)


#################
### DATA PREP ###
#################
### Remove edge trees and otherwise NA trees
#tdata <- tdata[tdata$Not.Edge %in% 1,] 
tdata <- tdata[!tdata$survival %in% NA,]
tdata <- droplevels(tdata)

### Z-TRANSFORM DATA
z.score <- function (data, scale=T, center=T) {
	xm <- mean (data, na.rm=TRUE)
	xsd <- sd(data, na.rm=TRUE)
  if(center==T & scale==T){ xtrans <- (data - xm)/( 2 * xsd) }
	if(center==F & scale==T){ xtrans <- (data)/( 2 * xsd) }
	if(center==T & scale==F){ xtrans <- (data - xm) }
	return(xtrans)
}

# ### Log-transform coefficients
# # log the NCI metric
# tdata$log.all.nci <- log(tdata$All.NCI)
# 
# # log the dbh
# tdata$log.dbh <- log(tdata$dbh)

# # log the tNCI and uNCI metrics
# tdata$log.tnci.hmax <- log(tdata$tnci.hmax)
# tdata$log.tnci.log.ldmc <- log(tdata$tnci.log.ldmc)
# tdata$log.tnci.log.lma <- log(tdata$tnci.log.lma)
# tdata$log.tnci.log.seed <- log(tdata$tnci.log.seed)
# tdata$log.tnci.wsg <- log(tdata$tnci.wsg)
# 
# tdata$log.unci.hmax <- log(tdata$unci.hmax + 1)
# tdata$log.unci.log.ldmc <- log(tdata$unci.log.ldmc + 1)
# tdata$log.unci.log.lma <- log(tdata$unci.log.lma + 1)
# tdata$log.unci.log.seed <- log(tdata$unci.log.seed + 1)
# tdata$log.unci.wsg <- log(tdata$unci.wsg + 1)


#save(tdata, file='Panama_AnalysisData_12.9.15.RDA')
# save(tdata, file='Panama_AnalysisData_6.14.16.RDA')
#save(tdata, file='Panama_AnalysisData_12.17.16.RDA')

save(tdata, file='Panama_AnalysisData_4.24.17.RDA')


head(traits)

cor(traits$LMA.mean, traits$WD.mean, use='p')

plot(traits$LMA.mean[traits$LMA.mean<200] , traits$WD.mean[traits$LMA.mean<200] )
cor.test(traits$LMA.mean[traits$LMA.mean<200] , traits$WD.mean[traits$LMA.mean<200] )
abline(lm(traits$WD.mean[traits$LMA.mean<200]  ~ traits$LMA.mean[traits$LMA.mean<200] ),col=2)
summary(lm(traits$WD.mean[traits$LMA.mean<200] ~ traits$LMA.mean[traits$LMA.mean<200]))

plot(log(traits$LMA.mean) , traits$WD.mean)
cor.test(log(traits$LMA.mean) , traits$WD.mean)
abline(lm(traits$WD.mean ~ log(traits$LMA.mean)),col=2)
summary(lm(traits$WD.mean ~ log(traits$LMA.mean)))


###############################################################
###   OLD WAY WAS TO CENTER WITHIN PLOTS WITHOUT SCALING   ###
###############################################################
# tdata <- tdata[order(tdata$plot, tdata$spcode, tdata$id, tdata$census),]
# tdata$log.nci.z <- unlist(tapply(tdata$log.nci, tdata$plot, scale, scale=F))
# tdata$log.tnci.z <- unlist(tapply(tdata$log.tnci, tdata$plot, scale, scale=F))
# tdata$log.dbh.z <- unlist(tapply(tdata$log.dbh, tdata$plot, scale, scale=F))
# tdata$growth.z <- unlist(tapply(tdata$growth, tdata$plot, scale, scale=F))





