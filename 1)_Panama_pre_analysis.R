library(sp)
library(reshape)

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")
setwd("K:/Bob/Panama/GIT/Panama_analysis/DATA")
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
bci4 <- bci4[ ! bci4$Status %in% c('dead','missing'),]
bci4 <- bci4[ ! is.na(bci4$DBH),]
bci4 <- droplevels(bci4)
bci5 <- bci5[ ! bci5$Status %in% c('missing'),]

# SET BORDER TO 20 METERS
bci4$Not.edge <- ifelse(bci4$PX > 20 & bci4$PX < 980 & bci4$PY > 20 & bci4$PY < 480, 1, 0)
bci5$Not.edge <- ifelse(bci5$PX > 20 & bci5$PX < 980 & bci5$PY > 20 & bci5$PY < 480, 1, 0)

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

bci4$survival <- ifelse(bci4$StatusFinal == 'alive', 1, 0)
bci4$days <- as.numeric(bci4$DateFinal - bci4$Date)
bci4$growth <- ((bci4$dbhFinal - bci4$dbhStart)/bci4$days) * 365
bci4$log.growth <- ((log(bci4$dbhFinal) - log(bci4$dbhStart))/(bci4$days)) * 365

bci <- bci4

# NOTE THAT I CHANGED 3 SPECIES CODES FROM THE ONLINE BCI SPCODE DATA 
# TO BE CONSISTENT WITH THE BCI TRAIT DATA...
spcodes <- read.csv("bci/spcodes.csv")
bci$spcode <- spcodes$spcode[match(bci$Latin, paste(spcodes$Genus, spcodes$SpeciesName))]

bci <- droplevels(bci)

bci <- bci[,c('Tag','Latin','spcode','PX','PY','Not.edge','plot', 'census', 'dbhStart', 'growth','survival','days')]
names(bci) <- c('tag','latin','spcode','x','y','Not.Edge','plot','census','dbh','growth','survival','days')
bci$id <- rownames(bci)

bci <- bci[!is.na(bci$survival),]
bci <- bci[!is.na(bci$x),]
bci <- bci[!is.na(bci$y),]
bci$tag <- NULL
bci$nci <- NA

# save(bci, file='bci/bci_preNCI.RDA')
# save(bci, file='bci/bci_preNCI_10.23.15.RDA')
# save(bci, file='bci/bci_preNCI_12.9.15.RDA')
load(file='bci/bci_preNCI_12.9.15.RDA')
head(bci)


################################
### ADD OTHER PANAMA PLOTS ###
################################
coc <- read.csv("cocoli/cocoli.csv")
cocsp <- read.table("cocoli/cocolisp.txt", header=T)

she <- read.csv("sherman/sherman.csv")
shesp <- read.table("sherman/shermansp.txt", header=T)

# plot(she$x, she$y, xlim=c(0,450), ylim=c(0,450), pch=16, cex=.2)
she$Not.Edge <- point.in.polygon(she$x, she$y, c(20,120,120,220,220,160,160,20), c(20,20,60,60,320,320,120,120))
she$Not.Edge <- ifelse(she$Not.Edge > 0, 1, 0)
# points(she$x, she$y, col= she$Not.Edge +1, pch=16, cex=.2)
# sort(table(she$spcode[she$Not.Edge ==1]))

# new <- point.in.polygon(she$x, she$y, c(140,240,240,140), c(340,340,440,440))
# col <- ifelse(new==1, 3, she$Not.Edge+1)
# points(she$x, she$y, col= col, pch=16, cex=.2)

# plot(coc$x, coc$y, xlim=c(0,200), ylim=c(0,300), pch=16, cex=.5)
coc$Not.Edge <- point.in.polygon(coc$x, coc$y, c(20,180,180,80,80,20), c(20,20,80,80,280,280))
coc$Not.Edge <- ifelse(coc$Not.Edge > 0, 1, 0)
# points(coc$x, coc$y, col= coc$Not.Edge +1, pch=16, cex=.5)
# sum(table(coc$spcode[coc$Not.Edge==1]))

she$plot <- "sherman"
coc$plot <- "cocoli"

coc <- coc[-8868,] # remove a duplicate tag

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

tmp$grow12 <- ((tmp$dbh2 - tmp$dbh1) / (tmp$date2 - tmp$date1))* 365
tmp$log.grow12 <- ((log(tmp$dbh2) - log(tmp$dbh1)) / (tmp$date2 - tmp$date1))* 365

valid1 <- ifelse(tmp$recr1=='A' & as.character(tmp$recr1) == as.character(tmp$recr2) & tmp$pom1 == tmp$pom2, 1, 0)
tmp$grow12 <- ifelse(valid1==1, tmp$grow12, NA)
tmp$log.grow12 <- ifelse(valid1==1, tmp$log.grow12, NA)

tmp$grow23 <- ((tmp$dbh3 - tmp$dbh2) / (tmp$date3 - tmp$date2))* 365
tmp$log.grow23 <- ((log(tmp$dbh3) - log(tmp$dbh2) / (tmp$date3 - tmp$date2)))* 365
valid2 <- ifelse(tmp$recr2=='A' & as.character(tmp$recr2) == as.character(tmp$recr3) & tmp$pom2 == tmp$pom3, 1, 0)
tmp$grow23 <- ifelse(valid2==1, tmp$grow23, NA)
tmp$log.grow23 <- ifelse(valid2==1, tmp$log.grow23, NA)

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
						"grow12","grow23","survive12","survive23","int12","int23")]

tmp2 <- reshape(tmp, varying=c(list(5:6), list(9:10), list(11:12), list(13:14)), direction='long')
tmp2 <- tmp2[order(tmp2$tag),]
names(tmp2)[7:11] <- c('census','dbh','growth','survival','days')
tmp2 <- tmp2[,c(1:12)]
tmp2$tag <- NULL

rownames(tmp2) <- NULL
tdata <- tmp2
tdata$nci <- NA

tdata <- tdata[ ! is.na(tdata$dbh),]

tdata$latin[tdata$plot=='cocoli'] <- paste(cocsp$genus, cocsp$species)[match(tdata$spcode, cocsp$spcode)][tdata$plot=='cocoli']
tdata$latin[tdata$plot=='sherman'] <- paste(shesp$genus, shesp$species)[match(tdata$spcode, shesp$spcode)][tdata$plot=='sherman']

tdata <- tdata[,c("latin","spcode","x","y","Not.Edge","plot","census","dbh","growth","survival","days","id","nci")]

tdata <- rbind(tdata, bci)

tdata <- tdata[tdata$dbh != 0 , ]

# save(tdata, file='panama_preNCI.RDA')
# save(tdata, file='panama_preNCI_10.23.15.RDA')
# save(tdata, file='panama_preNCI_12.9.15.RDA')
load('panama_preNCI_12.9.15.RDA')


##################################
### CHUNK 2 : ATTACH LFDP DATA ####
##################################
# load('panama_NCI_10.26.15.RDA')

luqsp <- read.csv("lfdp/LFDP_spcodes.csv", row.names=1)
load("lfdp/lfdp.RDA")
luq$latin <- luqsp$latin[match(luq$SPECIES, luqsp$spcode)]
luq <- as.data.frame(cbind(luq$latin, luq[,!colnames(luq) %in% 'latin']))
names(luq) <- names(tdata)
tdata <- rbind(tdata, luq)

##################################
### CHUNK 2b : IDENTIFY PALMS ####
##################################
bcisp <- read.csv("bci/spcodes.csv")
cocsp <- read.table("cocoli/cocolisp.txt", header=T)
shesp <- read.table("sherman/shermansp.txt", header=T)
bcisp <- bcisp[,c('spcode','Genus','SpeciesName','Family')]
names(bcisp) <- names(cocsp)
luqsp <- read.csv("lfdp/LFDP_spcodes.csv", row.names=1)
luqsp <- luqsp[,c('spcode','genus','species','family')]

bcisp$plot <- 'bci'
cocsp$plot <- 'cocoli'
shesp$plot <- 'sherman'
luqsp$plot <- 'lfdp'

spcodes <- rbind(bcisp, cocsp, shesp, luqsp)
spcodes$spcodeplot <- paste(spcodes$spcode, spcodes$plot)
spcodeplot <- paste(tdata$spcode, tdata$plot)

tdata$palm <-  spcodeplot %in% spcodes$spcodeplot[spcodes$family == "Arecaceae" ]

##################################
### CHUNK 3 : ATTACH TRAIT DATA ####
##################################

### EXCLUDE LFDP DATA FOR NOW... 
tdata <- tdata[!tdata$plot %in% 'lfdp',]
tdata <- droplevels(tdata)

#######################################################
#### THIS WAS TO USE THE PUBLIC WD DATA FROM CHAVE 2006
# wd <- read.csv("traits/ChaveWoodDensity.csv")[,c(1:4)]
# wd$binom <- tolower(paste(wd$GENUS, wd$SPECIES, sep='_'))
# tdata$wd <- wd$WSG[match(tdata$binom, wd$binom)]
#######################################################

bci <- read.csv("traits/BCITRAITS_20101220.csv")

# GET WSG 100C IF AVAILABLE, OTHERWISE WSG 60C AND THEN CHAVE WSG:
bci$WSG <- ifelse(!is.na(bci$SG100C_AVG), bci$SG100C_AVG, ifelse(!is.na(bci$SG60C_AVG), bci$SG60C_AVG, bci$WSG_CHAVE))

# SUBSET TO TRAITS OF PRIMARY INTEREST...
bci <- bci[,c('SP.', 'WSG','HEIGHT_AVG','SEED_DRY','LEAFAREA_AVI',
					'LEAFTHCK_AVI','LMALEAF_AVI','LDMC_AVI','AVG_LAMTUF',
					'RGR_10','RGR_50','RGR_100','MORT_10','MORT_100')]

### MATCH SP CODES WITH TRAITS..... 
traits <- bci[match(tdata$spcode, bci$SP.),]
tdata <- cbind(tdata, traits[,-1])

# setwd("/Users/Bob/Projects/Thesis/DATA/traits")
# luq <- read.csv("Site_Specific_mean_traits.csv", header=T, row.names=1)
# luq <- luq[luq$FOREST=="LUQ",]
# head(luq)
# tdata$WSG[tdata$plot=='lfdp'] <- luq$WD[match(tdata$spcode[tdata$plot=='lfdp'], luq$SPECIES)]

sum(!is.na(tdata$WSG))/nrow(tdata) # WD data for ~95% of observations here...

#######################################
###  START HERE WITH PROCESSED DATA ###
#######################################
# load("alldata_NCI.RDA")
#load("panama_NCI_Traits_tNCI_uNCI_10.26.15.RDA")
head(tdata)


Full.Neigh.Fun <- function(i, tdata){   
  # If stem is not on the edge
  if(tdata$Not.Edge[i] == 1){

    # Focal stem in row i
    foc.tree <- tdata[i,] 
    
    # Gets neighboring stems within 20 meters and from the same census year
    neighz <- tdata[(1:nrow(tdata))!=i 
                    & sqrt((tdata$x - foc.tree$x)^2 + (tdata$y - foc.tree$y)^2) <= 20
                    & tdata$census == foc.tree$census 
                    & tdata$plot == foc.tree$plot,]
    
    # If some of the neighboring stems have the same coordinates as focal tree, add a slight offset
    # Changes stems with same coordinates
    if(sum(((neighz$x == foc.tree$x) + (neighz$y == foc.tree$y)) == 2) > 0){
      neighz[((neighz$x == foc.tree$x) + (neighz$y == foc.tree$y)) == 2 ,3:4] <- neighz[((neighz$x == foc.tree$x) + (neighz$y == foc.tree$y)) == 2 ,3:4] + c(0.01, 0.01)
    }
    
    # Calculate NCI using DBH^2 and dist^-2
    All.NCI <- sum(na.omit(neighz$dbh^2/((neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2)))
    neighzbigger <- neighz[neighz$dbh >= foc.tree$dbh,]
    Big.NCI <- sum(na.omit(neighzbigger$dbh^2/((neighzbigger$x - foc.tree$x)^2 + (neighzbigger$y - foc.tree$y)^2)))

    wsg.diff <- foc.tree[,'WSG'] - neighz[,'WSG']
    lma.diff <- foc.tree[,'log.LMALEAF_AVI'] - neighz[,'log.LMALEAF_AVI']
    ldmc.diff <- foc.tree[,'log.LDMC_AVI'] - neighz[,'log.LDMC_AVI']
    ss.diff <- foc.tree[,'log.SEED_DRY'] - neighz[,'log.SEED_DRY']
    hmax.diff <- foc.tree[,'HEIGHT_AVG'] - neighz[,'HEIGHT_AVG']
    
    All.NCI.wsg.h <- sum(wsg.diff * na.omit(neighz$dbh^2/((neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2)), na.rm=T)
    All.NCI.wsg.a <- sum(abs(wsg.diff) * na.omit(neighz$dbh^2/((neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2)), na.rm=T)
    All.NCI.lma.h <- sum(lma.diff * na.omit(neighz$dbh^2/((neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2)), na.rm=T)
    All.NCI.lma.a <- sum(abs(lma.diff) * na.omit(neighz$dbh^2/((neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2)), na.rm=T)
    All.NCI.ldmc.h <- sum(ldmc.diff * na.omit(neighz$dbh^2/((neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2)), na.rm=T)
    All.NCI.ldmc.a <- sum(abs(ldmc.diff) * na.omit(neighz$dbh^2/((neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2)), na.rm=T)
    All.NCI.ss.h <- sum(ss.diff * na.omit(neighz$dbh^2/((neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2)), na.rm=T)
    All.NCI.ss.a <- sum(abs(ss.diff) * na.omit(neighz$dbh^2/((neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2)), na.rm=T)
    All.NCI.hmax.h <- sum(hmax.diff * na.omit(neighz$dbh^2/((neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2)), na.rm=T)
    All.NCI.hmax.a <- sum(abs(hmax.diff) * na.omit(neighz$dbh^2/((neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2)), na.rm=T)
    
    wsg.diff.bigger <- foc.tree[,'WSG'] - neighzbigger[,'WSG']
    lma.diff.bigger <- foc.tree[,'log.LMALEAF_AVI'] - neighzbigger[,'log.LMALEAF_AVI']
    ldmc.diff.bigger <- foc.tree[,'log.LDMC_AVI'] - neighzbigger[,'log.LDMC_AVI']
    ss.diff.bigger <- foc.tree[,'log.SEED_DRY'] - neighzbigger[,'log.SEED_DRY']
    hmax.diff.bigger <- foc.tree[,'HEIGHT_AVG'] - neighzbigger[,'HEIGHT_AVG']
    
    Big.NCI.wsg.h <- sum(wsg.diff.bigger * na.omit(neighzbigger$dbh^2/((neighzbigger$x - foc.tree$x)^2 + (neighzbigger$y - foc.tree$y)^2)), na.rm=T)
    Big.NCI.wsg.a <- sum(abs(wsg.diff.bigger) * na.omit(neighzbigger$dbh^2/((neighzbigger$x - foc.tree$x)^2 + (neighzbigger$y - foc.tree$y)^2)), na.rm=T)
    Big.NCI.lma.h <- sum(lma.diff.bigger * na.omit(neighzbigger$dbh^2/((neighzbigger$x - foc.tree$x)^2 + (neighzbigger$y - foc.tree$y)^2)), na.rm=T)
    Big.NCI.lma.a <- sum(abs(lma.diff.bigger) * na.omit(neighzbigger$dbh^2/((neighzbigger$x - foc.tree$x)^2 + (neighzbigger$y - foc.tree$y)^2)), na.rm=T)
    Big.NCI.ldmc.h <- sum(ldmc.diff.bigger * na.omit(neighzbigger$dbh^2/((neighzbigger$x - foc.tree$x)^2 + (neighzbigger$y - foc.tree$y)^2)), na.rm=T)
    Big.NCI.ldmc.a <- sum(abs(ldmc.diff.bigger) * na.omit(neighzbigger$dbh^2/((neighzbigger$x - foc.tree$x)^2 + (neighzbigger$y - foc.tree$y)^2)), na.rm=T)
    Big.NCI.ss.h <- sum(ss.diff.bigger * na.omit(neighzbigger$dbh^2/((neighzbigger$x - foc.tree$x)^2 + (neighzbigger$y - foc.tree$y)^2)), na.rm=T)
    Big.NCI.ss.a <- sum(abs(ss.diff.bigger) * na.omit(neighzbigger$dbh^2/((neighzbigger$x - foc.tree$x)^2 + (neighzbigger$y - foc.tree$y)^2)), na.rm=T)
    Big.NCI.hmax.h <- sum(hmax.diff.bigger * na.omit(neighzbigger$dbh^2/((neighzbigger$x - foc.tree$x)^2 + (neighzbigger$y - foc.tree$y)^2)), na.rm=T)
    Big.NCI.hmax.a <- sum(abs(hmax.diff.bigger) * na.omit(neighzbigger$dbh^2/((neighzbigger$x - foc.tree$x)^2 + (neighzbigger$y - foc.tree$y)^2)), na.rm=T)

    All.NCI.wsg.a <- ifelse(is.na(foc.tree[,'WSG']), NA, All.NCI.wsg.a)
    All.NCI.wsg.h <- ifelse(is.na(foc.tree[,'WSG']), NA, All.NCI.wsg.h)
    Big.NCI.wsg.a <- ifelse(is.na(foc.tree[,'WSG']), NA, Big.NCI.wsg.a)
    Big.NCI.wsg.h <- ifelse(is.na(foc.tree[,'WSG']), NA, Big.NCI.wsg.h)

    All.NCI.lma.a <- ifelse(is.na(foc.tree[,'log.LMALEAF_AVI']), NA, All.NCI.lma.a)
    All.NCI.lma.h <- ifelse(is.na(foc.tree[,'log.LMALEAF_AVI']), NA, All.NCI.lma.h)
    Big.NCI.lma.a <- ifelse(is.na(foc.tree[,'log.LMALEAF_AVI']), NA, Big.NCI.lma.a)
    Big.NCI.lma.h <- ifelse(is.na(foc.tree[,'log.LMALEAF_AVI']), NA, Big.NCI.lma.h)

    All.NCI.ldmc.a <- ifelse(is.na(foc.tree[,'log.LDMC_AVI']), NA, All.NCI.ldmc.a)
    All.NCI.ldmc.h <- ifelse(is.na(foc.tree[,'log.LDMC_AVI']), NA, All.NCI.ldmc.h)
    Big.NCI.ldmc.a <- ifelse(is.na(foc.tree[,'log.LDMC_AVI']), NA, Big.NCI.ldmc.a)
    Big.NCI.ldmc.h <- ifelse(is.na(foc.tree[,'log.LDMC_AVI']), NA, Big.NCI.ldmc.h)

    All.NCI.ss.a <- ifelse(is.na(foc.tree[,'log.SEED_DRY']), NA, All.NCI.ss.a)
    All.NCI.ss.h <- ifelse(is.na(foc.tree[,'log.SEED_DRY']), NA, All.NCI.ss.h)
    Big.NCI.ss.a <- ifelse(is.na(foc.tree[,'log.SEED_DRY']), NA, Big.NCI.ss.a)
    Big.NCI.ss.h <- ifelse(is.na(foc.tree[,'log.SEED_DRY']), NA, Big.NCI.ss.h)

    All.NCI.hmax.a <- ifelse(is.na(foc.tree[,'HEIGHT_AVG']), NA, All.NCI.hmax.a)
    All.NCI.hmax.h <- ifelse(is.na(foc.tree[,'HEIGHT_AVG']), NA, All.NCI.hmax.h)
    Big.NCI.hmax.a <- ifelse(is.na(foc.tree[,'HEIGHT_AVG']), NA, Big.NCI.hmax.a)
    Big.NCI.hmax.h <- ifelse(is.na(foc.tree[,'HEIGHT_AVG']), NA, Big.NCI.hmax.h)

return(
      as.data.frame(cbind(All.NCI, Big.NCI, 
          All.NCI.wsg.h, All.NCI.wsg.a,
          Big.NCI.wsg.h, Big.NCI.wsg.a,
          All.NCI.lma.h, All.NCI.lma.a,
          Big.NCI.lma.h, Big.NCI.lma.a,
          All.NCI.ldmc.h, All.NCI.ldmc.a,
          Big.NCI.ldmc.h, Big.NCI.ldmc.a,
          All.NCI.ss.h, All.NCI.ss.a,
          Big.NCI.ss.h, Big.NCI.ss.a,
          All.NCI.hmax.h, All.NCI.hmax.a,
          Big.NCI.hmax.h, Big.NCI.hmax.a
          ))
    )    
    # For edge stems   
  } else {
    rep(NA, 22)
  }  
}

# TEST THE FUNCTION
# samp <- sample(which(tdata$Not.Edge==T), 5)
# tmp <- do.call('rbind', lapply(samp, Full.Neigh.Fun))


library(parallel)
library(snow)

cl <- makeCluster(8, type='SOCK')
nci <- parLapply(cl, 1:nrow(tdata), Full.Neigh.Fun, tdata)
nci <- do.call('rbind', nci)
stopCluster(cl)


tdata <- cbind(tdata, nci)

# save(tdata, file="panama_NCI_Traits_11.14.15.RDA")


names(tdata)

cor(log(tdata$All.NCI+1), log(tdata$Big.NCI+1), use='p')
plot(log(tdata$All.NCI+1), log(tdata$Big.NCI+1))

tmp <- tdata[log(tdata$All.NCI+1) > 18,] 
dim(tmp)
table(tmp$plot)



###############################
### CURRENT CONSIDERATIONS: ###
###############################
load("panama_NCI_Traits_11.14.15.RDA")

### Remove LFDP data... (for now)
tdata <- tdata[tdata$plot != 'lfdp',]

### Remove palms (doing it this way to keep them for survival...)
tdata$growth[tdata$palm==T] <- NA

### What to do about growth outliers? 
### One option is the remove stems that grew more than a fixed amount (e.g. > 5 sd )
growth.include <- abs(tdata$growth) < sd(tdata$growth, na.rm=T) * 5
growth.include <- ifelse(is.na(growth.include), FALSE, growth.include)
tdata$Growth.Include <- growth.include

### Another is to follow Condit et al. 2004...
# 1. Negative growth must be smaller than  AND
# 2. Postive growth must be less than 75 mm / yr
tdata$Growth.Include.2 <- ((tdata$growth > 4 * (-(0.0062 * tdata$dbh + 0.904))) & (tdata$growth < 75))
tdata$Growth.Include.2 <- ifelse(is.na(tdata$Growth.Include.2), FALSE, tdata$Growth.Include.2)

# sum(!tdata$Growth.Include, na.rm=T)
# sum(!tdata$Growth.Include.2, na.rm=T)
# hist(tdata$growth[tdata$Growth.Include==T & tdata$Growth.Include.2==F], breaks=100)
# hist(tdata$growth[tdata$Growth.Include==F & tdata$Growth.Include.2==T], breaks=100)


#################
### DATA PREP ###
#################
### Remove edge trees and otherwise NA trees
tdata <- tdata[tdata$Not.Edge %in% 1,] 
tdata <- tdata[ ! tdata$survival %in% NA,]
tdata <- droplevels(tdata)

### Z-TRANSFORM DATA
z.score <- function (data) {
	xm<- mean (data, na.rm=TRUE)
	xsd<-sd(data, na.rm=TRUE)
	xtrans<-(data-xm)/(2*xsd)	
}

### Log-transform coefficients
# log the NCI metric
tdata$log.nci <- log(tdata$nci + 1)

# log the dbh
tdata$log.dbh <- log(tdata$dbh)

# log the tNCI and uNCI metrics
tdata$log.tnci.hmax <- log(tdata$tnci.hmax + 1)
tdata$log.unci.hmax <- log(tdata$unci.hmax + 1)
tdata$log.tnci.log.ldmc <- log(tdata$tnci.log.ldmc + 1)
tdata$log.unci.log.ldmc <- log(tdata$unci.log.ldmc + 1)
tdata$log.tnci.log.lma <- log(tdata$tnci.log.lma + 1)
tdata$log.unci.log.lma <- log(tdata$unci.log.lma + 1)
tdata$log.tnci.log.seed <- log(tdata$tnci.log.seed + 1)
tdata$log.unci.log.seed <- log(tdata$unci.log.seed + 1)
tdata$log.tnci.wsg <- log(tdata$tnci.wsg + 1)
tdata$log.unci.wsg <- log(tdata$unci.wsg + 1)


#########################################
####    Standardize coefficients within plots    ####
#########################################
tdata <- tdata[order(tdata$plot, tdata$spcode, tdata$id, tdata$census),]

tdata$growth.z <- unlist(tapply(tdata$growth, tdata$plot, scale, center=F))

tdata$log.nci.z <- unlist(tapply(tdata$log.nci, tdata$plot, z.score))

tdata$log.dbh.z <- unlist(tapply(tdata$log.dbh, tdata$plot, z.score))

tdata$log.tnci.wsg.z <- unlist(tapply(tdata$log.tnci.wsg, tdata$plot, z.score))
tdata$log.unci.wsg.z <- unlist(tapply(tdata$log.unci.wsg, tdata$plot, z.score))
tdata$log.tnci.log.ldmc.z <- unlist(tapply(tdata$log.tnci.log.ldmc, tdata$plot, z.score))
tdata$log.unci.log.ldmc.z <- unlist(tapply(tdata$log.unci.log.ldmc, tdata$plot, z.score))
tdata$log.tnci.log.lma.z <- unlist(tapply(tdata$log.tnci.log.lma, tdata$plot, z.score))
tdata$log.unci.log.lma.z <- unlist(tapply(tdata$log.unci.log.lma, tdata$plot, z.score))
tdata$log.tnci.log.seed.z <- unlist(tapply(tdata$log.tnci.log.seed, tdata$plot, z.score))
tdata$log.unci.log.seed.z <- unlist(tapply(tdata$log.unci.log.seed, tdata$plot, z.score))
tdata$log.tnci.hmax.z <- unlist(tapply(tdata$log.tnci.hmax, tdata$plot, z.score))
tdata$log.unci.hmax.z <- unlist(tapply(tdata$log.unci.hmax, tdata$plot, z.score))


###############################################################
###   OLD WAY WAS TO CENTER WITHIN PLOTS WITHOUT SCALING   ###
###############################################################
# tdata <- tdata[order(tdata$plot, tdata$spcode, tdata$id, tdata$census),]
# tdata$log.nci.z <- unlist(tapply(tdata$log.nci, tdata$plot, scale, scale=F))
# tdata$log.tnci.z <- unlist(tapply(tdata$log.tnci, tdata$plot, scale, scale=F))
# tdata$log.dbh.z <- unlist(tapply(tdata$log.dbh, tdata$plot, scale, scale=F))
# tdata$growth.z <- unlist(tapply(tdata$growth, tdata$plot, scale, scale=F))

save(tdata, file='Panama_AnalysisData_10.30.15.RDA')





