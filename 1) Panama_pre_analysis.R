library(sp)
library(reshape)


############################################
### CHUNK 1 : BUILD PANAMA DATA FROM RAW ####
############################################

### BCI ###
setwd("/Users/Bob/Projects/Postdoc/Panama/DATA/bci")
# bci4 <- read.csv("BCI_census4.csv")
# bci5 <- read.csv("BCI_census5.csv")
# save(bci4, file='bci4.RDA')
# save(bci5, file='bci5.RDA')
load('bci4.RDA')
load('bci5.RDA')

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

bci <- bci4

spcodes <- read.csv("spcodes.csv")
bci$spcode <- spcodes$spcode[match(bci$Latin, paste(spcodes$Genus, spcodes$SpeciesName))]
bci$Latin <- NULL


bci <- bci[,c('Tag','spcode','PX','PY','Not.edge','plot', 'census', 'dbhStart', 'growth','survival','days')]
names(bci) <- c('tag','spcode','x','y','Not.Edge','plot','census','dbh','growth','survival','days')
bci$id <- rownames(bci)

bci <- bci[!is.na(bci$survival),]
bci <- bci[!is.na(bci$x),]
bci <- bci[!is.na(bci$y),]
bci$tag <- NULL
bci$nci <- NA

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA/bci")
# save(bci, file='bci_preNCI.RDA')
load(file='bci_preNCI.RDA')
head(bci)


#######################
### ADD OTHER PANAMA PLOTS ###
#######################
setwd("/Users/Bob/Projects/Postdoc/Panama/DATA/cocoli")
coc <- read.csv("cocoli.csv")
head(coc)

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA/sherman")
she <- read.csv("sherman.csv")
head(she)

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
valid1 <- ifelse(tmp$recr1=='A' & as.character(tmp$recr1) == as.character(tmp$recr2) & tmp$pom1 == tmp$pom2, 1, 0)
tmp$grow12 <- ifelse(valid1==1, tmp$grow1, NA)

tmp$grow23 <- ((tmp$dbh3 - tmp$dbh2) / (tmp$date3 - tmp$date2))* 365
valid2 <- ifelse(tmp$recr2=='A' & as.character(tmp$recr2) == as.character(tmp$recr3) & tmp$pom2 == tmp$pom3, 1, 0)
tmp$grow23 <- ifelse(valid2==1, tmp$grow2, NA)

# SURVIVAL
tmp$survive12 <- ifelse(tmp$recr1=='A' & as.character(tmp$recr1) == as.character(tmp$recr2), 1, NA)
tmp$survive12 <- ifelse(tmp$recr1=='A' & tmp$recr2=='D' , 0, tmp$survive1)

tmp$survive23 <- ifelse(tmp$recr2=='A' & as.character(tmp$recr2) == as.character(tmp$recr3), 1, NA)
tmp$survive23 <- ifelse(tmp$recr2=='A' & tmp$recr3=='D' , 0, tmp$survive2)

tmp$int12 <- tmp$date2 - tmp$date1
tmp$int23 <- tmp$date3 - tmp$date2
head(tmp)


# Clean up
tmp <- tmp[,c("tag","spcode","x","y",'dbh1','dbh2',"Not.Edge","plot","grow12","grow23","survive12","survive23","int12","int23")]
head(tmp)

tmp2 <- reshape(tmp, varying=c(list(5:6), list(9:10), list(11:12), list(13:14)), direction='long')
tmp2 <- tmp2[order(tmp2$tag),]
names(tmp2)[7:11] <- c('census','dbh','growth','survival','days')
tmp2 <- tmp2[,c(1:12)]
tmp2$tag <- NULL

rownames(tmp2) <- NULL
tdata <- tmp2
tdata$nci <- NA

tdata <- tdata[ ! is.na(tdata$dbh),]

head(tdata)

tdata <- rbind(tdata, bci)

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")
# save(tdata, file='panama_preNCI.RDA')


### RUN THIS ... (POSSIBLY  ON OTHER MACHINE, TAKES OVERNIGHT)
Neigh.Fun <- function(i){ # a function to calculate trait NCI for the stem in row i
  if(tdata$Not.Edge[i] == 1){ #if stem is not on the edge
    foc.tree <- tdata[i,] #focal stem in row i
    neighz <- tdata[(1:nrow(tdata))!=i & sqrt((tdata$x - foc.tree$x)^2 + (tdata$y - foc.tree$y)^2) <= 20 & tdata$census == foc.tree$census & tdata$plot == foc.tree$plot,] #gets neighboring stems within 20 meters and from the same census year
    if(sum(((neighz$x == foc.tree$x) + (neighz$y == foc.tree$y)) == 2) > 0){ #if some of the neighboring stems have the same coordinates as focal tree, add a slight offset
      neighz[((neighz$x == foc.tree$x) + (neighz$y == foc.tree$y)) == 2 ,3:4] <- neighz[((neighz$x == foc.tree$x) + (neighz$y == foc.tree$y)) == 2 ,3:4] + c(0.25, 0.25) #changes stems with same coordinates
    }
    #Calculate NCI
    NCI <- sum(na.omit(neighz$dbh^2/((neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2))) #gives NCI using DBH^2 and dist^-2'
    NCI
  }else{
    NA #For edge stems   
  }	
}

tdata$nci <- unlist(lapply(1:nrow(tdata), Neigh.Fun))
# save(tdata, file='panama_NCI.RDA')

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")
load('panama_NCI.RDA')


##################################
### CHUNK 2 : ATTACH LFDP DATA ####
##################################
setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")
load("panama_NCI.RDA")

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA/lfdp")
load("lfdp.RDA")
head(luq)

names(luq) <- names(tdata)
tdata <- rbind(tdata, luq)

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")
# save(tdata, file='alldata_NCI.RDA')


##################################
### CHUNK 2b : IDENTIFY PALMS ####
##################################
setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")
load('alldata_NCI.RDA')

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA/bci")
bci_spcodes <- read.csv("spcodes.csv")
setwd("/Users/Bob/Projects/Postdoc/Panama/DATA/cocoli")
cocoli_spcodes <- read.table("cocolisp.txt", header=T)
setwd("/Users/Bob/Projects/Postdoc/Panama/DATA/sherman")
sherman_spcodes <- read.table("shermansp.txt", header=T)

bci_spcodes <- bci_spcodes[,c('spcode','Genus','SpeciesName','Family')]
names(bci_spcodes) <- names(cocoli_spcodes)

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA/lfdp")
lfdp_spcodes <- read.csv("LFDP_spcodes.csv")
lfdp_spcodes <- lfdp_spcodes[,c('spcode','spcode','spcode','Family')]
names(lfdp_spcodes) <- names(cocoli_spcodes)

bci_spcodes$plot <- 'bci'
cocoli_spcodes$plot <- 'cocoli'
sherman_spcodes$plot <- 'sherman'
lfdp_spcodes$plot <- 'luq'

spcodes <- rbind(bci_spcodes, cocoli_spcodes, sherman_spcodes, lfdp_spcodes)
spcodes$spcodeplot <- paste(spcodes$spcode, spcodes$plot)

spcodeplot <- paste(tdata$spcode, tdata$plot)

tdata$palm <-  spcodeplot %in% spcodes$spcodeplot[spcodes$family == "Arecaceae" ]

is.palm <- tdata$palm
##################################
### CHUNK 3 : ATTACH TRAIT DATA ####
##################################
setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")
load('alldata_NCI.RDA')

tdata$binom <- NA
tdata$wd <- NA

head(tdata)

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA/cocoli")
splist <- read.table("cocolisp.txt", header=T)
splist$binom <- tolower(paste(splist$genus, splist$species, sep='_'))
tdata$binom[tdata$plot=='cocoli'] <- splist$binom[match(tdata$spcode[tdata$plot=='cocoli'], splist$spcode)]

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA/sherman")
splist <- read.table("shermansp.txt", header=T)
splist$binom <- tolower(paste(splist$genus, splist$species, sep='_'))
tdata$binom[tdata$plot=='sherman'] <- splist$binom[match(tdata$spcode[tdata$plot=='sherman'], splist$spcode)]

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA/bci")
splist <- read.csv("spcodes.csv", header=T)
splist$binom <- tolower(paste(splist$Genus, splist$SpeciesName, sep='_'))
tdata$binom[tdata$plot=='bci'] <- splist$binom[match(tdata$spcode[tdata$plot=='bci'], splist$spcode)]

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA/bci")
splist <- read.csv("spcodes.csv", header=T)
splist$binom <- tolower(paste(splist$Genus, splist$SpeciesName, sep='_'))
tdata$binom[tdata$plot=='bci'] <- splist$binom[match(tdata$spcode[tdata$plot=='bci'], splist$spcode)]

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")
wd <- read.csv("WoodDensity.csv")[,c(1:4)]
wd$binom <- tolower(paste(wd$GENUS, wd$SPECIES, sep='_'))

tdata$wd <- wd$WSG[match(tdata$binom, wd$binom)]

setwd("/Users/Bob/Projects/Thesis/DATA/traits")
luq <- read.csv("Site_Specific_mean_traits.csv", header=T)
luq <- luq[luq$FOREST=="LUQ",]
tdata$wd[tdata$plot=='lfdp'] <- luq$WD[match(tdata$spcode[tdata$plot=='lfdp'], luq$SPECIES)]


sum(!is.na(tdata$wd))/nrow(tdata) # WD data for ~75% of observations here...

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")
#save(tdata, file='alldata_NCI.RDA')


####################################################################
####################################################################
####################################################################
##################################
### CHUNK 4 : GET TRAIT NCI 		    ####
##################################

# A function to calculate trait NCI for the stem in row i

Trait.Neigh.Fun <- function(i, trait){ 

# If stem is not on the edge
	if(tdata$Not.Edge[i] == 1){

# Focal stem in row i
    	foc.tree <- tdata[i,] 

# Gets neighboring stems within 20 meters and from the same census year
    	neighz <- tdata[(1:nrow(tdata))!=i & sqrt((tdata$x - foc.tree$x)^2 + (tdata$y - foc.tree$y)^2) <= 20 & tdata$census == foc.tree$census & tdata$plot == foc.tree$plot,] 

# If some of the neighboring stems have the same coordinates as focal tree, add a slight offset
# Changes stems with same coordinates
	    if(sum(((neighz$x == foc.tree$x) + (neighz$y == foc.tree$y)) == 2) > 0){ 
    		  neighz[((neighz$x == foc.tree$x) + (neighz$y == foc.tree$y)) == 2 ,3:4] <- neighz[((neighz$x == foc.tree$x) + (neighz$y == foc.tree$y)) == 2 ,3:4] + c(0.01, 0.01)
	    }
# Calculate NCI
# First calculate distance between focal stem and each other stem

	nci.abs.trait.diff <- abs(foc.tree[,trait] - neighz[,trait])

# Gives NCI using DBH^2 and dist^-2'
	nci.biomass <- na.omit(neighz$dbh^2/((neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2))

# ? # HOW TO DEAL WITH STEMS WITH MISSING TRAIT DATA?
    NCI.abs.trait <- sum(nci.abs.trait.diff * nci.biomass, na.rm=T)

# For edge stems   
  } else {
  	NA
  }	
}


tdata$tnci <- unlist(lapply(1:nrow(tdata), Trait.Neigh.Fun, 'wd'))

# save(tdata, file='data_10.13.15.RDA')
load("data_10.13.15.RDA")

head(tdata)

