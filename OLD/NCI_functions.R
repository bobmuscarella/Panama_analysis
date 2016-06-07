#master <- tdata
d
### RUN THIS ... (ON OTHER MACHINE, TAKES OVERNIGHT)
Neigh.Fun <- function(i){ # a function to calculate trait NCI for the stem in row i
  if(tdata$Not.Edge[i] == 1){ #if stem is not on the edge
    foc.tree <- tdata[i,] #focal stem in row i
    neighz <- tdata[(1:nrow(tdata))!=i 
                    & sqrt((tdata$x - foc.tree$x)^2 + (tdata$y - foc.tree$y)^2) <= 20 
                    & tdata$census == foc.tree$census 
                    & tdata$plot == foc.tree$plot,] #gets neighboring stems within 20 meters and from the same census year
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
## tdata$nci <- unlist(lapply(1:nrow(tdata), Neigh.Fun))

### GET THE COMPLETED FILE FROM OTHER COMPUTER AND SAVE NCI COLUMN TO MASTER
# load("postNCI_calculations/panama_NCI_notrait_10.23.15.RDA")
# master$nci <- tdata$nci
# tdata <- master
# rm(master)

# save(tdata, file='panama_NCI.RDA')
# save(tdata, file='panama_NCI_10.26.15.RDA')



#////////////////////\\\\\\\\\\\\\\\\\\\\\\\\
#///////////////////  \\\\\\\\\\\\\\\\\\\\\\\\
#\\\\\\\\\\\\\\\\\\\  ////////////////////////
#\\\\\\\\\\\\\\\\\\\\////////////////////////

######################################################
### CHUNK 4 : GET TRAIT NCI and UNKNOWN TRAIT NCI    ####
######################################################
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

Unk.Trait.Neigh.Fun <- function(i, trait){ 
  
  # If stem is not on the edge
  if(tdata$Not.Edge[i] == 1){
    
    # Focal stem in row i
    foc.tree <- tdata[i,] 
    
    # Gets neighboring stems within 20 meters and from the same census year
    neighz <- tdata[(1:nrow(tdata))!=i & sqrt((tdata$x - foc.tree$x)^2 + (tdata$y - foc.tree$y)^2) <= 20 & tdata$census == foc.tree$census & tdata$plot == foc.tree$plot & is.na(tdata[,trait]),] 
    
    # If some of the neighboring stems have the same coordinates as focal tree, add a slight offset
    # Changes stems with same coordinates
    if(sum(((neighz$x == foc.tree$x) + (neighz$y == foc.tree$y)) == 2) > 0){
      neighz[((neighz$x == foc.tree$x) + (neighz$y == foc.tree$y)) == 2 ,3:4] <- neighz[((neighz$x == foc.tree$x) + (neighz$y == foc.tree$y)) == 2 ,3:4] + c(0.01, 0.01)
    }
    # Calculate NCI
    # Gives NCI using DBH^2 and dist^-2'
    Unk.Trait.NCI <- sum(na.omit(neighz$dbh^2/((neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2)))
    Unk.Trait.NCI
    # For edge stems   
  } else {
    NA
  }	
}

### Log transform skewed traits *before* calculating tNCI...
# t <- 1
# par(mfrow=c(2,2))
# hist(bci[,names(traits)[-1][t]])
# hist(log(bci[,names(traits)[-1][t]]))
# plot(1,1,col=NA, axes=F)
# text(1,1,names(traits)[-1][t])

logtraits <- c('SEED_DRY','LEAFAREA_AVI','LEAFTHCK_AVI',
               'LMALEAF_AVI','LDMC_AVI','AVG_LAMTUF','RGR_10',
               'RGR_50','RGR_100','MORT_10','MORT_100')

for(i in 1:length(logtraits)) {
  newname <- paste('log', logtraits[i], sep='.')
  tdata[,newname] <- log(tdata[,logtraits[i]])
}

# Get tNCI for all traits of interest (biomass NCI weighted by absolute trait difference)
traits <- c('WSG','HEIGHT_AVG','log.SEED_DRY','log.LMALEAF_AVI','log.LDMC_AVI')

# traits <- c('WSG','HEIGHT_AVG','log.SEED_DRY','log.LEAFAREA_AVI','log.LEAFTHCK_AVI',
#					'log.LMALEAF_AVI','log.LDMC_AVI','log.AVG_LAMTUF','log.RGR_10',
#					'log.RGR_50','log.RGR_100','log.MORT_10','log.MORT_100')

####################
###   ARCHIVE...   #####
####################
# I USED THIS TEMPORARY SAVE TO RUN THE TRAIT NCI AT THE SAME TIME AS REGULAR NCI...
# save(tdata, file='panama_NCI_Traits_no.tNCI_10.26.15.RDA')
####################



# IN PRACTICE, I WILL SPLIT THIS UP INTO DIFFERENT R SESSIONS:
# Get tNCI for all traits of interest

for(1:length(traits)){
  trait <- traits[i]
  newname <- paste('tnci', traits[i], sep='.')
  tdata[,newname] <- unlist(lapply(1:nrow(tdata), Trait.Neigh.Fun, trait))
}

# Get uNCI for all traits of interest
for(1:length(traits)){
  trait <- traits[i]
  newname <- paste('unci', traits[i], sep='.')
  tdata[,newname] <- unlist(lapply(1:nrow(tdata), Unk.Trait.Neigh.Fun, trait))
}


##########################################################
###   APPEND tNCI and uNCI metrics from multiple R sessions...   #####
##########################################################
load('panama_NCI_Traits_no.tNCI_10.26.15.RDA')
master <- tdata
rm(tdata)
head(master, 21)

list.files('postNCI_calculations/')

### Get tNCI and uNCI by trait...

load("postNCI_calculations/panama_with_tnci.WSG_10.23.15.RDA")
tdata[19:21,]
master$tnci.wsg <- tdata$tnci.WSG
master$tnci.wsg[is.na(master$WSG)] <- NA
rm(tdata)
load("postNCI_calculations/panama_with_unci.WSG_10.23.15.RDA")
tdata[19:21,]
master$unci.wsg <- tdata$unci.WSG
master$unci.wsg[is.na(master$WSG)] <- NA
rm(tdata)


load("postNCI_calculations/panama_with_tnci.log.LDMC_AVI_10.23.15.RDA")
tdata[19:21,]
master$tnci.log.ldmc <- tdata$tnci.log.LDMC_AVI
master$tnci.log.ldmc[is.na(master$log.LDMC_AVI)] <- NA
rm(tdata)
load("postNCI_calculations/panama_with_unci.log.LDMC_AVI_10.23.15.RDA")
tdata[19:21,]
master$unci.log.ldmc <- tdata$unci.log.LDMC_AVI
master$unci.log.ldmc[is.na(master$log.LDMC_AVI)] <- NA
rm(tdata)


load("postNCI_calculations/panama_with_tnci.log.LMALEAF_AVI_10.23.15.RDA")
tdata[19:21,]
master$tnci.log.lma <- tdata$tnci.log.LMALEAF_AVI
master$tnci.log.lma[is.na(master$log.LMALEAF_AVI)] <- NA
rm(tdata)
load("postNCI_calculations/panama_with_unci.log.LMALEAF_AVI_10.23.15.RDA")
tdata[19:21,]
master$unci.log.lma <- tdata$unci.log.LMALEAF_AVI
master$unci.log.lma[is.na(master$log.LMALEAF_AVI)] <- NA
rm(tdata)


load("postNCI_calculations/panama_with_tnci.log.SEED_DRY_10.23.15.RDA")
tdata[19:21,]
master$tnci.log.seed <- tdata$tnci.log.SEED_DRY
master$tnci.log.seed[is.na(master$log.SEED_DRY)] <- NA
rm(tdata)
load("postNCI_calculations/panama_with_unci.log.SEED_DRY_10.23.15.RDA")
tdata[19:21,]
master$unci.log.seed <- tdata$unci.log.SEED_DRY
master$unci.log.seed[is.na(master$log.SEED_DRY)] <- NA
rm(tdata)


load("postNCI_calculations/panama_with_tnci.HEIGHT_AVG_10.23.15.RDA")
tdata[19:21,]
master$tnci.hmax <- tdata$tnci.HEIGHT_AVG
master$tnci.hmax[is.na(master$HEIGHT_AVG)] <- NA
rm(tdata)
load("postNCI_calculations/panama_with_unci.HEIGHT_AVG_10.23.15.RDA")
tdata[19:21,]
master$unci.hmax <- tdata$unci.HEIGHT_AVG
master$unci.hmax[is.na(master$HEIGHT_AVG)] <- NA
rm(tdata)


# check it out...
master[19:25,]

tdata <- master
# save(tdata, file="panama_NCI_Traits_tNCI_uNCI_10.26.15.RDA")





