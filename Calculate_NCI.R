library(data.table)
library(parallel)
library(snow)

setwd("K:/Bob/Panama/GIT/Panama_analysis/DATA")

#######################################
###  START HERE WITH PROCESSED DATA ###
#######################################
load("panama_traits_preNCI_6.1.16.RDA")

tdata <- as.data.table(tdata)

Neigh.Fun <- function(i, tdata, r=15){   
  # If stem is not on the edge
  if(tdata$Not.Edge[i] == 1){
    
    # Focal stem in row i
    foc.tree <- tdata[i,] 
    
    # Gets neighboring stems within r meters and from the same census year
    neighz <- tdata[(1:nrow(tdata))!=i 
                    & sqrt((tdata$x - foc.tree$x)^2 + (tdata$y - foc.tree$y)^2) <= r
                    & tdata$census == foc.tree$census 
                    & tdata$plot == foc.tree$plot,]
    
    neighz <- as.data.frame(neighz)
    dist2 <- (neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2    
    
    # If some of the neighbors have same coords as focal.tree, add an offset
    if(0 %in% dist2){
      neighz[dist2==0 ,3:4] <- neighz[dist2==0, 3:4] + 0.25
      neighz[dist2==0 ,3:4] <- neighz[dist2==0, 3:4] + 0.25
    }
    
    dist2 <- (neighz$x - foc.tree$x)^2 + (neighz$y - foc.tree$y)^2    
    diam2 <- (neighz$dbh/10)^2
    
    cneighz <- neighz[neighz$spcode %in% foc.tree$spcode,]
    cdist2 <- (cneighz$x - foc.tree$x)^2 + (cneighz$y - foc.tree$y)^2    
    cdiam2 <- (cneighz$dbh/10)^2
    
    # Calculate NCI using DBH^2 and dist^-2
    All.NCI <- sum(diam2/dist2)
    Con.NCI <- sum(cdiam2/cdist2)
    
    return(as.data.frame(cbind(All.NCI, Con.NCI)))    
    # For edge stems   
  } else {
    return(as.data.frame(cbind(All.NCI=NA, Con.NCI=NA)))
  }
}

# TEST THE FUNCTION
samp <- sample(which(tdata$Not.Edge==T), 10)
s <- Sys.time()
nci <- do.call('rbind', lapply(samp, Neigh.Fun, tdata))
etime <- Sys.time() - s

((nrow(tdata)/length(samp)) * etime)/60/60 # Should take ~4 hours in sequential...

## NOT PARALLEL
# nci <- do.call('rbind', lapply(1:nrow(tdata), Neigh.Fun, tdata))

### PARALLEL
cl <- makeCluster(8, type='SOCK')

### SUBSET TO TEST...
# nci <- parallel::parLapply(cl, 1:5000, Neigh.Fun, tdata)

### FULL DATASET
nci <- parLapply(cl, 1:nrow(tdata), Neigh.Fun, tdata)
nci <- do.call('rbind', nci)
stopCluster(cl)


save(nci, file='NCI_for_panama_traits_preNCI_6.1.16.RDA')










dim(tdata)
dim(nci)


# ADD 
setwd("K:/Bob/Panama/DATA")
load("panama_NCI_Traits_12.9.15.RDA")
oldnci <- tdata$nci


setwd("K:/Bob/Panama/DATA") 
tdata <- readRDS("Panama_AnalysisData_3.8.16.RDS")
newnci <- tdata$nci

oldnci[10:20]
newnci[10:20]

which( round(newnci,2) %in% round(46616.11,2))

head(oldnci[!is.na(oldnci)])



plot(log(tdata$Big.NCI+1), log(nci[,1]+1))
abline(0,1)

tdata <- cbind(tdata, nci)




tdata <- cbind(tdata, nci)

# save(tdata, file="panama_NCI_Traits_12.9.15.RDA")















