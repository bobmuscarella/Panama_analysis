#######################################
###  GROWTH ANALYSIS  -  SINGLE PLOT ###
#######################################
library(jagsUI)

### Z-TRANSFORM DATA
z.score <- function (data) {
  xm<- mean (data, na.rm=TRUE)
  xsd<-sd(data, na.rm=TRUE)
  xtrans<-(data-xm)/(2*xsd)  
  return(xtrans)
}

### Running on PC???
pc <- FALSE

#######################################
###  START HERE WITH PROCESSED DATA ###
#######################################
if(pc==T){ 
  setwd("K:/Bob/Panama/GIT/Panama_analysis/DATA") 
} else {
  setwd("/Users/Bob/Projects/Postdoc/Panama/GIT/Panama_analysis/DATA")
}
load("Panama_AnalysisData_12.9.15.RDA")

###########################
#### Prepare data for input  ####
###########################
### REMOVE GROWTH OUTLIERS
# tdata <- tdata[tdata$Growth.Include == TRUE,]
# tdata <- tdata[tdata$Growth.Include.2 == TRUE,]
d <- tdata[tdata$Growth.Include.3 == TRUE,]

### Change names for ease
names(d)[names(d)=='WSG'] <- 'wsg'
names(d)[names(d)=='log.LDMC_AVI'] <- 'log.ldmc'
names(d)[names(d)=='log.LMALEAF_AVI'] <- 'log.lma'
names(d)[names(d)=='log.SEED_DRY'] <- 'log.seed'
names(d)[names(d)=='HEIGHT_AVG'] <- 'hmax'

# Work with one trait at a time....
#  for (i in 1:length(traits)) {
i <- 1
traits <- c('wsg','log.ldmc','log.lma','log.seed','hmax')
trait <- traits[i]

# Remove species with NA for trait value
d <- d[!is.na (d[,trait]),]  					
d <- droplevels(d)


######################################################################################
#### Standardize and Center coefficients (within plots, within size classes 10cm) ####
######################################################################################
d$growth.z <- unlist(tapply(d$growth, d$plot, scale, center=F))

### Center / scale DBH within species, within size class
d$sc <- ifelse(d$dbh <= 100, 1, 2)
d <- d[order(d$plot, d$spcode, d$sc, d$id, d$census),]

### Center / scale other variableswithin size class
for (i in 1:2){
  d$log.dbh.z[d$sc %in% i] <- unlist(tapply(d$log.dbh[d$sc %in% i], d$plot[d$sc %in% i], z.score))
  d$log.nci.z[d$sc %in% i] <- unlist(tapply(d$log.nci[d$sc %in% i], d$plot[d$sc %in% i], z.score))
  d$log.tnci.wsg.z[d$sc %in% i] <- unlist(tapply(d$log.tnci.wsg[d$sc %in% i], d$plot[d$sc %in% i], z.score))
  d$log.tnci.log.ldmc.z[d$sc %in% i] <- unlist(tapply(d$log.tnci.log.ldmc[d$sc %in% i], d$plot[d$sc %in% i], z.score))
  d$log.tnci.log.lma.z[d$sc %in% i] <- unlist(tapply(d$log.tnci.log.lma[d$sc %in% i], d$plot[d$sc %in% i], z.score))
  d$log.tnci.log.seed.z[d$sc %in% i] <- unlist(tapply(d$log.tnci.log.seed[d$sc %in% i], d$plot[d$sc %in% i], z.score))
  d$log.tnci.hmax.z[d$sc %in% i] <- unlist(tapply(d$log.tnci.hmax[d$sc %in% i], d$plot[d$sc %in% i], z.score))
}


#################################
####    PREP Continues...    ####
#################################
d$plot <- substring(d$plot, 1, 3)

# CHOOSE A SIZE CLASS TO WORK WITH...
# for(sc in 1:2) {
sc <- 1

d <- d[d$sc %in% sc,]

d <- d[,c('spcode','plot','sc','census','growth.z','log.dbh.z','id','log.nci.z', 
          trait, 
          paste('log.tnci', trait, 'z', sep='.'))]

# SAMPLE DATA FOR EXPLORATORY...
# d <- d[sample(1:nrow(d), 5000),]
# d <- droplevels(d)
# rownames(d) <- NULL

# Make a speciesxplot column for correct indexing...
d$speciesxplot <- as.factor(paste(d$plot, d$spcode, sep='.'))	

# Drop factors for correct indexing
d <- droplevels(d)						

# Create an individual ID
d$indiv <- as.numeric(as.factor(d$id))

# Order for correct indexing
d <- d[order(d$indiv, d$census, d$spcode, d$plot),]


#################################
#### Organize the input data ####
#################################
data = list (
  ntree = nrow(d),
  nindiv = length(unique(d$id)),
  nspecies = length(levels(d$speciesxplot)),
  nplot = length(levels(as.factor(d$plot))),
  growth = as.numeric(d$growth.z),
  nci = as.numeric(d[,'log.nci.z']),
#  tnci = as.numeric(d[,paste('log.tnci.', trait, '.z', sep='')]),
  dbh = as.numeric(d$log.dbh.z),
  trait = z.score(tapply(d[,trait], d$speciesxplot, mean)),
  indiv = d$indiv,
  species = as.numeric(d$speciesxplot),
  plot = as.numeric(as.factor(substring(names(tapply(d[,trait], d$speciesxplot, mean)),1,3)))
)


##############################
#### Write the model file ####
##############################
if(pc==T){
  setwd("K:/Bob/Panama/MODELS") 
} else {
  setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
}

sink("growth_3level_NCI.bug")

cat(" model {
    
    for( i in 1:ntree ) {
    
    growth[i] ~ dnorm(mu[i], tau[1])
    
    mu[i] <- beta.1[species[i]]
             + beta.2[species[i]] * nci[i]
             + beta.3[species[i]] * dbh[i]
             + indiv.effect[indiv[i]]
    }
    
    for( j in 1:nspecies ) {
      beta.1[j] ~ dnorm(mu.beta.1[plot[j]] + beta.t.1[plot[j]] * trait[j], tau[2])
      beta.2[j] ~ dnorm(mu.beta.2[plot[j]] + beta.t.2[plot[j]] * trait[j], tau[3])
      beta.3[j] ~ dnorm(mu.beta.3[plot[j]], tau[4])
    }
    
    for( p in 1:nplot ) {
      mu.beta.1[p] ~ dnorm(mu.beta[1], tau[5])
      mu.beta.2[p] ~ dnorm(mu.beta[2], tau[6])
      mu.beta.3[p] ~ dnorm(mu.beta[3], tau[7])
      beta.t.1[p] ~ dnorm(beta.t[1], tau[8])
      beta.t.2[p] ~ dnorm(beta.t[2], tau[9])
    }

    for( i.a in 1:nindiv ) {
    indiv.effect[i.a] ~ dnorm(0, tau[10])
    }
    
    for( b in 1:2 ) {
    beta.t[b] ~ dnorm(0, 1E-4)
    }
    
    for( m in 1:3 ) {
    mu.beta[m] ~ dnorm(0, 1E-4)
    }
    
    for( t in 1:10 ) {
    tau[t] ~ dgamma(1E-3, 1E-3)
    }
    
    sigma <- 1 / sqrt(tau)
    
    }"
    , fill=TRUE)
sink()

################################################
### Set initial values, monitors, iterations and run model ###
################################################

# Set initial values
inits <- function (){
  list(
    beta.t = rnorm(2),
    mu.beta = rnorm(3),
    tau = rgamma(10, 1E3, 1E3))
}

if(pc==T){ 
  setwd("K:/Bob/Panama/MODELS") 
} else {
  setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
}

# Set monitors
params <- c(paste("beta.t.",1:2,sep=''), paste("mu.beta.",1:3,sep=''),"beta.t","mu.beta","sigma")

# Run model
adapt <- 1000
iter <- 5000
burn <- 4000
thin <- 5
chains <- 3

mod <- jagsUI::jags(data, inits, params, 
                    "growth_3level_NCI.bug", 
                    n.chains=chains, n.adapt=adapt, n.iter=iter, 
                    n.burnin=burn, n.thin=thin, parallel=F)

mod

plot(mod)

plot.coeffs(mod)

for(i in 1:3) {
  if(sum(unlist(mod$Rhat) > 1.1) > 0){
    mod <- update(mod, n.iter=2000, n.thin=4)
  }
}



if(pc==T){ 
  setwd("K:/Bob/Panama/RESULTS") 
} else {
  setwd("/Users/Bob/Projects/Postdoc/Panama/RESULTS")
}

file <- paste(focal.plot, trait, 'growth_2level_NCI_NCIXDBH.RDA', sep="_")
saveRDS(mod, file=file)

}


# YOU CAN USE e.g. THIS TO READ THE OUTPUT...
assign(paste(focal.plot,trait,sep='_'), readRDS(file))










plot.coeffs <- function(mod){
  x <- cbind(unlist(mod$q50),unlist(mod$q2.5),unlist(mod$q97.5))
  x <- x[-grep('deviance',rownames(x)),]
  x <- x[-grep('sigma',rownames(x)),]
  bg <- ifelse(sign(x[,2])==sign(x[,3]), 1, 0)
  plot(x[,1],ylim=c(min(x), max(x)), pch=21, bg=bg, axes=F, ylab='Std. Effect', xlab='', cex=1.5)
  arrows(1:nrow(x), x[,2], 1:nrow(x), x[,3], angle=90, code=3, len=0.05)
  abline(h=0,lty=2)
  axis(1, labels=rownames(x), at=1:nrow(x), las=2)
  axis(2); box()
  s <- which(!names(mod$sims.list) %in% c('sigma', 'deviance'))
  n <- lapply(mod$sims.list[s], function(x) apply(x, 2, quantile, c(0.05, 0.95)))
  b <- unlist(lapply(n, function(x)x[1,]))
  t <- unlist(lapply(n, function(x)x[2,]))
  segments(1:nrow(x), b, 1:nrow(x), t, lwd=4)
  points(x[,1], pch=21, bg=bg, cex=2)
}





# For trace/density plot of a single parameter...
plot(mod$samples[,'beta.t[1]'])
plot(mod)

par(mfrow=c(3,3), mar=c(5,4,1,1))

for(i in 1:2){
  mod <- list(cocmod, shemod)[[i]]
  x <- cbind(unlist(mod$q50),unlist(mod$q2.5),unlist(mod$q97.5))
  x <- x[-nrow(x),]
  x <- x[-grep('sigma',rownames(x)),]
  bg <- ifelse(sign(x[,2])==sign(x[,3]), 1, 0)
  plot(x[,1],ylim=c(min(x), max(x)), pch=21, bg=bg, axes=F, ylab='Std. Effect', xlab='', cex=1.5)
  arrows(1:nrow(x), x[,2], 1:nrow(x), x[,3], angle=90, code=3, len=0.05)
  abline(h=0,lty=2)
  axis(1, labels=rownames(x), at=1:nrow(x), las=2)
  axis(2); box()
  n <- lapply(mod$sims.list[1:2], function(x) apply(x, 2, quantile, c(0.05, 0.95)))
  b <- unlist(lapply(n, function(x)x[1,]))
  t <- unlist(lapply(n, function(x)x[2,]))
  segments(1:nrow(x), b, 1:nrow(x), t, lwd=4)
  points(x[,1], pch=21, bg=bg, cex=2)
}





if(pc==T){
  setwd("K:/Bob/Panama/RESULTS") 
} else {
  setwd("/Users/Bob/Projects/Postdoc/Panama/RESULTS")
}

print(paste('Finished', file, 'at', Sys.time()))
save(mod, file=file)

}
}






# Set initial values
inits <- function (){
  list(
    indiv.effect = rnorm(data$nindiv, 0, 1E-3),
    beta.t = rnorm(2),
    mu.beta = rnorm(4),
    tau = rgamma(6, 1E-3, 1E-3) + 1E-5)
}

if(pc==T){ 
  setwd("K:/Bob/Panama/MODELS") 
} else {
  setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
}

# Set monitors
params <- c(paste("beta",1:4,sep='.'),"beta.t","mu.beta","sigma")
#params <- c("beta.t","mu.beta","sigma")

file <- paste(trait, focal.plot, 'RDA', sep='.')
print(paste('Started', file, 'at', Sys.time()))

# Run model
mod <- jagsUI::jags(data, inits, params, 
                    "growth_2level_NCI_NCIXDBH.bug", 
                    n.chains=3, n.adapt=1000, n.iter=5000,
                    n.burnin=1000, n.thin=5, parallel=F)







###################################
# If the models need to be run longer.... 


setwd("/Users/Bob/Projects/Postdoc/Panama/RESULTS")
save(plotmod, file='testmod.RDA')
rm(plotmod)
load('testmod.RDA')
plotmod

mod2 <- update(plotmod, n.iter=100)



