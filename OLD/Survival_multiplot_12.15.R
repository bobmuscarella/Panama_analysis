##########################################
###  SURVIVAL ANALYSIS  -  MULTI PLOT ###
##########################################
library(jagsUI)

### Z-TRANSFORM DATA
z.score <- function (data, center=T) {
  xm<- mean (data, na.rm=TRUE)
  xsd<-sd(data, na.rm=TRUE)
  if(center==F){xm <- 0}
  xtrans<-(data-xm)/(2*xsd)  
}

### Running on PC???
pc <- T

#######################################
###  START HERE WITH PROCESSED DATA ###
#######################################
if(pc==T){ 
  setwd("K:/Bob/Panama/GIT/Panama_analysis/DATA") 
} else {
  setwd("/Users/Bob/Projects/Postdoc/Panama/GIT/Panama_analysis/DATA")
}

load("Panama_AnalysisData_12.9.15.RDA")
d <- tdata

###########################
#### Prepare data for input  ####
###########################
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
##### Ontogenetic size class: Center / scale DBH within species, within size class #####
# size.classes <- tapply(d$dbh, paste(d$spcode, d$plot), quantile, 0.5)
# d$size.class <- ifelse(d$dbh <= size.classes[match(paste(d$spcode, d$plot), names(size.classes))], 1, 2)
# sp.size.class <- paste(d$spcode, d$size.class, sep=".")
# for (i in 1:length(unique(sp.size.class))){
#   tmp <- sort(unique(sp.size.class))[i]
#   d$log.dbh.z[sp.size.class %in% tmp] <- z.score(d$log.dbh[sp.size.class %in% tmp])
# }
# d$log.dbh.z[is.na(d$log.dbh.z)] <- 0

###### Generic size class: Center / scale DBH within species, within size class #####
cutoff <- (-100)  # Setting cutoff to neg value will put everything in 1 size class.
d$size.class <- ifelse(d$dbh <= cutoff, 1, 2)

d <- d[order(d$plot, d$spcode, d$size.class, d$id, d$census),]

### Center / scale other variableswithin size class
for (i in 1:2){
  d$log.dbh.z[d$size.class %in% i] <- z.score(d$log.dbh[d$size.class %in% i])
  d$log.nci.z[d$size.class %in% i] <- z.score(d$log.nci[d$size.class %in% i])
  d$log.tnci.wsg.z[d$size.class %in% i] <- z.score(d$log.tnci.wsg[d$size.class %in% i])
  d$log.tnci.log.ldmc.z[d$size.class %in% i] <- z.score(d$log.tnci.log.ldmc[d$size.class %in% i])
  d$log.tnci.log.lma.z[d$size.class %in% i] <- z.score(d$log.tnci.log.lma[d$size.class %in% i])
  d$log.tnci.log.seed.z[d$size.class %in% i] <- z.score(d$log.tnci.log.seed[d$size.class %in% i])
  d$log.tnci.hmax.z[d$size.class %in% i] <- z.score(d$log.tnci.hmax[d$size.class %in% i])
}

### SET DIAM.z OF SINGLETON SPECIES (size classes) TO ZERO
d$log.dbh.z[is.na(d$log.dbh.z)] <- 0

#################################
####    PREP Continues...    ####
#################################
# CHOOSE A SIZE CLASS TO WORK WITH...
# for(sc in 1:2) {
#sc <- 1
#d <- d[d$size.class %in% sc,]

d <- d[,c('spcode','plot','census','survival','log.dbh.z','id','log.nci.z','days', 
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
  alive = as.numeric(d$survival),
  days = as.numeric(d$days),
  ncensus =  length(unique(paste(d$census, d$plot, sep=''))),
  census = as.numeric(as.factor(paste(d$census, d$plot, sep=''))),  
  nci = as.numeric(d[,'log.nci.z']),
  tnci = as.numeric(d[,paste('log.tnci.', trait, '.z', sep='')]),
  dbh = as.numeric(d$log.dbh.z),
  trait = z.score(tapply(d[,trait], d$speciesxplot, mean)),
  indiv = d$indiv,
  species = as.numeric(d$speciesxplot),
  nplot = length(levels(as.factor(d$plot))),
  plot = as.numeric(as.factor(substring(names(tapply(d[,trait], d$speciesxplot, mean)),1,3)))
)

### Add an indicator to set individual effect of non-rep indiv to zero
repindiv <- names(table(data$indiv))[table(data$indiv)>1]
data$indicator <- as.numeric(data$indiv %in% repindiv)



##############################
#### Write the model file ####
##############################
if(pc==T){
  setwd("K:/Bob/Panama/MODELS") 
} else {
  setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
}

sink("survive_3level_NCI_NCIXDBH.bug")

cat(" model {
    
    for( i in 1:ntree ) {
    
    alive[i] ~ dbern(t[i])
    
    t[i] <- pow(z[i], days[i]/365.25)
    
    logit(z[i]) <- beta.1[species[i]]
    + beta.2[species[i]] * nci[i]
    + beta.3[species[i]] * dbh[i]
    + beta.4[species[i]] * nci[i] * dbh[i]
    + indiv.effect[indiv[i]] * indicator[i]
    }
    
    for( j in 1:nspecies ) {
    beta.1[j] ~ dnorm(mu.beta.1[plot[j]] + beta.t.1[plot[j]] * trait[j], tau[1])
    beta.2[j] ~ dnorm(mu.beta.2[plot[j]] + beta.t.2[plot[j]] * trait[j], tau[2])
    beta.3[j] ~ dnorm(mu.beta.3[plot[j]] + beta.t.3[plot[j]] * trait[j], tau[3])
    beta.4[j] ~ dnorm(mu.beta.4[plot[j]] + beta.t.4[plot[j]] * trait[j], tau[4])
    }

    for( p in 1:nplot ) {
    mu.beta.1[p] ~ dnorm(0, 1E-4)
    mu.beta.2[p] ~ dnorm(0, 1E-4)
    mu.beta.3[p] ~ dnorm(0, 1E-4)
    mu.beta.4[p] ~ dnorm(0, 1E-4)
    beta.t.1[p] ~ dnorm(0, 1E-4)
    beta.t.2[p] ~ dnorm(0, 1E-4)
    beta.t.3[p] ~ dnorm(0, 1E-4)
    beta.t.4[p] ~ dnorm(0, 1E-4)
    }

    for( i.a in 1:nindiv ) {
    indiv.effect[i.a] ~ dnorm(0, tau[5])
    }
    
    for( t in 1:5 ) {
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
    beta.t.1 = rnorm(3),
    beta.t.2 = rnorm(3),
    beta.t.3 = rnorm(3),
    beta.t.4 = rnorm(3),
    mu.beta.1 = rnorm(3),
    mu.beta.2 = rnorm(3),
    mu.beta.3 = rnorm(3),
    mu.beta.4 = rnorm(3),
    tau = rgamma(5, 1E3, 1E3))
}

if(pc==T){ 
  setwd("K:/Bob/Panama/MODELS") 
} else {
  setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
}

# Set monitors
params <- c(paste("beta.t.",1:4,sep=''), 
            paste("mu.beta.",1:4,sep=''),
            "sigma")

# Run model
adapt <- 1000
iter <- 10000
burn <- 8000
thin <- 5
chains <- 3

Sys.time()
smod <- jagsUI::jags(data, inits, params, 
                     "survive_3level_NCI_NCIXDBH.bug", 
                     n.chains=chains, n.adapt=adapt, n.iter=iter, 
                     n.burnin=burn, n.thin=thin, parallel=T)

smod

plot.coeffs(smod, col=c(1,2,4))
abline(v=seq(3.5,(3.5*7), by=3))

smod <- update(smod, n.iter=1000)





for(i in 1:3) {
  if(sum(unlist(mod$Rhat) > 1.1) > 0){
    mod <- update(mod, n.iter=1000)
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





plot.coeffs <- function(mod, col=1, sigma=F){
  x <- cbind(unlist(mod$q50),unlist(mod$q2.5),unlist(mod$q97.5))
  x <- x[-grep('deviance',rownames(x)),]
  if(sigma==F) {x <- x[-grep('sigma',rownames(x)),]}
  bg <- ifelse(sign(x[,2])==sign(x[,3]), col, 0)
  plot(x[,1], ylim=c(min(x), max(x)), pch=21, col=col, bg=bg, axes=F, ylab='Std. Effect', xlab='', cex=1.5)
  arrows(1:nrow(x), x[,2], 1:nrow(x), x[,3], angle=90, code=3, len=0.05, col=col)
  abline(h=0,lty=2)
  axis(1, labels=rownames(x), at=1:nrow(x), las=2)
  axis(2); box()
  if(sigma==F) {s <- which(!names(mod$sims.list) %in% c('sigma', 'deviance'))}
  if(sigma==T) {s <- which(!names(mod$sims.list) %in% c('deviance'))}
  n <- lapply(mod$sims.list[s], function(x) apply(x, 2, quantile, c(0.05, 0.95)))
  b <- unlist(lapply(n, function(x)x[1,]))
  t <- unlist(lapply(n, function(x)x[2,]))
  segments(1:nrow(x), b, 1:nrow(x), t, lwd=4, col=col)
  points(x[,1], pch=21, bg=bg, cex=2, col=col)
}






















# For trace/density plot of a single parameter...
plot(initsmod$samples[,c('r','beta.t[1]')])
plot(mod)

par(mfrow=c(3,3), mar=c(5,4,1,1))

for(i in 1:2){
  mod <- list(cocmod, shemod)[[i]]
  x <- cbind(unlist(mod$q50),unlist(mod$q2.5),unlist(mod$q97.5))
  x <- x[-nrow(x),]
  x <- x[-grep('sigma',rownames(x)),]
  pch <- ifelse(sign(x[,2])==sign(x[,3]), 16, 1)
  plot(x[,1],ylim=c(min(x), max(x)), pch=pch, axes=F, ylab='Std. Effect', xlab='', cex=1.5)
  segments(1:nrow(x), x[,2], 1:nrow(x), x[,3])
  abline(h=0,lty=2)
  axis(1, labels=rownames(x), at=1:nrow(x), las=2)
  axis(2); box()
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



