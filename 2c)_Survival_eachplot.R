##########################################
###  SURVIVAL ANALYSIS  -  SINGLE PLOT ###
##########################################
library(jagsUI)

### Z-TRANSFORM DATA
z.score <- function (data) {
  xm<- mean (data, na.rm=TRUE)
  xsd<-sd(data, na.rm=TRUE)
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

d <- d[d$census==2,]

#################################
#### Prepare data for input  ####
#################################

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

################################################################
####    Standardize and Center coefficients within plots    ####
################################################################
d <- d[order(d$plot, d$spcode, d$id, d$census),]

d$growth.z <- unlist(tapply(d$growth, d$plot, scale, center=F))
d$log.nci.z <- unlist(tapply(d$log.nci, d$plot, z.score))
d$log.dbh.z <- unlist(tapply(d$log.dbh, d$plot, z.score))

d$log.tnci.wsg.z <- unlist(tapply(d$log.tnci.wsg, d$plot, z.score))
d$log.tnci.log.ldmc.z <- unlist(tapply(d$log.tnci.log.ldmc, d$plot, z.score))
d$log.tnci.log.lma.z <- unlist(tapply(d$log.tnci.log.lma, d$plot, z.score))
d$log.tnci.log.seed.z <- unlist(tapply(d$log.tnci.log.seed, d$plot, z.score))
d$log.tnci.hmax.z <- unlist(tapply(d$log.tnci.hmax, d$plot, z.score))

d$log.unci.wsg.z <- unlist(tapply(d$log.unci.wsg, d$plot, z.score))
d$log.unci.log.ldmc.z <- unlist(tapply(d$log.unci.log.ldmc, d$plot, z.score))
d$log.unci.log.lma.z <- unlist(tapply(d$log.unci.log.lma, d$plot, z.score))
d$log.unci.log.seed.z <- unlist(tapply(d$log.unci.log.seed, d$plot, z.score))
d$log.unci.hmax.z <- unlist(tapply(d$log.unci.hmax, d$plot, z.score))

d$plot <- substring(d$plot, 1, 3)

# CHOOSE A PLOT TO WORK WITH ONE PLOT AT A TIME...
plots <- c('bci','coc','she')

p <- 2
# for(p in 1:3) {

d <- d[,c('spcode','plot','census','survival','log.dbh.z','id','log.nci.z','days',
          trait, 
          paste('log.tnci', trait, 'z', sep='.'), 
          paste('log.unci', trait, 'z', sep='.'))]

focal.plot <- plots[p]
d <- d[d$plot %in% focal.plot,]

# SAMPLE DATA FOR EXPLORATORY...
# d <- d[sample(1:nrow(d), 1000),]
# d <- droplevels(d)
# rownames(d) <- NULL

# Make a speciesxplot column for correct indexing...
d$speciesxplot <- as.factor(paste(d$plot, d$spcode, sep='.'))	

# Drop factors for correct indexing
d <- droplevels(d)						

# Create an individual ID
d$indiv <- as.numeric(as.factor(d$id))				

# Order for correct indexing
d <- d[order(d$plot, d$spcode, d$id, d$census),]


#################################
#### Organize the input data ####
#################################
data = list (
  ntree = nrow(d),
  nindiv = length(unique(d$id)),
  nspecies = length(levels(d$speciesxplot)),
  alive = d$survival,
  nci = as.numeric(d[,'log.nci.z']),
  tnci = as.numeric(d[,paste('log.tnci.', trait, '.z', sep='')]),
  days = as.numeric(d$days),
  dbh = as.numeric(d$log.dbh.z),
  trait = z.score(tapply(d[,trait], d$speciesxplot, mean)),
  indiv = d$indiv,
  species = as.numeric(d$speciesxplot)
)


##############################
#### Write the model file ####
##############################
if(pc==T){ 
  setwd("K:/Bob/Panama/MODELS") 
} else {
  setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
}

sink("survive_2level_NCI_NCIXDBH.bug")

cat(" model {
    
    for( i in 1:ntree ) {
    
    alive[i] ~ dinterval(t[i], days[i])

    t[i] ~ dweib(r, mu[i])

    mu[i] <- exp(z[i])

    z[i] <- beta.1[species[i]]
    + beta.2[species[i]] * nci[i]
    + beta.3[species[i]] * dbh[i]
    + beta.4[species[i]] * nci[i] * dbh[i]
    + indiv.effect[indiv[i]]
    }
    
    for( j in 1:nspecies ) {
    beta.1[j] ~ dnorm(mu.beta[1] + beta.t[1] * trait[j], tau[1])
    beta.2[j] ~ dnorm(mu.beta[2] + beta.t[2] * trait[j], tau[2])
    beta.3[j] ~ dnorm(mu.beta[3], tau[3])
    beta.4[j] ~ dnorm(mu.beta[4] + beta.t[3] * trait[j], tau[4])
    }
    
    for( i.a in 1:nindiv ) {
    indiv.effect[i.a] ~ dnorm(0, tau[5])
    }
    
    for( b in 1:3 ) {
    beta.t[b] ~ dnorm(0, 1E-4)
    }
    
    for( m in 1:4 ) {
    mu.beta[m] ~ dnorm(0, 1E-4)
    }
    
    for( t in 1:5 ) {
    tau[t] ~ dgamma(1E-3, 1E-3)
    }
    
    sigma <- 1 / sqrt(tau)
    r ~ dgamma(1E-3, 1E-3)

    }"
, fill=TRUE)
sink()



sink("survive_2level_NCI_TNCI_NCIXDBH.bug")

cat(" model {
    
    for( i in 1:ntree ) {
    
    alive[i] ~ dinterval(t[i], days[i])

    t[i] ~ dweib(r, mu[i])

    mu[i] <- exp(z[i])

    z[i] <- beta.1[species[i]]
    + beta.2[species[i]] * nci[i]
    + beta.3[species[i]] * tnci[i]
    + beta.4[species[i]] * dbh[i]
    + beta.5[species[i]] * nci[i] * dbh[i]
    + indiv.effect[indiv[i]]
    }
    
    for( j in 1:nspecies ) {
    beta.1[j] ~ dnorm(mu.beta[1] + beta.t[1] * trait[j], tau[1])
    beta.2[j] ~ dnorm(mu.beta[2] + beta.t[2] * trait[j], tau[2])
    beta.3[j] ~ dnorm(mu.beta[3], tau[3])
    beta.4[j] ~ dnorm(mu.beta[4], tau[4])
    beta.5[j] ~ dnorm(mu.beta[5], tau[5])
    }
    
    for( i.a in 1:nindiv ) {
    indiv.effect[i.a] ~ dnorm(0, tau[6])
    }
    
    for(b in 1:2){
    beta.t[b] ~ dnorm(0, 1E-4)
    }
    
    for(m in 1:5){
    mu.beta[m] ~ dnorm(0, 1E-4)
    }
    
    for(t in 1:6){
    tau[t] ~ dgamma(1E-3, 1E-3)
    }
    
    r ~ dgamma(1E-3, 1E-3)
    sigma <- 1 / sqrt(tau)
    
    }"
, fill=TRUE)
sink()


################################################
### Set initial values, monitors, iterations and run model ###
################################################

# Set initial values
inits <- function (eps=0.1){
  list(
    beta.t = rnorm(3),
    mu.beta = rnorm(4),
    tau = rgamma(5, 1E3, 1E3),
    r = 2,
    t = with(data, days + ifelse(alive, eps, -eps)))
}

if(pc==T){ 
  setwd("K:/Bob/Panama/MODELS") 
} else {
  setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
}

# Set monitors
# params <- c(paste("beta",1:4,sep='.'),"beta.t","mu.beta","sigma")
params <- c("beta.t","mu.beta","sigma","r")

# Run model
adapt <- 2000
iter <- 20000
burn <- 15000
thin <- 5
chains <- 4

initsmod <- jagsUI::jags(data, inits, params, 
                    "survive_2level_NCI_NCIXDBH.bug", 
                    n.chains=chains, n.adapt=adapt, n.iter=iter, 
                    n.burnin=burn, n.thin=thin, parallel=T)

start <- inits()
start$beta.t <- initsmod$q50$beta.t
start$mu.beta <- initsmod$q50$mu.beta
start$r <- initsmod$q50$r
start$tau <- 1/(initsmod$q50$sigma^2)


adapt <- 500
iter <- 5000
burn <- 2500
thin <- 5
chains <- 4

mod <- jagsUI::jags(data, list(start,start,start,start), params, 
                    "survive_2level_NCI_NCIXDBH.bug", 
                    n.chains=chains, n.adapt=adapt, n.iter=iter, 
                    n.burnin=burn, n.thin=thin, parallel=T)





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



























# For trace/density plot of a single parameter...
plot(mod$samples[,'beta.t[1]'])
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



