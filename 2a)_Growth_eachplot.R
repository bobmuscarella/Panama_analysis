#######################################
###  GROWTH ANALYSIS  -  SINGLE PLOT ###
#######################################
library(jagsUI)

### Z-TRANSFORM DATA
z.score <- function (data) {
  xm<- mean (data, na.rm=TRUE)
  xsd<-sd(data, na.rm=TRUE)
  xtrans<-(data-xm)/(2*xsd)	
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
rownames(tdata) <- NULL

# Remove growth NA
tdata <- tdata[!is.na(tdata$growth),]
tdata <- tdata[tdata$Growth.Include.2 == TRUE,]

# Change names for ease
names(tdata)[names(tdata)=='WSG'] <- 'wsg'
names(tdata)[names(tdata)=='log.LDMC_AVI'] <- 'log.ldmc'
names(tdata)[names(tdata)=='log.LMALEAF_AVI'] <- 'log.lma'
names(tdata)[names(tdata)=='log.SEED_DRY'] <- 'log.seed'
names(tdata)[names(tdata)=='HEIGHT_AVG'] <- 'hmax'

# Work with one trait at a time....
traits <- c('wsg','log.ldmc','log.lma','log.seed','hmax')

tdata <- tdata[,c('spcode','plot','census','growth.z','log.dbh.z','id','log.nci.z', 
          traits, 
          paste('log.tnci', traits, 'z', sep='.'), 
          paste('log.unci', traits, 'z', sep='.'))]

tdata$plot <- substring(tdata$plot, 1, 3)
     
# CHOOSE A PLOT TO WORK WITH ONE PLOT AT A TIME...
plots <- c('bci','coc','she')

for(p in 1:3) {
  
focal.plot <- plots[p]
d <- tdata[tdata$plot %in% focal.plot,]

#	for (i in 1:length(traits)) {
i=1

trait <- traits[i]

# Remove species with NA for trait value
d <- d2[!is.na (d2[,trait]),]  					

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
  growth = as.numeric(d$growth.z),
  nci = as.numeric(d[,'log.nci.z']),
  #	tnci = as.numeric(d[,paste('log.tnci.', trait, '.z', sep='')]),
  #	unci = as.numeric(d[,paste('log.unci.', trait, '.z', sep='')]),
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

sink("growth_2level_NCI_NCIXDBH.bug")

cat(" model {
    
    for( i in 1:ntree ) {
    
    growth[i] ~ dnorm(mu[i], tau[1])
    mu[i] <- exp(z[i])
    
    z[i] <- beta.1[species[i]] 				          # grand mean of performance
    + beta.2[species[i]] * (nci[i]) 		        # biomass nci
    + beta.3[species[i]] * (dbh[i])			        # size-specific effect
    + beta.4[species[i]] * (nci[i]) * (dbh[i])	# interaction between size and crowding
    + indiv.effect[indiv[i]]  		              # to account for repeat sample of individuals
    }
    
    for( j in 1:nspecies ) {
    beta.1[j] ~ dnorm(mu.beta[1] + beta.t[1] * trait[j], tau[2])	# species-specific average growth
    beta.2[j] ~ dnorm(mu.beta[2] + beta.t[2] * trait[j], tau[3])	# speces-specific crowding effect
    beta.3[j] ~ dnorm(mu.beta[3], tau[4])				                  # plot x species-specific size effect
    beta.4[j] ~ dnorm(mu.beta[4] + beta.t[3] * trait[j], tau[5])	# plot x species-specific effect of dbh * nci interaction
    }
    
    for( i.a in 1:nindiv ) {
    indiv.effect[i.a] ~ dnorm(0, tau[6])
    }
    
    for( b in 1:3 ) {
    beta.t[b] ~ dnorm(0, 1E-4)
    }
    
    for( m in 1:4 ) {
    mu.beta[m] ~ dnorm(0, 1E-4)
    }
    
    for( t in 1:6 ) {
    tau[t] ~ dgamma(1E-3, 1E-3)
    }
    
    sigma <- 1 / sqrt(tau)
    
    }"
, fill=TRUE)
sink()



sink("growth_2level_NCI_TNCI_UNCI_NCIXDBH.bug")

cat(" model {
    
    for( i in 1:ntree ) {
    
    growth[i] ~ dnorm(mu[i], tau[1])
    
    mu[i] <- exp(z[i])
    
    z[i] <- indiv.effect[indiv[i]]  		# to account for repeat sample of individuals
    + beta.1[species[i]] 				# grand mean of performance
    + beta.2[species[i]] * (nci[i]) 		# biomass nci
    + beta.3[species[i]] * (tnci[i]) 		# absolute trait difference nci
    + beta.4[species[i]] * (dbh[i])			# size-specific effect
    + beta.5[species[i]] * (unci[i])		# biomass nci of neighbors with unknown traits
    + beta.6[species[i]] * (nci[i]) * (dbh[i])	# interaction between size and crowding
    }
    
    for( j in 1:nspecies ) {
    beta.1[j] ~ dnorm(mu.beta[1] + beta.t[1] * trait[j], tau[2])	# species-specific average growth
    beta.2[j] ~ dnorm(mu.beta[2] + beta.t[2] * trait[j], tau[3])	# speces-specific crowding effect
    beta.3[j] ~ dnorm(mu.beta[3], tau[4])				# plot x species-specific triat-hood effect
    beta.4[j] ~ dnorm(mu.beta[4], tau[5])				# plot x species-specific size effect
    beta.5[j] ~ dnorm(mu.beta[5], tau[6])				# plot x species-specific effect of unk trait neighbs
    beta.6[j] ~ dnorm(mu.beta[6] + beta.t[3] * trait[j], tau[7])				# plot x species-specific effect of dbh * nci interaction
    }
    
    for( i.a in 1:nindiv ) {
    indiv.effect[i.a] ~ dnorm(0, tau[8])
    }
    
    beta.t[1] ~ dnorm(0, 1E-4)
    beta.t[2] ~ dnorm(0, 1E-4)
    beta.t[3] ~ dnorm(0, 1E-4)
    
    mu.beta[1] ~ dnorm(0, 1E-4)
    mu.beta[2] ~ dnorm(0, 1E-4)
    mu.beta[3] ~ dnorm(0, 1E-4)
    mu.beta[4] ~ dnorm(0, 1E-4)
    mu.beta[5] ~ dnorm(0, 1E-4)
    mu.beta[6] ~ dnorm(0, 1E-4)
    
    tau[1] ~ dgamma(1E-3, 1E-3)
    tau[2] ~ dgamma(1E-3, 1E-3)
    tau[3] ~ dgamma(1E-3, 1E-3)
    tau[4] ~ dgamma(1E-3, 1E-3)
    tau[5] ~ dgamma(1E-3, 1E-3)
    tau[6] ~ dgamma(1E-3, 1E-3)
    tau[7] ~ dgamma(1E-3, 1E-3)
    tau[8] ~ dgamma(1E-3, 1E-3)
    
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
    beta.t = rnorm(3),
    mu.beta = rnorm(4),
    tau = rgamma(6, 1E-3, 1E-3) + 1E-5)
}

if(pc==T){ 
  setwd("K:/Bob/Panama/MODELS") 
} else {
  setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
}

# Set monitors
# params <- c(paste("beta",1:4,sep='.'),"beta.t","mu.beta","sigma")
params <- c("beta.t","mu.beta","sigma")

# Run model
adapt <- 2000
iter <- 5000
burn <- 3000
thin <- 5
chains <- 3

mod <- jagsUI::jags(data, inits, params, 
                   "growth_2level_NCI_NCIXDBH.bug", 
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


