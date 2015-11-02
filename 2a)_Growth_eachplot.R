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
	setwd("K:/Bob/Panama/DATA") 
	} else {
		setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")
		}

load("Panama_AnalysisData_10.30.15.RDA")
d <- tdata

###########################
#### Prepare data for input  ####
###########################
rownames(d) <- NULL

# Remove growth NA
d <- d[!is.na(d$growth),]
d <- d[d$Growth.Include.2 == TRUE,]

# Change names for ease
names(d)[names(d)=='WSG'] <- 'wsg'
names(d)[names(d)=='log.LDMC_AVI'] <- 'log.ldmc'
names(d)[names(d)=='log.LMALEAF_AVI'] <- 'log.lma'
names(d)[names(d)=='log.SEED_DRY'] <- 'log.seed'
names(d)[names(d)=='HEIGHT_AVG'] <- 'hmax'

# Work with one trait at a time....
traits <- c('wsg','log.ldmc','log.lma','log.seed','hmax')

d <- d[,c('spcode','plot','census','growth.z','log.dbh.z','id','log.nci.z', 
				traits, 
				paste('log.tnci', traits, 'z', sep='.'), 
				paste('log.unci', traits, 'z', sep='.'))]

# CHOOSE A PLOT TO WORK WITH ONE PLOT AT A TIME...
plots <- c('cocoli', 'bci', 'sherman')

#for(p in 1) {
p=1

  focal.plot <- plots[p]
	d <- d[d$plot %in% focal.plot,]

#	for (i in 1:length(traits)) {
i=1
		
    trait <- traits[i]
		
		# Remove species with NA for trait value
		d <- d[!is.na (d[,trait]),]  					


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
	tnci = as.numeric(d[,paste('log.tnci.', trait, '.z', sep='')]),
	unci = as.numeric(d[,paste('log.unci.', trait, '.z', sep='')]),
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

sink("growth_2level_NCI_TNCI_UNCI_NCIXDBH.bug")

cat(" model {
    
    for( i in 1:ntree ) {
    
	growth[i] ~ dnorm(mu[i], tau[1])
    
	mu[i] <- exp(z[i])
    
	z[i] <- indiv.effect[indiv[i]]			# to account for repeat sample of individuals
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
    
    z[i] <- beta.1[species[i]] 				  # grand mean of performance
    + beta.2[species[i]] * (nci[i]) 		# biomass nci
    + beta.3[species[i]] * (dbh[i])			# size-specific effect
    + beta.4[species[i]] * (nci[i]) * (dbh[i])	# interaction between size and crowding
    + indiv.effect[indiv[i]]  		      # to account for repeat sample of individuals
    }
    
    for( j in 1:nspecies ) {
    beta.1[j] ~ dnorm(mu.beta[1] + beta.t[1] * trait[j], tau[2])	# species-specific average growth
    beta.2[j] ~ dnorm(mu.beta[2] + beta.t[2] * trait[j], tau[3])	# speces-specific crowding effect
    beta.3[j] ~ dnorm(mu.beta[3], tau[4])				# plot x species-specific size effect
    beta.4[j] ~ dnorm(mu.beta[4] + beta.t[3] * trait[j], tau[5])				# plot x species-specific effect of dbh * nci interaction
    }

    for( i.a in 1:nindiv ) {
    indiv.effect[i.a] ~ dnorm(0, tau[6])
    }
    
    beta.t[1] ~ dnorm(0, 1E-4)
    beta.t[2] ~ dnorm(0, 1E-4)
    beta.t[3] ~ dnorm(0, 1E-4)
    
    mu.beta[1] ~ dnorm(0, 1E-4)
    mu.beta[2] ~ dnorm(0, 1E-4)
    mu.beta[3] ~ dnorm(0, 1E-4)
    mu.beta[4] ~ dnorm(0, 1E-4)
    
    tau[1] ~ dgamma(1E-3, 1E-3)
    tau[2] ~ dgamma(1E-3, 1E-3)
    tau[3] ~ dgamma(1E-3, 1E-3)
    tau[4] ~ dgamma(1E-3, 1E-3)
    tau[5] ~ dgamma(1E-3, 1E-3)
    tau[6] ~ dgamma(1E-3, 1E-3)
    
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
	mu.beta = rnorm(6),
	tau = rgamma(8, 1E-3, 1E-3) + 1E-5)
	}

if(pc==T){ 
	setwd("K:/Bob/Panama/MODELS") 
	} else {
		setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
		}


# Set monitors
# params <- c(paste("beta",1:4,sep='.'),"beta.t","mu.beta","sigma")
params <- c("beta.t","mu.beta","sigma")

file <- paste(trait, focal.plot, 'RDA', sep='.')
print(paste('Started', file, 'at', Sys.time()))

# Run model
mod <- jagsUI::jags(data, inits, params, 
                    "growth_2level_NCI_TNCI_UNCI_NCIXDBH.bug", 
                    n.chains=3, n.adapt=5000, n.iter=20000, 
                    n.burnin=10000, n.thin=50, parallel=F)


# For trace/density plot of a single parameter...
plot(mod$samples[,'beta.t[1]'])
plot(mod)



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
mod

mod <- update(mod, n.iter=10000)

mod

plot(mod)

str(mod)





###################################
# If the models need to be run longer.... 


setwd("/Users/Bob/Projects/Postdoc/Panama/RESULTS")
save(plotmod, file='testmod.RDA')
rm(plotmod)
load('testmod.RDA')
plotmod

mod2 <- update(plotmod, n.iter=100)

