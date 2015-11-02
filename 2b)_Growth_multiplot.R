#######################################
###  GROWTH ANALYSIS  -  MULTI-PLOT ###
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


#################################
#### Prepare data for input  ####
#################################
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
  nplot = length(unique(d$plot)),
  nspecies = length(levels(d$speciesxplot)),
  growth = as.numeric(d$growth.z),
  nci = as.numeric(d[,'log.nci.z']),
#  tnci = as.numeric(d[,paste('log.tnci.', trait, '.z', sep='')]),
#  unci = as.numeric(d[,paste('log.unci.', trait, '.z', sep='')]),
  dbh = as.numeric(d$log.dbh.z),
  trait = z.score(tapply(d[,trait], d$speciesxplot, mean)),
  indiv = d$indiv,
  species = as.numeric(d$speciesxplot),
  plot = as.numeric(as.factor(d$plot))
)


##############################
#### Write the model file ####
##############################
if(pc==T){ 
  setwd("K:/Bob/Panama/MODELS") 
    } else {
      setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
  }


sink("growth_3level_NCI_NCIXDBH.bug")
cat(" model {
    
    for( i in 1:ntree ) {
    
    growth[i] ~ dnorm(mu[i], tau[1])
    
    mu[i] <- exp(z[i])
    
    z[i] <- beta.1[species[i]]
    + beta.2[species[i]] * (nci[i])
    + beta.3[species[i]] * (dbh[i])
    + beta.4[species[i]] * (nci[i]) * (dbh[i])
    + indiv.effect[indiv[i]]
    }
    
    for( j in 1:nspecies ) {
    beta.1[j] ~ dnorm(mu.beta[1] + beta.t.1[plot[j]] * trait[j], tau[2])
    beta.2[j] ~ dnorm(mu.beta[2] + beta.t.2[plot[j]] * trait[j], tau[3])
    beta.3[j] ~ dnorm(beta.t.3[plot[j]], tau[4])
    beta.4[j] ~ dnorm(mu.beta[3] + beta.t.4[plot[j]] * trait[j], tau[5])

#    beta.3[j] ~ dnorm(mu.beta[3], tau[4])
#    beta.4[j] ~ dnorm(mu.beta[4] + beta.t.3[plot[j]] * trait[j], tau[5])
    }
    
    for( k in 1:nplot ) {
    beta.t.1[k] ~ dnorm(beta.t[1], tau[6])
    beta.t.2[k] ~ dnorm(beta.t[2], tau[7])
    beta.t.3[k] ~ dnorm(beta.t[3], tau[8])
    beta.t.4[k] ~ dnorm(beta.t[4], tau[9])
    }
    
    for( i.a in 1:nindiv ) {
    indiv.effect[i.a] ~ dnorm(0, tau[10])
    }
    
    beta.t[1] ~ dnorm(0, 1E-4)
    beta.t[2] ~ dnorm(0, 1E-4)
    beta.t[3] ~ dnorm(0, 1E-4)
    beta.t[4] ~ dnorm(0, 1E-4)
    
    mu.beta[1] ~ dnorm(0, 1E-4)
    mu.beta[2] ~ dnorm(0, 1E-4)
    mu.beta[3] ~ dnorm(0, 1E-4)
#    mu.beta[4] ~ dnorm(0, 1E-4)

    tau[1] ~ dgamma(1E-3, 1E-3)
    tau[2] ~ dgamma(1E-3, 1E-3)
    tau[3] ~ dgamma(1E-3, 1E-3)
    tau[4] ~ dgamma(1E-3, 1E-3)
    tau[5] ~ dgamma(1E-3, 1E-3)
    tau[6] ~ dgamma(1E-3, 1E-3)
    tau[7] ~ dgamma(1E-3, 1E-3)
    tau[8] ~ dgamma(1E-3, 1E-3)
    tau[9] ~ dgamma(1E-3, 1E-3)
    tau[10] ~ dgamma(1E-3, 1E-3)
    
    sigma <- 1 / sqrt(tau)
    }
    ", fill=TRUE)
sink()





sink("growth_3level_NCI_TNCI_UNCI_NCIXDBH.bug")
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
	beta.1[j] ~ dnorm(mu.beta[1] + beta.t.1[plot[j]] * trait[j], tau[2])	# species-specific average growth
	beta.2[j] ~ dnorm(mu.beta[2] + beta.t.2[plot[j]] * trait[j], tau[3])	# speces-specific crowding effect
	beta.3[j] ~ dnorm(beta.t.3[plot[j]], tau[4])				# plot x species-specific triat-hood effect
	beta.4[j] ~ dnorm(beta.t.4[plot[j]], tau[5])				# plot x species-specific size effect
	beta.5[j] ~ dnorm(beta.t.5[plot[j]], tau[6])				# plot x species-specific effect of unk trait neighbs
	beta.6[j] ~ dnorm(beta.t.6[plot[j]], tau[7])				# plot x species-specific effect of dbh * nci interaction
    }
    
    for( k in 1:nplot ) {
	beta.t.1[k] ~ dnorm(beta.t[1], tau[8])		# plot-specific trait effect on average growth
	beta.t.2[k] ~ dnorm(beta.t[2], tau[9])		# plot-specific trait effect on sensitivity to crowding
	beta.t.3[k] ~ dnorm(beta.t[3], tau[10])		# plot-specific trait effect on sensitivity to crowding
	beta.t.4[k] ~ dnorm(beta.t[4], tau[11])		# plot-specific trait effect on sensitivity to crowding
	beta.t.5[k] ~ dnorm(beta.t[5], tau[12])		# plot-specific trait effect on sensitivity to crowding
	beta.t.6[k] ~ dnorm(beta.t[6], tau[13])		# plot-specific trait effect on sensitivity to crowding
    }
    
    for( i.a in 1:nindiv ) {
	indiv.effect[i.a] ~ dnorm(0, tau[14])
    }
    
    beta.t[1] ~ dnorm(0, 1E-4)
	beta.t[2] ~ dnorm(0, 1E-4)
	beta.t[3] ~ dnorm(0, 1E-4)
	beta.t[4] ~ dnorm(0, 1E-4)
	beta.t[5] ~ dnorm(0, 1E-4)
	beta.t[6] ~ dnorm(0, 1E-4)
    
    mu.beta[1] ~ dnorm(0, 1E-4)
    mu.beta[2] ~ dnorm(0, 1E-4)
    
    tau[1] ~ dgamma(1E-3, 1E-3)
    tau[2] ~ dgamma(1E-3, 1E-3)
    tau[3] ~ dgamma(1E-3, 1E-3)
    tau[4] ~ dgamma(1E-3, 1E-3)
    tau[5] ~ dgamma(1E-3, 1E-3)
    tau[6] ~ dgamma(1E-3, 1E-3)
    tau[7] ~ dgamma(1E-3, 1E-3)
    tau[8] ~ dgamma(1E-3, 1E-3)
	tau[9] ~ dgamma(1E-3, 1E-3)
	tau[10] ~ dgamma(1E-3, 1E-3)
	tau[11] ~ dgamma(1E-3, 1E-3)
	tau[12] ~ dgamma(1E-3, 1E-3)
	tau[13] ~ dgamma(1E-3, 1E-3)
	tau[14] ~ dgamma(1E-3, 1E-3)
    
    sigma <- 1 / sqrt(tau)
    }
", fill=TRUE)
sink()


##########################
### Set initial values ###
##########################
if(pc==T){
  setwd("K:/Bob/Panama/MODELS") 
    } else {
      setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
  }


# Set initial
inits <- function (){
  list(
    beta.t = rnorm(4),
    mu.beta = rnorm(3),
    tau = rgamma(10, 0.1, 1E-3) + 1E-10 )
}


# Set monitors
# params <- c(paste('beta',1:4,sep='.'),"beta.t","mu.beta","sigma")
params <- c(paste('beta',1:4,sep='.'), paste('beta.t',1:4,sep='.'), "beta.t", "mu.beta", "sigma")



# Run model
file <- paste(trait, 'RDA', sep='.')
print(paste('Started', file, 'at', Sys.time()))

adapt <- 500
iter <- 1000
burnin <- 500
thin <- 1

mod <- jagsUI::jags(data, inits, params, 
                    "growth_3level_NCI_NCIXDBH.bug", 
                    n.chains=3, n.adapt=adapt, n.iter=iter, 
                    n.burnin=burnin, n.thin=thin, parallel=F)

mod

# Update model
paste("Start at:", Sys.time())
mod <- update(mod, n.iter=1000)
paste("Finish at:", Sys.time())



plot(mod$samples[,c('beta.t.1[1]','beta.t.1[2]','beta.t.1[3]')])







labs <- names(unlist(bcimod$q50))
plot(unlist(bcimod$q50)[-13], pch=21, bg=1, ylim=c(-3,3), axes=F, xlab='', ylab='effect size')
axis(2)
axis(1, seq(1.1,12.1,1), labels=labs[-13], las=2)
points(seq(1.1,12.1,1),unlist(cocmod$q50)[-13], pch=21, bg=2)
points(seq(1.2,12.2,1),unlist(shermod$q50)[-13], pch=21, bg=3)
segments(1:12, unlist(bcimod$q2.5)[-13], 1:12, unlist(bcimod$q97.5)[-13], col=1)
segments(seq(1.1,12.1,1), unlist(cocmod$q2.5)[-13], seq(1.1,12.1,1), unlist(cocmod$q97.5)[-13], col=2)
segments(seq(1.2,12.2,1), unlist(shermod$q2.5)[-13], seq(1.2,12.2,1), unlist(shermod$q97.5)[-13], col=3)
abline(h=0, lty=2)


labs <- names(unlist(fullmod$q50))[-21]
pch <- ifelse(unlist(fullmod$overlap0)==T, 1, 16)
plot(unlist(fullmod$q50)[-21], pch=pch, ylim=c(-3,3), axes=F, xlab='', ylab='effect size')
axis(2)
axis(1, seq(1,20,1), labels=labs, las=2)
segments(seq(1,20,1), unlist(fullmod$q2.5)[-21], seq(1,20,1), unlist(fullmod$q97.5)[-21])
abline(h=0, lty=2)








