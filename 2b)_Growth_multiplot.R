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


#################################
#### Prepare data for input  ####
#################################
rownames(d) <- NULL

### Remove growth NA
d <- d[!is.na(d$growth),]
d <- d[d$Growth.Include.3 == TRUE,]


### REMOVE SO MANY ZERO GROWTHS
# d <- d[d$growth != 0,]

### SAMPLE THE BCI DATA FOR SPEED
# bci <- d[d$plot %in% 'bci',]
# bci <- bci[sample(1:nrow(bci), 8000),]
# d <- d[!d$plot %in% 'bci',]
# d <- rbind(bci, d)
# table(d$plot)


### Change names for ease
names(d)[names(d)=='WSG'] <- 'wsg'
names(d)[names(d)=='log.LDMC_AVI'] <- 'log.ldmc'
names(d)[names(d)=='log.LMALEAF_AVI'] <- 'log.lma'
names(d)[names(d)=='log.SEED_DRY'] <- 'log.seed'
names(d)[names(d)=='HEIGHT_AVG'] <- 'hmax'

### Work with one trait at a time....
traits <- c('wsg','log.ldmc','log.lma','log.seed','hmax')

d <- d[,c('spcode','plot','census','growth.z','log.dbh.z','id','log.nci.z', 
          traits, 
          paste('log.tnci', traits, 'z', sep='.'), 
          paste('log.unci', traits, 'z', sep='.'))]


#	for (i in 1:length(traits)) {
i=1

trait <- traits[i]

### Remove species with NA for trait value
d <- d[!is.na (d[,trait]),]

# SAMPLE DATA FOR EXPLORATORY...
# d <- d[sample(1:nrow(d), 5000),]
# d <- droplevels(d)
# rownames(d) <- NULL

### Make a speciesxplot column for correct indexing...
d$speciesxplot <- as.factor(paste(d$plot, d$spcode, sep='.'))	

### Drop factors for correct indexing
d <- droplevels(d)						

### Create an individual ID
d$indiv <- as.numeric(as.factor(d$id))				

### Order for correct indexing
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

sink("growth_3level_NCI_NCIXDBH.bug")
cat(" model {
    
    for( i in 1:ntree ) {
    
    growth[i] ~ dnorm(mu[i], tau[1])
#    mu[i] <- (z[i])
    mu[i] <- exp(z[i])
    
    z[i] <- beta.1[species[i]]
    + beta.2[species[i]] * (nci[i])
    + beta.3[species[i]] * (dbh[i])
    + beta.4[species[i]] * (nci[i] * dbh[i])
    + indiv.effect[indiv[i]]
    }
    
    for( j in 1:nspecies ) {
    beta.1[j] ~ dnorm(mu.beta.1[plot[j]] + beta.t.1[plot[j]] * trait[j], tau[2])
    beta.2[j] ~ dnorm(mu.beta.2[plot[j]] + beta.t.2[plot[j]] * trait[j], tau[3])
    beta.3[j] ~ dnorm(mu.beta.3[plot[j]], tau[4])
    beta.4[j] ~ dnorm(mu.beta.4[plot[j]] + beta.t.3[plot[j]] * trait[j], tau[5])
    }
    
    for( k in 1:nplot ) {
    beta.t.1[k] ~ dnorm(beta.t[1], tau[6])
    beta.t.2[k] ~ dnorm(beta.t[2], tau[7])
    beta.t.3[k] ~ dnorm(beta.t[3], tau[8])
    mu.beta.1[k] ~ dnorm(mu.beta[1], tau[10])
    mu.beta.2[k] ~ dnorm(mu.beta[2], tau[11])
    mu.beta.3[k] ~ dnorm(mu.beta[3], tau[12])
    mu.beta.4[k] ~ dnorm(mu.beta[4], tau[12])
    }
    
    for( i.a in 1:nindiv ) {
    indiv.effect[i.a] ~ dnorm(0, tau[13])
    }
    
    for( b in 1:3 ) {
    beta.t[b] ~ dnorm(0, 1E-4)
    }

    for( m in 1:4 ) {
    mu.beta[m] ~ dnorm(0, 1E-4)
    }

    for ( t in 1:13 ) {
    tau[t] ~ dgamma(1E3, 1E3)
    }
    
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
    
    for( b in 1:6 ) {
    beta.t[b] ~ dnorm(0, 1E-4)
    }

    mu.beta[1] ~ dnorm(0, 1E-4)
    mu.beta[2] ~ dnorm(0, 1E-4)
    
    for( t in 1:14 ) {
    tau[t] ~ dgamma(1E3, 1E3)
    }

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
    beta.t = rnorm(3),
    mu.beta = rnorm(4),
    tau = rgamma(13, 1E3, 1E3))
}

# Set monitors
params <- c(paste('beta.t',1:4,sep='.'),"beta.t","mu.beta", "sigma")

# Run model
file <- paste('allplots', trait, 'RDA', sep='.')
print(paste('Started', file, 'at', Sys.time()))

adapt <- 1000
iter <- 5000
burnin <- 2500
thin <- 5
nchains <- 4

mod <- jagsUI::jags(data, inits, params,
                    "growth_3level_NCI_NCIXDBH.bug", 
                    n.chains=nchains, n.adapt=adapt, n.iter=iter, 
                    n.burnin=burnin, n.thin=thin, parallel=T)

print(paste('Finished', file, 'at', Sys.time()))

mod

# Update model
paste("Start at:", Sys.time())
mod <- update(mod, n.iter=5000)
paste("Finish at:", Sys.time())


dev <- length(unlist(mod$q50))
nvar <- length(unlist(mod$q2.5)[-dev])
labs <- names(unlist(mod$q50))[-dev]
pch <- ifelse(unlist(mod$overlap0)==T, 1, 16)
plot(unlist(mod$q50)[-dev], pch=pch, 
     ylim=c(-3,3), axes=F, xlab='', 
     ylab='effect size')
axis(2)
axis(1, seq(1,nvar,1), labels=labs, las=2)
segments(seq(1,nvar,1), 
         unlist(mod$q2.5)[-dev], 
         seq(1,nvar,1), 
         unlist(mod$q97.5)[-dev])
abline(h=0, lty=2)









plot(testmod$samples[,c('beta.t[1]','beta.t[2]','beta.t[3]')])







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









