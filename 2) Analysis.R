library(rjags)
library(jagsUI)

#########################################
####  START HERE WITH PROCESSED DATA ####
#########################################

# if running on the PC:
setwd("K:/Bob/Panama/DATA")

# if running on my Mac
setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")

# load("alldata_NCI.RDA")
load("data_10.13.15.RDA")


################################################################
####  CURRENT CONSIDERATIONS: 
################################################################

###############
#### 1) WHAT TO DO ABOUT GROWTH OUTLIERS?
## - For now, just leave them and see if it gives problems with convergence...
## - Alternatively, select a standard deviation and say it is probably measurement error?
# tdata$Growth.Include <- abs(tdata$growth) < (sd(tdata$growth, na.rm=T) * 15)
tdata$Growth.Include[is.na(tdata$growth)] <- F

###############
####  2) WHAT TO DO ABOUT PALMS FOR GROWTH ANALYSIS?
# - For now, remove them
tdata$Growth.Include[tdata$palm==T] <- F
tdata$Growth.Include[tdata$palm==F & !is.na(tdata$growth)] <- T

###############
####  3) WHAT TO DO ABOUT NULL SPECIES CODES (MULTIPLE SPECIES)?
# - For now, remove them
tdata <- tdata[tdata$spcode != 'NULL',]

###############
####  4) How to deal with missing trait data when doing trait NCI?
# See Lasky etal 2014 PNAS


################################################################
####  MOVING ON...
################################################################

#### REMOVE EDGE TREES AND SOME OTHERWISE NA TREES
tdata <- tdata[tdata$Not.Edge == 1,] 
tdata <- tdata[!is.na(tdata$nci),]
tdata <- tdata[!is.na(tdata$survival),]
tdata <- droplevels(tdata)

#### <<<REMOVE LFDP FOR NOW>>>
tdata <- tdata[tdata$plot != 'lfdp',]

#### Z-TRANSFORM DATA
z.score <- function (data) {
	xm<- mean (data, na.rm=TRUE)
	xsd<-sd(data, na.rm=TRUE)
	xtrans<-(data-xm)/(2*xsd)
}

## Log-transform coefficients
tdata$log.nci <- log(tdata$nci)
tdata$log.tnci <- log(tdata$tnci)
tdata$log.dbh <- log(tdata$dbh)

# ## *CANCEL* Standardize coefficients
# tdata <- tdata[order(tdata$plot, tdata$spcode, tdata$id, tdata$census),]
# tdata$log.nci.z <- unlist(tapply(tdata$log.nci, tdata$plot, z.score))
# tdata$log.dbh.z <- unlist(tapply(tdata$log.dbh, tdata$plot, z.score))
# tdata$growth.z <- unlist(tapply(tdata$growth, tdata$plot, scale, center=F))

## Center coefficients within plots (without scaling)
tdata <- tdata[order(tdata$plot, tdata$spcode, tdata$id, tdata$census),]
tdata$log.nci.z <- unlist(tapply(tdata$log.nci, tdata$plot, scale, scale=F))
tdata$log.tnci.z <- unlist(tapply(tdata$log.tnci, tdata$plot, scale, scale=F))
tdata$log.dbh.z <- unlist(tapply(tdata$log.dbh, tdata$plot, scale, scale=F))
tdata$growth.z <- unlist(tapply(tdata$growth, tdata$plot, scale, scale=F))

d <- tdata
rownames(d) <- NULL
head(d)

### Generate a species-plot column for correct indexing of model...
d$speciesxplot <- as.factor(paste(d$plot, d$spcode, sep='.'))



##############################
### TEST GROWTH MODELS WITH SUBSET DATASET
d <- d[!is.na(d$growth),]
#d <- d[sample(1:nrow(d), 5000),]
d <- d[order(d$plot, d$spcode, d$id, d$census),]
d <- droplevels(d)

d <- d[! is.na (d$wd),]  	# species with NA for trait value give problems when using traits...
d <- droplevels(d)

head(d)
unique(d$plot)

d$indiv <- as.numeric(as.factor(d$id))

##############################

			data = list (
					ntree = nrow(d),
					nindiv = length(unique(d$id)),
					nspecies = length(levels(d$speciesxplot)),
					nplot = length(levels(as.factor(d$plot))),
					growth = as.numeric(d$growth.z),
					nci = as.numeric(d$log.nci.z),
					nci.t = as.numeric(d$log.tnci.z),
#					nci.u.t = as.numeric(d$log.utnci.z),
					dbh = as.numeric(d$log.dbh.z),
					trait = as.numeric(tapply(d$wd, d$speciesxplot, mean)),
					indiv = d$indiv,
					species = as.numeric(d$speciesxplot),
					plot = as.numeric(tapply(as.numeric(as.factor(d$plot)), d$speciesxplot, mean))
				)

#### SET MODEL RUN LENGTHS
setwd("K:/Bob/Panama/MODELS")
setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")


#### Build the model
sink("growth_simple.bug")
cat(" model{
	for( i in 1:ntree ) {
		growth[i] ~ dnorm(mu[i], tau[4])
		mu[i] <-  exp(z[i])
		z[i]  <-  beta.1[species[i]] + (beta.2[species[i]] * dbh[i]) + (beta.3[species[i]] * nci[i])
		}

	for( j in 1:nspecies ) {
		beta.1[j] ~ dnorm(mu.beta[1], tau[1])
		beta.2[j] ~ dnorm(mu.beta[2], tau[2])
		beta.3[j] ~ dnorm(mu.beta[3], tau[3])
		}

	mu.beta[1] ~ dnorm(0, 1E-4)
	mu.beta[2] ~ dnorm(0, 1E-4)
	mu.beta[3] ~ dnorm(0, 1E-4)

	tau[1] ~ dgamma(1E-3, 1E-3)
	tau[2] ~ dgamma(1E-3, 1E-3)	
	tau[3] ~ dgamma(1E-3, 1E-3)
	tau[4] ~ dgamma(1E-3, 1E-3)	

	sigma <- 1 / sqrt(tau)
}

",fill=TRUE)
sink()


#### Specify the model
inits <- function (){
	list(
	mu.beta = rnorm(3),
	tau = rgamma(4, 1E-3, 1E-3))
#	tau = rep(1E-3, 4))
	}

start <- list(inits(), inits(), inits())


params <- c("beta.1","beta.2","beta.3","mu.beta","tau")
setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
mod <- jagsUI::jags(data, parameters.to.save=params, model.file="growth_simple.bug", n.iter=1000, n.chains=2, parallel=T)


### Initiate the model
n.adapt=100
n.update=1000
n.iter=1000

jm <- jags.model("growth_simple.bug", inits=inits, data=data, n.chains=3, n.adapt=n.adapt)






################################
#### Build the model
sink("growth_3level_notrait.bug")

cat(" model {

	for( i in 1:ntree ) {
		growth[i] ~ dnorm(mu[i], tau[7])
		mu[i] <-  exp(z[i])
		z[i]  <-  beta.1[species[i]] + (beta.2[species[i]] * nci[i]) + (beta.3[species[i]] * dbh[i])
		}

	for( j in 1:nspecies ) {
		beta.1[j] ~ dnorm(mu.beta.1[plot[j]], tau[1])	
		beta.2[j] ~ dnorm(mu.beta.2[plot[j]], tau[2])
		beta.3[j] ~ dnorm(mu.beta.3[plot[j]], tau[3])
		}

	for( k in 1:nplot ) {
		mu.beta.1[k] ~ dnorm(mu.beta[1], tau[4])
		mu.beta.2[k] ~ dnorm(mu.beta[2], tau[5])
		mu.beta.3[k] ~ dnorm(mu.beta[3], tau[6])
		}

	mu.beta[1] ~ dnorm(0, 1E-4)
	mu.beta[2] ~ dnorm(0, 1E-4)
	mu.beta[3] ~ dnorm(0, 1E-4)

	tau[1] ~ dgamma(1E-3, 1E-3)
	tau[2] ~ dgamma(1E-3, 1E-3)	
	tau[3] ~ dgamma(1E-3, 1E-3)
	tau[4] ~ dgamma(1E-3, 1E-3)	
	tau[5] ~ dgamma(1E-3, 1E-3)
	tau[6] ~ dgamma(1E-3, 1E-3)	
	tau[7] ~ dgamma(1E-3, 1E-3)	

	sigma <- 1 / sqrt(tau)
}
",fill=TRUE)
sink()

#### Specify the model
inits <- function (){
	list(
	mu.beta = rnorm(3),
	tau = rgamma(7, 1E-3, 1E-3) + 1E-10)
	}

start <- list(inits(), inits())


params <- c("beta.1","beta.2","beta.3","mu.beta","tau")
setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
mod <- jagsUI::jags(data, parameters.to.save=params, model.file="growth_3level_notrait.bug", n.iter=1000, n.chains=2, parallel=T, inits=start)





################################
#### Build the model
sink("growth_3level_trait.bug")

cat(" model {

	for( i in 1:ntree ) {

		growth[i] ~ dnorm(mu[i], tau[1])

		mu[i] <-  exp(z[i])

		z[i]  <-  beta.1[species[i]] 									# mean performance
					+ beta.2[species[i]] * (nci[i]) 				# biomass only nci
					+ beta.3[species[i]] * (nci.t[i]) 				# trait-based nci
					+ beta.4[species[i]] * (dbh[i])				# size-specific effect
#					+ beta.5[species[i]] * (nci.u.t[i])			# to account for species with unknown trait values...
					+ indiv.effect[indiv[i]]							# to account for repeat sample of individuals
		}

	for( j in 1:nspecies ) {
		beta.1[j] ~ dnorm(mu.beta[1] + beta.t.1[plot[j]] * trait[j], tau[2])		# species-specific average growth
		beta.2[j] ~ dnorm(mu.beta[2] + beta.t.2[plot[j]] * trait[j], tau[3])		# speces-specific crowding effect
		beta.3[j] ~ dnorm(mu.beta[3], tau[4])												# species-specific triat-hood effect
		beta.4[j] ~ dnorm(mu.beta[4], tau[5])												# species-specific size effect
#		beta.5[j] ~ dnorm(mu.beta[5], tau[X])												# sp-specific effect of unk trait neighbs
		}

	for( k in 1:nplot ) {
		beta.t.1[k] ~ dnorm(beta.t[1], tau[6])			# plot-specific trait effect on average growth
		beta.t.2[k] ~ dnorm(beta.t[2], tau[7])			# plot-specific trait effect on sensitivity to crowding
		}
		
	for( i.a in 1:nindiv ) {
		indiv.effect[i.a] ~ dnorm(0, tau[8])
		}

	beta.t[1] ~ dnorm(0, 1E-4)
	beta.t[2] ~ dnorm(0, 1E-4)

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
	tau[7] ~ dgamma(1E-3, 1E-3)	
	tau[8] ~ dgamma(1E-3, 1E-3)	
#	tau[9] ~ dgamma(1E-3, 1E-3)	

	sigma <- 1 / sqrt(tau)
}
", fill=TRUE)
sink()

#### Specify the model


inits <- function (){
	list(
	beta.t = rnorm(2),
	mu.beta = rnorm(4),
	tau = rgamma(8, 1E-3, 1E-3) + 1E-5)
	}

### USE jagsUI:
params <- c("beta.t.1","beta.t.2","beta.t","mu.beta","tau")

setwd("K:/Bob/Panama/MODELS")
setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")

mod <- jagsUI::jags(data, inits, params, 
                    "growth_3level_trait.bug", 
                    n.chains=3, n.adapt=1000, n.iter=5000, 
                    n.burnin=1000, n.thin=5, parallel=TRUE)

mod

paste("Start at:", Sys.time())
mod <- update(mod, n.iter=50000)
paste("Finish at:", Sys.time())







# Auto Run
mod <- jagsUI::autojags(data, inits, params, 
										"growth_3level_trait.bug", 
										n.chains=2, n.adapt=1000, 
										parallel=T, n.thin=5)


jagsUI::traceplot(mod)
plot(mod)


### Use original 'rjags' commands:
jm <- jags.model("growth_3level_trait.bug", data=data, n.chains=1, n.adapt=n.adapt)
update(jm, n.iter = n.update)
params <- c("beta.t.1","beta.t.2","beta.t","mu.beta","tau")
coda.results <- coda.samples(jm, variable.names=params, n.iter=n.iter, n.thin = 5)

plot(coda.results, ask=T)
gelman.diag(coda.results)
res <- summary(coda.results)

plot(data$trait, res[[2]][paste('beta.2[',1:nspecies,']',sep=''),3], bg=data$plot, pch=21, ylim=c(-.1,.1))
segments(data$trait, res[[2]][paste('beta.2[',1:nspecies,']',sep=''),1], data$trait, res[[2]][paste('beta.2[',1:nspecies,']',sep=''),5], col= data$plot, lwd=2)
abline(h=0)

############################################################
############################################################
### SURVIVAL 
############################################################
############################################################
d <- tdata
d <- d[d$Not.Edge==1,]
d <- d[order(d$plot, d$spcode, d$id, d$census),]

d$speciesxplot <- as.factor(paste(d$plot, d$spcode, sep='.'))
d <- d[!is.na(d$survival),]
d <- d[! is.na (d$wd),]
d$indiv <- as.numeric(as.factor(d$id))

rownames(d) <- NULL
d <- droplevels(d)

##############################
### TEST SURVIVAL MODELS WITH SUBSET DATASET
##############################
d <- d[sample(1:nrow(d), 5000),]
d$indiv <- as.numeric(as.factor(d$id))
d <- droplevels(d)

##############################
### BUILD DATA
##############################
			data = list (
					ntree = nrow(d),
					nindiv = length(unique(d$id)),
					nspecies = length(levels(d$speciesxplot)),
					nplot = length(levels(as.factor(d$plot))),
					alive = as.numeric(d$survival),
					time = as.numeric(d$days),
					nci = as.numeric(d$log.nci.z),
					nci.t = as.numeric(d$log.tnci.z),
#					nci.u.t = as.numeric(d$log.utnci.z),
					dbh = as.numeric(d$log.dbh.z),
					trait = as.numeric(tapply(d$wd, d$speciesxplot, mean)),
					indiv = d$indiv,
					species = as.numeric(d$speciesxplot),
					plot = as.numeric(tapply(as.numeric(as.factor(d$plot)), d$speciesxplot, mean))
				)


################################
#### Build the models
################################
setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
# setwd("K:\Users\Bob\Panama\MODELS")

sink("survival_3level_trait_tnci.bug")
cat(" model {

	for( i in 1:ntree ) {

		alive[i] ~ dinterval(t[i], time[i])

		t[i] ~ dweib(r, mu[i])

		mu[i] <- exp(z[i])

		z[i]  <- beta.1 [species[i]]
					+ beta.2 [species[i]] * (nci[i])
					+ beta.3 [species[i]] * (nci.t[i])
					+ beta.4 [species[i]] * (dbh[i])
#					+ beta.5 [species[i]] * (nci.u.t[i])
					+ indiv.effect [indiv[i]]
		}

	for( j in 1:nspecies ) {
		beta.1[j] ~ dnorm(mu.beta[1] + beta.t.1[plot[j]] * trait[j], tau[1])			# species-specific average survival
		beta.2[j] ~ dnorm(mu.beta[2] + beta.t.2[plot[j]] * trait[j], tau[2])			# speces-specific crowding effect
		beta.3[j] ~ dnorm(mu.beta[3], tau[3])													# species-specific size effect
		beta.4[j] ~ dnorm(mu.beta[4], tau[4])												# species-specific trait-mediate crowding
#		beta.5[j] ~ dnorm(mu.beta[5], tau[X])												# sp-specific unk trait-neighb crowding
		}

	for( k in 1:nplot ) {
		beta.t.1[k] ~ dnorm(beta.t[1], tau[5])			# plot-specific trait effect on average survival
		beta.t.2[k] ~ dnorm(beta.t[2], tau[6])			# plot-specific trait effect on sensitivity to crowding
		}

	for( i.a in 1:nindiv ) {
		indiv.effect[i.a] ~ dnorm(0, tau[7])
		}

	r ~ dgamma(1, 1E-3)

	beta.t[1] ~ dnorm(0, 1E-4)
	beta.t[2] ~ dnorm(0, 1E-4)

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
	tau[7] ~ dgamma(1E-3, 1E-3)	

	sigma <- 1 / sqrt(tau)
}
",fill=TRUE)
sink()
################################

inits <- function (){
	eps <- 0.1
	list(beta.t = rnorm(2),
	mu.beta = rnorm(4),
	tau = rgamma(7, 0.1, 1E-3) + 1E-10,  # for some reason, jagsUI won't allow tau.inits=0...
	r = 2,
	t = with(data, time + ifelse(alive, eps, -eps)))
	}

### USE jagsUI:
params <- c("beta.t.1","beta.t.2","beta.t","mu.beta","r","tau")
setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
# setwd("K:\Users\Bob\Panama\MODELS")
mod <- jagsUI::jags(data, inits, params, "survival_3level_trait_tnci.bug", n.chains=3, n.iter=1000)














#### Specify the model
n.adapt=100
n.update=500
n.iter=50
### Initiate the model
jm <- jags.model("survival_3level_trait_tnci.bug", data=data, n.chains=3, n.adapt=n.adapt, inits=start)

# Burnin the chain
 update(jm, n.iter = n.update)

# Generate coda object
params <- c("beta.t.1","beta.t.2","beta.t","mu.beta","r","tau")
coda.results <- coda.samples(model=jm, variable.names=params, n.iter=n.iter, thin = 5)

plot(coda.results, ask=T)

gelman.diag(coda.results)

res <- summary(coda.results)

plot(data$trait, res[[2]][paste('beta.2[',1:nspecies,']',sep=''),3], bg=data$plot, pch=21, ylim=c(-.1,.1))
segments(data$trait, res[[2]][paste('beta.2[',1:nspecies,']',sep=''),1], data$trait, res[[2]][paste('beta.2[',1:nspecies,']',sep=''),5], col= data$plot, lwd=2)
abline(h=0)



setwd("/Users/Bob/Projects/Postdoc/Panama/Results")
library(rjags)
load("coda_APR2_notrait_allplots.RDA")
plot(coda.results, ask=T)
load("summary_APR2_notrait_allplots.RDA")
res







library(RColorBrewer)
col <- brewer.pal(11, 'BrBG')[c(3,9,11)]
plot((1:3), cocoli.res[[2]][,3], ylim=c(-1.25,1.25), pch=21, bg=col[1], xlim=c(0.9,3.3), cex=2, axes=F, xlab='Parameter', ylab='Standardized Effect Size')
axis(1, at=c(1.1,2.1,3.1), labels=c('Avg. growth', 'DBH effect', 'NCI effect'), las=1, cex.axis=1)
axis(2)
segments((1:3), cocoli.res[[2]][,1], (1:3), cocoli.res[[2]][,5], col=col[1], lwd=4)
segments((1:3)+.1, bci.res[[2]][,1], (1:3)+.1, bci.res[[2]][,5], col=col[2], lwd=4)
points(1:3, cocoli.res[[2]][,3], pch=21, bg=col[1], cex=2)
points((1:3)+.1, bci.res[[2]][,3], pch=21, bg=col[2], cex=2)
segments((1:3)+.2, sherman.res[[2]][,1], (1:3)+.2, sherman.res[[2]][,5], col=col[3], lwd=4)
points((1:3)+.2, sherman.res[[2]][,3], pch=21, bg=col[3], cex=2)
abline(h=0, lty=3)
box()
legend('topleft', legend=c('cocoli (dry: 1950 mm/yr)', 'bci (moist: 2500 mm/yr)', 'sherman (moist: 2700-3000 mm/yr)'), pch=21, pt.bg=col, bty='n', pt.cex=2)









library(RColorBrewer)
col <- brewer.pal(11, 'BrBG')[c(1,3,9,11)]
col <- terrain.colors(8)[c(6,3,2,1)]


plot((1:3), res[[2]][c('mu.beta.1[2]','mu.beta.2[2]','mu.beta.3[2]'),3], ylim=c(-2,1.5), pch=21, bg=col[1], xlim=c(0.9,3.3), cex=1.5, axes=F, xlab='Parameter', ylab='Standardized Effect Size')
segments((1:3), res[[2]][c('mu.beta.1[2]','mu.beta.2[2]','mu.beta.3[2]'),1], (1:3), res[[2]][c('mu.beta.1[2]','mu.beta.2[2]','mu.beta.3[2]'),5], col=col[1], lwd=4)
points((1:3), res[[2]][c('mu.beta.1[2]','mu.beta.2[2]','mu.beta.3[2]'),3], pch=21, bg=col[1], cex=1.5)

segments((1:3)+.1, res[[2]][c('mu.beta.1[1]','mu.beta.2[1]','mu.beta.3[1]'),1], (1:3)+.1, res[[2]][c('mu.beta.1[1]','mu.beta.2[1]','mu.beta.3[1]'),5], col=col[2], lwd=4)
points((1:3)+.1, res[[2]][c('mu.beta.1[1]','mu.beta.2[1]','mu.beta.3[1]'),3], pch=21, bg=col[2], cex=1.5)

segments((1:3)+.2, res[[2]][c('mu.beta.1[4]','mu.beta.2[4]','mu.beta.3[4]'),1], (1:3)+.2, res[[2]][c('mu.beta.1[4]','mu.beta.2[4]','mu.beta.3[4]'),5], col=col[3], lwd=4)
points((1:3)+.2, res[[2]][c('mu.beta.1[4]','mu.beta.2[4]','mu.beta.3[4]'),3], pch=21, bg=col[3], cex=1.5)

segments((1:3)+.3, res[[2]][c('mu.beta.1[3]','mu.beta.2[3]','mu.beta.3[3]'),1], (1:3)+.3, res[[2]][c('mu.beta.1[3]','mu.beta.2[3]','mu.beta.3[3]'),5], col=col[4], lwd=4)
points((1:3)+.3, res[[2]][c('mu.beta.1[3]','mu.beta.2[3]','mu.beta.3[3]'),3], pch=21, bg=col[4], cex=1.5)

abline(h=0, lty=3)
box()

axis(1, at=1:3+.15, labels=c('Avg. growth', 'NCI effect', 'DBH effect'), las=1, cex.axis=1)
axis(2)

legend('topleft', legend=c('cocoli (dry: 1950 mm/yr)', 'bci (moist: 2500 mm/yr)', 'sherman (moist: 2850 mm/yr)','luquillo (wet: 3500 mm/yr)'), pch=21, pt.bg=col, bty='n', pt.cex=2)








setwd("/Users/Bob/Projects/Postdoc/Panama/Results")
library(rjags)
load("coda_APR4_trait_allplots.RDA")
#plot(coda.results, ask=T)
load("summary_APR4_trait_allplots.RDA")
res


library(RColorBrewer)
col <- brewer.pal(11, 'BrBG')[c(1,3,9,11)]
col <- terrain.colors(8)[c(6,3,2,1)]


plot((1:2), res[[2]][c('beta.t.1[2]','beta.t.2[2]'),3], ylim=c(-4,4), pch=21, bg=col[1], xlim=c(0.9,2.3), cex=1.5, axes=F, xlab='Parameter', ylab='Standardized Effect Size')
segments((1:2), res[[2]][c('beta.t.1[2]','beta.t.2[2]'),1], (1:2), res[[2]][c('beta.t.1[2]','beta.t.2[2]'),5], col=col[1], lwd=4)
points((1:2), res[[2]][c('beta.t.1[2]','beta.t.2[2]'),3], pch=21, bg=col[1], cex=1.5)

segments((1:2)+.1, res[[2]][c('beta.t.1[1]','beta.t.2[1]'),1], (1:2)+.1, res[[2]][c('beta.t.1[1]','beta.t.2[1]'),5], col=col[2], lwd=4)
points((1:2)+.1, res[[2]][c('beta.t.1[1]','beta.t.2[1]'),3], pch=21, bg=col[2], cex=1.5)

segments((1:2)+.2, res[[2]][c('beta.t.1[4]','beta.t.2[4]'),1], (1:2)+.2, res[[2]][c('beta.t.1[4]','beta.t.2[4]'),5], col=col[3], lwd=4)
points((1:2)+.2, res[[2]][c('beta.t.1[4]','beta.t.2[4]'),3], pch=21, bg=col[3], cex=1.5)

segments((1:2)+.3, res[[2]][c('beta.t.1[3]','beta.t.2[3]'),1], (1:2)+.3, res[[2]][c('beta.t.1[3]','beta.t.2[3]'),5], col=col[4], lwd=4)
points((1:2)+.3, res[[2]][c('beta.t.1[4]','beta.t.2[4]'),3], pch=21, bg=col[4], cex=1.5)

abline(h=0, lty=3)
box()

axis(1, at=1:2+.15, labels=c('Trait effect on avg. growth', 'Trait effect on NCI effect'), las=1, cex.axis=1)
axis(2)

legend('topleft', legend=c('cocoli (dry: 1950 mm/yr)', 'bci (moist: 2500 mm/yr)', 'sherman (moist: 2850 mm/yr)','luquillo (wet: 3500 mm/yr)'), pch=21, pt.bg=col, bty='n', pt.cex=2)




library(RColorBrewer)
col <- brewer.pal(11, 'BrBG')[c(1,3,9,11)]
col <- terrain.colors(8)[c(6,3,2,1)]


plot((1:2), res[[2]][c('beta.t.1[2]','beta.t.2[2]'),3], ylim=c(-4,4), pch=21, bg=col[1], xlim=c(0.9,3.3), cex=1.5, axes=F, xlab='Parameter', ylab='Standardized Effect Size')
segments((1:2), res[[2]][c('beta.t.1[2]','beta.t.2[2]'),1], (1:2), res[[2]][c('beta.t.1[2]','beta.t.2[2]'),5], col=col[1], lwd=4)
points((1:2), res[[2]][c('beta.t.1[2]','beta.t.2[2]'),3], pch=21, bg=col[1], cex=1.5)

segments((1:2)+.1, res[[2]][c('beta.t.1[1]','beta.t.2[1]'),1], (1:2)+.1, res[[2]][c('beta.t.1[1]','beta.t.2[1]'),5], col=col[2], lwd=4)
points((1:2)+.1, res[[2]][c('beta.t.1[1]','beta.t.2[1]'),3], pch=21, bg=col[2], cex=1.5)

segments((1:2)+.2, res[[2]][c('beta.t.1[4]','beta.t.2[4]'),1], (1:2)+.2, res[[2]][c('beta.t.1[4]','beta.t.2[4]'),5], col=col[3], lwd=4)
points((1:2)+.2, res[[2]][c('beta.t.1[4]','beta.t.2[4]'),3], pch=21, bg=col[3], cex=1.5)

segments((1:2)+.3, res[[2]][c('beta.t.1[3]','beta.t.2[3]'),1], (1:2)+.3, res[[2]][c('beta.t.1[3]','beta.t.2[3]'),5], col=col[4], lwd=4)
points((1:2)+.3, res[[2]][c('beta.t.1[4]','beta.t.2[4]'),3], pch=21, bg=col[4], cex=1.5)

segments(3+c(.1,.2,.3), res[[2]][c('mu.beta[1]','mu.beta[2]','mu.beta[3]'),1], 3+c(.1,.2,.3), res[[2]][c('mu.beta[1]','mu.beta[2]','mu.beta[3]'),5], col=1, lwd=4)
points(3+c(.1,.2,.3), res[[2]][c('mu.beta[1]','mu.beta[2]','mu.beta[3]'),3], pch=21, bg='white', cex=1.5)

abline(h=0, lty=3)
box()
axis(1, at=1:2+.15, labels=c('Trait effect on avg. growth', 'Trait effect on NCI effect'), las=1, cex.axis=1)
axis(2)

text(3+c(.1,.2,.3), -1, labels=c('intercept for avg. growth', 'intercept for NCI effect', 'intercept for DBh effect'), srt=90, adj=1)

legend('topleft', legend=c('cocoli (dry: 1950 mm/yr)', 'bci (moist: 2500 mm/yr)', 'sherman (moist: 2850 mm/yr)','luquillo (wet: 3500 mm/yr)'), pch=21, pt.bg=col, bty='n', pt.cex=2)




=======

	sigma <- 1 / sqrt(tau)
}

",fill=TRUE)
sink()


#### Specify the model
inits <- function (){
	list(
	mu.beta = rnorm(3),
	tau = rgamma(4, 1E-3, 1E-3))
#	tau = rep(1E-3, 4))
	}

start <- list(inits(), inits(), inits())


params <- c("beta.1","beta.2","beta.3","mu.beta","tau")
setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
mod <- jagsUI::jags(data, parameters.to.save=params, model.file="growth_simple.bug", n.iter=1000, n.chains=2, parallel=T)


### Initiate the model
n.adapt=100
n.update=1000
n.iter=1000

jm <- jags.model("growth_simple.bug", inits=inits, data=data, n.chains=3, n.adapt=n.adapt)






################################
#### Build the model
sink("growth_3level_notrait.bug")

cat(" model {

	for( i in 1:ntree ) {
		growth[i] ~ dnorm(mu[i], tau[7])
		mu[i] <-  exp(z[i])
		z[i]  <-  beta.1[species[i]] + (beta.2[species[i]] * nci[i]) + (beta.3[species[i]] * dbh[i])
		}

	for( j in 1:nspecies ) {
		beta.1[j] ~ dnorm(mu.beta.1[plot[j]], tau[1])	
		beta.2[j] ~ dnorm(mu.beta.2[plot[j]], tau[2])
		beta.3[j] ~ dnorm(mu.beta.3[plot[j]], tau[3])
		}

	for( k in 1:nplot ) {
		mu.beta.1[k] ~ dnorm(mu.beta[1], tau[4])
		mu.beta.2[k] ~ dnorm(mu.beta[2], tau[5])
		mu.beta.3[k] ~ dnorm(mu.beta[3], tau[6])
		}

	mu.beta[1] ~ dnorm(0, 1E-4)
	mu.beta[2] ~ dnorm(0, 1E-4)
	mu.beta[3] ~ dnorm(0, 1E-4)

	tau[1] ~ dgamma(1E-3, 1E-3)
	tau[2] ~ dgamma(1E-3, 1E-3)	
	tau[3] ~ dgamma(1E-3, 1E-3)
	tau[4] ~ dgamma(1E-3, 1E-3)	
	tau[5] ~ dgamma(1E-3, 1E-3)
	tau[6] ~ dgamma(1E-3, 1E-3)	
	tau[7] ~ dgamma(1E-3, 1E-3)	

	sigma <- 1 / sqrt(tau)
}
",fill=TRUE)
sink()

#### Specify the model
inits <- function (){
	list(
	mu.beta = rnorm(3),
	tau = rgamma(7, 1E-3, 1E-3) + 1E-10)
	}

start <- list(inits(), inits())


params <- c("beta.1","beta.2","beta.3","mu.beta","tau")
setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
mod <- jagsUI::jags(data, parameters.to.save=params, model.file="growth_3level_notrait.bug", n.iter=1000, n.chains=2, parallel=T, inits=start)





################################
#### Build the model
sink("growth_3level_trait.bug")

cat(" model {

	for( i in 1:ntree ) {

		growth[i] ~ dnorm(mu[i], tau[1])

		mu[i] <-  exp(z[i])

		z[i]  <-  beta.1[species[i]] 									# mean performance
					+ beta.2[species[i]] * (nci[i]) 				# biomass only nci
					+ beta.3[species[i]] * (nci.t[i]) 				# trait-based nci
					+ beta.4[species[i]] * (dbh[i])				# size-specific effect
#					+ beta.5[species[i]] * (nci.u.t[i])			# to account for species with unknown trait values...
					+ indiv.effect[indiv[i]]							# to account for repeat sample of individuals
		}

	for( j in 1:nspecies ) {
		beta.1[j] ~ dnorm(mu.beta[1] + beta.t.1[plot[j]] * trait[j], tau[2])		# species-specific average growth
		beta.2[j] ~ dnorm(mu.beta[2] + beta.t.2[plot[j]] * trait[j], tau[3])		# speces-specific crowding effect
		beta.3[j] ~ dnorm(mu.beta[3], tau[4])												# species-specific triat-hood effect
		beta.4[j] ~ dnorm(mu.beta[4], tau[5])												# species-specific size effect
#		beta.5[j] ~ dnorm(mu.beta[5], tau[X])												# sp-specific effect of unk trait neighbs
		}

	for( k in 1:nplot ) {
		beta.t.1[k] ~ dnorm(beta.t[1], tau[6])			# plot-specific trait effect on average growth
		beta.t.2[k] ~ dnorm(beta.t[2], tau[7])			# plot-specific trait effect on sensitivity to crowding
		}
		
	for( i.a in 1:nindiv ) {
		indiv.effect[i.a] ~ dnorm(0, tau[8])
		}

	beta.t[1] ~ dnorm(0, 1E-4)
	beta.t[2] ~ dnorm(0, 1E-4)

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
	tau[7] ~ dgamma(1E-3, 1E-3)	
	tau[8] ~ dgamma(1E-3, 1E-3)	
#	tau[9] ~ dgamma(1E-3, 1E-3)	

	sigma <- 1 / sqrt(tau)
}
", fill=TRUE)
sink()

#### Specify the model


inits <- function (){
	list(
	beta.t = rnorm(2),
	mu.beta = rnorm(4),
	tau = rgamma(8, 1E-3, 1E-3) + 1E-5)
	}

### USE jagsUI:
params <- c("beta.t.1","beta.t.2","beta.t","mu.beta","tau")
setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")

mod <- jagsUI::jags(data, inits, params, "growth_3level_trait.bug", n.chains=2, n.iter=1000, parallel=T)

update(mod, n.iter=5000)





# Auto Run
mod <- jagsUI::autojags(data, inits, params, 
										"growth_3level_trait.bug", 
										n.chains=2, n.adapt=1000, 
										parallel=T, n.thin=5)


jagsUI::traceplot(mod)
plot(mod)


### Use original 'rjags' commands:
jm <- jags.model("growth_3level_trait.bug", data=data, n.chains=1, n.adapt=n.adapt)
update(jm, n.iter = n.update)
params <- c("beta.t.1","beta.t.2","beta.t","mu.beta","tau")
coda.results <- coda.samples(jm, variable.names=params, n.iter=n.iter, n.thin = 5)

plot(coda.results, ask=T)
gelman.diag(coda.results)
res <- summary(coda.results)

plot(data$trait, res[[2]][paste('beta.2[',1:nspecies,']',sep=''),3], bg=data$plot, pch=21, ylim=c(-.1,.1))
segments(data$trait, res[[2]][paste('beta.2[',1:nspecies,']',sep=''),1], data$trait, res[[2]][paste('beta.2[',1:nspecies,']',sep=''),5], col= data$plot, lwd=2)
abline(h=0)

############################################################
############################################################
### SURVIVAL 
############################################################
############################################################

d <- tdata
d <- d[order(d$plot, d$spcode, d$id, d$census),]
rownames(d) <- NULL

### Generate a species-plot column for correct indexing of model...
d$speciesxplot <- as.factor(paste(d$plot, d$spcode, sep='.'))


##############################
### TEST SURVIVAL MODELS WITH SUBSET DATASET
##############################

d <- d[!is.na(d$survival),]
d <- d[sample(1:nrow(d), 5000),]
d <- d[order(d$plot, d$spcode, d$id, d$census),]
d <- droplevels(d)
d <- d[! is.na (d$wd),]
d <- droplevels(d)
d$indiv <- as.numeric(as.factor(d$id))

d <- droplevels(d)
table(d$survival)


			data = list (
					ntree = nrow(d),
					nindiv = length(unique(d$id)),
					nspecies = length(levels(d$speciesxplot)),
					nplot = length(levels(as.factor(d$plot))),
					alive = as.numeric(d$survival),
					time = as.numeric(d$days),
					nci = as.numeric(d$log.nci.z),
					nci.t = as.numeric(d$log.tnci.z),
#					nci.u.t = as.numeric(d$log.utnci.z),
					dbh = as.numeric(d$log.dbh.z),
					trait = as.numeric(tapply(d$wd, d$speciesxplot, mean)),
					indiv = d$indiv,
					species = as.numeric(d$speciesxplot),
					plot = as.numeric(tapply(as.numeric(as.factor(d$plot)), d$speciesxplot, mean))
				)

################################
#### Build the model

setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")

# setwd("K:\Users\Bob\Panama\MODELS")

sink("survival_3level_trait_tnci.bug")
cat(" model {

	for( i in 1:ntree ) {

		alive[i] ~ dinterval(t[i], time[i])

		t[i] ~ dweib(r, mu[i])

		mu[i] <- exp(z[i])

		z[i]  <- beta.1 [species[i]]
					+ beta.2 [species[i]] * (nci[i])
					+ beta.3 [species[i]] * (nci.t[i])
					+ beta.4 [species[i]] * (dbh[i])
#					+ beta.5 [species[i]] * (nci.u.t[i])
					+ indiv.effect [indiv[i]]
		}

	for( j in 1:nspecies ) {
		beta.1[j] ~ dnorm(mu.beta[1] + beta.t.1[plot[j]] * trait[j], tau[1])			# species-specific average survival
		beta.2[j] ~ dnorm(mu.beta[2] + beta.t.2[plot[j]] * trait[j], tau[2])			# speces-specific crowding effect
		beta.3[j] ~ dnorm(mu.beta[3], tau[3])													# species-specific size effect
		beta.4[j] ~ dnorm(mu.beta[4], tau[4])												# species-specific trait-mediate crowding
#		beta.5[j] ~ dnorm(mu.beta[5], tau[X])												# sp-specific unk trait-neighb crowding
		}

	for( k in 1:nplot ) {
		beta.t.1[k] ~ dnorm(beta.t[1], tau[5])			# plot-specific trait effect on average survival
		beta.t.2[k] ~ dnorm(beta.t[2], tau[6])			# plot-specific trait effect on sensitivity to crowding
		}

	for( i.a in 1:nindiv ) {
		indiv.effect[i.a] ~ dnorm(0, tau[7])
		}

	r ~ dgamma(1, 1E-3)

	beta.t[1] ~ dnorm(0, 1E-4)
	beta.t[2] ~ dnorm(0, 1E-4)

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
	tau[7] ~ dgamma(1E-3, 1E-3)	

	sigma <- 1 / sqrt(tau)
}
",fill=TRUE)
sink()
################################


inits <- function (){
	eps <- 0.1
	list(beta.t = rnorm(2),
	mu.beta = rnorm(4),
	tau = rgamma(7, 0.1, 1E-3) + 1E-10,
	r = 2,
	t = with(data, time + ifelse(alive, eps, -eps)))
	}
inits()


### USE jagsUI:
params <- c("beta.t.1","beta.t.2","beta.t","mu.beta","r","tau")
setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
mod <- jagsUI::jags(data, inits, params, "survival_3level_trait_tnci.bug", n.chains=3, n.iter=1000)











#### Specify the model
n.adapt=100
n.update=500
n.iter=50
### Initiate the model
jm <- jags.model("survival_3level_trait_tnci.bug", data=data, n.chains=3, n.adapt=n.adapt, inits=start)

# Burnin the chain
 update(jm, n.iter = n.update)

# Generate coda object
params <- c("beta.t.1","beta.t.2","beta.t","mu.beta","r","tau")
coda.results <- coda.samples(model=jm, variable.names=params, n.iter=n.iter, thin = 5)

plot(coda.results, ask=T)

gelman.diag(coda.results)

res <- summary(coda.results)


plot(data$trait, res[[2]][paste('beta.2[',1:nspecies,']',sep=''),3], bg=data$plot, pch=21, ylim=c(-.1,.1))
segments(data$trait, res[[2]][paste('beta.2[',1:nspecies,']',sep=''),1], data$trait, res[[2]][paste('beta.2[',1:nspecies,']',sep=''),5], col= data$plot, lwd=2)
abline(h=0)
