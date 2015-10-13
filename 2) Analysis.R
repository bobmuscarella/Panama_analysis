library(jagsUI)

#########################################
####  START HERE WITH PROCESSED DATA ####
#########################################

# if running on the PC:
setwd("K:\Users\Bob\Panama\DATA")

# if running on my Mac
setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")

# load("alldata_NCI.RDA")
load("data_10.13.15.RDA")




#############################
####  THINGS TO CONSIDER: 

####  1) How to deal with palms for growth analysis?

####  2) Whether / how to remove growth outliers?

####  3) Check species codes - some are NULL, etc.
# tdata <- tdata[tdata$spcode != 'NULL',]

####  4) Repeat measures of same individual in single plot...
# Incorporate individual effect...

####  5) How to deal with missing trait data when doing trait NCI?
# See Lasky etal 2014 PNAS

#############################

# REMOVE EDGE TREES
tdata <- tdata[tdata$Not.Edge == 1,] 
tdata <- droplevels(tdata)

# <<<REMOVE LFDP FOR NOW>>>
tdata <- tdata[tdata$plot != 'lfdp',] 
tdata <- droplevels(tdata)

# REMOVE QUESTIONABLE GROWTH OUTLIERS 
# (all obs > 15 sd of observed growth for now... need to adjust this later...)
tdata$Growth.Include <- (tdata$growth <= sd(tdata$growth, na.rm=T) * 15 & tdata$growth >= - sd(tdata$growth, na.rm=T) * 15)


#### Z-TRANSFORM DATA
z.score <- function (data) {
	xm<- mean (data, na.rm=TRUE)
	xsd<-sd(data, na.rm=TRUE)
	xtrans<-(data-xm)/(2*xsd)
}

tdata <- tdata[!is.na(tdata$nci),]
tdata <- tdata[!is.na(tdata$survival),]
tdata <- droplevels(tdata)

## Log-transform coefficients
tdata$log.nci <- log(tdata$nci)
tdata$log.tnci <- log(tdata$nci)
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

# Generate a species-plot column for correct indexing of model...
d$speciesxplot <- as.factor(paste(d$plot, d$spcode, sep='.'))



##############################
### TEST GROWTH MODELS WITH SUBSET DATASET
d <- d[!is.na(d$growth),]
d <- d[sample(1:nrow(d), 1000),]
d <- d[order(d$plot, d$spcode, d$id, d$census),]
d <- droplevels(d)

d <- d[! is.na (d$wd),]  	# species with NA for trait value give problems when using traits...
d <- droplevels(d)

head(d)
unique(d$plot)

##############################

ntree <- nrow(d)
nspecies <- length(levels(d$speciesxplot))
nplot <- length(levels(as.factor(d$plot)))

			data = list (
					ntree = ntree,
					nspecies = nspecies,
					nplot = nplot,
					growth = as.numeric(d$growth.z),
					survive = as.numeric(d$survival),
					nci = as.numeric(d$log.nci.z),
					tnci = as.numeric(d$log.tnci.z),
					dbh = as.numeric(d$log.dbh.z),
					trait = as.numeric(tapply(d$wd, d$speciesxplot, mean)),
					species = as.numeric(d$speciesxplot),
					plot = as.numeric(tapply(as.numeric(as.factor(d$plot)), d$speciesxplot, mean))
#					plot = as.numeric(as.factor(d$plot))
				)

#### SET MODEL RUN LENGTHS
setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
library(rjags)

#### Build the model
sink("growth_simple.bug")
cat(" model{
	for( i in 1:ntree ) {		growth[i] ~ dnorm(mu[i], tau[4])
		mu[i] <-  exp(z[i])

		z[i]  <-  beta.1[species[i]]
					+ beta.2[species[i]] * (dbh[i])
					+ beta.3[species[i]] * (nci[i])		}	for( j in 1:nspecies ) {		beta.1[j] ~ dnorm(mu.beta[1], tau[1])
		beta.2[j] ~ dnorm(mu.beta[2], tau[2])
		beta.3[j] ~ dnorm(mu.beta[3], tau[3])
		}	mu.beta[1] ~ dnorm(0, 1E-4)	mu.beta[2] ~ dnorm(0, 1E-4)	mu.beta[3] ~ dnorm(0, 1E-4)	tau[1] ~ dgamma(1E-3, 1E-3)	tau[2] ~ dgamma(1E-3, 1E-3)		tau[3] ~ dgamma(1E-3, 1E-3)
	tau[4] ~ dgamma(1E-3, 1E-3)	
	sigma <- 1 / sqrt(tau)}

",fill=TRUE)
sink()


################################
#### Build the model
sink("growth_3level_notrait.bug")

cat(" model {

	for( i in 1:ntree ) {
		growth[i] ~ dnorm(mu[i], tau[7])
		mu[i] <-  exp(z[i])
		z[i]  <-  beta.1[species[i]] + beta.2[species[i]] * (nci[i]) + beta.3[species[i]] * (dbh[i])
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

################################
#### Build the model
sink("growth_3level_trait.bug")

cat(" model {

	for( i in 1:ntree ) {

		growth[i] ~ dnorm(mu[i], tau[7])

		mu[i] <-  exp(z[i])

		z[i]  <-  beta.1[species[i]] 
					+ beta.2[species[i]] * (nci[i]) 
					+ beta.3[species[i]] * (tnci[i]) 
					+ beta.4[species[i]] * (dbh[i])
		}

	for( j in 1:nspecies ) {
		beta.1[j] ~ dnorm(mu.beta[1] + beta.t.1[plot[j]] * trait[j], tau[1])			# species-specific average growth
		beta.2[j] ~ dnorm(mu.beta[2] + beta.t.2[plot[j]] * trait[j], tau[2])			# speces-specific crowding effect
		beta.3[j] ~ dnorm(mu.beta[3], tau[3])												# species-specific triat-hood effect
		beta.4[j] ~ dnorm(mu.beta[4], tau[4])												# species-specific size effect
		}

	for( k in 1:nplot ) {
		beta.t.1[k] ~ dnorm(beta.t[1], tau[5])			# plot-specific trait effect on average growth
		beta.t.2[k] ~ dnorm(beta.t[2], tau[6])			# plot-specific trait effect on sensitivity to crowding
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

	sigma <- 1 / sqrt(tau)
}
", fill=TRUE)
sink()
################################

#### Specify the model
n.adapt=100
n.update=1000
n.iter=1000

inits <- function (){
	list(
	beta.t = rnorm(2),
	mu.beta = rnorm(4),
	tau = rgamma(7, 1E-3, 1E-3))
	}


### USE jagsUI:
library(jagsUI)
params <- c("beta.t.1","beta.t.2","beta.t","mu.beta","r","tau")
setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
mod <- jagsUI::jags(data, inits, params, "growth_3level_trait.bug", n.chains=1, n.iter=100)



### Initiate the model
jm <- jags.model("growth_3level_trait.bug", data=data, n.chains=1, n.adapt=n.adapt)



# Burnin the chain
update(jm, n.iter = n.update)

# Generate coda object
#params <- c("beta.t.1","beta.t.2","beta.t.3")
params <- c("beta.t.1","beta.t.2","beta.1","beta.2")
#params <- c("beta.t.1","beta.t.2")
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

##############################
### TEST SURVIVAL MODELS WITH SUBSET DATASET
d <- tdata[!is.na(tdata$wd),]
rownames(d) <- NULL
head(d)

# Generate a species-plot column for correct indexing of model...
d$speciesxplot <- as.factor(paste(d$plot, d$spcode, sep='.'))

d <- d[sample(1:nrow(d), 2000),]
d <- d[order(d$plot, d$spcode, d$id, d$census),]

# EXCLUDE SINGLETON SPECIES (TROUBLESHOOTING MODEL...)
##d <- d[d$spcode %in% names(table(d$spcode))[table(d$spcode) > 1], ]

unique(d$plot)


d <- droplevels(d)
table(d$survival)
##############################

ntree <- nrow(d)
nspecies <- length(levels(d$speciesxplot))
nplot <- length(levels(as.factor(d$plot)))

			data = list (
					ntree = ntree,
					nspecies = nspecies,
					nplot = nplot,
#					indiv = as.numeric(d$id),
					alive = as.numeric(d$survival),
					time = as.numeric(d$days),
#					time = as.numeric(scale(d$days)),
					nci = as.numeric(d$log.nci.z),
					tnci = as.numeric(d$log.tnci.z),
					dbh = as.numeric(d$log.dbh.z),
					trait = as.numeric(tapply(d$wd, d$speciesxplot, mean)),
					species = as.numeric(d$speciesxplot),
					plot = as.numeric(tapply(as.numeric(as.factor(d$plot)), d$speciesxplot, mean))
				)

################################
#### Build the model

setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
setwd("K:\Users\Bob\Panama\MODELS")

sink("survival_3level_trait_tnci.bug")
cat(" model {

	for( i in 1:ntree ) {

		alive[i] ~ dinterval(t[i], time[i])

		t[i] ~ dweib(r, mu[i])

		mu[i] <- exp(z[i])

		z[i]  <- beta.1[species[i]]
					+ beta.2[species[i]] * (nci[i])
					+ beta.3[species[i]] * (tnci[i])
					+ beta.4[species[i]] * (dbh[i])
		}

	for( j in 1:nspecies ) {
		beta.1[j] ~ dnorm(mu.beta[1] + beta.t.1[plot[j]] * trait[j], tau[1])			# species-specific average survival
		beta.2[j] ~ dnorm(mu.beta[2] + beta.t.2[plot[j]] * trait[j], tau[2])			# speces-specific crowding effect
		beta.3[j] ~ dnorm(mu.beta[3], tau[3])													# species-specific size effect
		beta.4[j] ~ dnorm(mu.beta[4], tau[4])												# species-specific trait-mediate crowding
		}

	for( k in 1:nplot ) {
		beta.t.1[k] ~ dnorm(beta.t[1], tau[5])			# plot-specific trait effect on average survival
		beta.t.2[k] ~ dnorm(beta.t[2], tau[6])			# plot-specific trait effect on sensitivity to crowding
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

	sigma <- 1 / sqrt(tau)
}
",fill=TRUE)
sink()
################################


#### Specify the model
n.adapt=100
n.update=500
n.iter=50


inits <- function (){
	eps <- 0.1
	list(beta.t = rnorm(2),
	mu.beta = rnorm(4),
	tau = rgamma(6, 0.1, 1E-3),
	r = 2,
	t = with(data, time + ifelse(alive, eps, -eps)))
	}
inits()


### Initiate the model
jm <- jags.model("survival_3level_trait_tnci.bug", data=data, n.chains=3, n.adapt=n.adapt, inits=start)




### USE jagsUI:
library(jagsUI)
params <- c("beta.t.1","beta.t.2","beta.t","mu.beta","r","tau")
setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
mod <- jagsUI::jags(data, inits, params, "survival_3level_trait_tnci.bug", n.chains=3, n.iter=100)


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




