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

### Remove growth NA
d <- d[!is.na(d$growth),]
d <- d[d$Growth.Include.2 == TRUE,]

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

sink("growth_3level_Ptaus_NCI_NCIXDBH.bug")
cat(" model {
    
    for( i in 1:ntree ) {
    
    growth[i] ~ dnorm(mu[i], tau1[1])
    mu[i] <- exp(z[i])
    
    z[i] <- beta.1[species[i]]
    + beta.2[species[i]] * (nci[i])
    + beta.3[species[i]] * (dbh[i])
    + beta.4[species[i]] * (nci[i] * dbh[i])
    + indiv.effect[indiv[i]]
    }
    
    for( j in 1:nspecies ) {
    beta.1[j] ~ dnorm(mu.beta.1[plot[j]] + beta.t.1[plot[j]] * trait[j], tau1[2])
    beta.2[j] ~ dnorm(mu.beta.2[plot[j]] + beta.t.2[plot[j]] * trait[j], tau1[3])
    beta.3[j] ~ dnorm(mu.beta.3[plot[j]], tau1[4])
    beta.4[j] ~ dnorm(mu.beta.4[plot[j]] + beta.t.3[plot[j]] * trait[j], tau1[5])
    }
    
    for( k in 1:nplot ) {
    beta.t.1[k] ~ dnorm(beta.t[1], tau2[k])
    beta.t.2[k] ~ dnorm(beta.t[2], tau3[k])
    beta.t.3[k] ~ dnorm(beta.t[3], tau4[k])
    mu.beta.1[k] ~ dnorm(mu.beta[1], tau5[k])
    mu.beta.2[k] ~ dnorm(mu.beta[2], tau6[k])
    mu.beta.3[k] ~ dnorm(mu.beta[3], tau7[k])
    mu.beta.4[k] ~ dnorm(mu.beta[4], tau8[k])
    }

    for( p in 1:nplot ) {
    tau2[p] ~ dgamma(1E-3, 1E-3)
    tau3[p] ~ dgamma(1E-3, 1E-3)
    tau4[p] ~ dgamma(1E-3, 1E-3)
    tau5[p] ~ dgamma(1E-3, 1E-3)
    tau6[p] ~ dgamma(1E-3, 1E-3)
    tau7[p] ~ dgamma(1E-3, 1E-3)
    tau8[p] ~ dgamma(1E-3, 1E-3)
    }
    
    for( i.a in 1:nindiv ) {
    indiv.effect[i.a] ~ dnorm(0, tau1[6])
    }

    for( b in 1:3 ) {
    beta.t[b] ~ dnorm(0, 1E-4)
    }

    for( m in 1:4 ) {
    mu.beta[m] ~ dnorm(0, 1E-4)
    }

    for ( t in 1:6 ) {
    tau1[t] ~ dgamma(1E-3, 1E-3)
    }
    
    sigma1 <- 1 / sqrt(tau1)
    sigma2 <- 1 / sqrt(tau2)
    sigma3 <- 1 / sqrt(tau3)
    sigma4 <- 1 / sqrt(tau4)
    sigma5 <- 1 / sqrt(tau5)
    sigma6 <- 1 / sqrt(tau6)
    sigma7 <- 1 / sqrt(tau7)
    sigma8 <- 1 / sqrt(tau8)
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
    tau1 = rgamma(6, 1E3, 1E3),
    tau2 = rgamma(3, 1E3, 1E3),
    tau3 = rgamma(3, 1E3, 1E3),
    tau4 = rgamma(3, 1E3, 1E3),
    tau5 = rgamma(3, 1E3, 1E3),
    tau6 = rgamma(3, 1E3, 1E3),
    tau7 = rgamma(3, 1E3, 1E3),
    tau8 = rgamma(3, 1E3, 1E3))
}

# Set monitors
params <- c(paste('beta.t',1:3,sep='.'), paste('mu.beta',1:4,sep='.'),"beta.t","mu.beta", paste('sigma',1:8,sep=''))

# Run model
adapt <- 2500
iter <- 10000
burnin <- 5000
thin <- 5
chains <- 3

mod <- jagsUI::jags(data, inits, params,
                    "growth_3level_Ptaus_NCI_NCIXDBH.bug", 
                    n.chains=chains, n.adapt=adapt, n.iter=iter, 
                    n.burnin=burnin, n.thin=thin, parallel=F)

mod
plot(mod)



x <- lapply(mod$samples, function(x)as.vector(x[,'beta.t[1]']))
plot(x[[1]], type='l')
points(x[[2]], type='l', col=2)
points(x[[3]], type='l', col=3)


plot(rep(1:length(x[[1]]), times=3), c(x[[1]],x[[2]],x[[3]]), col=1, type='l')

plot(mod$samples[,'beta.t[1]'], ylim=c(-1,1))

mod

# Update model
paste("Start at:", Sys.time())
mod <- update(mod, n.iter=1000)
paste("Finish at:", Sys.time())





pdf(file='growth_WSG_3L_Ptaus_NCI_NCIXDBH.bug')
x <- cbind(unlist(mod$q50),unlist(mod$q2.5),unlist(mod$q97.5))
x <- x[-nrow(x),]
x <- x[-grep('sigma', rownames(x)),]
col <- c(rep(2:4, times=7), rep(1, 100))
pch <- ifelse(sign(x[,2])==sign(x[,3]), 16, 1)
plot(x[,1],ylim=c(min(x), max(x)), pch=pch, axes=F, ylab='Std. Effect', xlab='', col=col, cex=2)
segments(1:nrow(x), x[,2], 1:nrow(x), x[,3], col=col)
abline(h=0,lty=2)
axis(1, labels=rownames(x), at=1:nrow(x), las=2)
axis(2); box()
legend('bottomleft', pch=16, col=c(2:4,1), legend=c('bci','cocoli','sherman','pooled plots'), lty=1)
mtext("Model: growth_WSG_3L_Ptaus_NCI_NCIXDBH",3)
dev.off()





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


labs <- names(unlist(fullmod$q50))[-21]
pch <- ifelse(unlist(fullmod$overlap0)==T, 1, 16)
plot(unlist(fullmod$q50)[-21], pch=pch, ylim=c(-3,3), axes=F, xlab='', ylab='effect size')
axis(2)
axis(1, seq(1,20,1), labels=labs, las=2)
segments(seq(1,20,1), unlist(fullmod$q2.5)[-21], seq(1,20,1), unlist(fullmod$q97.5)[-21])
abline(h=0, lty=2)









