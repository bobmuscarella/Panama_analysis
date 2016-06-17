#######################################
###  GROWTH ANALYSIS
###  SINGLE-PLOT
###  ONLY NCI, DBH BY SIZE CLASS
###  Results in "K:/Bob/Panama/RESULTS/..."
#######################################
library(jagsUI)
library(rjags)
library(sp)

### Running on PC???
pc <- T

#######################################
###  START HERE WITH PROCESSED DATA ###
#######################################
if(pc==T){
  setwd("K:/Bob/Panama/DATA") 
  load("Panama_AnalysisData_6.14.16.RDA") # tdata
  load("panama_ITV_traits_6.7.16.RDA") # traits
}


###########################
# PREP FOR TESTING:
if(pc==F){
  setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")
  load("Panama_AnalysisData_6.7.16.RDA") # tdata
  load("panama_ITV_traits_6.7.16.RDA") # traits
}

tdata <- tdata[!is.na(tdata$growth) & !is.na(tdata$All.NCI) & tdata$Not.Edge==1,]
tdata <- droplevels(tdata)

tdata$log.all.nci <- log(tdata$All.NCI + 1)
tdata$log.dbh <- log(tdata$dbh)
traits$log.LMA.mean <- log(traits$LMA.mean)
###########################


###########################
#### Prepare data for input  ####
###########################
# If you want to drop species with NA for both traits
#d <- tdata[tdata$spplot %in% traits$sp[!is.na(traits$WD.mean)] | tdata$spplot %in% traits$sp[!is.na(traits$log.LMA.mean)],]

# Or not
d <- tdata

d <- droplevels(d)


##################################################
#### Start loop to model each plot separately ####
##################################################
#p <- 1
for(p in c(1,3)){

dp <- d[d$plot==p,]

######################################################################################
#### Standardize and Center coefficients (within plots, within size classes 10cm) ####
######################################################################################
###### Generic size class: Center / scale DBH within species, within size class #####
cutoff <- (100)  # Setting cutoff to neg value will put everything in 1 size class.
dp$size.class <- ifelse(dp$dbh <= cutoff, 1, 2)

#size <- 1
    for(size in 1:2) {
dps <- dp[dp$size.class %in% size,]

dps$sd5.growth  <- (abs(dps$growth) < (sd(dps$growth)*5))

#dps <- dps[dps$sd5.growth,]
dps <- dps[dps$sd5.growth & dps$Growth.Include.3,]

### Center / scale other variables within size class
dps$log.dbh.z <- as.vector(scale(dps$log.dbh))
dps$log.all.nci.z <- as.vector(scale(dps$log.all.nci))
dps$growth.z <- as.vector(scale(dps$growth, center=F))

dps <- dps[,c('spplot','plot','census','growth.z','growth','dbh',
              'log.dbh.z','id','log.all.nci.z')]

# Create an individual ID
dps <- droplevels(dps)
dps$indiv <- as.numeric(as.factor(dps$id))

# Order for correct indexing
dps <- dps[order(dps$spplot, dps$indiv, dps$census),]

### GET AND SCALE SPECIES MEAN TRAIT VALUES AND AND INTRASPECIFIC SD'S...
wd.mean <- as.vector(traits$WD.mean[match(unique(dps$spplot), traits$sp)])
wd.sd <- as.vector(traits$WD.sd[match(unique(dps$spplot), traits$sp)])
wd.mean.z <- as.vector(scale(wd.mean))
wd.sp.sd <- sd(wd.mean, na.rm=T)
wd.sd.z <- wd.sd / wd.sp.sd
wd.sd.z[is.na(wd.sd.z)] <- mean(wd.sd.z, na.rm=T)
wd.tau.z <- 1/(wd.sd.z^2)
wd.mean.z[is.na(wd.mean.z)] <- rep(0, length(wd.mean.z[is.na(wd.mean.z)]))

lma.mean <- as.vector(traits$log.LMA.mean[match(unique(dps$spplot), traits$sp)])
lma.sd <- as.vector(traits$log.LMA.sd[match(unique(dps$spplot), traits$sp)])
lma.mean.z <- as.vector(scale(lma.mean))
lma.sp.sd <- sd(lma.mean, na.rm=T)
lma.sd.z <- lma.sd / lma.sp.sd
lma.sd.z[is.na(lma.sd.z)] <- mean(lma.sd.z, na.rm=T)
lma.tau.z <- 1/(lma.sd.z^2)
lma.mean.z[is.na(lma.mean.z)] <- rep(0, length(lma.mean.z[is.na(lma.mean.z)]))

wd.sp.tau <- 1/(wd.sp.sd^2)
lma.sp.tau <- 1/(lma.sp.sd^2)
tcor <- cor(wd.mean.z, lma.mean.z)
omega <- matrix(nrow=2, ncol=2, data=c(wd.sp.tau, tcor, tcor, lma.sp.tau))
tmeans.z <- cbind(wd.mean.z, lma.mean.z)
colnames(tmeans.z) <- NULL      
ttaus.z <- cbind(wd.tau.z, lma.tau.z)
colnames(ttaus.z) <- NULL

#################################
#### Organize the input data ####
#################################
data = list (
  ntree = nrow(dps),
  nindiv = length(unique(dps$id)),
  indiv = dps$indiv,
  nspecies = length(unique(dps$spplot)),
  species = as.numeric(as.factor(dps$spplot)),
  growth = as.numeric(dps$growth.z),
  allnci = as.numeric(dps[,'log.all.nci.z']),
  dbh = as.numeric(dps$log.dbh.z),
  tmeans.z = tmeans.z,
  ttaus.z = ttaus.z,
  omega = omega
)

### Add an indicator to set individual effect of non-rep indiv to zero
if(p!=2){
  repindiv <- names(table(data$indiv))[table(data$indiv) > 1]
  data$indicator <- as.numeric(data$indiv %in% repindiv)
}



##############################
#### Write the model file? ####
##############################
write.mods <- F
if(write.mods==T){
setwd("K:/Bob/Panama/GIT/Panama_Analysis/MODELS") 

sink("growth_6.2.16_bci.bug")
cat(" model {
    
    for( i in 1:ntree ) {
    
    growth[i] ~ dnorm(mu[i], tau[1])
    
    mu[i] <- beta.1[species[i]]
    + beta.2[species[i]] * allnci[i]
    + beta.3[species[i]] * dbh[i]
    }
    
    for( j in 1:nspecies ) {
    ### MULTIVARIATE TRAITS ###
    beta.1[j] ~ dnorm(mu.beta[1] + (beta.wd[1] * t.pred[j,1]) + (beta.lma[1] * t.pred[j,2]), tau[2])
    beta.2[j] ~ dnorm(mu.beta[2] + (beta.wd[2] * t.pred[j,1]) + (beta.lma[2] * t.pred[j,2]), tau[3])
    beta.3[j] ~ dnorm(mu.beta[3], tau[4])
    
    ### TURN ON ITV (v3) ###
    t.pred[j,1:2] ~ dmnorm(pred.tmeans[j,], omega[,])
          for (N in 1:2){
          pred.tmeans[j,N] ~ dnorm(tmeans.z[j,N], pred.tau[j,N])
          pred.tau[j,N] ~ dgamma(sh[j,N], ra[j,N])
          ra[j,N] <- (ttaus.z[j,N] + sqrt(ttaus.z[j,N]^2 + 4*sd[N]^2))/(2*sd[N]^2)
          sh[j,N] <- (1 + ttaus.z[j,N]*ra[j,N])
          }
    }

    ### PRIORS ####
    for( m in 1:3 ) {
    mu.beta[m] ~ dnorm(0, 1E-3)
    }
    
    for( b in 1:2 ) {
    beta.wd[b] ~ dnorm(0, 1E-3)
    beta.lma[b] ~ dnorm(0, 1E-3)
    sd[b] ~ dunif(0, 100)
    }

    for( t in 1:5 ) {
    tau[t] ~ dgamma(1E-3, 1E-3)
    }
    
    sigma <- 1 / sqrt(tau)
    pred.sigma <- 1 / sqrt(pred.tau)
    
    }"
    , fill=TRUE)
sink()



sink("growth_6.2.16_coc.she.bug")
cat(" model {
    
    for( i in 1:ntree ) {
    
    growth[i] ~ dnorm(mu[i], tau[1])
    
    mu[i] <- beta.1[species[i]]
    + beta.2[species[i]] * allnci[i]
    + beta.3[species[i]] * dbh[i]
    + indiv.effect[indiv[i]] * indicator[i]
    }
    
    for( j in 1:nspecies ) {
    ### MULTIVARIATE TRAITS ###
    beta.1[j] ~ dnorm(mu.beta[1] + (beta.wd[1] * t.pred[j,1]) + (beta.lma[1] * t.pred[j,2]), tau[2])
    beta.2[j] ~ dnorm(mu.beta[2] + (beta.wd[2] * t.pred[j,1]) + (beta.lma[2] * t.pred[j,2]), tau[3])
    beta.3[j] ~ dnorm(mu.beta[3], tau[4])
    
    ### TURN OFF ITV
    #           wd.pred[j] <- wd.mean.z[j]
    #           lma.pred[j] <- lma.mean.z[j]
    
    ### TURN ON ITV (v1) ###
    #              for (N in 1:2){ t.pred[j,N] ~ dnorm(tmeans.z[j,N], ttaus.z[j,N]) }
    
    ### TURN ON ITV (v2) ###
    #             t.pred[j,1:2] ~ dmnorm(tmeans.z[j,], omega[,])
    
    ### TURN ON ITV (v3) ###
    t.pred[j,1:2] ~ dmnorm(pred.tmeans[j,], omega[,])
    for (N in 1:2){
    pred.tmeans[j,N] ~ dnorm(tmeans.z[j,N], pred.tau[j,N])
    pred.tau[j,N] ~ dgamma(sh[j,N], ra[j,N])
    ra[j,N] <- (ttaus.z[j,N] + sqrt(ttaus.z[j,N]^2 + 4*sd[N]^2))/(2*sd[N]^2)
    sh[j,N] <- (1 + ttaus.z[j,N]*ra[j,N])
    }
    }
    
    ### PRIORS ####
    for( m in 1:3 ) {
    mu.beta[m] ~ dnorm(0, 1E-3)
    }
    
    for( b in 1:2 ) {
    beta.wd[b] ~ dnorm(0, 1E-3)
    beta.lma[b] ~ dnorm(0, 1E-3)
    sd[b] ~ dgamma(1E-3, 1E-3) #dunif(0, 100)
    }
    
    for( i.a in 1:nindiv ) {
    indiv.effect[i.a] ~ dnorm(0, tau[5])
    }
    
    for( t in 1:5 ) {
    tau[t] ~ dgamma(1E-3, 1E-3)
    }
    
    sigma <- 1 / sqrt(tau)
    pred.sigma <- 1 / sqrt(pred.tau)
    
    }"
    , fill=TRUE)
sink()
}

################################################
### Set initial values, monitors, iterations and run model ###
################################################

if(p==2){  inits <- function (){
    list(
      beta.wd = rnorm(2),
      beta.lma = rnorm(2),
      mu.beta = rnorm(3),
      sd =  rgamma(2, 1E3, 1E3),#runif(2, 0, 100),
      tau = rgamma(4, 1E3, 1E3))
  }}
if(p!=2){  inits <- function (){
  list(
    beta.wd = rnorm(2),
    beta.lma = rnorm(2),
    mu.beta = rnorm(3),
    sd = rgamma(2, 1E3, 1E3),#runif(2, 0, 100),
    tau = rgamma(5, 1E3, 1E3))
}}

# Set monitors & run model
params <- c('beta.wd','beta.lma','mu.beta','sigma')

adapt <- 5000
iter <- 250000
burn <- 200000
thin <- 50
chains <- 3

setwd("K:/Bob/Panama/GIT/Panama_Analysis/MODELS") 
modfile <- ifelse(p!=2, "growth_6.2.16_coc.she.bug", "growth_6.2.16_bci.bug")
warning(paste("Now working on:", paste(ifelse(p==1,'Cocoli',ifelse(p==2,'BCI','Sherman')), ifelse(size==1,'< 10cm','> 10cm'),sep=" ")))

mod <- jagsUI::jags(data, inits, params, modfile, 
                    n.chains=chains, n.adapt=adapt, 
                    n.iter=iter, n.burnin=burn, n.thin=thin, 
                    parallel=T, store.data=T)

# mod
# plot(mod)
# plot.params.2(mod)
# plot(unlist(mod$samples[,'beta.wd[1]']), unlist(mod$samples[,'beta.wd[2]']))
# plot(unlist(mod$samples[,'beta.lma[1]']), unlist(mod$samples[,'beta.lma[2]']))
# plot(unlist(mod$samples[,'beta.wd[1]']), unlist(mod$samples[,'beta.lma[1]']))
# plot(unlist(mod$samples[,'beta.wd[2]']), unlist(mod$samples[,'beta.lma[2]']))
# plot(unlist(mod$samples[,'mu.beta[1]']), unlist(mod$samples[,'mu.beta[2]']))

for(reps in 1:10){
  if(max(unlist(mod$Rhat)) > 1.1){
    print(paste('Doing update #', reps))
    mod <- update(mod, n.iter=10000)
  }
}

setwd("K:/Bob/Panama/RESULTS/_6.14.16/growth") 

file <- paste(ifelse(p==1,'coc',ifelse(p==2,'bci','she')), ifelse(size==1,'sm','lg'), 'Rdata',sep=".")
saveRDS(mod, file=file)

}
}





plot(unlist(mod$samples[,'beta.wd[1]']), unlist(mod$samples[,'beta.wd[2]']))
plot(unlist(mod$samples[,'beta.lma[1]']), unlist(mod$samples[,'beta.lma[2]']))
plot(unlist(mod$samples[,'mu.beta[1]']), unlist(mod$samples[,'mu.beta[2]']))

plot(mod$q50['t.pred'][[1]][,1], mod$q50['t.pred'][[1]][,2])




plot.params <- function(mod){
  x <- cbind(unlist(mod$q2.5), unlist(mod$q50), unlist(mod$q97.5))
  x <- x[-grep('deviance', rownames(x)),]
  x <- x[-grep('sigma', rownames(x)),]
  bg <- ifelse(sign(x[,1]) == sign(x[,3]), 1, 'white')
  plot(x[,2], axes=F, ylim=range(x), xlab='', pch=21, bg=bg, cex=2)
  arrows(1:nrow(x), x[,1], 1:nrow(x), x[,3], len=0.1, code=3, angle=90)
  abline(h=0, lty=2)
  points(x[,2], pch=21, bg=bg, cex=2)
  axis(1, labels=rownames(x), at=1:nrow(x), las=2)
  axis(2); box()
}


plot.params.2 <- function(mod){
  x <- cbind(unlist(mod$q2.5), unlist(mod$q50), unlist(mod$q97.5))
  x <- x[-grep('deviance', rownames(x)),]
  x <- x[-grep('sigma', rownames(x)),]
  # x <- x[-grep('pred', rownames(x)),]
  #   x <- x[-grep('wd.mu', rownames(x)),]
  #   x <- x[-grep('lma.mu', rownames(x)),]
  bg <- ifelse(sign(x[,1]) == sign(x[,3]), 1, 'white')
  plot(x[,2], axes=F, ylim=range(x), xlab='', pch=21, bg=bg, cex=2)
  arrows(1:nrow(x), x[,1], 1:nrow(x), x[,3], len=0.1, code=3, angle=90)
  abline(h=0, lty=2)
  points(x[,2], pch=21, bg=bg, cex=2)
  axis(1, labels=rownames(x), at=1:nrow(x), las=2)
  axis(2); box()
}

