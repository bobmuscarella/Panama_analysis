#######################################
###  GROWTH ANALYSIS
###  SINGLE-PLOT
###  ONLY NCI, DBH BY SIZE CLASS
###  Results in "K:/Bob/Panama/RESULTS/..."
#######################################
library(rjags)
load.module('mix')
set.factory("mix::TemperedMix", 'sampler', FALSE)
load.module('glm')
library(coda)
library(runjags)
# library(jagsUI)

### Running on PC???
pc <- T

#######################################
###  START HERE WITH PROCESSED DATA ###
#######################################
if(pc==T){
  setwd("K:/Bob/Panama/DATA") 
  load("Panama_AnalysisData_12.17.16.RDA") # tdata
  load("panama_ITV_traits_6.7.16.RDA") # traits
}


###########################
# PREP FOR TESTING:
if(pc==F){
  setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")
  load("Panama_AnalysisData_6.14.16.RDA") # tdata
  load("Panama_AnalysisData_12.17.16.RDA") # tdata
  load("panama_ITV_traits_6.7.16.RDA") # traits
}

tdata <- tdata[!is.na(tdata$growth) & !is.na(tdata$All.NCI) & tdata$Not.Edge==1,]
tdata <- droplevels(tdata)

tdata$log.all.nci <- log(tdata$All.NCI)
tdata$log.dbh <- log(tdata$dbh)
traits$log.LMA.mean <- log(traits$LMA.mean)
###########################


###########################
#### Prepare data for input  ####
###########################
# Get species with WD data from Panama
wdsp <- names(rowSums(table(tdata$spplot, tdata$wd.source)[,1:2]))[rowSums(table(tdata$spplot, tdata$wd.source)[,1:2])>0]
# Get species with LMA data from Panama
lmasp <- names(rowSums(table(tdata$spplot, tdata$lma.source)[,1:3]))[rowSums(table(tdata$spplot, tdata$lma.source)[,1:3])>0]

foctraitsp <- unique(wdsp, lmasp)
tdata <- tdata[tdata$spplot %in% foctraitsp,]

# If you want to drop species with NA for both traits
#d <- tdata[tdata$spplot %in% traits$sp[!is.na(traits$WD.mean)] | tdata$spplot %in% traits$sp[!is.na(traits$log.LMA.mean)],]

# If you want to drop species with NA for either trait
d <- tdata[tdata$spplot %in% traits$sp[!is.na(traits$WD.mean)] & tdata$spplot %in% traits$sp[!is.na(traits$log.LMA.mean)],]

# Or not
#d <- tdata

d <- droplevels(d)



### Center / scale variables across full dataset
d$log.dbh.zall <- as.vector(scale(d$log.dbh))
d$log.all.nci.zall <- as.vector(scale(d$log.all.nci))

d$log.all.size.nci.zall[d$dbh<100] <- as.vector(scale(d$log.all.nci[d$dbh<100]))
d$log.all.size.nci.zall[d$dbh>=100] <- as.vector(scale(d$log.all.nci[d$dbh>=100]))


##################################################
#### Start loop to model each plot separately ####
##################################################
# p <- 1
for(p in c(1)) {

dp <- d[d$plot==p,]

######################################################################################
#### Standardize and Center coefficients (within plots, within size classes 10cm) ####
######################################################################################
###### Generic size class: Center / scale DBH within species, within size class #####
cutoff <- (100)  # Setting cutoff to neg value will put everything in 1 size class.
dp$size.class <- ifelse(dp$dbh <= cutoff, 1, 2)


# size <- 2
for(size in 1:2) {
dps <- dp[dp$size.class %in% size,]


dps$sd5.growth  <- (abs(dps$growth) < (sd(dps$growth)*5))
#dps <- dps[dps$Growth.Include.3,]
dps <- dps[dps$sd5.growth & dps$Growth.Include.3,]

dps$sd5.RGR  <- (abs(dps$RGR) < (sd(dps$RGR)*5))
dps <- dps[dps$sd5.RGR & dps$Growth.Include.3,]

### Center / scale other variables within size class
dps$log.dbh.z <- as.vector(scale(dps$log.dbh))
dps$log.all.nci.z <- as.vector(scale(dps$log.all.nci))
dps$growth.z <- as.vector(scale(dps$growth, center=F))
dps$RGR.z <- as.vector(scale(dps$growth, center=F))

dps <- dps[,c('spplot','plot','census','growth.z','RGR.z','growth','dbh',
              'log.dbh.zall','log.all.nci.zall', 'log.dbh.z',
              'id','log.all.nci.z','days','log.all.size.nci.zall')]

### Create an individual ID
dps <- droplevels(dps)
dps$indiv <- as.numeric(as.factor(dps$id))

### Order for correct indexing
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

wd.sd[is.na(wd.sd)] <- mean(wd.sd[!is.na(wd.sd)])
lma.sd[is.na(lma.sd)] <- mean(lma.sd[!is.na(lma.sd)])
wd.tau <- 1/(wd.sd^2)
lma.tau <- 1/(lma.sd^2)


tmeans <- cbind(wd.mean, lma.mean)
colnames(tmeans) <- NULL      
tmeans.z <- cbind(wd.mean.z, lma.mean.z)
colnames(tmeans.z) <- NULL      
ttaus.z <- cbind(wd.tau.z, lma.tau.z)
colnames(ttaus.z) <- NULL

omega <- matrix(nrow=2, ncol=2, data=c(wd.sp.tau, tcor, tcor, lma.sp.tau))

omegas <- list()
for(i in 1:length(wd.mean.z)){
  omegas[[i]] <- matrix(ncol=2,nrow=2,data=c(wd.tau[i],#wd.tau.z[i],
                                             tcor, tcor,
                                             lma.tau[i]))#lma.tau.z[i]))
}
omegas <- array(unlist(omegas), dim=c(2,2,length(omegas)))



#################################
#### Organize the input data ####
#################################
data = list (
  N = nrow(dps),
  tree = dps$indiv,
  n.tree = length(unique(dps$indiv)),
  n.sp = length(unique(dps$spplot)),
  sp = as.numeric(as.factor(dps$spplot)),
<<<<<<< HEAD
#  obs.growth = as.numeric(dps$growth),
  obs.growth = as.numeric(dps$RGR.z),
=======
  obs.growth = as.numeric(dps$growth),#.z),
  #  obs.growth = dps$RGR,#as.numeric(dps$growth),#.z),
>>>>>>> origin/master
  days = dps$days/365,
#  log.nci = as.numeric(dps[,'log.all.nci.zall']),
#  log.nci = as.numeric(dps[,'log.all.nci.z']),
  log.nci = as.numeric(dps[,'log.all.size.nci.zall']),
  dbh = as.vector(dps$dbh),
#  log.dbh = dps$log.dbh.zall,
  log.dbh = dps$log.dbh.z,
  tmeans.z = tmeans.z,
  omegas = omegas
)

### Add an indicator to set individual effect of non-rep indiv to zero
if(p!=2){
  repindiv <- names(table(data$tree))[table(data$tree) > 1]
  data$indicator <- as.numeric(data$tree %in% repindiv)
}

##############################
#### Write the model file? ####
##############################
setwd("K:/Bob/Panama/GIT/Panama_Analysis/MODELS") 

sink("Growth_Model_wITV_1.31.17_error.bug")
  cat(" model {
      
      for (i in 1:N){

      obs.growth[i] ~ dnormmix(mu[1:2, i], m.tau[1:2, i], f)
      mu[1,i] <- true.growth[i]
      mu[2,i] <- true.growth[i]
      sd1[i] <- (0.927 + 0.0038 * (dbh[i] + 15)) * 1.414 / days[i]
      sd2[i] <- 25.6 * 1.414 / days[i]
      m.tau[1,i] <- pow(sd1[i], -2)
      m.tau[2,i] <- pow(sd2[i], -2)
      
      true.growth[i] <- exp(true.log.growth[i])

      true.log.growth[i] ~ dnorm(predict.log.growth[i], tau[1])

      predict.log.growth[i] <- b0[sp[i]] 
              + b1[sp[i]] * log.nci[i] 
              + b2[sp[i]] * log.dbh[i] 
              + indiv.effect[tree[i]] * indicator
      }
      
      for( j in 1:n.sp ) {
        b0[j] ~ dnorm(mu.beta[1] + (beta.wd[1] * t.pred[j,1]) + (beta.lma[1] * t.pred[j,2]), tau[2])
        b1[j] ~ dnorm(mu.beta[2] + (beta.wd[2] * t.pred[j,1]) + (beta.lma[2] * t.pred[j,2]), tau[3])
        b2[j] ~ dnorm(mu.beta[3], tau[4])
        t.pred[j,1:2] ~ dmnorm(tmeans.z[j,], omegas[1:2,1:2,j])
      }

      ### prior and random effect ##########
      f[1] <- 0.9724
      f[2] <- 0.0276

      for( i.a in 1:n.tree ) {
      indiv.effect[i.a] ~ dnorm(0, tau[5])
      }

      for( t in 1:5 ) {
      tau[t] ~ dgamma(1E-3, 1E-3)
      }

      for( m in 1:3 ) {
      mu.beta[m] ~ dnorm(0, 1E-3)
      }

      for( b in 1:2 ) {
        beta.wd[b] ~ dnorm(0, 1E-3)
        beta.lma[b] ~ dnorm(0, 1E-3)
      }

}"
    , fill=TRUE)
  sink()
  
  

sink("Growth_Model_indiv_wITV_1.31.17_error.bug")
  cat(" model {
      
      for (i in 1:N){
      
      obs.growth[i] ~ dnormmix(mu[1:2, i], m.tau[1:2, i], f)
      mu[1,i] <- true.growth[i]
      mu[2,i] <- true.growth[i]
      sd1[i] <- (0.927 + 0.0038 * (dbh[i] + 15)) * 1.414 / days[i]
      sd2[i] <- 25.6 * 1.414 / days[i]
      m.tau[1,i] <- pow(sd1[i], -2)
      m.tau[2,i] <- pow(sd2[i], -2)
      
      true.growth[i] <- exp(true.log.growth[i])
      
      true.log.growth[i] ~ dnorm(predict.log.growth[i], tau[1])
      
      predict.log.growth[i] <- b0[sp[i]] 
      + b1[sp[i]] * log.nci[i] 
      + b2[sp[i]] * log.dbh[i] 
#      + indiv.effect[tree[i]] * indicator
      }
      
      for( j in 1:n.tree ) {
      b0[j] ~ dnorm(mu.beta[1] + (beta.wd[1] * t.pred[j,1]) + (beta.lma[1] * t.pred[j,2]), tau[2])
      b1[j] ~ dnorm(mu.beta[2] + (beta.wd[2] * t.pred[j,1]) + (beta.lma[2] * t.pred[j,2]), tau[3])
      b2[j] ~ dnorm(mu.beta[3], tau[4])
      t.pred[j,1:2] ~ dmnorm(tmeans.z[sp[j],], omegas[1:2,1:2,sp[j]])
      }
      
      ### prior and random effect ##########
      f[1] <- 0.9724
      f[2] <- 0.0276

      # for( i.a in 1:n.tree ) {
      # indiv.effect[i.a] ~ dnorm(0, tau[5])
      # }

      for( t in 1:4 ) {
      # for( t in 1:5 ) {
      tau[t] ~ dgamma(1E-3, 1E-3)
      }
      
      for( m in 1:3 ) {
      mu.beta[m] ~ dnorm(0, 1E-3)
      }
      
      for( b in 1:2 ) {
      beta.wd[b] ~ dnorm(0, 1E-3)
      beta.lma[b] ~ dnorm(0, 1E-3)
      }
      
}"
    , fill=TRUE)
  sink()
  

  
  sink("Growth_Model_noITV_2.9.17_error.bug")
  cat(" model {
      
      for (i in 1:N){
      
      obs.growth[i] ~ dnormmix(mu[1:2, i], m.tau[1:2, i], f)
      mu[1,i] <- true.growth[i]
      mu[2,i] <- true.growth[i]
      sd1[i] <- (0.927 + 0.0038 * (dbh[i] + 15)) * 1.414 / days[i]
      sd2[i] <- 25.6 * 1.414 / days[i]
      m.tau[1,i] <- pow(sd1[i], -2)
      m.tau[2,i] <- pow(sd2[i], -2)
      
      true.growth[i] <- exp(true.log.growth[i])
      
      true.log.growth[i] ~ dnorm(predict.log.growth[i], tau[1])
      
      predict.log.growth[i] <- b0[sp[i]] 
      + b1[sp[i]] * log.nci[i] 
      + b2[sp[i]] * log.dbh[i]
      + indiv.effect[tree[i]] * indicator[i]
      }
      
      for( j in 1:n.sp ) {
      b0[j] ~ dnorm(mu.beta[1] + (beta.wd[1] * tmeans.z[j,1]) + (beta.lma[1] * tmeans.z[j,2]), tau[2])
      b1[j] ~ dnorm(mu.beta[2] + (beta.wd[2] * tmeans.z[j,1]) + (beta.lma[2] * tmeans.z[j,2]), tau[3])
      b2[j] ~ dnorm(mu.beta[3], tau[4])
      }
      
      ### prior and random effect ##########
      f[1] <- 0.9724
      f[2] <- 0.0276
      
      for( i.a in 1:n.tree ) {
      indiv.effect[i.a] ~ dnorm(0, tau[5])
      }
      
      for( t in 1:5 ) {
      tau[t] ~ dgamma(1E-3, 1E-3)
      }
      
      for( m in 1:3 ) {
      mu.beta[m] ~ dnorm(0, 1E-3)
      }
      
      for( b in 1:2 ) {
      beta.wd[b] ~ dnorm(0, 1E-3)
      beta.lma[b] ~ dnorm(0, 1E-3)
      }
      
}"
    , fill=TRUE)
  sink()  

  
sink("Growth_Model_noITV_2.9.17_noerror.bug")
  cat(" model {
      
      for (i in 1:N){

      # obs.growth[i] ~ dnormmix(mu[1:2, i], m.tau[1:2, i], f)
      # mu[1,i] <- true.growth[i]
      # mu[2,i] <- true.growth[i]
      # sd1[i] <- (0.927 + 0.0038 * (dbh[i] + 15)) * 1.414 / days[i]
      # sd2[i] <- 25.6 * 1.414 / days[i]
      # m.tau[1,i] <- pow(sd1[i], -2)
      # m.tau[2,i] <- pow(sd2[i], -2)
      # 
      # true.growth[i] <- exp(true.log.growth[i])
      # 
      # true.log.growth[i] ~ dnorm(predict.log.growth[i], tau[1])
      
      obs.growth[i] ~ dnorm(predict.growth[i], tau[1])

      predict.growth[i] <- exp(predict.log.growth[i])

      predict.log.growth[i] <- b0[sp[i]] 
      + b1[sp[i]] * log.nci[i] 
      + b2[sp[i]] * log.dbh[i]
#      + indiv.effect[tree[i]] * indicator[i]
      }
      
      for( j in 1:n.sp ) {
      b0[j] ~ dnorm(mu.beta[1] + (beta.wd[1] * tmeans.z[j,1]) + (beta.lma[1] * tmeans.z[j,2]), tau[2])
      b1[j] ~ dnorm(mu.beta[2] + (beta.wd[2] * tmeans.z[j,1]) + (beta.lma[2] * tmeans.z[j,2]), tau[3])
      b2[j] ~ dnorm(mu.beta[3] + (beta.wd[2] * tmeans.z[j,1]) + (beta.lma[2] * tmeans.z[j,2]), tau[4])
      }
      
      ### prior and random effect ##########
      # f[1] <- 0.9724
      # f[2] <- 0.0276
      
      # for( i.a in 1:n.tree ) {
      # indiv.effect[i.a] ~ dnorm(0, tau[5])
      # }
      
      for( t in 1:4 ) {
      # for( t in 1:5 ) {
      tau[t] ~ dgamma(1E-3, 1E-3)
      }
      
      for( m in 1:3 ) {
      mu.beta[m] ~ dnorm(0, 1E-3)
      }
      
      for( b in 1:2 ) {
      beta.wd[b] ~ dnorm(0, 1E-3)
      beta.lma[b] ~ dnorm(0, 1E-3)
      }
      
}"
<<<<<<< HEAD
    , fill=TRUE)
  sink()  
  
  
  
  sink("Growth_Model_noITV_noerror_RGR.bug")
  cat(" model {
      
      for (i in 1:N){
      
      obs.growth[i] ~ dnorm(predict.growth[i], tau[1])
      
      predict.growth[i] <- exp(predict.log.growth[i])
      
      predict.log.growth[i] <- b0[sp[i]] 
      + b1[sp[i]] * log.nci[i] 
      + b2[sp[i]] * log.dbh[i]
      + indiv.effect[tree[i]] * indicator[i]
      }
=======
  , fill=TRUE)
sink()





##### Write the model with intraspecific variation #####
sink("Growth_Model_noITV_1.31.17.bug")
cat(" model {
    
    for (i in 1:N){
    
        obs.growth[i] ~ dnorm(mu[i], tau[1])
        mu[i] <- exp(z[i])
        z[i] <- b0[sp[i]]
                + b1[sp[i]] * log.nci[i] 
                + b2[sp[i]] * log.dbh[i] 
                + i.tree[tree[i]]
    }
    
    for (j in 1:n.sp){
    b0[j] ~ dnorm(mu.beta[1] + (beta.wd[1] * tmeans.z[j,1]) + (beta.lma[1] * tmeans.z[j,2]), tau[2])
    b1[j] ~ dnorm(mu.beta[2] + (beta.wd[2] * tmeans.z[j,1]) + (beta.lma[2] * tmeans.z[j,2]), tau[3])
    b2[j] ~ dnorm(mu.beta[3], tau[4])
    }
    
    ### prior and random effect ##########
    for (ind in 1:n.tree){
    i.tree[ind] ~ dnorm(0, tau[5])
    }
    
    for (t in 1:5){
    tau[t] ~ dunif(0.0001, 10000)
    }

    for (m in 1:3){
    mu.beta[m] ~ dnorm(0, 1.0E-6)
    }

    for (b in 1:2){
    beta.wd[b] ~ dnorm(0, 1.0E-6)
    beta.lma[b] ~ dnorm(0, 1.0E-6)
    }

}"
  , fill=TRUE)
sink()








}
>>>>>>> origin/master

      for( j in 1:n.sp ) {
      b0[j] ~ dnorm(mu.beta[1] + (beta.wd[1] * tmeans.z[j,1]) + (beta.lma[1] * tmeans.z[j,2]), tau[2])
      b1[j] ~ dnorm(mu.beta[2] + (beta.wd[2] * tmeans.z[j,1]) + (beta.lma[2] * tmeans.z[j,2]), tau[3])
      b2[j] ~ dnorm(mu.beta[3], tau[4])
      }
      
      ### prior and random effect ##########
      for( i.a in 1:n.tree ) {
      indiv.effect[i.a] ~ dnorm(0, tau[5])
      }
      
      for( t in 1:5 ) {
      tau[t] ~ dgamma(1E-3, 1E-3)
      }
      
      for( m in 1:3 ) {
      mu.beta[m] ~ dnorm(0, 1E-3)
      }
      
      for( b in 1:2 ) {
      beta.wd[b] ~ dnorm(0, 1E-3)
      beta.lma[b] ~ dnorm(0, 1E-3)
      }
      
}"
    , fill=TRUE)
  sink()  
  
################################################
### Set initial values, monitors, iterations and run model ###
################################################
params <- c('mu.beta','beta.wd','beta.lma')

<<<<<<< HEAD
warning(paste("Now working on:", paste(ifelse(p==1,'Cocoli',ifelse(p==2,'BCI','Sherman')), ifelse(size==1,'< 10cm','> 10cm'),sep=" ")), immediate. = T)
  
mod <- run.jags(model='K:/Bob/Panama/GIT/Panama_Analysis/MODELS/Growth_Model_noITV_2.9.17_error.bug', 
                monitor=params, data=data, n.chains=3, 
                burnin=5000, sample=1000, adapt=1000, modules=c('glm'), 
                thin=3, method='parallel')



mod <- run.jags(model='K:/Bob/Panama/GIT/Panama_Analysis/MODELS/Growth_Model_noITV_2.9.17_error.bug',
                monitor=params, data=data, n.chains=3,
                burnin=2000, sample=750, adapt=500, modules=c('glm','mix'),
                factories='mix::TemperedMix sampler off', thin=3, method='parallel')

for(i in 1:10){
  if(any(mod$psrf$psrf[,1] > 1.0999)){
    warning(paste("Now working on update",i), immediate. = T)
    mod <- extend.jags(mod, burnin=5000, sample=1000, adapt=0,
                      thin=3, method='parallel', combine=F)
  }
}

setwd("K:/Bob/Panama/RESULTS/_2.14.17/growth") 
file <- paste(ifelse(p==1,'coc',ifelse(p==2,'bci','she')), ifelse(size==1,'sm','lg'), 'Rdata',sep=".")
saveRDS(mod, file=file)
}
}





=======
mod <- run.jags(model='Growth_Model_noITV_1.31.17.bug', monitor=params, data=data,
                n.chains=3, burnin=2500, sample=250, adapt=1000, modules=c('glm'), thin=3)

#mod <- run.jags(model='Growth_Model_noITV_1.31.17.bug', monitor=params, data=data,
#                n.chains=3, burnin=500, sample=250, adapt=100, modules=c('mix','glm'),
#                factories='mix::TemperedMix sampler off', thin=3)#, method='parallel')

mod.ext <- extend.jags(mod, burnin=1000, sample=500, adapt=1000, thin=3, combine=F)

mod.ext

mod.ext <- extend.jags(mod, burnin=10000, sample=500, adapt=0, 
                       thin=3, combine=F, method='parallel')
>>>>>>> origin/master

mod.ext

plot(mod.ext$mcmc)

plot(mod.ext)


<<<<<<< HEAD
### WORK ON BCI LARGE, RESTART MODEL WITH INITIAL VALUES SET 
setwd("K:/Bob/Panama/RESULTS/_2.14.17/growth") 

mod <- readRDS(list.files()[4])

mod

inits <- function(){
  list(
    mu.beta=c(rnorm(1,4,0.001),rnorm(1,-0.5,0.001),rnorm(1,-0.4,0.001)),
    beta.wd=c(rnorm(1,3.3,0.001),rnorm(1,-0.6,0.001)),
    beta.lma=c(rnorm(1,0.6,0.001),rnorm(1,-2.7,0.001))
  )
}

inits()

setwd("K:/Bob/Panama/GIT/Panama_Analysis/MODELS") 
mod2 <- run.jags(model='survival_2.2.17_bci2.bug', monitor=params, data=mod$data, 
                 n.chains=3, inits=inits, burnin=100, sample=100, adapt=100, 
                 modules=c('glm','mix'), factories='mix::TemperedMix sampler off', 
                 thin=3, method='parallel')

mod2 <- extend.jags(mod2, burnin=500, sample=500, adapt=0,
                   thin=3, method='parallel', combine=F)
=======
>>>>>>> origin/master


mod2 <- extend.jags(mod2, burnin=2000, sample=500, adapt=0,
                    thin=3, method='parallel', combine=F)



mod2$dic



extract.runjags(mod.ext, 'dic')




############################################
############### OLD WAY ####################
############################################

# if(p==2){  inits <- function (){
#     list(
#       beta.wd = rnorm(2),
#       beta.lma = rnorm(2),
#       mu.beta = rnorm(3),
#       sd =  rgamma(2, 1E3, 1E3),#runif(2, 0, 100),
#       tau = rgamma(4, 1E3, 1E3))
#   }}
# if(p!=2){  inits <- function (){
#   list(
#     beta.wd = rnorm(2),
#     beta.lma = rnorm(2),
#     mu.beta = rnorm(3),
#     sd = rgamma(2, 1E3, 1E3),#runif(2, 0, 100),
#     tau = rgamma(5, 1E3, 1E3))
# }}


# Set monitors & run model
params <- c('beta.wd','beta.lma','mu.beta','sigma')#,'R2','y.rep')

adapt <- 1000
iter <- 250000
burn <- 200000
thin <- 50
chains <- 3

setwd("K:/Bob/Panama/GIT/Panama_Analysis/MODELS") 
modfile <- ifelse(p!=2, "growth_12.17.16_coc.she.bug", "growth_12.17.16_bci.bug")
warning(paste("Now working on:", paste(ifelse(p==1,'Cocoli',ifelse(p==2,'BCI','Sherman')), ifelse(size==1,'< 10cm','> 10cm'),sep=" ")))

mod <- jagsUI::jags(data, inits, params, modfile, 
                    n.chains=chains, n.adapt=adapt, 
                    n.iter=iter, n.burnin=burn, n.thin=thin, 
                    parallel=F, store.data=T)

# mod
# plot(mod)
# plot.params.2(mod)
# plot(unlist(mod$samples[,'beta.wd[1]']), unlist(mod$samples[,'beta.wd[2]']))
# plot(unlist(mod$samples[,'beta.lma[1]']), unlist(mod$samples[,'beta.lma[2]']))
# plot(unlist(mod$samples[,'beta.wd[1]']), unlist(mod$samples[,'beta.lma[1]']))
# plot(unlist(mod$samples[,'beta.wd[2]']), unlist(mod$samples[,'beta.lma[2]']))
# plot(unlist(mod$samples[,'mu.beta[1]']), unlist(mod$samples[,'mu.beta[2]']))


x <- unlist(mod$q50)
x <- x[grep('y.rep', names(x))]
plot(x, data$growth)
summary(lm(data$growth ~ x))



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

