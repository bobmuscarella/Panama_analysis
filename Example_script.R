##### Load the packages #####
library(rjags)
load.module('mix')
set.factory("mix::TemperedMix", 'sampler', FALSE)
load.module('glm')
library(coda)
library(jagsUI)
library(runjags)


##### Load the data #####
# data <- readRDS('Example_data.RDS')

# data$days <- data$days/365

data <- list (
  N = nrow(dps),
  tree = dps$indiv,
  n.tree = length(unique(dps$indiv)),
  n.sp = length(unique(dps$spplot)),
  sp = as.numeric(as.factor(dps$spplot)),
  obs.growth = as.numeric(dps$growth),#.z),
  days = dps$days/365,
  log.nci = as.numeric(dps[,'log.all.nci.z']),
  dbh = as.vector(dps$dbh),
  log.dbh = dps$log.dbh.z,
  tmeans.z = tmeans.z,
  ttaus.z = ttaus.z,
  # trait.corr = trait.corr
  # wd.z = tmeans.z[,1],
  # lma.z = tmeans.z[,2],
  # wd.tau.z = ttaus.z[,1],
  # lma.tau.z = ttaus.z[,2],
  omegas = omegas#,
  # omega = omega
)



##### Look at the data #####
names(data)
# N = the total number of growth observations
# tree = the individual stem ID
# n.tree = the number of individual stems
# n.sp = the total number of species
# sp = the species ID for each stem
# obs.growth = the observed absolute growth rate (mm/yr)
# days = the number of days between measurements for each observation
# log.nci = the log neighborhood crowding index
# dbh = stem diameter (mm)
# log.dbh = log stem diameter

lapply(data, head)


##### Write the model with intraspecific variation #####
sink("Test_Growth_Model.bug")
cat(" model {
    
    for (i in 1:N){
    predict.log.growth[i] <- b0[tree[i]] 
                              + b1[tree[i]] * log.dbh[i] 
                              + b2[tree[i]] * log.nci[i] 
                              + i.tree[tree[i]]
    
    true.log.growth[i] ~ dnorm(predict.log.growth[i], process.tau)
    true.growth[i] <- exp(true.log.growth[i])
    
    obs.growth[i] ~ dnormmix(mu[1:2, i], tau[1:2, i], f)
    mu[1,i] <- true.growth[i]
    mu[2,i] <- true.growth[i]
    
    sd1[i] <- (0.927 + 0.0038 * (dbh[i] - 45)) * 1.414 / days[i]
    tau[1,i] <- pow(sd1[i], -2)
    sd2[i] <- 25.6 * 1.414 / days[i]
    tau[2,i] <- pow(sd2[i], -2)
    }
    
    for (j in 1:n.tree){
    
    b0[j] <- b0.overall + (b0.wd * t.pred[j,1]) + (b0.lma * t.pred[j,2]) + b0.sp[sp[j]]
    b1[j] <- b1.overall + (b1.wd * t.pred[j,1]) + (b1.lma * t.pred[j,2]) + b1.sp[sp[j]]
    b2[j] <- b2.overall + b2.sp[sp[j]]
    
    t.pred[j,1:2] ~ dmnorm(tmeans.z[sp[j],], omegas[1:2, 1:2, sp[j]])
    }

    ### prior and random effect ##########
    f[1] <- 0.9724
    f[2] <- 0.0276

    for(sp in 1:n.sp){
      b0.sp[sp] ~ dnorm(0, 10000)
      b1.sp[sp] ~ dnorm(0, 10000)
      b2.sp[sp] ~ dnorm(0, 10000)
    }

    for (i in 1:n.tree){
      i.tree[i] ~ dnorm(0, i.tree.prec)
    }
    
    process.tau ~ dunif(0.0001, 10000)
    i.tree.prec ~ dunif(0.0001, 10000)
    
    b0.overall ~ dnorm(0, 1.0E-6)
    b1.overall ~ dnorm(0, 1.0E-6)
    b2.overall ~ dnorm(0, 1.0E-6)

    b0.wd ~ dnorm(0, b0.wd.prec)
    b1.wd ~ dnorm(0, b1.wd.prec)
    b0.lma ~ dnorm(0, b0.lma.prec)
    b1.lma ~ dnorm(0, b1.lma.prec)

    b0.wd.prec ~ dunif(0.0001, 10000)
    b0.lma.prec ~ dunif(0.0001, 10000)
    b1.wd.prec ~ dunif(0.0001, 10000)
    b1.lma.prec ~ dunif(0.0001, 10000)

    }"
    , fill=TRUE)
sink()


##### Write the model with intraspecific variation #####
sink("Test_Growth_Model_noITV.bug")
cat(" model {
    
    for (i in 1:N){
    predict.log.growth[i] <- b0[sp[i]] 
    + b1[sp[i]] * log.nci[i] 
    + b2[sp[i]] * log.dbh[i] 
    + i.tree[tree[i]]
    
    true.log.growth[i] ~ dnorm(predict.log.growth[i], process.tau)
    true.growth[i] <- exp(true.log.growth[i])
    
    obs.growth[i] ~ dnormmix(mu[1:2, i], tau[1:2, i], f)
    mu[1,i] <- true.growth[i]
    mu[2,i] <- true.growth[i]
    
    sd1[i] <- (0.927 + 0.0038 * (dbh[i] + 15)) * 1.414 / days[i]
#    sd1[i] <- (0.927 + 0.0038 * (dbh[i] - 45)) * 1.414 / days[i]
    tau[1,i] <- pow(sd1[i], -2)
    sd2[i] <- 25.6 * 1.414 / days[i]
    tau[2,i] <- pow(sd2[i], -2)
    }
    
    for (j in 1:n.sp){
    b0[j] ~ dnorm(b0.overall + (b0.wd * tmeans.z[j,1]) + (b0.lma * tmeans.z[j,2]), b0.sp.prec)
    b1[j] ~ dnorm(b1.overall + (b1.wd * tmeans.z[j,1]) + (b1.lma * tmeans.z[j,2]), b1.sp.prec)
    b2[j] ~ dnorm(b2.overall, b2.sp.prec)

    # b0[j] <- b0.overall + (b0.wd * tmeans.z[j,1]) + (b0.lma * tmeans.z[j,2]) + b0.sp[j]
    # b1[j] <- b1.overall + (b1.wd * tmeans.z[j,1]) + (b1.lma * tmeans.z[j,2]) + b1.sp[j]
    # b2[j] <- b2.overall + b2.sp[j]
    # 
    # b0.sp[j] ~ dnorm(0, b0.sp.prec)
    # b1.sp[j] ~ dnorm(0, b1.sp.prec)
    # b2.sp[j] ~ dnorm(0, b2.sp.prec)
    }
    
    ### prior and random effect ##########
    f[1] <- 0.9724
    f[2] <- 0.0276
    
    for (i in 1:n.tree){
    i.tree[i] ~ dnorm(0, i.tree.prec)
    }
    
    process.tau ~ dunif(0.0001, 10000)
    i.tree.prec ~ dunif(0.0001, 10000)

    b0.sp.prec ~ dunif(0.0001, 10000)
    b1.sp.prec ~ dunif(0.0001, 10000)
    b2.sp.prec ~ dunif(0.0001, 10000)
    
    b0.overall ~ dnorm(0, 1.0E-6)
    b1.overall ~ dnorm(0, 1.0E-6)
    b2.overall ~ dnorm(0, 1.0E-6)

    b0.wd ~ dnorm(0, 1.0E-6)
    b1.wd ~ dnorm(0, 1.0E-6)
    b0.lma ~ dnorm(0, 1.0E-6)
    b1.lma ~ dnorm(0, 1.0E-6)

    }"
    , fill=TRUE)
sink()


##### Run the model #####

params <- c('b0.overall','b1.overall','b2.overall','b0.wd','b0.lma','b1.wd','b1.lma')

j.mod <- jags.model(file="Test_Growth_Model.bug", data=data, n.chains=3, n.adapt=1000)
update(j.mod, n.iter=1000)
samp <- jags.samples(j.mod, params, n.iter=900, n.thin=3)


summary(samp$b0.wd, quantile, c(.025,0.5,.975))$stat





# samp

mod <- jagsUI::jags(data, inits=NULL, params, 'Test_Growth_Model_noITV.bug',
                    n.chains=3, n.adapt=1000,
                    n.iter=5000, n.burnin=2500, n.thin=1,
                    parallel=F, store.data=F, modules=c('mix','glm'))

mod <- update(mod, n.iter=5000, n.burnin=2500, n.thin=2, modules=c('mix','glm'))


par(mfrow=c(3,3))
for(i in 1:8){
plot(as.vector(mod$samples[[1]][,8]), as.vector(mod$samples[[1]][,i]))
}


mod <- run.jags(model='Test_Growth_Model2.bug', monitor=params, data=data, 
                n.chains=3, burnin=500, sample=250, adapt=100, modules=c('mix','glm'), 
                factories='mix::TemperedMix sampler off', thin=3)#, method='parallel')

plot(mod, c('trace'), vars=c(params), layout=c(3,3))

mod.ext <- extend.jags(mod, sample=600, thin=3, burnin=1000)
mod.ext

mod.ext <- extend.jags(mod.ext, sample=600, thin=3, burnin=1000)
mod.ext

plot(mod.ext, c('trace'), vars=c(params), layout=c(3,3))



