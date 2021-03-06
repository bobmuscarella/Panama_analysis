#######################################
###  GROWTH ANALYSIS
###  SINGLE-PLOT
###  ONLY NCI, DBH BY SIZE CLASS
###  Results in "K:/Bob/Panama/RESULTS/_12.17.15/growth/nci_dbh_sizes"
#######################################

library(jagsUI)

### Z-TRANSFORM DATA
z.score <- function (data, center=T, scale=T) {
  xm <- ifelse(center==T, mean (data, na.rm=TRUE), 0)
  xsd <- ifelse(scale==T, sd(data, na.rm=TRUE), 2)
  xtrans <- (data - xm) / (2 * xsd)
  return(xtrans)
}

### Running on PC???
pc <- T

#######################################
###  START HERE WITH PROCESSED DATA ###
#######################################
setwd("K:/Bob/Panama/GIT/Panama_analysis/DATA") 
load("Panama_AnalysisData_12.9.15.RDA")

###########################
#### Prepare data for input  ####
###########################
tdata <- tdata[tdata$Growth.Include.3 == TRUE,]
tdata <- tdata[tdata$growth > (-75),]

names(tdata)[names(tdata)=='WSG'] <- 'wsg'
names(tdata)[names(tdata)=='log.LDMC_AVI'] <- 'log.ldmc'
names(tdata)[names(tdata)=='log.LMALEAF_AVI'] <- 'log.lma'
names(tdata)[names(tdata)=='log.SEED_DRY'] <- 'log.seed'
names(tdata)[names(tdata)=='HEIGHT_AVG'] <- 'hmax'

tdata$plot <- ifelse(tdata$plot=='bci', 2, ifelse(tdata$plot=='cocoli', 1, 3))

#####################
#### Start Loop  ####
#####################

traits <- c('wsg','log.ldmc','log.lma','log.seed','hmax')

for (trt in 1:length(traits)) {
  
  trait <- traits[trt]
  
  # Remove species with NA for trait value
  d <- tdata[!is.na (tdata[,trait]),]  					
  d <- droplevels(d)
  
  ######################################################################################
  #### Standardize and Center coefficients (within plots, within size classes 10cm) ####
  ######################################################################################
  ###### Generic size class: Center / scale DBH within species, within size class #####
  cutoff <- (100)  # Setting cutoff to neg value will put everything in 1 size class.
  d$size.class <- ifelse(d$dbh <= cutoff, 1, 2)
  
  ### Center / scale other variableswithin size class
  for (i in 1:2){
    d$log.dbh.z[d$size.class %in% i] <- z.score(d$log.dbh[d$size.class %in% i])
    d$log.nci.z[d$size.class %in% i] <- z.score(d$log.nci[d$size.class %in% i])
    d$growth.z[d$size.class %in% i] <- z.score(d$growth[d$size.class %in% i], center=F)
  }
  
  for(p in 1:3){
    dp <- d[d$plot==p,]
    
    for(size in 1:2) {
      dps <- dp[dp$size.class %in% size,]
      
      dps <- dps[,c('spcode','plot','size.class','census','growth.z','log.dbh.z','id','log.nci.z', 
                    trait)]
      
      # Drop factors for correct indexing
      dps <- droplevels(dps)
      
      # Create an individual ID
      dps$indiv <- as.numeric(as.factor(dps$id))
      
      # Order for correct indexing
      dps <- dps[order(dps$indiv, dps$census, dps$spcode, dps$plot),]
      
      #################################
      #### Organize the input data ####
      #################################
      data = list (
        ntree = nrow(dps),
        nindiv = length(unique(dps$id)),
        nspecies = length(unique(dps$spcode)),
        growth = as.numeric(dps$growth.z),
        nci = as.numeric(dps[,'log.nci.z']),
        ncensus =  length(unique(paste(dps$census, dps$plot, sep=''))),
        census = as.numeric(as.factor(paste(d$census, dps$plot, sep=''))),  
        dbh = as.numeric(dps$log.dbh.z),
        trait = z.score(tapply(dps[,trait], dps$spcode, mean)),
        indiv = dps$indiv,
        species = as.numeric(dps$spcode)
      )
      
      ### Add an indicator to set individual effect of non-rep indiv to zero
      repindiv <- names(table(data$indiv))[table(data$indiv) > 1]
      data$indicator <- as.numeric(data$indiv %in% repindiv)
      
      ##############################
      #### Write the model file ####
      ##############################
      setwd("K:/Bob/Panama/MODELS") 
      
      sink("growth_2level_bci.bug")
      cat(" model {
          
          for( i in 1:ntree ) {
          
          growth[i] ~ dnorm(mu[i], tau[1])
          
          mu[i] <- beta.1[species[i]]
          + beta.2[species[i]] * nci[i]
          + beta.3[species[i]] * dbh[i]
          }
          
          for( j in 1:nspecies ) {
          beta.1[j] ~ dnorm(mu.beta[1] + beta.t[1] * trait[j], tau[2])
          beta.2[j] ~ dnorm(mu.beta[2] + beta.t[2] * trait[j], tau[3])
          beta.3[j] ~ dnorm(mu.beta[3], tau[4])
          }
          
          for( m in 1:3 ) {
          mu.beta[m] ~ dnorm(0, 1E-3)
          }
          
          for( b in 1:2 ) {
          beta.t[b] ~ dnorm(0, 1E-3)
          }
          
          for( t in 1:4 ) {
          tau[t] ~ dgamma(1E-3, 1E-3)
          }
          
          sigma <- 1 / sqrt(tau)
          
    }"
    , fill=TRUE)
      sink()
      
      
      sink("growth_2level_coc.she.bug")
      cat(" model {
          
          for( i in 1:ntree ) {
          
          growth[i] ~ dnorm(mu[i], tau[1])
          
          mu[i] <- beta.1[species[i]]
          + beta.2[species[i]] * nci[i]
          + beta.3[species[i]] * dbh[i]
          + indiv.effect[indiv[i]] * indicator[i]
          }
          
          for( j in 1:nspecies ) {
          beta.1[j] ~ dnorm(mu.beta[1] + beta.t[1] * trait[j], tau[2])
          beta.2[j] ~ dnorm(mu.beta[2] + beta.t[2] * trait[j], tau[3])
          beta.3[j] ~ dnorm(mu.beta[3], tau[4])
          }
          
          for( m in 1:3 ) {
          mu.beta[m] ~ dnorm(0, 1E-3)
          }
          
          for( b in 1:2 ) {
          beta.t[b] ~ dnorm(0, 1E-3)
          }
          
          for( i.a in 1:nindiv ) {
          indiv.effect[i.a] ~ dnorm(0, tau[5])
          }
          
          for( t in 1:5 ) {
          tau[t] ~ dgamma(1E-3, 1E-3)
          }
          
          sigma <- 1 / sqrt(tau)
          
  }"
    , fill=TRUE)
      sink()
      ################################################
      ### Set initial values, monitors, iterations and run model ###
      ################################################
      
      if(p!=2){
        inits <- function (){
          list(
            beta.t = rnorm(2),
            mu.beta = rnorm(3),    
            tau = rgamma(5, 1E3, 1E3))
        }
      }
      
      if(p==2){
        inits <- function (){
          list(
            beta.t = rnorm(2),
            mu.beta = rnorm(3),    
            tau = rgamma(4, 1E3, 1E3))
        }
      }
      
      if(pc==T){ 
        setwd("K:/Bob/Panama/MODELS") 
      } else {
        setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
      }
      
      # Set monitors
      params <- c('beta.t','mu.beta','sigma')
      
      # Run model
      adapt <- 2000
      iter <- 25000
      burn <- 22000
      thin <- 5
      chains <- 3
      
      modfile <- ifelse(p!=2, "growth_2level_coc.she.bug", "growth_2level_bci.bug")
      
      print(paste("Now working on:", paste(ifelse(p==1,'coc',ifelse(p==2,'bci','she')), trait, ifelse(size==1,'sm','lg'),sep=".")))
      print(paste('trait =',trait))
      print(paste('p =',p))
      print(paste('size =',size))
      
      
      mod <- jagsUI::jags(data, inits, params, modfile, 
                          n.chains=chains, n.adapt=adapt, n.iter=iter, 
                          n.burnin=burn, n.thin=thin, parallel=T)
      
      for(reps in 1:10){
        if(max(unlist(mod$Rhat)) > 1.1){
          print(paste('Doing update #', reps))
          mod <- update(mod, n.iter=5000)
        }
      }
      
      setwd("K:/Bob/Panama/RESULTS/_12.17.15/growth/nci_dbh_sizes") 
      
      file <- paste(trait, ifelse(p==1,'coc',ifelse(p==2,'bci','she')), ifelse(size==1,'sm','lg'), 'Rdata',sep=".")
      saveRDS(mod, file=file)
      
      }
    }
  }

