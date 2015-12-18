#######################################
###  SURVIVAL ANALYSIS
###  ALL-PLOTS
###  ONLY NCI, DBH BY SIZE CLASS
###  Results in "K:/Bob/Panama/RESULTS/_12.17.15/survival/3level/"
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
  }
  
    dp <- d
      
    for(size in 1:2) {
      dps <- dp[dp$size.class %in% size,]
      
      dps <- dps[,c('spcode','plot','size.class','census','survival','days','log.dbh.z','id','log.nci.z',
                    trait)]
      
      # Drop factors for correct indexing
      dps <- droplevels(dps)
      
      dps$speciesxplot <- paste(dps$plot, dps$spcode, sep='.')
      
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
        nspecies = length(unique(dps$speciesxplot)),
        survival = as.numeric(dps$survival),
        days = as.numeric(dps$days),
        nci = as.numeric(dps[,'log.nci.z']),
        ncensus = length(unique(paste(dps$census, dps$plot, sep=''))),
        census = as.numeric(as.factor(paste(dps$census, dps$plot, sep=''))),  
        dbh = as.numeric(dps$log.dbh.z),
        trait = z.score(tapply(dps[,trait], dps$speciesxplot, mean)),
        indiv = dps$indiv,
        species = as.numeric(as.factor(dps$speciesxplot)),
        nplot = length(levels(as.factor(dps$plot))),
        plot = as.numeric(substring(names(tapply(dps[,trait], dps$speciesxplot, mean)),1,1))
      )
      
      ### Add an indicator to set individual effect of non-rep indiv to zero
      repindiv <- names(table(data$indiv))[table(data$indiv) > 1]
      data$indicator <- as.numeric(data$indiv %in% repindiv)
      
      ##############################
      #### Write the model file ####
      ##############################
      setwd("K:/Bob/Panama/MODELS") 
      
      sink("survival_3level_NCI.bug")
      
      cat(" model {
          
          for( i in 1:ntree ) {

          alive[i] ~ dbern(t[i])
          
          t[i] <- pow(z[i], days[i]/365.25)
          
          logit(z[i]) <- beta.1[species[i]]
          + beta.2[species[i]] * nci[i]
          + beta.3[species[i]] * dbh[i]
          + indiv.effect[indiv[i]] * indicator[i]
          + census.effect[census[i]]
          }
          
          for( j in 1:nspecies ) {
          beta.1[j] ~ dnorm(mu.beta.1[plot[j]] + beta.t.1[plot[j]] * trait[j], tau[1])
          beta.2[j] ~ dnorm(mu.beta.2[plot[j]] + beta.t.2[plot[j]] * trait[j], tau[2])
          beta.3[j] ~ dnorm(mu.beta.3[plot[j]], tau[3])
          }
          
          for( m in 1:nplot ) {
          mu.beta.1[m] ~ dnorm(0, 1E-3)
          mu.beta.2[m] ~ dnorm(0, 1E-3)
          mu.beta.3[m] ~ dnorm(0, 1E-3)
          }

          for( p in 1:nplot ) {
          beta.t.1[p] ~ dnorm(0, 1E-3)
          beta.t.2[p] ~ dnorm(0, 1E-3)
          }
          
          for( c.a in 1:ncensus ) {
          census.effect[c.a] ~ dnorm(0, tau[4])
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

      inits <- function (){
        list(
          beta.t.1 = rnorm(3),
          beta.t.1 = rnorm(3),
          mu.beta.1 = rnorm(3),    
          mu.beta.2 = rnorm(3),    
          mu.beta.3 = rnorm(3),    
          tau = rgamma(5, 1E3, 1E3))
      }
      
      setwd("K:/Bob/Panama/MODELS") 

      # Set monitors
      params <- c('beta.t.1','beta.t.2','mu.beta.1','mu.beta.2','mu.beta.3','sigma')
      
      # Run model
      adapt <- 2000
      iter <- 25000
      burn <- 22000
      thin <- 5
      chains <- 3
      
      modfile <- "survival_3level_NCI.bug"
      
      print(paste("Now working on:", trait, ifelse(size==1,'sm','lg'),sep="."))
      
      mod <- jagsUI::jags(data, inits, params, modfile, 
                          n.chains=chains, n.adapt=adapt, n.iter=iter, 
                          n.burnin=burn, n.thin=thin, parallel=T)
      
      for(reps in 1:10){
        if(max(unlist(mod$Rhat)) > 1.1){
          print(paste('Doing update #', reps))
          mod <- update(mod, n.iter=5000)
        }
      }
      
      setwd("K:/Bob/Panama/RESULTS/_12.17.15/survival/3level/") 
      
      file <- paste('allplots', trait, ifelse(size==1,'sm','lg'), 'Rdata',sep=".")
      saveRDS(mod, file=file)
      
      }
    }

