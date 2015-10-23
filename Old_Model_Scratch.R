# THIS SCRATCH PAD HAS CODE FROM OLD (SIMPLER) MODELS USED IN DEVELOPMENT OF THE CURRENT MODELS...

setwd("K:/Bob/Panama/MODELS")
setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")

################################
#### 2 LEVEL WITHOUT TRAITS
################################
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
	}

params <- c("beta.1","beta.2","beta.3","mu.beta","tau")
mod <- jagsUI::jags(data, parameters.to.save=params, 
                    model.file="growth_simple.bug", n.iter=1000, 
                    n.chains=2, parallel=T)


################################
#### 3 LEVEL WITHOUT TRAITS
################################
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

params <- c("beta.1","beta.2","beta.3","mu.beta","tau")
setwd("/Users/Bob/Projects/Postdoc/Panama/MODELS")
mod <- jagsUI::jags(data, parameters.to.save=params, 
                    model.file="growth_3level_notrait.bug", 
                    n.iter=1000, n.chains=2, parallel=T, inits=inits())
