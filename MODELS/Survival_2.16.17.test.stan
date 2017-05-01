
    data {
    int<lower=0> N; 			# number of observations
    int<lower=0> nsp; 			# number of species
    int<lower=0> sp[N]; 			# species of each observation
    int<lower=0> nind; 			# number of unique individuals
    int<lower=0> ind[N]; 			# individual ID for each survival observation
    int<lower=0> ncov; 			# number of covariates
    
    vector[N] days; 			# days (years) of interval period
    vector[nsp] meanWD; 			# trait values for each species
    vector[nsp] meanLMA; 			# trait values for each species
    
    matrix[N,ncov] COV; 			# covariate values for each obs. 1st col is all equal to ‘1’ for intercept, 2nd column is x2, …
    int<lower=0,upper=1> alive[N]; 	# survival observations
    }
    
    parameters {
    real t[4]; 				# second level regression of trait effects
    real b[nsp,ncov]; 			# species specific parameters for each covariate effect, including intercept
    real int1[ncov]; 			# hyper distribution means for species specific covariate effects
    real int2[nind]; 			# individual random effects
    real <lower=0.00001> sigmab[ncov]; 	# standard deviations of hyper-distributions for species-specific covariate effects
    real <lower=0.00001> sigmaind; 	# standard deviation of individual random effects
    }
    
    model {
    sigmaind ~ gamma(100,100); 				# prior for individual tree random effects standard deviation
    
    for(i in 1:2) {
    t[i] ~ normal(0,100); 			# priors for 2nd level regression of trait effects
    }
    
    for(i in 1:ncov) {
    int1[i] ~ normal(0,100); 		# priors for hyper distribution means for sp-specific covariate effects
    sigmab[i] ~ gamma(100,100); 		# priors for standard deviations of hyper distributions for sp-specific covariate effects
    }
    
    for(j in 1:nsp) {
    b[j,1] ~ normal(int1[1] + (t[1] * meanWD[j]) + (t[3] * meanLMA[j]), sigmab[1]); 	# sp specific intercept, including second level regression giving trait effect
    b[j,2] ~ normal(int1[2] + (t[2] * meanWD[j]) + (t[4] * meanLMA[j]), sigmab[2]); 	# sp specific covariate effect of root collar diameter
    b[j,3] ~ normal(int1[3], sigmab[3]); 							# sp specific effect of NCIS, i.e. trait similarity with neighbors
    }
    
    for(i in 1:nind) {
    int2[i] ~ normal(0,sigmaind); 		# individual random effects
    }
    
    for(i in 1:N) {
    alive[i] ~ bernoulli(pow(inv_logit(b[sp[i],1]*COV[i,1] 
    + b[sp[i],2]*COV[i,2] 
    + b[sp[i],3]*COV[i,3] 
    + int2[ind[i]]), days[i]));
    }

    }

    
