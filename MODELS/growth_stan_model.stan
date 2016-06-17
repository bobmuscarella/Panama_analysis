data {
    int<lower=0> ntree;
    int<lower=0> nindiv;
    int<lower=0> nspecies;
    int<lower=0> indiv[ntree];
    int<lower=0> species[ntree];
    int<lower=0> indicator[ntree];

    vector[ntree] growth;
    vector[ntree] allnci;
    vector[ntree] dbh;
    vector[nspecies] wdmeanz;
    vector[nspecies] wdtauz;
}


parameters {
    real beta1[nspecies];
    real beta2[nspecies];
    real beta3[nspecies];
    real mubeta[3];
    real betat[2];
    real tpred[nspecies];
    real <lower=0.00001> tau[5];
    real indiveffect[nindiv];
}

transformed parameters {
   real mu[ntree];
	for(i in 1:ntree){
		mu[i] <- beta1[species[i]] 
			+ beta2[species[i]] * allnci[i]
			+ beta3[species[i]] * dbh[i]
			+ indiveffect[indiv[i]] * indicator[i];
	}
}

model {

	for(i in 1:ntree) {
		growth[i] ~ normal(mu[i], tau[1]);

	}
	
	for(j in 1:nspecies) {
		beta1[j] ~ normal(mubeta[1] + (betat[1] * tpred[j]), tau[2]);
		beta2[j] ~ normal(mubeta[2] + (betat[2] * tpred[j]), tau[3]);
		beta3[j] ~ normal(mubeta[3], tau[4]);

		tpred[j] ~ normal(wdmeanz[j], wdtauz[j]);
		}
	}

## Priors
	for(i in 1:3){
		mubeta[i] ~ normal(0,100);
	}

	for(b in 1:2){
		betat[b] ~ normal(0, 1E-3);
	}

	for(t in 1:5){
		tau[t] ~ gamma(100,100);
	}

	for(i.a in 1:nindiv){
		indiveffect[i] ~ normal(0, tau[5]);
	}
}