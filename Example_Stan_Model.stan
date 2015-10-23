# THIS IS THE STAN MODEL FOR SURVIVAL FROM JESSE
# I'VE ADDED A SECTION TO SIMULATE DATA TO TEST IT OUT...

stanmod2 <- '
data {
int<lower=0> nspecies;
int<lower=0> ntree_die;
int<lower=0> ntree_live;
int<lower=1,upper=nspecies> species_die[ntree_die];
int<lower=1,upper=nspecies> species_live[ntree_live];
real<lower=0> time_die[ntree_die];
real<lower=0> time_live[ntree_live];
int<lower=0> nind;
int<lower=1,upper=nind> ind_die[ntree_die];
int<lower=1,upper=nind> ind_live[ntree_live];
int<lower=0> ncens;
real precip[ncens];
int<lower=1,upper=ncens> census_die[ntree_die];
int<lower=1,upper=ncens> census_live[ntree_live];
real coll_live[ntree_live];
real coll_die[ntree_die];
int<lower=0> nplot;
int<lower=1,upper=nplot> plot_die[ntree_die];
int<lower=1,upper=nplot> plot_live[ntree_live];
real cons_live[ntree_live];
real cons_die[ntree_die];
real light_live[ntree_live];
real light_die[ntree_die];

}

parameters {
real<lower=0> r;
real beta[nspecies];
real mu;
real<lower=0,upper=1000> sigma1;
real ind_effect[nind];
real<lower=0,upper=1000> sigma2;
real beta_coll;
real beta_plot[nplot];
real<lower=0,upper=1000> sigma3;
real beta_cons[nspecies];
real mu_cons;
real<lower=0,upper=1000> sigma4;
real beta_light[nspecies];
real mu_light;
real<lower=0,upper=1000> sigma5;
real mu_prec;
real beta_prec[nspecies];
real<lower=0,upper=1000> sigma6;
real mu_prec_cons;
real beta_prec_cons[nspecies];
real<lower=0,upper=1000> sigma7;

}

model {
r ~ exponential(0.001); #priors
mu ~ normal(0, 100);
mu_prec ~ normal(0, 100);
beta_coll ~ normal(0, 100);
mu_cons ~ normal(0, 100);
mu_light ~ normal(0, 100);
mu_prec_cons ~ normal(0, 100);

for (n in 1:ntree_live) { #liklehoods
increment_log_prob(weibull_ccdf_log(time_live[n], r, exp(-(
		beta[species_live[n]] 
	+ ind_effect[ind_live[n]] 
	+ beta_prec[species_live[n]] * precip[census_live[n]] 
	+ beta_coll * coll_live[n] 
	+ beta_plot[plot_live[n]] 
	+ beta_cons[species_live[n]] * cons_live[n] 
	+ beta_light[species_live[n]] * light_live[n]
	+ beta_prec_cons[species_live[n]] * precip[census_live[n]] * cons_live[n]	
	)/r))); #for those that survive
}
for (n in 1:ntree_die) {
increment_log_prob(weibull_cdf_log(time_die[n], r, exp(-(
	beta[species_die[n]] 
	+ ind_effect[ind_die[n]] 
	+ beta_prec[species_die[n]] * precip[census_die[n]] 
	+ beta_coll * coll_die[n] 
	+ beta_plot[plot_die[n]] 
	+ beta_cons[species_die[n]] * cons_die[n] 
	+ beta_light[species_die[n]] * light_die[n] 
	+ beta_prec_cons[species_die[n]] * precip[census_die[n]] * cons_die[n]	
	)/r))); #for those that die i essentially want the cdf
}

beta ~ normal(mu, sigma1); #hyperdistributions
beta_plot ~ normal(0, sigma3);
beta_cons ~ normal(mu_cons, sigma4); #hyperdistributions
beta_light ~ normal(mu_light, sigma5); #hyperdistributions
beta_prec ~ normal(mu_prec, sigma6); #hyperdistributions
beta_prec_cons ~ normal(mu_prec_cons, sigma7); #hyperdistributions

ind_effect ~ normal(0, sigma2);

}

'

## TEST WITH SOME SIMULATED DATA...
ndie <- 100
nlive <- 200

d <- list(nspecies=2, 
			ntree_die=ndie, 
			ntree_live=nlive, 
			species_die=sample(1:2, ndie, replace=T),
			species_live=sample(1:2, nlive, replace=T),
			time_die=abs(rnorm(ndie)),
			time_live=abs(rnorm(nlive)+1),
			nind=ndie+nlive,
			ind_die=1:ndie,
			ind_live=(ndie+1):(nlive+ndie),
			ncens=2,
			precip=c(0.1,0.2),
			census_die=sample(1:2, ndie, replace=T),
			census_live=sample(1:2, nlive, replace=T),
			coll_die=abs(rnorm(ndie)+.1),
			coll_live=abs(rnorm(nlive)),
			nplot=3,
			plot_die=sample(1:3, ndie, replace=T),
			plot_live=sample(1:3, nlive, replace=T),
			cons_die=rnorm(ndie),
			cons_live=rnorm(nlive),
			light_die=rnorm(ndie),
			light_live=rnorm(nlive)
			)

inits <- function() { list ( r=2 ) }
inits()

fit <- stan(model_code = stanmod2, data = d, iter = 1000, chains = 2, thin = 50, init=inits)


