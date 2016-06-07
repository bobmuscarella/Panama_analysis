
# Generate dbhs
dbh1 <- rlnorm(1000, 2)
dbh2 <- dbh1 + dbh1*abs(rnorm(1000))

plot <- sample(1:5, 1000, replace=T)

# Sum DBH of stems larger than focal stem
bl <- vector(length=length(dbh1))
for(i in 1:length(dbh1)){
  tdbh <- dbh1[-i][plot==plot[i]]
  bl[i] <- sum(tdbh[tdbh > dbh[i]])
}

# Sum DBH of all other stems 
bt <- vector(length=length(dbh1))
for(i in 1:length(dbh1)){
  tdbh <- dbh1[plot==plot[i]]
  bt[i] <- sum(tdbh)
}

plot(bt, bl, cex=dbh/4, ylim=c(0,1600),xlim=c(0,1600))
plot(dbh1, dbh2-dbh1)

# GROWTH
G <- dbh2-dbh1

# FUNCTION FOR ABOVE-GROUND COMPETITION
(lambda[1] * dbh1^alpha) / (1 + ((lambda[1]/lambda[2]) * exp(1) ^ (lambda[3]*bl)))

alpha ~ dnorm(0,1E-3)
lambda[1] ~ dnorm(0,1E-3)
lambda[2] ~ dnorm(0,1E-3)
lambda[3] ~ dnorm(0,1E-3)


# FUNCTION FOR BELOW-GROUND COMPETITION
(lambda[1] * dbh1^alpha) / (1 + (lambda[4]*bt))

alpha ~ dnorm(0,1E-3)
lambda[1] ~ dnorm(0,1E-3)
lambda[2] ~ dnorm(0,1E-3)
lambda[3] ~ dnorm(0,1E-3)
lambda[4] ~ dnorm(0,1E-3)


# FUNCTION FOR COMBINED ABOVE AND BELOW-GROUND COMPETITION
(lambda[1] * dbh1^alpha) / ((1 + (lambda[1]/lambda[2])*(exp(1)^(lambda[3]*bl)))*(1 + (lambda[4]*bt)))

alpha ~ dnorm(0,1E-3)
lambda[1] ~ dnorm(0,1E-3)
lambda[2] ~ dnorm(0,1E-3)
lambda[3] ~ dnorm(0,1E-3)
lambda[4] ~ dnorm(0,1E-3)








nci <- seq(1, 100, 1)
dbh1 <- rep(-1, 100)
dbh2 <- rep(1, 100)

y1 <-  (-0.25 * nci) + (0.5 * dbh1) + (0.1 * nci * dbh1)
y2 <-  (-0.25 * nci) + (0.5 * dbh2) + (0.1 * nci * dbh2)
plot(nci, y1, ylim=c(-50,100), type='l')
points(nci, y2, col=4, type='l')
abline(h=0, lty=2)
legend('topleft', legend=c('Small DBH','Large DBH'), col=c(1,4), lty=1)





library(maptools)
setwd("/Users/Bob/Projects/Thesis/GIS Work/data")
x <- readShapePoints("Sites")
coordinates(x)
