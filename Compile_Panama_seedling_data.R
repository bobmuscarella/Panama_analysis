library(gdata)
setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")

#####################################
###   READ AND PROCESS 2013 DATA   ###
#####################################

out2013 <- data.frame()
files0 <- list.files('panama_seedlings')[1]
files1 <- list.files(paste('panama_seedlings', files0, sep='/'))

for (f1 in 1:length(files1)){

	files2 <- list.files(paste('panama_seedlings', files0, files1[f1], sep='/'))	

		for (f2 in 1:length(files2)){				

			tmp <- read.xls(paste("panama_seedlings", files0, files1[f1], files2[f2], sep="/"))
			x <- sort(as.character(tmp[(nrow(tmp)-3):nrow(tmp), ncol(tmp)]))
			date <- x[x != '_'][1]								
			tmp <- tmp[ !(as.character(tmp[,1]) == "" | is.na(as.character(tmp[,1]))) , ]
			colnames(tmp) <- NULL
			tmp$date <- date
			tmp$year <- files0
			tmp$file <- paste(files0, files1[f1], files2[f2], sep="/")
			tmp$status <- 'new'
			out2013 <- rbind(out2013, tmp)
			print(paste("panama_seedlings", files0, files1[f1], files2[f2], sep="/"))
	}
}
colnames(out2013) <- c('Q20','P5','TAG','SPP','ALT',
										'DAP','TALLO','CODIGOS','NOTAS',
										'DATE','YEAR','FILE','STATUS')

Alt.2013 <- rep(NA, nrow(out2013))
out2013 <- as.data.frame(cbind(out2013[1:4], Alt.2013, out2013[,5:ncol(out2013)]))



#####################################
###   READ AND PROCESS 2014 DATA   ###
#####################################

out2014 <- data.frame()

files0 <- list.files('panama_seedlings')[2]
files1 <- list.files(paste('panama_seedlings', files0, sep='/'))

for (f1 in 1:length(files1)){

	files2 <- list.files(paste('panama_seedlings', files0, files1[f1], sep='/'))

	for (f2 in 1:length(files2)){

		tmp1 <- read.xls(paste("panama_seedlings", files0, files1[f1], files2[f2], sep="/"), sheet=1)[,1:10]
		tmp1$date <- as.character(tmp1 [ , ncol(tmp1)] [tmp1[ , (ncol(tmp1) -1 )] == 'FECHA:'][1])
		tmp1$year <- files0
		tmp1 <- tmp1[ unlist(lapply(as.character(tmp1[,1]), nchar)) == 5, ]				
		tmp1 <- tmp1[ tmp1[,2] != "" , ]
		tmp1$file <- paste(files0, files1[f1], files2[f2], sep="/")
		tmp1$status <- 'old'
		head(tmp1)
		names(tmp1) <- c('Q20','P5','TAG','SPP','Alt.2013',
		'ALT','DAP','TALLO','CODIGOS','NOTAS','DATE','YEAR','FILE','STATUS')
		
		tmp2 <- read.xls(paste("panama_seedlings", files0, files1[f1], files2[f2], sep="/"), sheet=2)[,1:9]
		tmp2$date <- as.character(tmp2 [ , ncol(tmp2)] [tmp2[ , (ncol(tmp2) -1 )] == 'FECHA:'][1])
		tmp2$year <- files0
		tmp2 <- tmp2[ unlist(lapply(as.character(tmp2[,1]), nchar)) == 5, ]				
		tmp2 <- tmp2[ tmp2[,2] != "" , ]
		tmp2$file <- paste(files0, files1[f1], files2[f2], sep="/")
		tmp2$status <- 'new'

		names(tmp2) <- c('Q20','P5','TAG','SPP','ALT','DAP','TALLO',
		'CODIGOS','NOTAS','DATE','YEAR','FILE','STATUS')

		Alt.2013 <- rep(NA, nrow(tmp2))
		tmp2 <- cbind(tmp2[,1:4], Alt.2013, tmp2[,5:ncol(tmp2)])

		rownames(tmp1) <- rownames(tmp2) <- NULL

		mat <- as.data.frame(rbind(tmp1, tmp2))
		out2014 <- rbind(out2014, mat)
		print(paste("panama_seedlings", files0, files1[f1], files2[f2], sep="/"))
		
	}

}

################################
###   COMBINE 2013 / 2014 DATA   ###
################################
out <- rbind(out2013, out2014)

sites <- c('Buena', 'Charco', 'Metro', 'Oleo', 'Pacifico', 'Santa', 'Sherman', 'Soberania')
sitenames <- c('Bu', 'Ch', 'Me', 'Ol', 'Pa', 'Sr', 'Sh', 'So')

out$SITE <- NA
for (i in 1:length(sites)){
	out$SITE[grep(sites[i], out$FILE)] <- sitenames[i]
}

out$uID <- paste(out$SITE, out$Q20, out$P5, out$TAG, out$SPP, sep='.')
out <- out[out$SPP != '_',]
out <- out[out$SPP != '',]

out <- droplevels(out)

sort(unique(out$DATE))
out$DATE[out$DATE %in% "13/Dic/13"] <- "2013-12-13"
out$DATE[out$DATE %in% "13/Dic/2013"] <- "2013-12-13"
out$DATE[out$DATE %in% "14-Oct_2013"] <- "2013-10-14"
out$DATE[out$DATE %in% "16/Dic/13"] <- "2013-12-16"
out$DATE[out$DATE %in% "17-Oct_2013"] <- "2013-10-17"
out$DATE[out$DATE %in% "17/Dic/13"] <- "2013-12-17"
out$DATE[out$DATE %in% "17/Dic/2013"] <- "2013-12-17"
out$DATE[out$DATE %in% "18/Dic/2013"] <- "2013-12-18"
out$DATE[out$DATE %in% "19/Dic/2013"] <- "2013-12-19"
out$DATE[out$DATE %in% "28-Oct_2013"] <- "2013-10-28"
out <- droplevels(out)

out$Date <- as.Date(out$DATE)
sort(unique(out$DATE))

################################
###   WORKING FRONT....


# Subset uIDs that occurred in both years (not necessarily alive)
test <- out[out$uID %in% names(table(out$uID))[table(out$uID)==2],]
test <- test[order(test$uID, test$YEAR),]
head(test)

# Subset uIDs that were only observed once
singleobs <- names(table(out$uID))[ table(out$uID) == 1 ]
singleobsyear <- out$YEAR[out$uID %in% singleobs]

# How many total deaths were there?
sum(table(out$CODIGOS)[names(table(out$CODIGOS)) %in% c('D','T','N','D,P','T,P','N,P')])

# Note as alive=T if not one of these death codes
out$alive <- !out$CODIGOS %in% c('D','T','N','D,P','T,P','N,P')

# Unelegant way to remove a couple stems with 3 observations (appear to be errors)
threeobs <- names(table(out$uID))[ table(out$uID) > 2 ]
out[out$uID %in% threeobs,]
out <- out[!rownames(out) %in% c(7362, 8026), ]
subtag <- c('a','b','a','b')
out$uID[out$uID == "Sr.04,03.31.189.LICAHY"] <- paste('Sr.04,03.31.189.', subtag, '.LICAHY', sep='')


# Post hoc removal of stems that have errors
out <- out[!out$uID %in% "Sr.04,03.31.193.PSYCB2",]
out <- out[!out$uID %in% "Ol.02,04.42.118.CROTBI",]
out <- out[!out$uID %in% "Ol.02,04.42.121.CROTBI",]
out <- out[!out$uID %in% "Ol.02,04.42.102.CROTBI",]
out <- out[!out$uID %in% "Ol.02,04.42.165.CROTBI",]

head(out)
head(out[out$alive==F,])
table(out$uID, out$YEAR)

out <- droplevels(out)

res <- matrix(nrow=length(unique(out$uID)), ncol=14)

for(i in 1:length(unique(out$uID))){

	ind <- unique(out$uID)[i]
	tmp <- out[out$uID == ind,]
	tmp$ALT <- as.numeric(as.character(tmp$ALT))
	tmp$Alt.2013 <- as.numeric(as.character(tmp$Alt.2013))  
  
	if(length(unique(tmp$SPP)) == 1){
	  SPP <- as.character(unique(tmp$SPP))
	} else {
	  SPP <- "SPP_CONFLICT"
	}
  
  if(length(tmp$ALT[tmp$YEAR==2013]) > 0){
		ALT2013 <- as.numeric(as.character(tmp$ALT[tmp$YEAR==2013]))
	} else {
		ALT2013 <- NA
	}

	if(length(tmp$ALT[tmp$YEAR==2014]) > 0){
		ALT2014 <- as.numeric(as.character(tmp$ALT[tmp$YEAR==2014]))
	} else {
		ALT2014 <- NA
	}

	if(length(tmp$DATE[tmp$YEAR==2013]) > 0){
		DATE2013 <- as.Date(tmp$DATE[tmp$YEAR==2013], format="%Y-%m-%d")
	} else {
		DATE2013 <- NA
	}

	if(length(tmp$DATE[tmp$YEAR==2014]) > 0){
		DATE2014 <- as.Date(tmp$DATE[tmp$YEAR==2014], format="%Y-%m-%d")
	} else {
		DATE2014 <- NA
	}
	
	if(!is.na(DATE2013) & !is.na(DATE2014)){
		INT <- as.numeric(DATE2014 - DATE2013)
	} else {
		INT <- NA
	}
	
	GROWTH <- ALT2014 - ALT2013
	SURVIVE <- !is.na(ALT2013) & !is.na(ALT2014)
	RECRUIT <- is.na(ALT2013) & !is.na(ALT2014)
	DIE <- !is.na(ALT2013) & is.na(ALT2014)
	STATUS <- c('s','r','d')[c(SURVIVE, RECRUIT, DIE)]
	
	STATUS <- ifelse(length(STATUS) > 0, STATUS, NA)

	SITE <- unique(tmp$SITE)
	Q20 <- unique(tmp$Q20)
	P5 <- unique(tmp$P5)
	TAG <- unique(tmp$TAG)
  CODES <- paste(tmp$CODIGOS[tmp$CODIGOS != "_"], collapse=', ')

	res[ i , ] <- c(ind, SPP, ALT2013, ALT2014, GROWTH, DATE2013, DATE2014, INT, STATUS, SITE, Q20, P5, TAG, CODES)
	print( paste(i, 'out of', length(unique(out$uID))))
}

res <- as.data.frame(res)

res <- res[!is.na(res[,1]),]

names(res) <- c('uID','spcode','ht13','ht14','growth','date13','date14','int','status','site','q20','p5','tag','codes')


par(mfrow=c(2,2))


res$growth <- as.numeric(as.character(res$growth))
boxplot(res$growth[!is.na(res$growth) & !rownames(res) %in% grep('R', res$codes)] ~ res$site[!is.na(res$growth) & !rownames(res) %in% grep('R', res$codes)])





tmp <- res[!rownames(res) %in% grep('R', res$codes),]
plot(tapply(tmp$growth, tmp$site, median, na.rm=T), ylim=c(-1000,1000))
q10 <- tapply(tmp$growth, tmp$site, quantile, na.rm=T, prob=0.2)
q90 <- tapply(tmp$growth, tmp$site, quantile, na.rm=T, prob=0.8)
segments(1:8, MIN, 1:8, MAX)




tmp <- res[!is.na(res$growth) & !rownames(res) %in% grep('R', res$codes),]
table(tmp$site, tmp$spcode)


### SOME PROBABLE SPECIES CODE ERRORS: ###
### REALLY I NEED A TABLE OF REAL SPECIES CODES... ###
# "ardigl" NEEDS TO BE "ARDIGL"
# "BIGNONACEAE" AS A SPCODE???
# SHOULD "CASTLE" BE "CASTEL"?


head(res)




pdf(file='test.pdf')

par(mfrow=c(2,2))

plot(table(res$status, res$site)[1,], ylim=c(0,1200), type='b', col=2, pch=21, bg=2, lwd=2, axes=F, xlab='Site', ylab='Number of Individuals')
points(table(res$status, res$site)[3,], type='b', col=3, pch=21, bg=3, lwd=2)
points(table(res$status, res$site)[2,], type='b', col=4, pch=21, bg=4, lwd=2)
axis(1, at=1:length(unique(res$site)), labels=unique(res$site))
axis(2)
legend("topleft", pch=16, col=c(3,2,4), legend=c('Survive (2013-14)','Die (2013-14)','Recruit (2014)'), bty='n', inset=0.1, lty=1)

n.uID <- tapply(res$uID, res$site, function(x) length(unique(x)))
plot(n.uID, ylim=c(0,2000), pch=21, bg='grey', ylab='Total number of unique individuals', axes=F, xlab='Site')
axis(1, at=1:length(unique(res$site)), labels=unique(res$site))
axis(2)

comm <- table(res$site, res$spcode)
library(vegan)
raresp <- rarefy(x=comm, sample=min(rowSums(comm>0)))
plot(raresp, ylab='Rarefied Species Richness', xlab='Site', axes=F, pch=21, bg='grey', ylim=c(10,50))
axis(1, at=1:length(unique(res$site)), labels=unique(res$site))
axis(2)

boxplot(res$growth[!rownames(res) %in% grep('R', res$codes)] ~ res$site[!rownames(res) %in% grep('R', res$codes)], ylim=c(-200,500), notch=T, ylab='Growth (Alt)', xlab='Site', col='grey')
abline(h=0, lty=2, col=2, lwd=2)

dev.off()





tmp <- res[!rownames(res) %in% grep('R', res$codes),]

plot(table(tmp$status, tmp$site)[1,], ylim=c(0,1200), type='b', col=2, pch=21, bg=2, lwd=2, axes=F, xlab='Site', ylab='Number')
points(table(tmp$status, tmp$site)[3,], type='b', col=3, pch=21, bg=3, lwd=2)
points(table(tmp$status, tmp$site)[2,], type='b', col=4, pch=21, bg=4, lwd=2)
axis(1, at=1:length(unique(tmp$site)), labels=unique(tmp$site))
axis(2)
legend("topleft", pch=16, col=c(3,2,4), legend=c('Survive','Die','Recruit'), bty='n', inset=0.1, lty=1)






