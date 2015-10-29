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


# Fix some typos found along the way...
out2013$Q20[out2013$Q20 == "02,OO"] <- "02,00"
out2014$Q20[out2014$Q20 == "02,OO"] <- "02,00"




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

# Get rid of lines without species
out <- out[out$SPP != '_',]
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

out$DATE <- as.Date(out$DATE)
out <- droplevels(out)

sort(unique(out$DATE))

######################
######################
######################
### Check out the species codes...
sort(unique(as.character(out$SPP)))

### Manually correct spcode errors:
out$SPP[out$SPP=='ardigl'] <- 'ARDIGL'

### Find stems with > 2 observations overall (manual error check)
x <- names(table(out$uID))[ table(out$uID) > 2 ]
out[order(out$uID),][out$uID[order(out$uID)] %in% x, ]
out <- out[!rownames(out) %in% c(7362, 8026), ]
subtag <- c('a','b','a','b')
out$uID[out$uID == "Sr.04,03.31.189.LICAHY"] <- paste('Sr.04,03.31.189.', subtag, '.LICAHY', sep='')

### Find (and remove) stems with > 2 observations in 1 year (errors?)
x <- rownames(table(out$uID, out$YEAR))[rowSums(table(out$uID, out$YEAR) > 1) > 0]
out[order(out$uID),][out$uID[order(out$uID)] %in% x, ]
out <- out[!out$uID %in% x, ]

### Find stems with blank SPP
### Replace SPP and uID if changed from blank to non-blank between 2013-2014
tmp <- out[out$SPP=='',]
for(i in 1:length(tmp$uID)){
  tmp2 <- out[grep(tmp$uID[i], out$uID),]
  if(length(unique(tmp2$uID)) > 1){
    new.uID <- unique(as.character(tmp2$uID[tmp2$SPP!='']))[1]
    new.SPP <- unique(as.character(tmp2$SPP[tmp2$SPP!='']))[1]
    out$SPP[out$uID==tmp$uID[i]] <- new.SPP
    out$uID[out$uID==tmp$uID[i]] <- new.uID
  }
}

### Identify stems with more than three observations (errors?)
uid <- rownames(table(out$uID, out$YEAR))
uid <- uid[rowSums(table(out$uID, out$YEAR) > 1) > 0]
out[out$uID %in% uid,] # explore the wierdos...

# Manually remove these stems...
out <- out[!(out$uID == "Sr.04,02.33.184.PSYCPO" & out$TAG == 1),]
out <- out[!out$uID %in% c("Sr.04,04.31.172.MICOEL","Sr.04,02.33.80.MABEOC","Sr.01,04.41.73.MYRCGA"),]


### Convert to wide format
wide <- reshape(out, v.names=c("ALT","Alt.2013", "SPP","DAP","CODIGOS",
                "NOTAS","DATE","STATUS"), idvar="uID", drop=c("FILE","TALLO"),
                direction = "wide", timevar='YEAR')

### Convert numeric columns back to numeric...
num.cols <- c("ALT.2013","Alt.2013.2013","DAP.2013","ALT.2014","Alt.2013.2014","DAP.2014")
  for(i in 1:length(num.cols)){
    wide[,num.cols[i]] <- as.numeric(as.character(wide[,num.cols[i]]))  
  }

wide$GROWTH <- wide$ALT.2014 - wide$ALT.2013
wide$DAYS <- as.numeric(wide$DATE.2014 - wide$DATE.2013)
wide$RECRUIT <- ifelse(is.na(wide$STATUS.2013), 2014, NA)
wide$SURVIVE <- ifelse(!is.na(wide$ALT.2013) & !is.na(wide$ALT.2014), 1, ifelse(wide$STATUS.2014 %in% 'new', NA, 0))

data <- wide[,c('uID','SITE','SPP.2014','ALT.2013','ALT.2014','STATUS.2014','GROWTH','DAYS','SURVIVE','CODIGOS.2013','CODIGOS.2014')]
names(data) <- c('uID','SITE','SPCODE','HT13','HT14','STATUS14','GROWTH','DAYS','ALIVE','code13','code14')
data <- droplevels(data)

save(data, file='panama_seedlings/Panama_gradient_seedlings.RDA')


######################
######################
######################



pdf(file='panama_seedlings/Preliminary_plots.pdf')

par(mfrow=c(2,2), oma=c(1,1,0,0), mar=c(4,4,2,2))

plot(table(data$ALIVE, data$SITE)[1,], ylim=c(0,1800), type='b', col=2, pch=21, bg=2, lwd=2, axes=F, xlab='Site', ylab='Number of Individuals')
points(table(data$ALIVE, data$SITE)[2,], type='b', col=3, pch=21, bg=3, lwd=2)
points(tapply(data$ALIVE, data$SITE, function(x) sum(is.na(x))), type='b', col=4, pch=21, bg=4, lwd=2)
axis(1, at=1:length(unique(data$SITE)), labels=unique(data$SITE))
axis(2)
legend("topleft", pch=16, col=c(3,2,4), legend=c('Survive (2013-14)','Die (2013-14)','Recruit (2014)'), bty='n', inset=0.1, lty=1)

comm <- table(data$SITE, data$SPCODE)
library(vegan)
raresp <- rarefy(x=comm, sample=min(rowSums(comm>0)))
plot(raresp, ylab='Rarefied Richness', xlab='Site', axes=F, pch=21, bg='grey', ylim=c(10,50))
axis(1, at=1:length(unique(data$SITE)), labels=unique(data$SITE))
axis(2)

boxplot(data$GROWTH[!rownames(data) %in% grep('R', data$code14)] ~ data$SITE[!rownames(data) %in% grep('R', data$code14)], ylim=c(-300,500), notch=T, ylab='Ht Growth (mm)', xlab='Site', col='grey', axes=F)
abline(h=0, lty=2, col=2, lwd=2)
axis(1, at=1:length(unique(data$SITE)), labels=unique(data$SITE), cex.axis=.75)
axis(2)
box()

dev.off()
