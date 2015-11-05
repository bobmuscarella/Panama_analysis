library(reshape)

###########################################
### COMPILE PANAMA SMALL PLOT TREE DATA ###
###########################################
setwd("/Users/Bob/Projects/Postdoc/Panama/DATA")

sites <- list.files('small plots')
res <- data.frame()

for (i in 1:length(sites)){
  subfiles <- list.files(paste('small plots/',sites[i],sep=''))
  subfiles2 <- subfiles[grep('.txt', subfiles)] 
    for (f in 1:(length(subfiles2)-1)){
      tab <- read.table(paste('small plots',sites[i],subfiles2[f],sep='/'), header=T, sep='\t')[,-1]
      tab$site <- sites[i]
#      tmp <- unlist(strsplit(subfiles2[f], '[.]'))[1]
#      tab$census <- substring(tmp, nchar(tmp), nchar(tmp))
      res <- rbind(res, tab)
  }  
}
res$Date <- as.Date(res$Date)
res$year <- format(res$Date,'%Y')

### WORK ONLY WITH MAIN STEMS
res <- res[res$Stem %in% 'main',]

### TRY WORKING ONLY WITH CENSUS 1 and 2 FOR NOW...
res <- res[res$Census %in% c(1,2),]

### Make a unique ID
res$uID <- paste(res$site, res$TreeID, res$StemID, sep='_')
res <- droplevels(res)
wide <- reshape(res, v.names=c("Latin","Quadrat","PX","PY","TreeID","StemID","year",
        "DBH","HOM","Date","Codes","Status","site","year"), idvar=c("uID"), 
        drop=c("SubSpecies","Tag","StemTag","Stem"), direction = "wide", 
        timevar="Census")

wide$int <- as.numeric(wide$Date.2 - wide$Date.1)

### WORK ONLY WITH STEMS ALIVE AT TIME 1
wide <- wide[wide$Status.1=='alive',]

### WORK ONLY WITH STEMS WITH INTERVAL <= 5 years
(wide$int/365)[(wide$int/365)<=6]
wide[(wide$int/365)]

test <- wide[!is.na(wide$int) & (wide$int/365) <= 6,]
dim(test)
head(test)

table(test$site.1)


sum(!is.na(test$Status.2))

hist(wide$DBH.1[!is.na(wide$Status.2)], col=2)
hist(wide$DBH.1[is.na(wide$Status.2)], add=T, col=rgb(0,0,1,.5))
hist(wide$DBH.1[wide$Status.2=='alive'], add=T, col=rgb(1,1,1,.5))

# SEE THE PLOT LAYOUT
unique(res$site)
plot <- unique(res$site)[8]
plot(res$PX[res$site==site], res$PY[res$site==site], col=res$Census[res$site==site])



