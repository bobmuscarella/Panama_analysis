### COMPILE DATA FROM RECENSUSED 1-HA PLOTS IN PANAMA

setwd("/Users/Bob/Projects/Postdoc/Panama/DATA/small plots")

out <- data.frame()

for(i in 1:length(list.files())){
  plot <- list.files()[i]
  f0 <- list.files(plot)
  f1 <- (1:length(f0))[-grep('tax', f0)]
  files <- f1[-grep('zip', f0[f1])]
  infiles <- paste(list.files()[i], f0[files], sep="/")
    for (j in 1:length(infiles)){
      tmp <- read.table(infiles[j], sep='\t', header=T)[,-1]
      tmp$plot <- plot
      tmp$census <- substring(strsplit(infiles[j], "_")[[1]][2],2,2)
      tmp$year <- as.numeric(format(as.Date(tmp$Date),'%Y'))
      out <- rbind(out, tmp)
  }
}
  
head(out)
  
table(out$plot, out$year)

unique(out$plot)
plot <- unique(out$plot)[8]
plot(out$PX[out$plot==plot], out$PY[out$plot==plot], col=out$census[out$plot==plot])

hist(out$DBH, breaks=1000, xlim=c(0,200))


