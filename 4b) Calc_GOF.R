setwd("/Users/Bob/Projects/Postdoc/Panama/RESULTS/_6.14.16")

get.y.pred <- function(m){
  ie <- ifelse(is.null(m$q50$indiv.effect), 0, m$q50$indiv.effect)
  y.pred <- (m$q50$beta.1[m$data$species]) +
    (m$q50$beta.2[m$data$species] * m$data$allnci) +
    (m$q50$beta.1[m$data$species] * m$data$dbh) + ie
  return(y.pred)
}


ntfiles <- list.files('growth')[grepl('splevel_notrait_growth', list.files('growth'))]
trtfiles <- list.files('growth')[!grepl('splevel_notrait_growth', list.files('growth')) & grepl('splevel', list.files('growth'))]


res <- matrix(ncol=4, nrow=12)

for (i in seq_along(ntfiles)){

ntmod <- readRDS(paste('growth',ntfiles[i],sep='/'))
trtmod <- readRDS(paste('growth',trtfiles[i],sep='/'))

# Get predicted value for each
y.pred.nt <- get.y.pred(ntmod)
y.pred.trt <- get.y.pred(trtmod)

# Get R2
y.obs <- ntmod$data$growth
ntr2 <- summary(lm(y.obs ~ y.pred.nt))$r.squared
trtr2 <- summary(lm(y.obs ~ y.pred.trt))$r.squared
ntdev <- ntmod$q50$deviance
trtdev <- trtmod$q50$deviance
deltaDIC_nt <- (ntdev - min(c(ntdev, trtdev)))
deltaDIC_trt <- (trtdev - min(c(ntdev, trtdev)))


#name <- 'growth', substring(trtfiles[i], nchar(trtfiles[i])-11, nchar(trtfiles[i])-6), 
res[i,] <- c(ntr2, trtr2, deltaDIC_nt, deltaDIC_trt)

}

res <- as.data.frame(res)
res[,1:2] <- round(res[,1:2], 3)
res[,3:4] <- round(res[,3:4], 1)

res$metric <- rep(c('Growth','Survival'), each=6)
res$plot <- rep(rep(c('BCI','Cocoli','Sherman'), each=2), times=2)
res$size <- rep(c('Large','Small'), times=6)
colnames(res)[1:4] <- c('R2_nt','R2_trt','Deviance_nt','Deviance_trt')

index <- paste(ifelse(res$plot=='Cocoli', 1, ifelse(res$plot=='BCI', 2, 3)), 
               ifelse(res$size=='Small', 1, 2), sep='.')

res[order(res$metric, index),]







