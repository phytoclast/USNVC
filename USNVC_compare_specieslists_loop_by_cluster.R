
library(stringr)
library(BiodiversityR)
library(cluster)
library(ape)
library(dendextend)
library(dplyr)
library(dynamicTreeCut)
library(rpart)
library(rpart.plot)
library(goeveg)
library(proxy)
library(isopam)
library(optpart)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

betasim <- function(p){
  d <- matrix(1, nrow = nrow(p), ncol = nrow(p))
  rownames(d) <- rownames(p)
  colnames(d) <- rownames(p)
  for(j in 1:nrow(p)){
    for(k in 1:nrow(p)){
      d[j,k] <- 1-sum((p[j,]*p[k,])^0.5)/sqrt(sum(p[j,])*sum(p[k,]))
    }}
  d<-as.dist(d)
}
makeplot <- function(amethod,d,t,k){
  filename <- paste0('output/vegplots_',amethod,'.png')
  t <- as.hclust(t)
  #make cuts and reformat dendrogram
  ngroups=k
  groups <- cutree(t, k = k)
  
  soilplot <- names(groups)
  clust <- unname(groups)
  groupdf <- as.data.frame(cbind(soilplot, clust))
  groupdf$clust <- (as.numeric(as.character(groupdf$clust)))
  maxcluster <- max(groupdf$clust)
  numberzeros <- nrow(groupdf[(groupdf$clust == 0),])
  whichrecords <- which(groupdf$clust == 0)
  if (nrow(groupdf[groupdf$clust == 0,]) != 0){
    for (i in 1:numberzeros){ #assign all zero clusters to unique cluster number.
      groupdf[whichrecords[i],]$clust <- maxcluster+i}}
  
  newlabels <- t$labels
  newlabels <- as.data.frame(newlabels)
  newlabels$row <- row(newlabels)
  newlabels <- merge(newlabels, groupdf, by.x='newlabels', by.y ='soilplot')
  newlabels$newlabels <- paste(newlabels$clust, newlabels$newlabels)
  newlabels <- newlabels[order(newlabels$row),1]
  newtree <- t
  newtree$labels <- newlabels
  
  dend1 <- color_branches(as.hclust(newtree), k = k)
  dend1 <- color_labels(dend1, k = k)
  
  #output file
  
  w <- 800
  h <- nrow(d)*12+80
  u <- 12
  png(filename=filename,width = w, height = h, units = "px", pointsize = u)
  
  par(mar = c(2,0,1,13))
  plot(dend1, horiz = TRUE, main=paste('floristic simularity', amethod,' method of plots'), font=1, cex=0.84)
  dev.off()
  
}
ecoregion <-  read.delim('data/d_usfs_ecoregion2007.txt', encoding = 'UTF-8', na.strings = '', stringsAsFactors=FALSE)
vegecoregion <-  read.delim('data/UnitXEcoregionUsfs2007.txt', encoding = 'UTF-8', na.strings = '', stringsAsFactors=FALSE)
vegecoregion <- merge(vegecoregion, ecoregion, by='usfs_ecoregion_2007_id')

unit <-  read.delim('data/unit.txt', encoding = 'UTF-8', na.strings = '', stringsAsFactors=FALSE)
states <-read.delim('data/d_subnation.txt')
vegstates <-read.delim('data/UnitXSubnation.txt')
vegstates <- merge(vegstates, states, by='subnation_id')

USNVClist <- readRDS('data/USNVClist.RDS')
USNVClist <- subset(USNVClist, !grepl('\\.', acctaxon))
plotdata <- readRDS('data/Com.Sp.mean.RDS')
plotdata$soilplot <- str_replace_all(plotdata$soilplot, '\\)', '.')
plotdata$soilplot <- str_replace_all(plotdata$soilplot, '\\(', '.')
#plotdata <- subset(plotdata, !Observation_Label %in% c('S12062901', 'W14090801'))
if(F){plotdata$over <- 100*(1-10^(apply(log10(1-(plotdata[,c('Subcanopy', 'Tree')]/100.001)), MARGIN = 1, FUN='sum')))
plotdata.summary <- aggregate(
  list(Trees = log10(1-(plotdata[,c('over')]/100.001)), 
       Shrubs = log10(1-(plotdata[,c('Shrub')]/100.001)), 
       Field = log10(1-(plotdata[,c('Field')]/100.001))),
       by= list(Observation_ID = plotdata$Observation_ID), FUN='sum')
plotdata.summary[,c('Trees','Shrubs','Field')] <- 100*(1-10^plotdata.summary[,c('Trees','Shrubs','Field')])
Woodland <-  subset(plotdata.summary, Trees >= 25)$Observation_ID
Shrubland <-  subset(plotdata.summary, Trees < 10 & Shrubs >= 25)$Observation_ID
Grassland <-  subset(plotdata.summary, Trees < 10 & Shrubs < 10 & Field >= 25)$Observation_ID

plotdata[plotdata$Observation_ID %in% Woodland,c('Field', 'Shrub')] <- 
  plotdata[plotdata$Observation_ID %in% Woodland,c('Field', 'Shrub')]/100

plotdata[plotdata$Observation_ID %in% Shrubland,c('Field', 'Subcanopy', 'Tree')] <- 
  plotdata[plotdata$Observation_ID %in% Shrubland,c('Field', 'Subcanopy', 'Tree')]/100

plotdata[plotdata$Observation_ID %in% Grassland,c('Shrub', 'Subcanopy', 'Tree')] <- 
  plotdata[plotdata$Observation_ID %in% Grassland,c('Shrub', 'Subcanopy', 'Tree')]/100

plotdata$sqrttotal <- (100*(1-10^(apply(log10(1-(plotdata[,c('Field', 'Shrub', 'Subcanopy', 'Tree')]/100.001)), MARGIN = 1, FUN='sum'))))^0.5
}
plotdata$logtotal <- (log10(100*(1-10^(apply(log10(1-(plotdata[,c('Field', 'Shrub', 'Subcanopy', 'Tree')]/100.001)), MARGIN = 1, FUN='sum'))))+2)

plotmatrix <- makecommunitydataset(plotdata, row = 'soilplot', column = 'Species', value = 'logtotal', drop = TRUE)



if (F){
  k=8
  c = 2
  pv <- 0.01
  #flex
  d <- vegdist(plotmatrix, method='bray', binary=FALSE, na.rm=T)
  t <- flexbeta(d, beta =  -0.25)
  s <- stride(seq= 2:k,arg2= t)
  o <- optimclass(comm=plotmatrix, stride=s, pval = pv, counts = c)
  oo <- data.frame(cbind(clusters = o$clusters, flex.sp = o$sig.spc, flex.cl = o$sig.clust))
  #ward
  d <- vegdist(plotmatrix, method='bray', binary=FALSE, na.rm=T)
  t <- agnes(d, method='ward')
  s <- stride(seq= 2:k,arg2= as.hclust(t))
  o <- optimclass(comm=plotmatrix, stride=s, pval = pv, counts = c)
  oo <- cbind(oo, ward.sp = o$sig.spc, ward.cl = o$sig.clust)
  #upgma
  d <- vegdist(plotmatrix, method='bray', binary=FALSE, na.rm=T)
  t <- agnes(d, method='average')
  s <- stride(seq= 2:k,arg2= as.hclust(t))
  o <- optimclass(comm=plotmatrix, stride=s, pval = pv, counts = c)
  oo <- cbind(oo, upgma.sp = o$sig.spc, upgma.cl = o$sig.clust)
  #jacc
  d <- vegdist(plotmatrix, method='jaccard', binary=FALSE, na.rm=T)
  t <- agnes(d, method='average')
  s <- stride(seq= 2:k,arg2= as.hclust(t))
  o <- optimclass(comm=plotmatrix, stride=s, pval = pv, counts = c)
  oo <- cbind(oo, jacc.sp = o$sig.spc, jacc.cl = o$sig.clust)
  #diana
  d <- vegdist(plotmatrix, method='bray', binary=FALSE, na.rm=T)
  t <- diana(d)
  s <- stride(seq= 2:k,arg2= as.hclust(t))
  o <- optimclass(comm=plotmatrix, stride=s, pval = pv, counts = c)
  oo <- cbind(oo, diana.sp = o$sig.spc, diana.cl = o$sig.clust)
  #flexbin
  d <- vegdist(plotmatrix, method='bray', binary=TRUE, na.rm=T)
  t <- flexbeta(d, beta =  -0.25)
  s <- stride(seq= 2:k,arg2= t)
  o <- optimclass(comm=plotmatrix, stride=s, pval = pv, counts = c)
  oo <- cbind(oo, flexbin.sp = o$sig.spc, flexbin.cl = o$sig.clust)
  #wardbin
  d <- vegdist(plotmatrix, method='bray', binary=TRUE, na.rm=T)
  t <- agnes(d, method='ward')
  s <- stride(seq= 2:k,arg2= as.hclust(t))
  o <- optimclass(comm=plotmatrix, stride=s, pval = pv, counts = c)
  oo <- cbind(oo, wardbin.sp = o$sig.spc, wardbin.cl = o$sig.clust)
  #upgmabin
  d <- vegdist(plotmatrix, method='bray', binary=TRUE, na.rm=T)
  t <- agnes(d, method='average')
  s <- stride(seq= 2:k,arg2= as.hclust(t))
  o <- optimclass(comm=plotmatrix, stride=s, pval = pv, counts = c)
  oo <- cbind(oo, upgmabin.sp = o$sig.spc, upgmabin.cl = o$sig.clust)
  #dianabin
  d <- vegdist(plotmatrix, method='bray', binary=TRUE, na.rm=T)
  t <- diana(d)
  s <- stride(seq= 2:k,arg2= as.hclust(t))
  o <- optimclass(comm=plotmatrix, stride=s, pval = pv, counts = c)
  oo <- cbind(oo, dianabin.sp = o$sig.spc, dianabin.cl = o$sig.clust)
  oo <- oo[,c(1, (1:9)*2+1, (1:9)*2)]
}
if (T){
  k=8
  c = 2
  pv <- 0.01
  #ward ----
  amethod <- 'bray-ward' 
  d <- vegdist(plotmatrix, method='bray', binary=FALSE, na.rm=T)
  t <- agnes(d, method='ward')
  s <- stride(seq= 2:k,arg2= as.hclust(t))
  sil<-NA
  for(j in 2:k){
    sil0 <- (t %>% cutree(k=j) %>% silhouette(d)) 
    siltable <-  cbind(plot = rownames(plotmatrix), sil0) %>% as.data.frame()
    silneg <- subset(siltable, sil_width <0) %>% nrow()
    silmean <- siltable$sil_width %>% as.numeric() %>% mean()
    sil0 <- as.data.frame(cbind(silmean=silmean, silneg=silneg))
    if(is.na(sil)[1]){sil <- sil0}else{sil <- rbind(sil,sil0)}
  }
  o <- optimclass(comm=plotmatrix, stride=s, pval = pv, counts = c)
  oo <- data.frame(cbind(clusters = o$clusters, method = amethod, ind.sp = o$sig.spc, ind.clust = o$sig.clust, sil.mean = sil$silmean, sil.neg = sil$silneg))
  makeplot(amethod,d,t,k)
  #upgma ----
  amethod <- 'bray-upgma' 
  d <- vegdist(plotmatrix, method='bray', binary=FALSE, na.rm=T)
  t <- agnes(d, method='average')
  s <- stride(seq= 2:k,arg2= as.hclust(t))
  sil<-NA
  for(j in 2:k){
    sil0 <- (t %>% cutree(k=j) %>% silhouette(d)) 
    siltable <-  cbind(plot = rownames(plotmatrix), sil0) %>% as.data.frame()
    silneg <- subset(siltable, sil_width <0) %>% nrow()
    silmean <- siltable$sil_width %>% as.numeric() %>% mean()
    sil0 <- as.data.frame(cbind(silmean=silmean, silneg=silneg))
    if(is.na(sil)[1]){sil <- sil0}else{sil <- rbind(sil,sil0)}
  }
  o <- optimclass(comm=plotmatrix, stride=s, pval = pv, counts = c)
  oo0 <- data.frame(cbind(clusters = o$clusters, method = amethod, ind.sp = o$sig.spc, ind.clust = o$sig.clust, sil.mean = sil$silmean, sil.neg = sil$silneg))
  oo <- rbind(oo,oo0)
  makeplot(amethod,d,t,k)
  #diana ----
  amethod <- 'bray-diana' 
  d <- vegdist(plotmatrix, method='bray', binary=FALSE, na.rm=T)
  t <- diana(d)
  s <- stride(seq= 2:k,arg2= as.hclust(t))
  sil<-NA
  for(j in 2:k){
    sil0 <- (t %>% cutree(k=j) %>% silhouette(d)) 
    siltable <-  cbind(plot = rownames(plotmatrix), sil0) %>% as.data.frame()
    silneg <- subset(siltable, sil_width <0) %>% nrow()
    silmean <- siltable$sil_width %>% as.numeric() %>% mean()
    sil0 <- as.data.frame(cbind(silmean=silmean, silneg=silneg))
    if(is.na(sil)[1]){sil <- sil0}else{sil <- rbind(sil,sil0)}
  }
  o <- optimclass(comm=plotmatrix, stride=s, pval = pv, counts = c)
  oo0 <- data.frame(cbind(clusters = o$clusters, method = amethod, ind.sp = o$sig.spc, ind.clust = o$sig.clust, sil.mean = sil$silmean, sil.neg = sil$silneg))
  oo <- rbind(oo,oo0)
  makeplot(amethod,d,t,k)
  #flexp00 ----
  amethod <- 'bray-flexp00' 
  d <- vegdist(plotmatrix, method='bray', binary=FALSE, na.rm=T)
  t <- flexbeta(d, beta =  0.00)
  s <- stride(seq= 2:k,arg2= t)
  sil<-NA
  for(j in 2:k){
    sil0 <- (t %>% cutree(k=j) %>% silhouette(d)) 
    siltable <-  cbind(plot = rownames(plotmatrix), sil0) %>% as.data.frame()
    silneg <- subset(siltable, sil_width <0) %>% nrow()
    silmean <- siltable$sil_width %>% as.numeric() %>% mean()
    sil0 <- as.data.frame(cbind(silmean=silmean, silneg=silneg))
    if(is.na(sil)[1]){sil <- sil0}else{sil <- rbind(sil,sil0)}
  }
  o <- optimclass(comm=plotmatrix, stride=s, pval = pv, counts = c)
  oo0 <- data.frame(cbind(clusters = o$clusters, method = amethod, ind.sp = o$sig.spc, ind.clust = o$sig.clust, sil.mean = sil$silmean, sil.neg = sil$silneg))
  oo <- rbind(oo,oo0)  
  makeplot(amethod,d,t,k)
  #flexn10 ----
  amethod <- 'bray-flexn10' 
  d <- vegdist(plotmatrix, method='bray', binary=FALSE, na.rm=T)
  t <- flexbeta(d, beta =  -0.10)
  s <- stride(seq= 2:k,arg2= t)
  sil<-NA
  for(j in 2:k){
    sil0 <- (t %>% cutree(k=j) %>% silhouette(d)) 
    siltable <-  cbind(plot = rownames(plotmatrix), sil0) %>% as.data.frame()
    silneg <- subset(siltable, sil_width <0) %>% nrow()
    silmean <- siltable$sil_width %>% as.numeric() %>% mean()
    sil0 <- as.data.frame(cbind(silmean=silmean, silneg=silneg))
    if(is.na(sil)[1]){sil <- sil0}else{sil <- rbind(sil,sil0)}
  }
  o <- optimclass(comm=plotmatrix, stride=s, pval = pv, counts = c)
  oo0 <- data.frame(cbind(clusters = o$clusters, method = amethod, ind.sp = o$sig.spc, ind.clust = o$sig.clust, sil.mean = sil$silmean, sil.neg = sil$silneg))
  oo <- rbind(oo,oo0) 
  makeplot(amethod,d,t,k)
  #flexn25 ----
  amethod <- 'bray-flexn25' 
  d <- vegdist(plotmatrix, method='bray', binary=FALSE, na.rm=T)
  t <- flexbeta(d, beta =  -0.25)
  s <- stride(seq= 2:k,arg2= t)
  sil<-NA
  for(j in 2:k){
    sil0 <- (t %>% cutree(k=j) %>% silhouette(d)) 
    siltable <-  cbind(plot = rownames(plotmatrix), sil0) %>% as.data.frame()
    silneg <- subset(siltable, sil_width <0) %>% nrow()
    silmean <- siltable$sil_width %>% as.numeric() %>% mean()
    sil0 <- as.data.frame(cbind(silmean=silmean, silneg=silneg))
    if(is.na(sil)[1]){sil <- sil0}else{sil <- rbind(sil,sil0)}
  }
  o <- optimclass(comm=plotmatrix, stride=s, pval = pv, counts = c)
  oo0 <- data.frame(cbind(clusters = o$clusters, method = amethod, ind.sp = o$sig.spc, ind.clust = o$sig.clust, sil.mean = sil$silmean, sil.neg = sil$silneg))
  oo <- rbind(oo,oo0)
  makeplot(amethod,d,t,k)
  #flexn30 ----
  amethod <- 'bray-flexn30' 
  d <- vegdist(plotmatrix, method='bray', binary=FALSE, na.rm=T)
  t <- flexbeta(d, beta =  -0.3)
  s <- stride(seq= 2:k,arg2= t)
  sil<-NA
  for(j in 2:k){
    sil0 <- (t %>% cutree(k=j) %>% silhouette(d)) 
    siltable <-  cbind(plot = rownames(plotmatrix), sil0) %>% as.data.frame()
    silneg <- subset(siltable, sil_width <0) %>% nrow()
    silmean <- siltable$sil_width %>% as.numeric() %>% mean()
    sil0 <- as.data.frame(cbind(silmean=silmean, silneg=silneg))
    if(is.na(sil)[1]){sil <- sil0}else{sil <- rbind(sil,sil0)}
  }
  o <- optimclass(comm=plotmatrix, stride=s, pval = pv, counts = c)
  oo0 <- data.frame(cbind(clusters = o$clusters, method = amethod, ind.sp = o$sig.spc, ind.clust = o$sig.clust, sil.mean = sil$silmean, sil.neg = sil$silneg))
  oo <- rbind(oo,oo0)
  makeplot(amethod,d,t,k)
  #flexn35 ----
  amethod <- 'bray-flexn35' 
  d <- vegdist(plotmatrix, method='bray', binary=FALSE, na.rm=T)
  t <- flexbeta(d, beta =  -0.35)
  s <- stride(seq= 2:k,arg2= t)
  sil<-NA
  for(j in 2:k){
    sil0 <- (t %>% cutree(k=j) %>% silhouette(d)) 
    siltable <-  cbind(plot = rownames(plotmatrix), sil0) %>% as.data.frame()
    silneg <- subset(siltable, sil_width <0) %>% nrow()
    silmean <- siltable$sil_width %>% as.numeric() %>% mean()
    sil0 <- as.data.frame(cbind(silmean=silmean, silneg=silneg))
    if(is.na(sil)[1]){sil <- sil0}else{sil <- rbind(sil,sil0)}
  }
  o <- optimclass(comm=plotmatrix, stride=s, pval = pv, counts = c)
  oo0 <- data.frame(cbind(clusters = o$clusters, method = amethod, ind.sp = o$sig.spc, ind.clust = o$sig.clust, sil.mean = sil$silmean, sil.neg = sil$silneg))
  oo <- rbind(oo,oo0)
  makeplot(amethod,d,t,k)
  oo[,c(1,3:6)] <- apply(oo[,c(1,3:6)], MARGIN = 2, FUN = 'as.numeric')
  oo <- oo[order(oo$clusters, -oo$sil.mean),]
}
if (F){
  amethod <- 'kulczynski-isopam' 
  k=4
  d <- vegdist(plotmatrix, method='kulczynski', binary=FALSE, na.rm=T)
  pamtree <- isopam(plotmatrix, distance = 'kulczynski', stopat = c(1,7))
  t <- pamtree$dendro
  makeplot(amethod,d,t,k)
  pamtab <- isotab(pamtree, level = 3)
  pamtab <- pamtab$tab
}
if (F){
  amethod <- 'bray-agnes' 
  k=8
  d <- vegdist(plotmatrix, method='bray', binary=FALSE, na.rm=T)
  t <- agnes(d, method='average')
  makeplot(amethod,d,t,k)
}
if (F){
  amethod <- 'bray-diana' 
  k=8
  d <- vegdist(plotmatrix, method='bray', binary=FALSE, na.rm=T)
  t <- diana(d)
  makeplot(amethod,d,t,k)
}

d <- vegdist(plotmatrix, method='bray', binary=FALSE, na.rm=T)
t <- flexbeta(d, beta =  -0.25)

sil0 <- (t %>% cutree(k=8) %>% silhouette(d)) 
siltable <-  cbind(plot = rownames(plotmatrix), sil0) %>% as.data.frame()
siltable[,c(2:4)] <- apply(siltable[,c(2:4)], MARGIN = 2, FUN = 'as.numeric')
silneg <- subset(siltable, sil_width <0) %>% nrow()
silmean <- siltable$sil_width %>% as.numeric() %>% mean()
#continue ----
groups <- cutree(t, k = k)
#Get species importance by cluster
soilplot <- names(groups)
clust <- unname(groups)
groupdf <- as.data.frame(cbind(soilplot, clust))

plotgroup <- merge(groupdf, plotdata, by='soilplot')
plotgroupsum <- aggregate(plotgroup$Total, by=list(clust = plotgroup$clust, Species = plotgroup$Species), FUN='sum')
names(plotgroupsum)[names(plotgroupsum)=='x'] <-'sum'
plotgroupmax <- aggregate(plotgroupsum$sum, by=list(clust = plotgroupsum$clust), FUN='max')
names(plotgroupmax)[names(plotgroupmax)=='x'] <-'max'
plotgroupsum <- merge(plotgroupsum, plotgroupmax, by='clust')
plotgroupsum$Imp <- plotgroupsum$sum/plotgroupsum$max
rm(plotgroupmax, plotgroup)

states <- c('MI','IN','OH')
core <- c('MI')
ecoregion <- c('222')
level <- 'Association'
level <- unique(unit[unit$hierarchylevel %in% level,'element_global_id'])
states <- unique(vegstates[vegstates$subnation_code %in% states,'element_global_id'])
core <- unique(vegstates[vegstates$subnation_code %in% core,'element_global_id'])
ecoregion <- unique(vegecoregion[vegecoregion$usfs_ecoregion_2007_concat_cd %in% ecoregion,'element_global_id'])

plotassociations <- as.data.frame(lapply(as.data.frame(cbind(clust = 'x', 'element_global_id'=0, 'scientificname'='x')), as.character), stringsAsFactors=FALSE)

for(i in 1:k){
  g <- subset(plotgroupsum, clust %in% i)
  g$Imp <- g$Imp^1 #aggregate data should not be square rooted
  
  gtotal <- sum(g$Imp)
  vegtotal <- aggregate(USNVClist$x, by=list(element_global_id = USNVClist$element_global_id, scientificname = USNVClist$scientificname), FUN='sum')
  names(vegtotal)[names(vegtotal)=='x'] <-'vegtotal'
  gmerge <- merge(g, USNVClist, by.x = 'Species', by.y = 'acctaxon')
  gmerge$intersect <- (gmerge$Imp*gmerge$x)^0.5
  gintersect <- aggregate(gmerge$intersect, by=list(element_global_id = gmerge$element_global_id, scientificname = gmerge$scientificname), FUN='sum')
  names(gintersect)[names(gintersect)=='x'] <-'intersect'
  g <- merge(gintersect, vegtotal, by=c('element_global_id', 'scientificname'))
  g$affinity <- g$intersect/(g$vegtotal*gtotal)^0.5*100
  g$state <- 'no'
  if(nrow(g[g$element_global_id %in% states,])>0){
    g[g$element_global_id %in% states,]$state <- 'near'}
  if(nrow(g[g$element_global_id %in% core,])>0){
    g[g$element_global_id %in% core,]$state <- 'yes'}
  g$ecoregion <- 'no'
  if(nrow(g[g$element_global_id %in% ecoregion,])>0){
    g[g$element_global_id %in% ecoregion,]$ecoregion <- 'yes'}
  g$level <- 'no'
  g[g$element_global_id %in% level,]$level <- 'yes'
  g <- subset(g, level == 'yes')
  g$best <- g$affinity
  g[g$ecoregion == 'no',]$best <- g[g$ecoregion == 'no',]$best*0.75
  g[g$element_global_id %in% states,]$best <- g[g$element_global_id %in% states,]$best*1/0.75
  g[g$element_global_id %in% core,]$best <- g[g$element_global_id %in% core,]$best*1/0.75
  g$best <- g$best/max(g$best)*100
  g <- g[,!colnames(g)%in% c('intersect','vegtotal')]
  rm(gmerge, vegtotal, gintersect)
  g <- subset(g, best >= 0 & level == 'yes')# & state == 'yes')
  g <- g[order(g$best, decreasing = TRUE),]
  plotassociations1 <- as.data.frame(lapply(as.data.frame(cbind(clust = i,g[1,1:2])), as.character), stringsAsFactors=FALSE)
  
  plotassociations <- rbind(plotassociations,plotassociations1)
}
rm(plotassociations1)
plotassociations <- plotassociations[-1,]
dmatrix <-  as.data.frame(as.matrix(d))
write.csv(plotassociations, 'output/allclustassociations.csv', row.names = F)