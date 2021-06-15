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
plotdata$over <- 100*(1-10^(apply(log10(1-(plotdata[,c('Subcanopy', 'Tree')]/100.001)), MARGIN = 1, FUN='sum')))
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

plotmatrix <- makecommunitydataset(plotdata, row = 'soilplot', column = 'Species', value = 'sqrttotal', drop = TRUE)



if (T){
  amethod <- 'bray-ward' 
  k=9
  d <- vegdist(plotmatrix, method='bray', binary=FALSE, na.rm=T)
  t <- agnes(d, method='ward')
  makeplot(amethod,d,t,k)
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
write.csv(plotassociations, 'output/allclustassociations.csv', row.names = F)