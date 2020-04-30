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
  
  #make cuts and reformat dendrogram
  ngroups=k
  groups <- cutree(t, k = ngroups)
  
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
  
  newlabels <- t$order.lab
  newlabels <- as.data.frame(newlabels)
  newlabels$row <- row(newlabels)
  newlabels <- merge(newlabels, groupdf, by.x='newlabels', by.y ='soilplot')
  newlabels$newlabels <- paste(newlabels$clust, newlabels$newlabels)
  newlabels <- newlabels[order(newlabels$row),1]
  newtree <- t
  newtree$order.lab <- newlabels
  
  dend1 <- color_branches(as.hclust(newtree), k = ngroups)
  dend1 <- color_labels(dend1, k = ngroups)
  
  #output file
  
  w <- 800
  h <- nrow(d)*12+80
  u <- 12
  png(filename=filename,width = w, height = h, units = "px", pointsize = u)
  
  par(mar = c(2,0,1,13))
  plot(dend1, horiz = TRUE, main=paste('floristic simularity', amethod,' method of plots'), font=1, cex=0.84)
  dev.off()
  
}

states <-read.delim('data/d_subnation.txt')
vegstates <-read.delim('data/UnitXSubnation.txt')
vegstates <- merge(vegstates, states, by='subnation_id')

USNVClist <- readRDS('data/USNVClist.RDS')
USNVClist <- subset(USNVClist, !grepl('\\.', acctaxon))
plotdata <- readRDS('data/plotdata.RDS')




plotmatrix <- makecommunitydataset(plotdata, row = 'soilplot', column = 'Species', value = 'sqrttotal', drop = TRUE)



if (T){
  amethod <- 'bray-ward' 
  k=8
  d <- vegdist(plotmatrix, method='bray', binary=FALSE, na.rm=T)
  t <- agnes(d, method='ward')
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
cluster <- '6'

states <- c('MI', 'IN')
states <- unique(vegstates[vegstates$subnation_code %in% states,'element_global_id'])
g <- subset(plotgroupsum, clust==cluster)

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
g[g$element_global_id %in% states,]$state <- 'yes'
g$best <- g$affinity/max(g$affinity)*100
g <- g[,!colnames(g)%in% c('intersect','vegtotal')]
rm(gmerge, vegtotal, gintersect)
g <- subset(g, best >= 50)
g <- g[order(g$best, decreasing = TRUE),]

