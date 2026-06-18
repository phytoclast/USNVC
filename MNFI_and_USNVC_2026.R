library(vegnasis)
library(xlsx)
#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

mnfi <- read.xlsx("data/mnficommunities.xlsx",1)
mnfi <- mnfi |> mutate(community = gsub('\\W+',' ',community),
                       community = trimws(community))

gsub('\\u00A0','@','Bog ')
common <- subset(mnfi, !is.na(plants))
rare <- subset(mnfi, !is.na(rare))

common <- common |> mutate(taxa = stringr::str_split_i(plants, '\\(|\\)', 2),
                           multi = ifelse(grepl(',|\\Wand\\W',taxa), T,F),
                           adund = 'common')
multi0 <- subset(common, multi) |> mutate(genus = extractTaxon(taxa,report = "genus"),
                                          pre = gsub('\\Wand\\W|\\W', '@', taxa),
                                          pre = gsub('@+', '@', pre))
#non-breaking white space '\\u00A0'
for(i in 1:8){#i=3
  k=i*2

  multi00 <- multi0 |> mutate(pre = str_split_i(pre, '@', k),
                              taxa = paste(genus,pre))
multi00 <- subset(multi00, !is.na(pre) & !pre %in% 'others')
multi00 <- multi00[,colnames(common)]
if(i==1){multi <- multi00}else{multi <- rbind(multi, multi00)}
}
nonmulti <- subset(common, !multi & !is.na(taxa))

common <- rbind(multi, nonmulti)
rare <- rare |> mutate(plants = rare,
                       rare = NULL,
                       taxa = stringr::str_split_i(plants, '\\(|\\)', 1),
                       adund = 'rare',
                       )
comcol <- base::intersect(colnames(common), colnames(rare))
mnfi <- rbind(common[,comcol], rare[,comcol])


####Talley species by MLRA

library(sf)
library(terra)
library(vegnasis)
library(climatools)
fips <- read_sf("C:/GIS/Base/americas3.shp")
mlra <- read_sf("C:/GIS/Ecoregion/MLRA_2018.shp")


thisproj <- climatools::setProjection(prj = "equalarea.azimuthal", lat=45, lon=-100)
thesestates <- fips |> subset(STATECODE %in% c('US.AK','US.ME','US.PR','US.HI')) |> sf::st_transform(thisproj)

r <- rast(resolution=c(5000,5000), crs=thisproj, extent=ext(thesestates))
mlra <- mlra |> st_transform(thisproj)
mlra.r <- rasterize(mlra, y=r, field='MLRA')
plot(mlra.r)
fips <- fips |> st_transform(thisproj)
fips.r <- rasterize(fips, y=r, field='DATAOUT')
plot(fips.r)
mlrafips <- c(mlra.r,fips.r)
mlrafips <- terra::crosstab(mlrafips, long=TRUE)
bfips <- readRDS("C:/scripts/bonap2020/data/final/b.fips.RDS")
bfips <- bfips |> subset(Status %in% 'N')
bmlra <- bfips |> left_join(mlrafips, by=join_by(cty.data==DATAOUT))
bmlra <- subset(bmlra, !is.na(ac.binomial) & !is.na(MLRA), select=c(ac.binomial, MLRA, cty.data, n)) |> unique()
bmlra <- bmlra |> group_by(ac.binomial, MLRA) |> summarise(n = sum(n))
bmlra <- bmlra |> group_by(MLRA) |> mutate(mlratotal = max(n)) |> ungroup() |> group_by(ac.binomial) |> mutate(spptotal = sum(n)) |> ungroup()
bmlra <- bmlra |> mutate(pmlra = 100*n/mlratotal, pspp = 100*n/spptotal, rv = ifelse(pspp >= 15 | pmlra >= 20, 1,0)) 
missing <- subset(bmlra, rv %in% 1)
missing <- subset(bmlra, !ac.binomial %in% missing$ac.binomial)



####Talley species by USNVC
#species to choose from
spp <- vegnasis::syns3
sppacc <- c(spp$bonap, spp$kew, spp$wplants, spp$usda) |> unique()
#USNVC association names
ass <- read.csv('data/USNVC v3.0.3 2026-03-25/unit.csv')
#USNVC species mentions
desc <- read.csv('data/USNVC v3.0.3 2026-03-25/unitDescription.csv')

timA <- Sys.time()
spbyass <- NULL
spbydesc <- NULL
for(i in 1:length(sppacc)){
thistax <- sppacc[i]
#thistax <- 'Pinus strobus'

thisass <- ass[grepl(thistax, ass$scientificName),]$ELEMENT_GLOBAL_ID
if(length(thisass)>0){
  spbyass0 <- data.frame(taxon=thistax, ELEMENT_GLOBAL_ID=thisass)
  if(is.null(spbyass)){spbyass <- spbyass0}else{spbyass <- rbind(spbyass,spbyass0)}}

thisdesc <- desc[grepl(thistax, desc$Floristics),]$ELEMENT_GLOBAL_ID
if(length(thisdesc)>0){
  spbydesc0 <- data.frame(taxon=thistax, ELEMENT_GLOBAL_ID=thisdesc)
  if(is.null(spbydesc)){spbydesc <- spbydesc0}else{spbydesc <- rbind(spbydesc,spbydesc0)}}
}
Sys.time() - timA#Time difference of 2.221369 hours
spbyass$indicator <- 1
spbydesc$indicator <- 0

usnvcspp <- rbind(spbyass,spbydesc) |> unique() 
# colnames(usnvcspp) <- c("taxon","ELEMENT_GLOBAL_ID","indicator" )
write.csv(usnvcspp, 'usnvcspp.csv', row.names = F, na='')

#get species hidden in parentheses
usnvcspp <- read.csv('usnvcspp.csv')
asspar <- subset(ass, grepl('[A-Z]{1}[a-z]+\\s\\([a-z]', scientificName), select=c(ELEMENT_GLOBAL_ID, scientificName))


asspar <- asspar |> mutate(genus = stringr::str_extract(scientificName, '[A-Z]{1}[a-z]+\\s\\('), 
                           genus = stringr::str_remove(genus,'\\s\\('),
                           epithets = stringr::str_extract(scientificName, '[A-Z]{1}[a-z]+\\s\\(.+\\)'),
                           epithets = stringr::str_remove(epithets, '[A-Z]{1}[a-z]+\\s\\('),
                           epithets = stringr::str_split_i(epithets, '\\)', 1))
allasspar <- NULL
for(i in 1:10){
asspar0 <- asspar |> mutate(epithet = stringr::str_split_i(epithets, ',\\s*', i))
if(is.null(allasspar)){allasspar = asspar0}else{allasspar = rbind(allasspar, asspar0)}
}
allasspar <- allasspar |> subset(!is.na(epithet)) |> mutate(taxon = paste(genus, epithet), indicator = 1) 
usnvcspp2 <- rbind(usnvcspp,allasspar[,c("taxon", "ELEMENT_GLOBAL_ID","indicator")]) |> unique()  
usnvcspp2 <- usnvcspp2 |> group_by(taxon, ELEMENT_GLOBAL_ID) |> summarise(indicator=sum(indicator))


usnvcspp2 <- ass[,c('classificationCode','databaseCode','PARENT_ID','ELEMENT_GLOBAL_ID','hierarchyLevel','colloquialName', 'scientificName')] |> left_join(usnvcspp2) |> arrange(classificationCode, scientificName, -indicator, taxon) |> subset(hierarchyLevel %in% 'Association')
write.csv(usnvcspp2, 'usnvcspp2.csv', row.names = F, na='')
