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
# usnvcspp2 <- usnvcspp2 |> group_by(taxon, ELEMENT_GLOBAL_ID) |> summarise(indicator=sum(indicator))
# usnvcspp2 <- ass[,c('classificationCode','databaseCode','PARENT_ID','ELEMENT_GLOBAL_ID','hierarchyLevel','colloquialName', 'scientificName')] |> left_join(usnvcspp2) |> arrange(classificationCode, scientificName, -indicator, taxon) |> subset(hierarchyLevel %in% 'Association')
write.csv(usnvcspp2, 'usnvcspp2.csv', row.names = F, na='')


#Extract any bryophytes or lichens ----

usda <- read.csv('data/plantlst.txt')
usda <- usda |> mutate(taxon = extractTaxon(Scientific.Name.with.Author), auth=extractTaxon(Scientific.Name.with.Author, 'author'), acc=Symbol, syn=ifelse(Synonym.Symbol %in% '',Symbol,Synonym.Symbol))
usda.backbone <- subset(usda, acc==syn, select=c("Common.Name","Family", "taxon","auth","acc")) |> mutate(acgenus = extractTaxon(harmonize.taxa(taxon, fix = TRUE), 'genus'))
colnames(usda.backbone)
apg <- vegnasis::apg
familylink <- vegnasis::familylink
apglink <- apg |> left_join(familylink)
usda.backbone <- usda.backbone |> merge(apglink, by.x = 'acgenus', by.y = 'genus', all.x = TRUE)
usdanva <- subset(usda.backbone, is.na(phylum))

bry <- read.delim('data/0000116-260623161305970.csv')
lic <- read.delim('data/0000193-260623161305970.csv')
alg <- read.delim('data/0000322-260623161305970.csv')
lic2 <- read.delim('data/0000964-260623161305970.csv')
alg2 <- read.delim('data/0000982-260623161305970.csv')

nva <-rbind(bry,lic,alg,lic2,alg2)
nva <- nva |> mutate(actaxon = extractTaxon(acceptedScientificName),acauth = extractTaxon(acceptedScientificName, 'author'), 
                     acgenus = extractTaxon(actaxon, 'genus'), taxon = extractTaxon(scientificName),auth = extractTaxon(scientificName, 'author'),
                     genus = extractTaxon(taxon, 'genus'))



nvataxonomy <- subset(nva, (taxonRank %in% 'GENUS' | (numberOfOccurrences >= 100 & taxonRank %in% 'SPECIES'))  & taxonomicStatus %in% c('PROVISIONALLY_ACCEPTED','ACCEPTED') & (kingdom %in% c("Bacteria","Fungi") | class %in% c("Phaeophyceae","Xanthophyceae") | phylum %in% c("Bryophyta","Marchantiophyta","Anthocerotophyta","Rhodophyta","Charophyta","Chlorophyta")), select=c(kingdom, phylum, class, order, family, genus)) |> unique() |> arrange(kingdom, phylum, class, order, family, genus) |> group_by(genus) |> mutate(n = length(genus)) |> ungroup()
correct = data.frame(
  genus = c('Volvox','Verdigellas','Variolaria','Valdonia','Urospora','Trismegistia','Timmiella','Thrombium','Thorea','Tetmemorus','Tephromela','Teilingia','Syringoderma','Symphyodon','Striaria','Stigonema','Stictyosiphon','Stichococcus','Stereocaulon','Staurodesmus','Staurastrum','Sporastatia','Splachnobryum', 'Speerschneidera', 'Sematophyllum','Scouleria','Sclerococcum','Schizochlamys','Schimmelmannia','Schaereria','Saelania', 'Ropalospora', 'Rivularia','Rhadinocladia','Pyrenastrum', 'Pycnora','Pterospermella','Pterogonidium', 'Pterobryella', 'Pseudolithophyllum', 'Pseudolithoderma','Pseudohypnella','Protoderma','Porina','Polycoccum','Pleurotaenium','Pleopsidium','Plectonema','Platygramme','Pilotrichella','Phymatoceros','Phaeostrophion','Phaeographina','Pellia','Pallavicinia','Orthodontium','Opegrapha','Nowellia','Nothoceros','Neohodgsonia','Myriotrichia','Mycomicrothelia','Mycobilimbia', 'Muellerella','Microzonia','Microthamnion','Microphiale','Micromitrium','Microcoleus','Micrasterias','Mesochaete','Mastigocoleus','Malcolmiella', 'Malbranchea','Luisierella','Lopadium', 'Lobothallia', 'Lithoderma', 'Limbella', 'Lichenostigma', 'Leucomium', 'Leptotheca', 'Leptopterigynandrum', 'Leptobryum','Lecidea', 'Leathesia','Klebsormidium', 'Kirchneriella','Jamesoniella','Hypopterygium', 'Hypodontium', 'Hypocenomyce','Hypnella','Hymenodon','Hydropogonella','Hydrocoleum','Hyalotheca','Hummia','Hookeriopsis','Holodontium','Heterocladium','Helminthocarpon','Haplogloia','Halorhipis','Halecania','Grateloupia','Graphina','Gracilariopsis','Glyphomitrium','Glyphium','Gloeocystis','Geminella','Garovaglia','Fimbriaria','Euptychium','Euastrum','Epibryon','Drummondia','Dolichospermum','Distichium','Dirinaria','Diphyscium','Diorygma','Dimerella','Dictyosiphon','Delamarea','Dactylospora','Cyrtopus','Cyathophorum','Crucigenia','Crocynia','Cosmocladium','Cosmarium','Compsonema','Coenogonium','Circinaria','Chrysoblastella','Chroodactylon','Chondria','Characium','Chaetosphaeridium','Chaetomitrium','Ceramium','Carbacanthographis','Capsosiphon','Calyptrozyma','Calothrix','Brachytrichia','Bissetia','Bilimbia','Baculifera','Audouinella','Asterella','Aspicilia', 'Aspergillus','Arthopyrenia', 'Ankyra','Aneura','Analipus','Phaeoceros','Hymenoloma','Halochlorococcum','Tolypothrix','Punctaria','Micrasterias','Leathesia','Hypnea','Cornicularia','Coilodesme','Loxospora', 'Aspicilia','Soranthera','Schizothrix','Amandinea','Dirinaria','Opegrapha','Porina','Rhizofabronia','Tephromela','Ulvella','Wilsoniella'), 
  
  family = c('Volvocaceae','Palmophyllaceae','Pertusariaceae','Seligeriaceae', 'Ulotrichaceae','Pylaisiadelphaceae','Timmiellaceae','Protothelenellaceae','Thoreaceae','Desmidiaceae','Tephromelataceae','Desmidiaceae','Syringodermataceae','Symphyodontaceae','Chordariaceae','Stigonemataceae','Chordariaceae','Prasiolaceae','Stereocaulaceae','Desmidiaceae','Desmidiaceae','Sporastatiaceae','Splachnobryaceae','Leprocaulaceae','Sematophyllaceae', 'Scouleriaceae',	'Sclerococcaceae','Schizochlamydaceae','Schimmelmanniaceae','Schaereriaceae','Saelaniaceae', 'Ropalosporaceae', 'Nostocaceae','Chordariaceae','Pyrenulaceae','Pycnoraceae','Pterospermataceae','Pylaisiadelphaceae','Pterobryellaceae', 'Corallinaceae', 'Lithodermataceae', 'Hypnaceae','Chaetophoraceae','Porinaceae','Polycoccaceae','Desmidiaceae','Acarosporaceae','Oscillatoriaceae','Graphidaceae','Lembophyllaceae','Phymatocerotaceae','Phaeostrophionaceae','Graphidaceae','Pelliaceae','Pallaviciniaceae','Orthodontiaceae','Opegraphaceae','Cephaloziaceae','Dendrocerotaceae','Neohodgsoniaceae','Chordariaceae','Trypetheliaceae','Ramalinaceae','Verrucariaceae','Syringodermataceae','Microthamniaceae','Coenogoniaceae','Micromitriaceae','Microcoleaceae','Desmidiaceae','Aulacomniaceae','Nostochopsidaceae','Pilocarpaceae','Onygenaceae','Timmiellaceae','Lopadiaceae', 'Megasporaceae', 'Lithodermataceae', 'Amblystegiaceae','Phaeococcomycetaceae', 'Leucomiaceae', 'Orthodontiaceae', 'Thuidiaceae', 'Bryaceae', 'Lecideaceae', 'Chordariaceae','Klebsormidiaceae', 'Selenastraceae','Adelanthaceae', 'Hypopterygiaceae', 'Hypodontiaceae','Ophioparmaceae','Pilotrichaceae','Orthodontiaceae','Sematophyllaceae','Microcoleaceae','Desmidiaceae','Chordariaceae','Pilotrichaceae','Hypodontiaceae','Lembophyllaceae','Graphidaceae','Chordariaceae','Chordariaceae','Leprocaulaceae','Halymeniaceae','Graphidaceae','Gracilariaceae','Rhabdoweisiaceae','Mytilinidiaceae','Radiococcaceae','Chlorellaceae','Ptychomniaceae','Rhodomelaceae','Ptychomniaceae','Desmidiaceae','Epibryaceae','Drummondiaceae','Aphanizomenonaceae','Distichiaceae','Caliciaceae','Diphysciaceae','Graphidaceae','Coenogoniaceae','Chordariaceae','Chordariaceae','Dactylosporaceae','Cyrtopoaceae','Hypopterygiaceae','Scenedesmaceae','Ramalinaceae','Desmidiaceae','Desmidiaceae','Scytosiphonaceae','Coenogoniaceae','Megasporaceae','Chrysoblastellaceae','Stylonemataceae','Rhodomelaceae','Characiaceae','Chaetosphaeridiaceae','Symphyodontaceae','Ceramiaceae','Graphidaceae','Ulotrichaceae','Helotiaceae','Rivulariaceae','Scytonemataceae','Neckeraceae','Ramalinaceae','Caliciaceae','Acrochaetiaceae', 'Aytoniaceae','Megasporaceae', 'Aspergillaceae', 'Trypetheliaceae', 'Sphaeropleaceae','Aneuraceae','Scytosiphonaceae','Notothyladaceae','Hymenolomataceae','Oltmannsiellopsidaceae','Tolypothrichaceae','Chordariaceae','Desmidiaceae','Chordariaceae','Cystocloniaceae','Parmeliaceae','Chordariaceae','Sarrameanaceae', 'Megasporaceae','Chordariaceae','Trichocoleusaceae','Caliciaceae','Caliciaceae','Opegraphaceae','Porinaceae','Rhizofabroniaceae','Tephromelataceae','Ulvellaceae','Ditrichaceae'), correct=1)
#	Analipus -> Ralfsiaceae; Mastigocoleus -> Hapalosiphonaceae; Rivularia ->	Rivulariaceae

nvataxonomy <- left_join(nvataxonomy, correct) |> unique()
nvataxonomy <- nvataxonomy |> subset(n %in% 1 | correct %in% 1, select=-c(correct,n))
nvataxonomy <- nvataxonomy |> mutate(family = case_when(genus %in% 'Analipus' ~ 'Ralfsiaceae',
                                                        genus %in% 'Mastigocoleus' ~ 'Hapalosiphonaceae',
                                                        genus %in% 'Rivularia' ~ 'Rivulariaceae',
                                                        TRUE ~ family
))
nvataxonomy <- nvataxonomy |> mutate(order = case_when(family %in% 'Ralfsiaceae' ~ 'Ralfsiales',
                                                        TRUE ~ order
))


usdalich <- subset(usda, grepl(tolower('lichen'), Common.Name)) |> mutate(genus = extractTaxon(taxon, 'genus'))
nvataxonomy <- nvataxonomy |> mutate(usdalichen = ifelse(genus %in% usdalich$genus,1,0), usdalichen2 = ifelse(family %in% usdalich$Family,1,0))
nvataxonomy <- nvataxonomy |> mutate(type =  case_when(phylum %in% 'Cyanobacteria' ~ 'Cyanobacteria',
                                                       class %in% 'Phaeophyceae' ~ 'brown algae',
                                                       class %in% 'Xanthophyceae' ~ 'yellow-green algae',
                                                       order %in% 'Arthoniales' ~ 'lichen',
                                                       order %in% 'Vezdaeales' ~ 'lichen',
                                                       order %in% 'Monoblastiales' ~ 'lichen',
                                                       order %in% 'Strigulales' ~ 'lichen',
                                                       #class %in% 'Eurotiomycetes' ~ 'lichen',
                                                       class %in% 'Lecanoromycetes' ~ 'lichen',
                                                       class %in% 'Lichinomycetes' ~ 'lichen',
                                                       class %in% 'Candelariomycetes' ~ 'lichen',
                                                       class %in% 'Coniocybomycetes' ~ 'lichen',
                                                       family %in% 'Trypetheliaceae' ~ 'lichen',
                                                       family %in% 'Pyrenulaceae' ~ 'lichen',
                                                       family %in% 'Verrucariaceae' ~ 'lichen',
                                                       phylum %in% 'Anthocerotophyta' ~ 'bryophyte',
                                                       phylum %in% 'Bryophyta' ~ 'bryophyte',
                                                       phylum %in% 'Chlorophyta' ~ 'green algae',
                                                       phylum %in% 'Marchantiophyta' ~ 'bryophyte',
                                                       phylum %in% 'Chlorophyta' ~ 'green algae',
                                                       phylum %in% 'Rhodophyta' ~ 'red algae',
                                                       phylum %in% 'Charophyta' ~ 'green algae',
                                                       usdalichen %in% 1 | usdalichen2 %in% 1 ~ '2lich',
                                                       TRUE ~ NA)) |> subset(!is.na(type) & !type %in% '2lich', select=-c(usdalichen,usdalichen2))

nvaselect <- nva #|> subset((taxonomicStatus %in% c('PROVISIONALLY_ACCEPTED','ACCEPTED') | numberOfOccurrences >= 500 | taxon %in% usda$taxon) & acgenus %in% nvataxonomy$genus & taxonRank %in% c('SPECIES', 'VARIETY', 'SUBSPECIES'), select=c("acgenus","actaxon","acauth","taxon","auth" ))

write.csv(nvataxonomy, 'data/nvataxonomy.csv', row.names = F, na='')
write.csv(nvaselect, 'data/nva.csv', row.names = F, na='')

#######Talley nonvascular species by USNVC
library(vegnasis)
library(xlsx)
#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#species to choose from
nva <- read.csv('data/nva.csv')
nvataxonomy <- read.csv('data/nvataxonomy.csv')
sppacc <- subset(nva, !grepl('\\.',taxon)) 
sppacc <- unique(c(sppacc$actaxon, sppacc$taxon, nvataxonomy$genus))
#USNVC association names
ass <- read.csv('data/USNVC v3.0.3 2026-03-25/unit.csv')
#USNVC species mentions
desc <- read.csv('data/USNVC v3.0.3 2026-03-25/unitDescription.csv')

timA <- Sys.time()
spbyass <- NULL
spbydesc <- NULL
for(i in 1:length(sppacc)){#i=1
  thistax <- sppacc[i]
  #thistax <- 'Nostoc'
  
  thisass <- ass[grepl(thistax, ass$scientificName),]$ELEMENT_GLOBAL_ID
  if(length(thisass)>0){
    spbyass0 <- data.frame(taxon=paste0(thistax,'\\W|',thistax,'$'), ELEMENT_GLOBAL_ID=thisass)
    if(is.null(spbyass)){spbyass <- spbyass0}else{spbyass <- rbind(spbyass,spbyass0)}}
  
  thisdesc <- desc[grepl(thistax, desc$Floristics),]$ELEMENT_GLOBAL_ID
  if(length(thisdesc)>0){
    spbydesc0 <- data.frame(taxon=paste0(thistax,'\\W|',thistax,'$'), ELEMENT_GLOBAL_ID=thisdesc)
    if(is.null(spbydesc)){spbydesc <- spbydesc0}else{spbydesc <- rbind(spbydesc,spbydesc0)}}
}
Sys.time() - timA#Time difference of 1.324091 hours
spbyass$indicator <- 1
spbydesc$indicator <- 0

usnvcnva <- rbind(spbyass,spbydesc) |> unique() 
# colnames(usnvcnva) <- c("taxon","ELEMENT_GLOBAL_ID","indicator" )
write.csv(usnvcnva, 'usnvcnva.csv', row.names = F, na='')



#sweep up vascular genera ----

spp <- vegnasis::syns3
sppacc <- c(spp$bonap, spp$kew, spp$wplants, spp$usda) |> unique()
sppacc <- extractTaxon(sppacc, 'genus') |> unique()
#USNVC association names
ass <- read.csv('data/USNVC v3.0.3 2026-03-25/unit.csv')
#USNVC species mentions
desc <- read.csv('data/USNVC v3.0.3 2026-03-25/unitDescription.csv')

timA <- Sys.time()
spbyass <- NULL
spbydesc <- NULL
for(i in 1:length(sppacc)){#i=1
  thistax <- sppacc[i]
  #thistax <- 'Pinus strobus'
  
  thisass <- ass[grepl(paste0(thistax,'\\W|',thistax,'$'), ass$scientificName),]$ELEMENT_GLOBAL_ID
  if(length(thisass)>0){
    spbyass0 <- data.frame(taxon=thistax, ELEMENT_GLOBAL_ID=thisass)
    if(is.null(spbyass)){spbyass <- spbyass0}else{spbyass <- rbind(spbyass,spbyass0)}}
  
  thisdesc <- desc[grepl(paste0(thistax,'\\W|',thistax,'$'), desc$Floristics),]$ELEMENT_GLOBAL_ID
  if(length(thisdesc)>0){
    spbydesc0 <- data.frame(taxon=thistax, ELEMENT_GLOBAL_ID=thisdesc)
    if(is.null(spbydesc)){spbydesc <- spbydesc0}else{spbydesc <- rbind(spbydesc,spbydesc0)}}
}
Sys.time() - timA#Time difference of 11.32324 mins
spbyass$indicator <- 1
spbydesc$indicator <- 0

usnvcspp <- rbind(spbyass,spbydesc) |> unique() 
write.csv(usnvcspp, 'usnvcgen.csv', row.names = F, na='')




#Combine lists and clean them up ----
library(vegnasis)
library(xlsx)
#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#vasculars
spp <- vegnasis::syns3
familylink <- vegnasis::familylink
apg <- vegnasis::apg
#nonvasculars
nva <- read.csv('data/nva.csv')
nvataxonomy <- read.csv('data/nvataxonomy.csv')
#USNVC association names
ass <- read.csv('data/USNVC v3.0.3 2026-03-25/unit.csv')
#USNVC species mentions
desc <- read.csv('data/USNVC v3.0.3 2026-03-25/unitDescription.csv')

#fix misapplied genera in respective vascular and non-vascular lists
familylink <- subset(familylink, !genus %in% c('Botrydium','Hedwigia','Hookeria','Ephemerum','Carteria'))
nvataxonomy0 <-  subset(nvataxonomy, genus %in% familylink$genus)
spp0  <-  spp |> mutate(genus=extractTaxon(kew, 'genus')) |> subset(genus %in% nvataxonomy$genus)
nvataxonomy <- subset(nvataxonomy, !genus %in% familylink$genus)

usnvcspp <- read.csv('usnvcspp2.csv')
# usnvcspp <- usnvcspp |> mutate(taxon=harmonize.taxa(taxon, fix=TRUE, sensu = 'kew'), taxon=extractTaxon(taxon, 'binomial'))
usnvcgen <- read.csv('usnvcgen.csv')
usnvcnva <- read.csv('usnvcnva.csv')

usnvcspp <- rbind(usnvcspp[,colnames(usnvcnva)], usnvcnva, usnvcgen) |> unique()
usnvcspp <- usnvcspp |> mutate(genus = extractTaxon(taxon, 'genus'))
usnvcspp <- usnvcspp |> group_by(genus, ELEMENT_GLOBAL_ID, indicator) |> mutate(n=length(taxon)) |> ungroup()
usnvcspp <- usnvcspp |> mutate(flag=ifelse((genus != taxon | n == 1) & nchar(taxon) >= 1,0,1))
usnvcspp <- usnvcspp |> subset(!flag %in% 1, select=-c(genus, n, flag))




usnvcspp3 <- ass[,c('classificationCode','databaseCode','PARENT_ID','ELEMENT_GLOBAL_ID','hierarchyLevel','colloquialName', 'scientificName')] |> left_join(usnvcspp) |> arrange(classificationCode, scientificName, -indicator, taxon) |> subset(hierarchyLevel %in% 'Association')
