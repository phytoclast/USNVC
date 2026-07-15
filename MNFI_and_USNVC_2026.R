library(vegnasis)
library(xlsx)

library(WorldFlora)

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

mnfi <- mnfi |> mutate(taxa = trimws(gsub('\\u00A0|\\xA0',' ', taxa)))
x <- rbind(
  c('Oxalis acetosella','Oxalis acetosella ssp. montana'),
  c('Phragmites australis subsp. americanus','Phragmites australis ssp. americanus'),
  c('Potamogeton pectinata','Potamogeton pectinatus'),
  c('Arnica cordiformis','Arnica cordifolia'),
  c('Calystegia spithamea','Calystegia spithamaea'),
  c('Drepanocladus polygamus','Campylium polygamum'),
  c('Drepanocladus trifarius','Calliergon trifarium'),
  c('Eriophorum viridi-carinatum','Eriophorum viridicarinatum'),
  c('Pterospora andromeda','Pterospora andromedea'),
  c('Lycopodiella margueriteae','Lycopodiella margueritiae'),
  c('Symphyotrichum novae','Symphyotrichum novae-angliae'),
  c('Potamogeton spp','Potamogeton spp.'))
mnfi <- mnfi |> mutate(taxa = gsub('\\ssp\\.|\\sspp\\.','', taxa))
for(i in 1:nrow(x)){
mnfi <- mnfi |> mutate(taxa = gsub(paste0('^',x[i,1],'$'),x[i,2], taxa))
}
saveRDS(mnfi, 'mnfi.RDS')

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
saveRDS(bmlra, 'bmlra.RDS')


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
library(vegnasis)
library(xlsx)

#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
# nvae <- subset(nva, grepl('ë',taxon))

# nva.all <- subset(nva, taxonomicStatus %in% c('SYNONYM','ACCEPTED') & taxonRank %in% c('SPECIES', 'VARIETY', 'SUBSPECIES'), select = c(actaxon, acauth, taxon, auth))

nva.all <- nva  |> subset(taxonRank %in% c('SPECIES', 'VARIETY', 'SUBSPECIES'), select=c(numberOfOccurrences,acceptedScientificName, scientificName))

nva.all2 <- data.frame(numberOfOccurrences = nva.all$numberOfOccurrences, acceptedScientificName=nva.all$acceptedScientificName,  scientificName=nva.all$acceptedScientificName)

nva.all <- rbind(nva.all, nva.all2) |> unique() |> mutate(actaxon = extractTaxon(acceptedScientificName),acauth = extractTaxon(acceptedScientificName, 'author'), genus = extractTaxon(actaxon, 'genus'), taxon = extractTaxon(scientificName),auth = extractTaxon(scientificName, 'author')) |> mutate(taxon=gsub('ë','e',taxon), actaxon=gsub('ë','e',actaxon))

nva.all <- nva.all |> mutate(isusda = ifelse(scientificName %in% usda$Scientific.Name.with.Author,1,0), ac = ifelse(taxon == actaxon,1,0)) |> group_by(taxon, ac) |> mutate(nhom = length(taxon), maxocc = max(numberOfOccurrences)) |> group_by(taxon) |> mutate(maxac = max(ac), maxusda = max(isusda)) |> ungroup() |>
  mutate(keep = case_when(ac == 0 & numberOfOccurrences < 5 ~ 0,
                          maxusda == isusda & maxusda == 1 & maxac == 0 ~ 1,
                          maxocc == numberOfOccurrences & maxac == ac & (maxusda == 0 | maxac==1) ~ 1,
                          TRUE ~ 0)) |> subset(keep == 1, select = c(numberOfOccurrences, actaxon, acauth, taxon, auth))

nva.all <- nva.all |> mutate(nauth =  nchar(auth), sameauth = ifelse(taxon==actaxon & auth==acauth,1,0)) |> group_by(taxon) |> mutate(n=length(taxon)) |> group_by(taxon, actaxon) |> mutate(n2=length(taxon), maxsameauth = max(sameauth)) |> group_by(taxon) |> mutate(maxn2 = max(n2)) |> group_by(taxon, actaxon)  |> mutate(maxnchar = max(nauth)) |> ungroup() |> mutate(keep = (maxn2==n2 & maxnchar==nauth & maxsameauth == 0) | (maxsameauth ==1 & sameauth==1))

nva.all <- nva.all |> subset(keep, select = c(numberOfOccurrences, actaxon,acauth,taxon,auth)) |> group_by(taxon) |> mutate(n=length(taxon)) |> ungroup()

nva.all <- nva.all |> mutate(keep=ifelse(n>1,0,1))
nvalist <- unique(nva.all[nva.all$n > 1,]$taxon)
for (i in 1:length(nvalist)){#i=3
  thistaxon <-  nvalist[i]#'Aphanolejeunea diaphana'
  usdadist <- usda |> subset(taxon %in% thistaxon)
  nvanarrow <- nva.all |> subset(taxon %in% thistaxon)
  if(nrow(usdadist) > 0){
    nvanarrow <- nvanarrow |> mutate(d = stringdist::stringdist(auth, usdadist$auth), mind = min(d))
    nvanarrow <- nvanarrow |> subset(mind==d)
    nva.all <- nva.all |> mutate(keep = case_when(taxon %in% thistaxon & auth %in% nvanarrow$auth ~ -1,
                                                  TRUE ~ keep))
  }
}
nva.all <- nva.all |> group_by(taxon) |> mutate(rnk = rank(paste(keep,auth)))
nva.all <- nva.all |> subset(rnk == 1,  select = c(actaxon,acauth,taxon,auth))
#URLs: https://www.inaturalist.org/taxa/inaturalist-taxonomy.dwca.zip
inattax <- read.csv('data/taxa.csv')
inatgen <- subset(inattax, !is.na(genus) & !is.na(family) & kingdom %in% c("Plantae","Chromista","Fungi","Bacteria") & !genus %in% '' & !family %in% '', select=c("kingdom","phylum","class","order","family","genus")) |> unique()
inatspp <- subset(inattax, !is.na(genus) & !is.na(family) & kingdom %in% c("Plantae","Chromista","Fungi","Bacteria") & !genus %in% '' & !family %in% '', select=c("kingdom","phylum","class","order","family", "scientificName")) |> unique()
nva.all <- nva.all |> mutate(genus = extractTaxon(actaxon, 'genus'))
nva.all <- nva.all |>   left_join(inatgen)

nva.misstax <- subset(nva.all, is.na(family), select=-c(kingdom,phylum,class,order,family)) |> mutate(genus=extractTaxon(actaxon, 'genus'))
#nva.misstax <- nva.misstax |> left_join(inatspp, by=join_by(taxon==scientificName))
#nva.misstax2 <- subset(nva.misstax, is.na(family), select=-c(kingdom,phylum,class,order,family))

nvaac <- subset(nva, select=c(numberOfOccurrences, kingdom, phylum, class, order, family, acceptedScientificName)) |> mutate(genus=extractTaxon(acceptedScientificName, 'genus'), actaxon = extractTaxon(acceptedScientificName)) |> group_by(kingdom, phylum, class, order, family, genus) |> summarise(n=sum(numberOfOccurrences)) |> ungroup()

nvataxonomy <- subset(nvaac, genus %in% nva.misstax$genus, select=c(n, kingdom, phylum, class, order, family, genus)) |> arrange(kingdom, phylum, class, order, family, genus) |> unique() |> group_by(genus) |> mutate(n2 = length(genus)) |> ungroup()

correct0 = data.frame(
  genus = c('Spiloma', 'Naevia', 'Lithothamnion' ,'Elharveya','Chaetomitrium','Volvox','Verdigellas','Variolaria','Valdonia','Urospora','Trismegistia','Timmiella','Thrombium','Thorea','Tetmemorus','Tephromela','Teilingia','Syringoderma','Symphyodon','Striaria','Stigonema','Stictyosiphon','Stichococcus','Stereocaulon','Staurodesmus','Staurastrum','Sporastatia','Splachnobryum', 'Speerschneidera', 'Sematophyllum','Scouleria','Sclerococcum','Schizochlamys','Schimmelmannia','Schaereria','Saelania', 'Ropalospora', 'Rivularia','Rhadinocladia','Pyrenastrum', 'Pycnora','Pterospermella','Pterogonidium', 'Pterobryella', 'Pseudolithophyllum', 'Pseudolithoderma','Pseudohypnella','Protoderma','Porina','Polycoccum','Pleurotaenium','Pleopsidium','Plectonema','Platygramme','Pilotrichella','Phymatoceros','Phaeostrophion','Phaeographina','Pellia','Pallavicinia','Orthodontium','Opegrapha','Nowellia','Nothoceros','Neohodgsonia','Myriotrichia','Mycomicrothelia','Mycobilimbia', 'Muellerella','Microzonia','Microthamnion','Microphiale','Micromitrium','Microcoleus','Micrasterias','Mesochaete','Mastigocoleus','Malcolmiella', 'Malbranchea','Luisierella','Lopadium', 'Lobothallia', 'Lithoderma', 'Limbella', 'Lichenostigma', 'Leucomium', 'Leptotheca', 'Leptopterigynandrum', 'Leptobryum','Lecidea', 'Leathesia','Klebsormidium', 'Kirchneriella','Jamesoniella','Hypopterygium', 'Hypodontium', 'Hypocenomyce','Hypnella','Hymenodon','Hydropogonella','Hydrocoleum','Hyalotheca','Hummia','Hookeriopsis','Holodontium','Heterocladium','Helminthocarpon','Haplogloia','Halorhipis','Halecania','Grateloupia','Graphina','Gracilariopsis','Glyphomitrium','Glyphium','Gloeocystis','Geminella','Garovaglia','Fimbriaria','Euptychium','Euastrum','Epibryon','Drummondia','Dolichospermum','Distichium','Dirinaria','Diphyscium','Diorygma','Dimerella','Dictyosiphon','Delamarea','Dactylospora','Cyrtopus','Cyathophorum','Crucigenia','Crocynia','Cosmocladium','Cosmarium','Compsonema','Coenogonium','Circinaria','Chrysoblastella','Chroodactylon','Chondria','Characium','Chaetosphaeridium','Chaetomitrium','Ceramium','Carbacanthographis','Capsosiphon','Calyptrozyma','Calothrix','Brachytrichia','Bissetia','Bilimbia','Baculifera','Audouinella','Asterella','Aspicilia', 'Aspergillus','Arthopyrenia', 'Ankyra','Aneura','Analipus','Phaeoceros','Hymenoloma','Halochlorococcum','Tolypothrix','Punctaria','Micrasterias','Leathesia','Hypnea','Cornicularia','Coilodesme','Loxospora', 'Aspicilia','Soranthera','Schizothrix','Amandinea','Dirinaria','Opegrapha','Porina','Rhizofabronia','Tephromela','Ulvella','Wilsoniella'),

  family = c('Xylographaceae', 'Arthoniaceae', 'Lithophyllaceae', 'Hypnaceae' ,'Symphyodontaceae','Volvocaceae','Palmophyllaceae','Pertusariaceae','Seligeriaceae', 'Ulotrichaceae','Pylaisiadelphaceae','Timmiellaceae','Protothelenellaceae','Thoreaceae','Desmidiaceae','Tephromelataceae','Desmidiaceae','Syringodermataceae','Symphyodontaceae','Chordariaceae','Stigonemataceae','Chordariaceae','Prasiolaceae','Stereocaulaceae','Desmidiaceae','Desmidiaceae','Sporastatiaceae','Splachnobryaceae','Leprocaulaceae','Sematophyllaceae', 'Scouleriaceae',	'Sclerococcaceae','Schizochlamydaceae','Schimmelmanniaceae','Schaereriaceae','Saelaniaceae', 'Ropalosporaceae', 'Nostocaceae','Chordariaceae','Pyrenulaceae','Pycnoraceae','Pterospermataceae','Pylaisiadelphaceae','Pterobryellaceae', 'Corallinaceae', 'Lithodermataceae', 'Hypnaceae','Chaetophoraceae','Porinaceae','Polycoccaceae','Desmidiaceae','Acarosporaceae','Oscillatoriaceae','Graphidaceae','Lembophyllaceae','Phymatocerotaceae','Phaeostrophionaceae','Graphidaceae','Pelliaceae','Pallaviciniaceae','Orthodontiaceae','Opegraphaceae','Cephaloziaceae','Dendrocerotaceae','Neohodgsoniaceae','Chordariaceae','Trypetheliaceae','Ramalinaceae','Verrucariaceae','Syringodermataceae','Microthamniaceae','Coenogoniaceae','Micromitriaceae','Microcoleaceae','Desmidiaceae','Aulacomniaceae','Nostochopsidaceae','Pilocarpaceae','Onygenaceae','Timmiellaceae','Lopadiaceae', 'Megasporaceae', 'Lithodermataceae', 'Amblystegiaceae','Phaeococcomycetaceae', 'Leucomiaceae', 'Orthodontiaceae', 'Thuidiaceae', 'Bryaceae', 'Lecideaceae', 'Chordariaceae','Klebsormidiaceae', 'Selenastraceae','Adelanthaceae', 'Hypopterygiaceae', 'Hypodontiaceae','Ophioparmaceae','Pilotrichaceae','Orthodontiaceae','Sematophyllaceae','Microcoleaceae','Desmidiaceae','Chordariaceae','Pilotrichaceae','Hypodontiaceae','Lembophyllaceae','Graphidaceae','Chordariaceae','Chordariaceae','Leprocaulaceae','Halymeniaceae','Graphidaceae','Gracilariaceae','Rhabdoweisiaceae','Mytilinidiaceae','Radiococcaceae','Chlorellaceae','Ptychomniaceae','Rhodomelaceae','Ptychomniaceae','Desmidiaceae','Epibryaceae','Drummondiaceae','Aphanizomenonaceae','Distichiaceae','Caliciaceae','Diphysciaceae','Graphidaceae','Coenogoniaceae','Chordariaceae','Chordariaceae','Dactylosporaceae','Cyrtopoaceae','Hypopterygiaceae','Scenedesmaceae','Ramalinaceae','Desmidiaceae','Desmidiaceae','Scytosiphonaceae','Coenogoniaceae','Megasporaceae','Chrysoblastellaceae','Stylonemataceae','Rhodomelaceae','Characiaceae','Chaetosphaeridiaceae','Symphyodontaceae','Ceramiaceae','Graphidaceae','Ulotrichaceae','Helotiaceae','Rivulariaceae','Scytonemataceae','Neckeraceae','Ramalinaceae','Caliciaceae','Acrochaetiaceae', 'Aytoniaceae','Megasporaceae', 'Aspergillaceae', 'Trypetheliaceae', 'Sphaeropleaceae','Aneuraceae','Scytosiphonaceae','Notothyladaceae','Hymenolomataceae','Oltmannsiellopsidaceae','Tolypothrichaceae','Chordariaceae','Desmidiaceae','Chordariaceae','Cystocloniaceae','Parmeliaceae','Chordariaceae','Sarrameanaceae', 'Megasporaceae','Chordariaceae','Trichocoleusaceae','Caliciaceae','Caliciaceae','Opegraphaceae','Porinaceae','Rhizofabroniaceae','Tephromelataceae','Ulvellaceae','Ditrichaceae'), correct=1)

nvataxonomy <- left_join(nvataxonomy, correct0) |> unique() |> mutate(correct = ifelse(is.na(correct),0,correct))
nvataxonomy <- nvataxonomy |> group_by(genus) |> mutate(maxocc=max(n), maxcor = max(correct)) |> ungroup()
nvataxonomy <- nvataxonomy |> subset((family != '' & n2 %in% 1 & nchar(order) > 1) | correct %in% 1 | (family != '' & nchar(order) > 1 & n > 10 & maxocc == n & maxcor == 0), select=-c(maxocc, maxcor,correct,n,n2))#| !genus %in% correct0$genus)#, select=-c(correct,n))
nva.misstax <- nva.misstax |> left_join(nvataxonomy)
nva.misstax <- subset(nva.misstax, !is.na(family))
nva.all0 <- rbind(nva.all, nva.misstax)
nva.all <- subset(nva.all0, !is.na(family)) |> mutate(genus=extractTaxon(actaxon, 'genus')) |> unique() |> arrange(kingdom, phylum, class, order, family, genus, actaxon, taxon)

nvataxonomy <- subset(nva.all, select=c(kingdom, phylum, class, order, family, genus)) |> unique()
# ordermissing <- subset(nvataxonomy, order =='')
# correct = data.frame(rbind(c('Aphanopsidaceae','Lecanorales','Lecanoromycetes','Ascomycota'),
#                            c('Characeae','Charales','Charophyceae','Charophyta'),
#                            c('Epigloeaceae','Lecanorales','Lecanoromycetes','Ascomycota'),
#                            c('Aphanopsidaceae','Lecanorales','Lecanoromycetes','Ascomycota'),

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
                                                       family %in% 'Pyrenidiaceae' ~ 'lichen',
                                                       family %in% 'Strangosporaceae' ~ 'lichen',
                                                       family %in% 'Harpidiaceae' ~ 'lichen',
                                                       family %in% 'Aphanopsidaceae' ~ 'lichen',
                                                       phylum %in% 'Anthocerotophyta' ~ 'bryophyte',
                                                       phylum %in% 'Bryophyta' ~ 'bryophyte',
                                                       phylum %in% 'Chlorophyta' ~ 'green algae',
                                                       phylum %in% 'Marchantiophyta' ~ 'bryophyte',
                                                       phylum %in% 'Chlorophyta' ~ 'green algae',
                                                       phylum %in% 'Rhodophyta' ~ 'red algae',
                                                       phylum %in% 'Charophyta' ~ 'green algae',
                                                       genus %in% 'Pyrenothrix' ~ 'lichen',
                                                       usdalichen %in% 1 | usdalichen2 %in% 1 ~ '2lich',
                                                       TRUE ~ 'unk')) |> subset(type != 'unk' & !type %in% '2lich', select=-c(usdalichen,usdalichen2))
nvataxonomy <- nvataxonomy |> group_by(genus) |> mutate(n = length(genus)) |> ungroup()
nva.all.tax <- nva.all |> left_join(nvataxonomy[,c('genus', 'type')]) |> subset(!is.na(type)) |> unique()

write.csv(nvataxonomy, 'data/nvataxonomy.csv', row.names = F, na='')
write.csv(nva.all.tax, 'data/nva.csv', row.names = F, na='')

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
  #thistax <- 'Vaucheria'

  thisass <- ass[grepl(paste0(thistax,'\\W|',thistax,'$'), ass$scientificName),]$ELEMENT_GLOBAL_ID
  if(length(thisass)>0){
    spbyass0 <- data.frame(taxon=thistax, ELEMENT_GLOBAL_ID=thisass)
    if(is.null(spbyass)){spbyass <- spbyass0}else{spbyass <- rbind(spbyass,spbyass0)}}

  thisdesc <- desc[grepl(paste0(thistax,'\\W|',thistax,'$'), desc$Floristics),]$ELEMENT_GLOBAL_ID
  if(length(thisdesc)>0){
    spbydesc0 <- data.frame(taxon=thistax, ELEMENT_GLOBAL_ID=thisdesc)
    if(is.null(spbydesc)){spbydesc <- spbydesc0}else{spbydesc <- rbind(spbydesc,spbydesc0)}}
}
Sys.time() - timA#Time difference of 2.653698 hours
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
#vascular forms
habs <- vegnasis::taxon.habits
#nonvasculars
nva <- read.csv('data/nva.csv')
nvataxonomy <- read.csv('data/nvataxonomy.csv')

#fix misapplied genera in respective vascular and non-vascular lists
familylink <- subset(familylink, !genus %in% c('Botrydium','Hedwigia','Hookeria','Ephemerum','Carteria'))
nvataxonomy0 <-  subset(nvataxonomy, genus %in% familylink$genus)
spp0  <-  spp |> mutate(genus=extractTaxon(kew, 'genus')) |> subset(genus %in% nvataxonomy$genus)
nvataxonomy <- subset(nvataxonomy, !genus %in% familylink$genus)



usda <- read.csv('data/plantlst.txt')
usda <- usda |> mutate(taxon = extractTaxon(Scientific.Name.with.Author), auth=extractTaxon(Scientific.Name.with.Author, 'author'),genus=extractTaxon(Scientific.Name.with.Author, 'genus'), acc=Symbol, syn=ifelse(Synonym.Symbol %in% '',Symbol,Synonym.Symbol))
#weed out homonyms
usda <- usda |> subset(!grepl('^non\\s|\\snon\\s|auct|sensu\\s',auth))

#Narrow to only non-vasculars then create a core set of taxa to include
missing.lichen.fams <- usda |> subset((grepl('lichen|moss|liverwort',Common.Name)|genus %in% nva$genus) & !Family %in% nvataxonomy$family & !Family %in% apg$family, select=Family) |> unique()


#write.csv(missing.lichen.fams, 'data/missing.lichen.fams3.csv', row.names = F)
missing.lichen.fams <- read.csv('data/missing.lichen.fams.csv')
missing.lichen.fams <- missing.lichen.fams |> subset(!type %in% 'na')
#list of accepteds
usda.backbone <- subset(usda, (Family %in% nvataxonomy$family | Family %in% missing.lichen.fams$Family) & grepl(' ', taxon), select=c("Common.Name","Family", "taxon","auth","acc"))
usda.acc <- data.frame(acc = usda.backbone$acc, actaxon=usda.backbone$taxon, acauth=usda.backbone$auth, taxon=usda.backbone$taxon, auth=usda.backbone$auth)
#list of synonyms
usda.syns <- usda |> subset(syn!=acc, select=c("taxon","auth","acc"))
usda.df1 <- data.frame(acc = usda.backbone$acc, actaxon=usda.backbone$taxon, acauth=usda.backbone$auth)
usda.syns <- usda.df1 |> left_join(usda.syns) |> subset(!is.na(taxon))

usda.df <- usda.acc |> rbind(usda.syns)
#need to remove redundant names...
usda.df <- usda.df |> group_by(taxon) |> mutate(n=length(taxon)) |> ungroup()
redund <- usda.df |> subset(n > 1) |> arrange(taxon); # write.csv(redund, 'data/redund2.csv', row.names = F, na='', fileEncoding = "latin1")
redund <- read.csv('data/redund2.csv',fileEncoding = "latin1")
usda.df <- usda.df |> left_join(redund[,c('taxon', 'auth', 'keep')])
usda.df <- usda.df |> subset(n == 1 | keep==1, select=-c(n,keep))
nva.new <-  usda.df |> subset(!taxon %in% nva$taxon)|> group_by(taxon) |> mutate(n=length(taxon)) |> ungroup() |> arrange(taxon)
#---combine nva with usda ----
#nva taxa
nva1 <- data.frame(taxon = nva$taxon, auth = nva$auth, gbif = nva$actaxon) |> unique()#|> group_by(taxon) |> mutate(n=length(taxon)) |> ungroup()
nva1a <- data.frame(taxon=usda.df$taxon, usda=usda.df$actaxon) |> unique() #|> group_by(taxon) |> mutate(n=length(taxon)) |> ungroup()
nva1 <- nva1 |> left_join(nva1a)
nva2 <- data.frame(taxon = nva.new$taxon, auth = nva.new$auth, gbif = NA, usda = nva.new$actaxon)
nva1 <- nva1 |> rbind(nva2)
nva1 <- nva1 |> mutate(gbif=case_when(taxon %in% c('Polychidium intricatulum','Polychidium umhausense') ~ 'Ricasolia amplissima',
                                      TRUE ~ gbif))

#rejoin taxonomy
nva <- read.csv('data/nva.csv')
nvataxonomy <- nva |> subset(!phylum %in% c('Foraminifera','Ciliophora') & !family %in% c('Chaetomiaceae'), select = c("kingdom","phylum","class","order","family","genus","type")) |> unique()
nva1.taxonomy <- nva1 |> mutate(genus1 = extractTaxon(gbif, 'genus'), genus2=extractTaxon(usda, 'genus'), genus=ifelse(is.na(gbif), genus2, genus1)) |> unique()
nva1.taxonomy <- nva1.taxonomy |> left_join(nvataxonomy)
inattax <- read.csv('data/taxa.csv')
inatgen <- subset(inattax, !is.na(genus) & !is.na(family) & kingdom %in% c("Plantae","Chromista","Fungi","Bacteria") & !genus %in% '' & !family %in% '', select=c("kingdom","phylum","class","order","family","genus")) |> unique()
inatfam <- subset(inattax, !is.na(genus) & !is.na(family) & kingdom %in% c("Plantae","Chromista","Fungi","Bacteria") & !genus %in% '' & !family %in% '', select=c("kingdom","phylum","class","order","family")) |> unique()
colnames(inatgen) <- c("kingdom1","phylum1","class1","order1","family1","genus")
nva1.taxonomy <- nva1.taxonomy |> left_join(inatgen)
nva1.taxonomy <- nva1.taxonomy |> mutate(kingdom = ifelse(is.na(kingdom), kingdom1, kingdom),
                                         phylum = ifelse(is.na(phylum), phylum1, phylum),
                                         class = ifelse(is.na(class), class1, class),
                                         order = ifelse(is.na(order), order1, order),
                                         family = ifelse(is.na(family), family1, family)) |> subset(select= -c(kingdom1,phylum1,class1,order1,family1))
missingfams <- subset(nva1.taxonomy, is.na(family))
missingfams <- nva1.taxonomy |> subset(genus2 %in% missingfams$genus & !is.na(family), select = c("kingdom","phylum","class","order","family","genus2")) |> unique()
colnames(missingfams) <- c("kingdom1","phylum1","class1","order1","family1","genus")
nva1.taxonomy <- nva1.taxonomy |> left_join(missingfams)
nva1.taxonomy <- nva1.taxonomy |> mutate(kingdom = ifelse(is.na(kingdom), kingdom1, kingdom),
                                         phylum = ifelse(is.na(phylum), phylum1, phylum),
                                         class = ifelse(is.na(class), class1, class),
                                         order = ifelse(is.na(order), order1, order),
                                         family = ifelse(is.na(family), family1, family)) |> subset(select= -c(kingdom1,phylum1,class1,order1,family1))
missingfams <- subset(nva1.taxonomy, is.na(family)) |> unique()
correct= data.frame(rbind(c('Anamylospora',	'Baeomycetaceae'),
                          c('Solenospora',	'Catillariaceae'),
                          c('Eopyrenula',	'Dacampiaceae'),
                          c('Hafellnera',	'Schaereriaceae'),
                          c('Lauderlindsaya',	'Verrucariaceae'),
                          c('Pyrenocollema',	'Xanthopyreniaceae')))
colnames(correct) <- c('genus', 'family')
correct <- correct |> left_join(inatfam)
colnames(correct) <- c("genus","family1","kingdom1","phylum1","class1","order1")
nva1.taxonomy <- nva1.taxonomy |> left_join(correct)
nva1.taxonomy <- nva1.taxonomy |> mutate(kingdom = ifelse(is.na(kingdom), kingdom1, kingdom),
                                         phylum = ifelse(is.na(phylum), phylum1, phylum),
                                         class = ifelse(is.na(class), class1, class),
                                         order = ifelse(is.na(order), order1, order),
                                         family = ifelse(is.na(family), family1, family)) |> subset(select= -c(kingdom1,phylum1,class1,order1,family1))
nva1.taxonomy <- nva1.taxonomy |> mutate(type =  case_when(phylum %in% 'Cyanobacteria' ~ 'Cyanobacteria',
                                                       class %in% 'Phaeophyceae' ~ 'brown algae',
                                                       class %in% 'Xanthophyceae' ~ 'yellow-green algae',
                                                       order %in% 'Arthoniales' ~ 'lichen',
                                                       order %in% 'Vezdaeales' ~ 'lichen',
                                                       order %in% 'Monoblastiales' ~ 'lichen',
                                                       order %in% 'Strigulales' ~ 'lichen',
                                                       order %in% 'Thelocarpales' ~ 'lichen',

                                                       #class %in% 'Eurotiomycetes' ~ 'lichen',
                                                       class %in% 'Lecanoromycetes' ~ 'lichen',
                                                       class %in% 'Lichinomycetes' & !order %in% 'Sareales' ~ 'lichen',
                                                       class %in% 'Candelariomycetes' ~ 'lichen',
                                                       class %in% 'Coniocybomycetes' ~ 'lichen',
                                                       family %in% 'Trypetheliaceae' ~ 'lichen',
                                                       family %in% 'Pyrenulaceae' ~ 'lichen',
                                                       family %in% 'Verrucariaceae' ~ 'lichen',
                                                       family %in% 'Pyrenidiaceae' ~ 'lichen',
                                                       family %in% 'Strangosporaceae' ~ 'lichen',
                                                       family %in% 'Harpidiaceae' ~ 'lichen',
                                                       family %in% 'Aphanopsidaceae' ~ 'lichen',
                                                       phylum %in% 'Anthocerotophyta' ~ 'bryophyte',
                                                       phylum %in% 'Bryophyta' ~ 'bryophyte',
                                                       phylum %in% 'Chlorophyta' ~ 'green algae',
                                                       phylum %in% 'Marchantiophyta' ~ 'bryophyte',
                                                       phylum %in% 'Chlorophyta' ~ 'green algae',
                                                       phylum %in% 'Rhodophyta' ~ 'red algae',
                                                       phylum %in% 'Charophyta' ~ 'green algae',
                                                       genus %in% 'Pyrenothrix' ~ 'lichen',
                                                       TRUE ~ 'unk')) |> subset(type != 'unk')
nva1.taxonomy <- nva1.taxonomy |> unique() |> group_by(taxon) |> mutate(n = length(genus)) |> ungroup()
nva2.taxonomy <- nva1.taxonomy |> select(c(kingdom,phylum,class,order,family, type))|> unique() |> group_by(family, type) |> mutate(n = length(family)) |> ungroup()
colnames(inatfam) <- c("kingdom1","phylum1", "class1","order1","family")
familytaxonomy <- nva2.taxonomy |> unique() |> left_join(inatfam)
familytaxonomy <- familytaxonomy |> mutate(wrong = ifelse(order != order1 | class != class1, 1,0))
familytaxonomy <- familytaxonomy |> subset(n == 1 | wrong ==0)
familytaxonomy <- familytaxonomy |> mutate(kingdom = ifelse(wrong == 1 & !is.na(wrong), kingdom1, kingdom),
       phylum = ifelse(wrong == 1 & !is.na(wrong), phylum1, phylum),
       class = ifelse(wrong == 1 & !is.na(wrong), class1, class),
       order = ifelse(wrong == 1 & !is.na(wrong), order1, order))
familytaxonomy <- familytaxonomy |> select(c(kingdom,phylum,class,order,family, type)) |> arrange(kingdom,phylum,class,order,family)
genustaxonomy <- nva1.taxonomy |> select(c(family, genus, type))|> unique() #|> group_by(family, genus) |> mutate(n = length(family)) |> ungroup()
nva <- nva1.taxonomy |> select(c(taxon, auth, gbif, usda))|> unique()
write.csv(familytaxonomy, 'nvafamilytaxonomy.csv', na='', row.names = F)
write.csv(genustaxonomy, 'nvagenustaxonomy.csv', na='', row.names = F)
write.csv(nva, 'nvanomenclature.csv', na='', row.names = F)


# nvaextra <- nva1 |> mutiate(usda1=usda) |> subset(!is.na(usda) & !is.na(gbif), select=c(gbif, usda1)) |> unique()
# nvaextra <- nvaextra |> group_by(gbif) |> mutate(n=length(gbif)) |> ungroup()
# nva.fill <- nva1 |> left_join(nvaextra)
# nva.fill <- nva.fill |> mutate(usda=ifelse(is.na(usda), usda1,usda))|> select(-c(usda1)) |> unique()
# nva.fill <- nva.fill |> group_by(taxon, gbif, usda) |> mutate(n=length(taxon)) |> ungroup()




#join with USNVC
library(vegnasis)
#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#USNVC association names
ass <- read.csv('data/USNVC v3.0.3 2026-03-25/unit.csv')
#USNVC species mentions
desc <- read.csv('data/USNVC v3.0.3 2026-03-25/unitDescription.csv')

bio <- subset(ass, hierarchyLevel %in% 'Biome', select=c('PARENT_ID','ELEMENT_GLOBAL_ID','databaseCode','colloquialName'))
sbi <- subset(ass, hierarchyLevel %in% 'Subbiome', select=c('PARENT_ID','ELEMENT_GLOBAL_ID','databaseCode','colloquialName'))
ebi <- subset(ass, hierarchyLevel %in% 'Ecobiome', select=c('PARENT_ID','ELEMENT_GLOBAL_ID','databaseCode','colloquialName'))
dvn <- subset(ass, hierarchyLevel %in% 'Division', select=c('PARENT_ID','ELEMENT_GLOBAL_ID','databaseCode','colloquialName'))
mgp <- subset(ass, hierarchyLevel %in% 'Macrogroup', select=c('PARENT_ID','ELEMENT_GLOBAL_ID','databaseCode','colloquialName'))
grp <- subset(ass, hierarchyLevel %in% 'Group', select=c('PARENT_ID','ELEMENT_GLOBAL_ID','databaseCode','colloquialName'))
ali <- subset(ass, hierarchyLevel %in% 'Alliance', select=c('PARENT_ID','ELEMENT_GLOBAL_ID','databaseCode', 'scientificName'))
ass1 <- subset(ass, hierarchyLevel %in% 'Association', select=c('PARENT_ID','ELEMENT_GLOBAL_ID','databaseCode', 'scientificName'))
bio <- bio |> mutate(databaseCode = as.factor(databaseCode))
sbi <- sbi |> mutate(databaseCode = as.factor(databaseCode))
ebi <- ebi |> mutate(databaseCode = as.factor(databaseCode))
dvn <- dvn |> mutate(databaseCode = as.factor(databaseCode))
mgp <- mgp |> mutate(databaseCode = as.factor(databaseCode))
grp <- grp |> mutate(databaseCode = as.factor(databaseCode))
ali <- ali |> mutate(databaseCode = as.factor(databaseCode))
ass1 <- ass1 |> mutate(databaseCode = as.factor(databaseCode))
colnames(bio) <- c('realmID','biomeID','biomeCode','biome')
colnames(sbi) <- c('biomeID','subbiomeID','subbiomeCode','subbiome')
colnames(ebi) <- c('subbiomeID','ecobiomeID','ecobiomeCode','ecobiome')
colnames(dvn) <- c('ecobiomeID','divisionID','divisionCode','division')
colnames(mgp) <- c('divisionID','macrogroupID','macrogroupCode','macrogroup')
colnames(grp) <- c('macrogroupID','groupID','groupCode','group')
colnames(ali) <- c('groupID','allianceID','allianceCode','alliance')
colnames(ass1) <- c('allianceID','associationID','associationCode','association')
ass1 <- bio |> left_join(sbi) |> left_join(ebi)|> left_join(dvn)|> left_join(mgp)|> left_join(grp)|> left_join(ali)|> left_join(ass1)
ass1 <- ass1 |> mutate(ELEMENT_GLOBAL_ID =  case_when(!is.na(associationID) ~ associationID,
                                                      !is.na(allianceID) ~ allianceID,
                                                      !is.na(groupID) ~ groupID,
                                                      !is.na(macrogroupID) ~ macrogroupID,
                                                      !is.na(divisionID) ~ divisionID,
                                                      !is.na(ecobiomeID) ~ ecobiomeID,
                                                      !is.na(subbiomeID) ~ subbiomeID,
                                                      !is.na(biomeID) ~ biomeID))
ass1 <- subset(ass1, select = -c(biomeID,subbiomeID,subbiomeID,ecobiomeID,divisionID,macrogroupID,groupID,allianceID))

usfseco <- desc |> subset(select = c(ELEMENT_GLOBAL_ID, Subnations, usfsEcoregions2007))
ass1 <- ass1 |> left_join(usfseco, by=join_by(ELEMENT_GLOBAL_ID==ELEMENT_GLOBAL_ID))

usnvcspp <- read.csv('usnvcspp2.csv')
# usnvcspp <- usnvcspp |> mutate(taxon=harmonize.taxa(taxon, fix=TRUE, sensu = 'kew'), taxon=extractTaxon(taxon, 'binomial'))
usnvcgen <- read.csv('usnvcgen.csv')
usnvcnva <- read.csv('usnvcnva.csv')

usnvcspp <- rbind(usnvcspp[,colnames(usnvcnva)], usnvcnva, usnvcgen) |> unique()
usnvcspp <- usnvcspp |> mutate(genus = extractTaxon(taxon, 'genus'))
usnvcspp <- usnvcspp |> group_by(genus, ELEMENT_GLOBAL_ID, indicator) |> mutate(n=length(taxon)) |> ungroup()
usnvcspp <- usnvcspp |> mutate(flag=ifelse((genus != taxon | n == 1) & nchar(taxon) >= 1,0,1))
usnvcspp <- usnvcspp |> subset(!flag %in% 1, select=-c(genus, n, flag))
# nvaforms <- nvataxonomy |> subset(select=c(genus, type)) |> unique()

usnvcspp <- usnvcspp |> mutate(type = vegnasis::fill.type(taxon), habit=get.habit.name(get.habit.code(taxon)))
usnvcspp$type <- factor(usnvcspp$type, levels=c("tree","shrub/vine","forb","grass/grasslike","moss","lichen","microbiotic crust")) 

usnvcspp3 <- ass[,c('classificationCode','databaseCode','PARENT_ID','ELEMENT_GLOBAL_ID','hierarchyLevel','colloquialName', 'scientificName')] |> left_join(usnvcspp) |> arrange(classificationCode, scientificName,  type, desc(indicator), taxon)  |> left_join(usfseco)#|> subset(hierarchyLevel %in% 'Association')

gennul <- subset(usnvcspp3, type %in% 'NA' | is.na(type), select = taxon) |> unique() 
usnvcspp3 <- subset(usnvcspp3, !(type %in% 'NA' | is.na(type))) |> unique()
write.csv(usnvcspp3, 'Plants_by_Plant_Association.csv', row.names = F, na='')
write.csv(ass1, 'usnvc3hierarchy.csv', row.names = F, na='')
saveRDS(usnvcspp3, 'usnvcspp3.RDS')
saveRDS(ass1, 'ass1.RDS')



#Query communities ----
library(vegnasis)
#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#load data ----
bmlra <- readRDS('bmlra.RDS')
mnfi <- readRDS('mnfi.RDS')
ass1 <- readRDS('ass1.RDS')
usnvcspp3 <- readRDS('usnvcspp3.RDS')

#query inputs ----
region <- 'Northeast' #native region for species
thisMLRA <- '96' #native MLRA for species
states <- c('MI') #narrow list of communities
ecoreg <- c('212H') #narrow list of communities
species <- c('Tsuga canadensis')
ElementCode <- c('CEGL005163') #Lookup specific community
thisMNFI <- c('Mesic Northern Forest') #Lookup MNFI community

#run queries ----
q1 <- paste0('\\b',states,'\\b', collapse  = '|')
q2 <- paste0('\\b',ecoreg,'.\\b', collapse  = '|')
q3 <- paste0('^',species,'$', collapse  = '|')
qspp <- usnvcspp3 |> subset(grepl(q3, taxon))
mlraspp <- bmlra |> subset(MLRA %in% thisMLRA)
q4 <- paste0('^',ElementCode,'$', collapse  = '|')

qass <- ass1 |> subset(((grepl(q1, Subnations) & grepl(q2, usfsEcoregions2007))| 
                         (grepl(q1, Subnations) & nchar(usfsEcoregions2007)<2)| 
                         (grepl(q2, usfsEcoregions2007) & nchar(Subnations)<2)) & ELEMENT_GLOBAL_ID %in% qspp$ELEMENT_GLOBAL_ID)

qmnfi <- mnfi |> subset(community %in% thisMNFI, select=c(community, taxa, adund)) |>   
  mutate(accepted =  harmonize.taxa(taxa, TRUE), 
         type = fill.type(accepted),
         type = factor(type, levels = c("tree","shrub/vine","forb","grass/grasslike","moss","lichen","microbiotic crust")),
         habit = get.habit.name(get.habit.code(accepted)),
         nativity = fill.nativity(accepted, region = region),
         mlranative = case_when(accepted %in% mlraspp$ac.binomial ~ 1,                                                                                                                     type %in% c('moss', 'lichen', 'microbiotic crust') ~ 1,                                                                                                                           !grepl('\\s',taxa) ~ 1,
                                TRUE ~ 0),
         acceptedusda =  harmonize.taxa(accepted, TRUE, 'usda'),
         usdasym = fill.usda.symbols(acceptedusda),
         acceptedusda2 =  harmonize.taxa(taxa, TRUE, 'usda'),
         usdasym2 = fill.usda.symbols(acceptedusda2),
         usdasym = ifelse(is.na(usdasym),usdasym2,usdasym),
         acceptedusda=NULL,usdasym2=NULL,acceptedusda2=NULL) |> arrange(community, type, taxa)

splist <- usnvcspp3 |> subset(grepl(q4,databaseCode), select=c(databaseCode, indicator, taxon, type, habit))  |> 
  mutate(accepted =  harmonize.taxa(taxon, TRUE), 
         nativity = fill.nativity(accepted, region = region),
         mlranative = case_when(accepted %in% mlraspp$ac.binomial ~ 1,                                                                                                                     type %in% c('moss', 'lichen', 'microbiotic crust') ~ 1,                                                                                                                           !grepl('\\s',taxon) ~ 1,
                                TRUE ~ 0),
         accepted =  harmonize.taxa(taxon, TRUE, 'usda'),
         usdasym = fill.usda.symbols(accepted),
         accepted=NULL)
qmnfi
splist
