library(stringr)
syn <- readRDS('data/syn.RDS')
unit <-  read.delim('data/unit.txt', encoding = 'UTF-8', na.strings = '', stringsAsFactors=FALSE)
unitDescription <-  read.delim('data/unitDescription.txt', encoding = 'UTF-8', na.strings = '', stringsAsFactors=FALSE, quote = "")
unitDescription <- merge(unit[,c('element_global_id', 'scientificname')], unitDescription, by='element_global_id')

listofacc <- as.character(unique(syn$acctaxon))
veg <- as.data.frame(cbind(element_global_id=0, scientificname='x',acctaxon ='x')) 
#i <- which(listofacc == 'Abies balsamea')
#listofacc[i]
for (i in 1:length(listofacc)){
search <- syn[syn$acctaxon %in% listofacc[i],]$syntaxon

veg1 <- unitDescription[
  grepl(paste(search, collapse = "|"),unitDescription$typeconcept)|
    grepl(paste(search, collapse = "|"),unitDescription$floristics)
  ,c('element_global_id', 'scientificname')]
if(nrow(veg1) > 0){
veg$element_global_id <- as.integer(veg$element_global_id)

veg1$acctaxon <- listofacc[i]
veg <- rbind(veg, veg1)}
}

saveRDS(veg, 'output/vegDesc.RDS')

