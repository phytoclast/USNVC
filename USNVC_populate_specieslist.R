library(stringr)
syn <- readRDS('data/syn.RDS')
unit <-  read.delim('data/unit.txt', encoding = 'UTF-8', na.strings = '', stringsAsFactors=FALSE)
listofacc <- as.character(unique(syn$acctaxon))
veg <- as.data.frame(cbind(element_global_id=0, scientificname='x',acctaxon ='x')) 
#i <- which(listofacc == 'Athyrium asplenioides')
#listofacc[i]
for (i in 1:length(listofacc)){
search <- syn[syn$acctaxon %in% listofacc[i],]$syntaxon

veg1 <- unit[
  grepl(paste(search, collapse = "|"),gsub("\\(", "", unit$scientificname))|
    grepl(paste(search, collapse = "|"),gsub(" \\(.{0,15},", "", unit$scientificname))|
    grepl(paste(search, collapse = "|"),gsub(" \\(.{0,15},.{0,15},", "", unit$scientificname))
  ,c('element_global_id', 'scientificname')]
if(nrow(veg1) > 0){
veg$element_global_id <- as.integer(veg$element_global_id)

veg1$acctaxon <- listofacc[i]
veg <- rbind(veg, veg1)}
}

#saveRDS(veg, 'output/veg.RDS')

