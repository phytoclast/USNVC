library(stringr)
obsspp <- read.delim("data/Observed_Species.txt", encoding = 'UTF-8', na.strings = '')
obsspp[obsspp$AcTaxon == 'Phalaris' & !is.na(obsspp$AcTaxon) ,]$AcTaxon <- 'Phalaris arundinacea'
obsspp <- obsspp[!grepl("\\?", obsspp$AcTaxon) | is.na(obsspp$AcTaxon) ,]
obsspp <- subset(obsspp, !is.na(specific_epithet) | AcTaxon %in% c('Sphagnum', 'Chara'))
obsspp <- subset(obsspp, !substr(AcTaxon,1,1) %in% '-'& !AcTaxon %in% '' & !is.na(Habit))

splist <- as.data.frame(as.factor(unique(as.character(obsspp$AcTaxon))))

checklist <- read.delim("data/checklist.txt", encoding = 'UTF-8', na.strings = '')
accepted <- subset(checklist, grepl(' ',Scientific.Name) & Acsy == 'AC', select = c('Newseq','Scientific.Name', 'Author'))
colnames(accepted) <- c('id', 'acctaxon', 'accauth')
syn <- subset(checklist, grepl(' ',Scientific.Name) , select = c('Newac','Scientific.Name', 'Author'))
colnames(syn) <- c('id', 'syntaxon', 'synauth')

syn <- merge(accepted, syn, by='id')

dup <- aggregate(syn$syntaxon, by=list(syntaxon = syn$syntaxon), FUN='length')
dup <- (subset(dup, x>1))[,1]

syn <- subset(syn, !((grepl('non ', synauth)|grepl('auct', synauth)) & (syntaxon != acctaxon) & syntaxon %in% dup))

dup2 <- aggregate(syn$syntaxon, by=list(syntaxon = syn$syntaxon), FUN='length')
dup2 <- (subset(dup2, x>1))[,1]

dup3 <- subset(syn, syntaxon == acctaxon & syntaxon %in% dup2)[,'syntaxon']
syn <- subset(syn, !(syntaxon != acctaxon & syntaxon %in% dup3))
