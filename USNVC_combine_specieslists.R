vegDesc <- readRDS('output/vegDesc.RDS')
veg <- readRDS('output/veg.RDS')
veg$amount <- 1
vegDesc$amount <- 0.1
USNVClist <- rbind(veg,vegDesc)
USNVClist <- subset(USNVClist, element_global_id !=1)
USNVClist <- aggregate(USNVClist[,'amount'], by=list(element_global_id = USNVClist$element_global_id,
                                                     scientificname = USNVClist$scientificname,
                                                     acctaxon = USNVClist$acctaxon
                                                     ), FUN='max')
saveRDS(USNVClist, 'data/USNVClist.RDS')