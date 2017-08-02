args = commandArgs(TRUE)

greatin = args[1]
greatout = args[2]

greatnow = read.table(greatin,skip=5,sep="\t")
clean_1 = greatnow[which(greatnow[,8] > 1.5),]
clean_2 = clean_1[which(clean_1[,19] > 5),]
clean_3 = clean_2[which(clean_2[,7] < 0.05 & clean_2[,16] < 0.05),]
clean_4 = clean_3[,1:22]
clean_5 = clean_4[order(clean_4[,5]),]
colnames(clean_5) = c("Ontology","ID","Desc","BinomRank","BinomP","BinomBonfP","BinomFdrQ","RegionFoldEnrich","ExpRegions","ObsRegions","GenomeFrac","SetCov","HyperRank","HyperP","HyperBonfP","HyperFdrQ","GeneFoldEnrich","ExpGenes","ObsGenes","TotalGenes","GeneSetCov","TermCov")
write.table(clean_5,greatout,row.names=F,quote=F,sep="\t")