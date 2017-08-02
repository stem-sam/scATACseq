#module load gmp/5.0.2 mpfr/3.1.0 mpc/0.8.2 gcc/4.9.1 R/3.2.1
library(Matrix)
library(limma)
library(RColorBrewer)
library(irlba,lib.loc="/net/shendure/vol1/home/cusanovi/R/x86_64-unknown-linux-gnu-library/3.1/")
library(data.table)
library(Rtsne)
library(densityClust)
library(methods)
library(viridis)
# gz set from 'TRUE' to 'FALSE'
gz = FALSE
upperandlower = FALSE
cellcutoff=0.1
uppercellcutoff=1-cellcutoff
sitecutoff=30000
siteperccutoff=0.05
#inmat = "/net/shendure/vol10/projects/mouse_atlas/nobackup/tissues/Spleen_62016/SC_atac_combined.Spleen_62016_P2.true.nodups.peaks.dhsmatrix.txt"
inmat = "/net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/matrix/scatac_hEB.q10.sort.true.nodups_peaks.dhsmatrix.txt"
#cellsout = "/net/shendure/vol10/projects/mouse_atlas/nobackup/tissues/Spleen_62016/SC_atac_combined.Spleen_62016_P2.true.nodups.peaks.densitypeaks.cellclusters."
cellsout = "/net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/tSNE/scatac.summits.densitypeaks.cellclusters."
#tsneout = "/net/shendure/vol10/projects/mouse_atlas/nobackup/tissues/Spleen_62016/SC_atac_combined.Spleen_62016_P2.true.nodups.peaks.lesswrong.tsne.0.05cells.0.01cellsites."
if(upperandlower){
  tsneout = paste0("/net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/tSNE/scatac.summits.tsne.",cellcutoff,"cellsupperandlower.",siteperccutoff,"cellsites")
}else{
  tsneout = paste0("/net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/tSNE/scatac.peaks.lesswrong.tsne.",cellcutoff,"cells.",siteperccutoff,"cellsites")
}
#promwindows = read.table("/net/shendure/vol10/projects/mouse_atlas/nobackup/master_combined_peaks_update.within.2_5kb.of.tss.whitelist.bed")
promwindows = read.table("/net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/sort/macs/master_combined_peaks.within.2.5kb.of.tss.whitelist.bed")
countreport = read.table("/net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/fastq/scatac_hEB.q10.sort.true.nodups.nodups.cells.duplicate_report.txt",header=T)
oldclusters = read.table("/net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/AlphabetExperiments/sc_drone_dal/SCatac_dal.true.nodups.5kbwindows.0.1cells.0.05cellsites.6components.3clusters.cellclusters.txt")
#oldclusters = read.table("/net/shendure/vol10/projects/mouse_atlas/nobackup/tissues/Spleen_62016/SC_atac_combined.Spleen_62016_P2.true.nodups.5kbwindows.cellclusters.5percnew.txt")
#sexcolors = read.table("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/SCatac_tha.10to12.true.nodups.xycalls.txt")
#sexcounts = read.table("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/SCatac_tha.10to12.true.nodups.xycounts.txt",header=T)
numclust = length(unique(oldclusters[,2]))


#Colors from https://www.r-bloggers.com/the-paul-tol-21-color-salute/
tol2qualitative=c("#4477AA", "#CC6677")
tol3qualitative=c("#4477AA", "#DDCC77", "#CC6677")
tol4qualitative=c("#4477AA", "#117733", "#DDCC77", "#CC6677")
tol5qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677")
tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol8qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499")
tol9qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol11qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol12qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")
tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
tol15rainbow=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")
tol18rainbow=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
dc25rainbow = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","lightgray","darkgray","black","orange")

if(gz == TRUE){
	bigmat = fread(paste0("zcat ",inmat),header=T)
} else {
	bigmat = fread(inmat,header=T)
}
bigmat = as.data.frame(bigmat)

print("Reformatting data...")
mastersites = bigmat[,4]
sortmat = as.matrix(bigmat[,-c(1:4)])
nonzeroes = rowSums(sortmat) > 0

sortmat.nonzero = sortmat[nonzeroes,]
mastersites.nonzero = mastersites[nonzeroes]
currnonzs = unlist(apply(sortmat.nonzero,2,function(x) which(x > 0)))
namers = gsub('[[:digit:]]+', '', names(currnonzs))
namenums = c(1:length(unique(namers)))[unlist(as.factor(namers))]
mastermat = cbind(currnonzs,namenums)
namelist = unique(namers)

rm(bigmat)
rm(sortmat)
rm(sortmat.nonzero)
gc()

bigmat.bin = sparseMatrix(i=mastermat[,1],j=mastermat[,2],x=rep(1,times=dim(mastermat)[1]))
rownames(bigmat.bin) = mastersites.nonzero
colnames(bigmat.bin) = namelist
sites = strsplit2(mastersites.nonzero,"_")
label = paste0(sites[,1],":",sites[,2],"-",sites[,3])
annot = data.frame(row.names=label, cell_spec=mastersites.nonzero)

num_cells_ncounted = rowSums(bigmat.bin)

#annot.ncounts = rownames(annot)[num_cells_ncounted >= num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[sitecutoff]]]
#ncounts = bigmat.bin[num_cells_ncounted >= num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[sitecutoff]],]
annot.ncounts = rownames(annot)[num_cells_ncounted >= dim(bigmat.bin)[2]*siteperccutoff]
ncounts = bigmat.bin[num_cells_ncounted >= dim(bigmat.bin)[2]*siteperccutoff,]
new_counts = colSums(ncounts)
if(upperandlower){
  ncounts = ncounts[,new_counts >= quantile(new_counts,probs=cellcutoff) & new_counts <= quantile(new_counts,probs=uppercellcutoff)]
}else{
  ncounts = ncounts[,new_counts >= quantile(new_counts,probs=cellcutoff)]
}

annot.ncounts = annot.ncounts[rowSums(ncounts) > 0]
ncounts = ncounts[rowSums(ncounts) > 0,]
# sexsites = c(grep("chrY",rownames(ncounts)),grep("chrX",rownames(ncounts)))
# ncounts.nosex = ncounts[-sexsites,]
# annot.ncounts.nosex = annot.ncounts[-sexsites]

print("Normalizing data...")
nfreqs <- t(t(ncounts) / colSums(ncounts))
tf_idf_counts <- nfreqs * log(1 + ncol(ncounts) / rowSums(ncounts))
dim(tf_idf_counts)

print("Running tSNE...")
set.seed(0)
SVDtsne = irlba(tf_idf_counts, 50, 50)
d_diagtsne = matrix(0, nrow=length(SVDtsne$d), ncol=length(SVDtsne$d))
diag(d_diagtsne) <- SVDtsne$d
SVDtsne_vd = t(d_diagtsne %*% t(SVDtsne$v))

for (n in c(5, 50, 100)){
  set.seed(0)
  tsnetfidf = Rtsne(SVDtsne_vd, perplexity = n, pca=F)
  tsnedist = dist(tsnetfidf$Y)
  cell_pal <- brewer.pal(numclust, "Paired")
  dclust = densityClust(tsnedist,gaussian=T)
  rhoer=50
  deltar = 2.5
  dclust = findClusters(dclust, rho = rhoer, delta = deltar, verbose = verbose)
  pdf(paste0(tsneout,".perplexity",n,".rho",rhoer,".delta",deltar,".pdf"))
  plot(tsnetfidf$Y,pch=20,main="No Clusters")
  plot(tsnetfidf$Y,pch=20,col=cell_pal[oldclusters[match(colnames(tf_idf_counts),oldclusters[,1]),2]],main="Old Clusters")
  #plot(tsnetfidf$Y,pch=20,col=c("blue","red")[sexcolors[match(colnames(tf_idf_counts),sexcolors[,1]),2]],main="Sex Coloring")
  plot(tsnetfidf$Y,pch=20,col=dc25rainbow[as.factor(dclust$clusters)],main="New Peak Density Clusters")
  plot(tsnetfidf$Y,pch=20,col=dc25rainbow[as.factor(dclust$clusters)],main=paste0("New Peak Density Clusters\n",dim(tf_idf_counts)[1]," sites\n",dim(tf_idf_counts)[2]," cells"))
  text(tsnetfidf$Y[dclust$peaks,1],tsnetfidf$Y[dclust$peaks,2],labels=dclust$clusters[dclust$peaks],cex=3)
  plot(dclust$rho,dclust$delta,pch=20)
  points(dclust$rho[dclust$peaks],dclust$delta[dclust$peaks],col="red",pch=20)
  text(dclust$rho[dclust$peaks]-2,dclust$delta[dclust$peaks]+2,labels=dclust$clusters[dclust$peaks])
  abline(v=rhoer)
  abline(h=deltar)
  dev.off()

}

