#module load gmp/5.0.2 mpfr/3.1.0 mpc/0.8.2 gcc/4.9.1 R/3.2.1
library(Matrix)
library(limma)
library(RColorBrewer)
#library(stringr)
library(proxy)
library(irlba,lib.loc="/net/shendure/vol1/home/cusanovi/R/x86_64-unknown-linux-gnu-library/3.1/")
library(gplots)
library(data.table)
#library(dbscan)
library(Rtsne)
library(densityClust,lib.loc="/net/shendure/vol1/home/cusanovi/R/x86_64-unknown-linux-gnu-library/3.1/")
library(methods)
#library(monocle)

#Function to perform LSA
my_lsa <- function (x, dims = 2) 
{
  n <- nrow(x)
  p <- ncol(x)
  SVD = irlba(x, dims, dims)
  if (is.function(dims)) {
    dims = dims(SVD$d)
  }
  if (dims < 2) 
    dims = 2
  if (any(SVD$d <= sqrt(.Machine$double.eps))) {
    warning("[lsa] - there are singular values which are zero.")
  }
  space = NULL
  space$tk = SVD$u[, 1:dims]
  space$dk = SVD$v[, 1:dims]
  space$sk = SVD$d[1:dims]
  rownames(space$tk) = rownames(x)
  rownames(space$dk) = colnames(x)
  class(space) = "LSAspace"
  return(space)
}

gz = FALSE
numclust = 5
numcomponents = 6
cellcutoff=0.1
sitecutoff=20000
siteperccutoff=0.05
outheader = "/net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/AlphabetExperiments/sc_drone_dal/SCatac_dal.true.nodups.5kbwindows."
inmat = "/net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/fastq/scatac_hEB.q10.sort.true.nodups_5kb.dhsmatrix.txt"
outhists = paste0(outheader,cellcutoff,"cells.",siteperccutoff,"cellsites.",numcomponents,"components.",numclust,"clusters.binaryhists.pdf")
outheatmap = paste0(outheader,cellcutoff,"cells.",siteperccutoff,"cellsites.",numcomponents,"components.",numclust,"clusters.ward.heatmap.jpg")
sitesout = paste0(outheader,cellcutoff,"cells.",siteperccutoff,"cellsites.",numcomponents,"components.",numclust,"clusters.siteclusters.txt")
cellsout = paste0(outheader,cellcutoff,"cells.",siteperccutoff,"cellsites.",numcomponents,"components.",numclust,"clusters.cellclusters.txt")
tsneout = paste0(outheader,cellcutoff,"cells.",siteperccutoff,"cellsites.",numcomponents,"components.",numclust,"clusters.tsne.pdf")
compplotsout = paste0(outheader,cellcutoff,"cells.",siteperccutoff,"cellsites.",numcomponents,"components.",numclust,"clusters.componentplots.pdf")
promwindows = read.table("/net/shendure/vol10/projects/scATAC/nobackup/genomes/annotations/hg19/hg19.hs37d5.5kb.windows.within.2_5kb.of.tss.whitelist.bed")

if(gz){
  outprefix = strsplit(inmat,"[.]dhsmatrix[.]txt[.]gz")[[1]]
}else{
  outprefix = strsplit(inmat,"[.]dhsmatrix[.]txt")[[1]]
}
rdsfile = paste0(outprefix,".dhsmatrix.rds")
sparserdsfile = paste0(outprefix,".dhsmatrix.binary.sparse.rds")

#Read in matrix
if(file.exists(rdsfile)){
  bigmat = readRDS(rdsfile)
}else{
  if(gz == TRUE){
    bigmat = fread(paste0("zcat ",inmat),header=T)
  } else {
    bigmat = fread(inmat,header=T)
  }
  bigmat = as.data.frame(bigmat)
  saveRDS(bigmat,rdsfile)
}
print("Formatting data...")
sortmat = as.matrix(bigmat[,-c(1:4)])
#Make a list of sites
sites = strsplit2(bigmat[,4],"_")
label = paste0(sites[,1],":",sites[,2],"-",sites[,3])

if(file.exists(sparserdsfile)){
  bigmat.bin = readRDS(sparserdsfile)
  namelist = readRDS(paste0(outprefix,".dhsmatrix.binary.namelist.rds"))
}else{
  currnonzs = unlist(apply(sortmat,2,function(x) which(x > 0)))
  namers = gsub('[[:digit:]]+', '', names(currnonzs))
  namenums = c(1:length(unique(namers)))[unlist(as.factor(namers))]
  mastermat = cbind(currnonzs,namenums)
  namelist = unique(namers)
  bigmat.bin = sparseMatrix(i=mastermat[,1],j=mastermat[,2],x=rep(1,times=dim(mastermat)[1]))
  saveRDS(bigmat.bin,sparserdsfile)
  saveRDS(namelist,paste0(outprefix,".dhsmatrix.binary.namelist.rds"))
  rm(currnonzs)
  rm(mastermat)
  rm(namers)
  rm(namenums)
  gc()
}


num_cells_ncounted = rowSums(bigmat.bin)
num_sites_ncounted = colSums(bigmat.bin)

pdf(outhists)
hist(log10(num_cells_ncounted),main="Number of Cells Each Site is Observed In",breaks=50)
abline(v=log10(dim(bigmat.bin)[2]*siteperccutoff),lwd=2,col="indianred")
hist(log10(num_sites_ncounted),main="Number of Sites Each Cell Uses",breaks=50)
abline(v=log10(quantile(num_sites_ncounted,probs=cellcutoff)),lwd=2,col="indianred")
#dev.off()

annot.ncounts = label[num_cells_ncounted >= num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[sitecutoff]]]
#ncounts = bigmat.bin[num_cells_ncounted >= dim(bigmat.bin)[2]*siteperccutoff,]
ncounts = bigmat.bin[num_cells_ncounted >= num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[sitecutoff]],]

#
#nquants = sortmat.nonzero[num_cells_ncounted >= num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[sitecutoff]],]
#
new_counts = colSums(ncounts)
ncounts = ncounts[,new_counts >= quantile(new_counts,probs=cellcutoff)]
#
#nquants = nquants[,new_counts >= quantile(new_counts,probs=cellcutoff)]
#
hist(log10(new_counts),main="Number of Sites Each Cell Uses",breaks=50)
abline(v=log10(quantile(new_counts,probs=cellcutoff)),lwd=2,col="indianred")
dev.off()

annot.ncounts = annot.ncounts[rowSums(ncounts) > 0]
#
#nquants = nquants[rowSums(ncounts) > 0,]
#
ncounts = ncounts[rowSums(ncounts) > 0,]
ncounts = Matrix(ncounts,sparse=T)

nfreqs <- t(t(ncounts) / colSums(ncounts))

tf_idf_counts <- nfreqs * log(1 + ncol(ncounts) / rowSums(ncounts))
#nfreqs <- t(t(nquants) / colSums(nquants))
#tf_idf_counts <- nfreqs * log(1 + ncol(nquants) / rowSums(nquants))
set.seed(0)
LLL <- my_lsa(tf_idf_counts,numcomponents)

sk_diag <- matrix(0, nrow=length(LLL$sk), ncol=length(LLL$sk))
diag(sk_diag) <- LLL$sk
sk_diag[1,1] = 0
LLL_d <- t(sk_diag %*% t(LLL$dk))
hclust_cells <- hclust(proxy::dist(LLL_d, method="cosine"), method="ward.D2")
#hclust_cells <- hclust(proxy::dist(LLL_d), method="ward.D2")
TTT_d <- t(sk_diag %*% t(LLL$tk))
hclust_genes <- hclust(proxy::dist(TTT_d, method="cosine"), method="ward.D2")
#hclust_genes <- hclust(proxy::dist(TTT_d), method="ward.D2")

genes_tree_cut <- cutree(hclust_genes, numclust)
gene_pal <- brewer.pal(length(unique(genes_tree_cut)), "Paired")
cells_tree_cut <- cutree(hclust_cells, numclust)
cell_pal <- brewer.pal(length(unique(cells_tree_cut)), "Paired")

LSI_out <-  t(t(sk_diag %*% t(LLL$dk)) %*% t(LLL$tk))
LSI_out <- t(scale(t(LSI_out)))

scale_max <- 1.5
scale_min <- -1.5
LSI_out[LSI_out > scale_max] <- scale_max
LSI_out[LSI_out < scale_min] <- scale_min

hmcols <- colorpanel(100, "steelblue", "white", "tomato")
jpeg(outheatmap, width=4, height=5, units="in", res=600,type="cairo")
heatmap.2(LSI_out, 
          col=hmcols,
          ColSideColors=cell_pal[as.factor(cells_tree_cut)],
          RowSideColors=gene_pal[as.factor(genes_tree_cut)],
          Rowv = as.dendrogram(hclust_genes), Colv = as.dendrogram(hclust_cells),
          labRow=FALSE, labCol=FALSE, trace="none",  scale="none",
          useRaster=TRUE)
dev.off()

pdf(compplotsout)
par(mfrow=c(3,3))
plot(LLL$dk[,1],log10(colSums(ncounts)),pch=20,col="dodgerblue2",xlab="PC1",ylab="Read Depth")
for(i in seq(1,dim(LLL$dk)[2],2)){
    plot(LLL$dk[,i],LLL$dk[,i+1],pch=20,col="mediumseagreen",xlab=i,ylab=i+1)
}
dev.off()

set.seed(0)
tsnetfidf = Rtsne(LLL_d,pca=F)

pdf(tsneout)
plot(tsnetfidf$Y,pch=20,col=cell_pal[as.factor(cells_tree_cut)],main="Standard Components")

set.seed(0)
LLLtsne <- my_lsa(tf_idf_counts,50)
sk_diagtsne <- matrix(0, nrow=length(LLLtsne$sk), ncol=length(LLLtsne$sk))
diag(sk_diagtsne) <- LLLtsne$sk
LLL_dtsne <- t(sk_diagtsne %*% t(LLLtsne$dk))
set.seed(0)
tsnetfidf = Rtsne(LLL_dtsne,pca=F)
plot(tsnetfidf$Y,pch=20,col=cell_pal[as.factor(cells_tree_cut)],main="50 Components")
dev.off()

sitebeds = strsplit2(annot.ncounts,"-|:")
site_clusters = cbind(annot.ncounts,genes_tree_cut)
site_clusters = cbind(sitebeds,site_clusters)
colnames(site_clusters) = c("Chrom","Start","End","Site","Cluster")
cell_clusters = cbind(colnames(ncounts),cells_tree_cut)
colnames(cell_clusters) = c("Cell","Cluster")
write.table(site_clusters,sitesout,row.names=F,col.names=F,sep="\t",quote=F)
write.table(cell_clusters,cellsout,row.names=F,col.names=F,sep="\t",quote=F)

promlabel = paste0(promwindows[,1],":",promwindows[,2],"-",promwindows[,3])
prommatches = intersect(promlabel,label)

promind = match(prommatches,label)
prom.bin = bigmat.bin[promind,]
#
#prom.quant = sortmat.nonzero[promind,]
#

label.prom = label[promind]

num_cells_ncounted.prom = rowSums(prom.bin)
num_sites_ncounted.prom = colSums(prom.bin)

pdf(paste0(outhists,".prom.pdf"))
hist(log10(num_cells_ncounted.prom),main="Number of Cells Each Site is Observed In",breaks=50)
abline(v=log10(dim(prom.bin)[2]*siteperccutoff),lwd=2,col="indianred")
hist(log10(num_sites_ncounted.prom),main="Number of Sites Each Cell Uses",breaks=50)
abline(v=log10(quantile(num_sites_ncounted.prom,probs=cellcutoff)),lwd=2,col="indianred")
#dev.off()

# if(length(num_cells_ncounted.prom) > sitecutoff){
# 	annot.ncounts.prom = label.prom[num_cells_ncounted.prom >= dim(prom.bin)[2]*siteperccutoff]
# 	ncounts.prom = prom.bin[num_cells_ncounted.prom >= dim(prom.bin)[2]*siteperccutoff,]
# }else{
# 	annot.ncounts.prom = label.prom
# 	ncounts.prom = prom.bin
# }
annot.ncounts.prom = label.prom[num_cells_ncounted.prom >= dim(prom.bin)[2]*siteperccutoff]
ncounts.prom = prom.bin[num_cells_ncounted.prom >= dim(prom.bin)[2]*siteperccutoff,]

#
#nquants.prom = prom.bin[num_cells_ncounted.prom >= num_cells_ncounted.prom[order(num_cells_ncounted.prom,decreasing=T)[sitecutoff]],]
#
new_counts = colSums(ncounts.prom)
ncounts.prom = ncounts.prom[,new_counts >= quantile(new_counts,probs=cellcutoff)]
#
#nquants.prom = nquants.prom[,new_counts >= quantile(new_counts,probs=cellcutoff)]
#
hist(log10(new_counts),main="Number of Sites Each Cell Uses",breaks=50)
abline(v=log10(quantile(new_counts,probs=cellcutoff)),lwd=2,col="indianred")
dev.off()

annot.ncounts.prom = annot.ncounts.prom[rowSums(ncounts.prom) > 0]
#
#nquants.prom = nquants.prom[rowSums(ncounts.prom) > 0,]
#
ncounts.prom = ncounts.prom[rowSums(ncounts.prom) > 0,]
ncounts.prom = Matrix(ncounts.prom,sparse=T)

nfreqs.prom <- t(t(ncounts.prom) / colSums(ncounts.prom))
tf_idf_counts.prom <- nfreqs.prom * log(1 + ncol(ncounts.prom) / rowSums(ncounts.prom))
#
#nfreqs.prom <- t(t(nquants.prom) / colSums(nquants.prom))
#tf_idf_counts.prom <- nfreqs.prom * log(1 + ncol(nquants.prom) / rowSums(nquants.prom))
#
set.seed(0)
LLL.prom <- my_lsa(tf_idf_counts.prom,numcomponents)

sk_diag.prom <- matrix(0, nrow=length(LLL.prom$sk), ncol=length(LLL.prom$sk))
diag(sk_diag.prom) <- LLL.prom$sk
sk_diag.prom[1,1] = 0
LLL_d.prom <- t(sk_diag.prom %*% t(LLL.prom$dk))
hclust_cells.prom <- hclust(proxy::dist(LLL_d.prom, method="cosine"), method="ward.D2")
#hclust_cells.prom <- hclust(proxy::dist(LLL_d.prom), method="ward.D2")
TTT_d.prom <- t(sk_diag.prom %*% t(LLL.prom$tk))
hclust_genes.prom <- hclust(proxy::dist(TTT_d.prom, method="cosine"), method="ward.D2")
#hclust_genes.prom <- hclust(proxy::dist(TTT_d.prom), method="ward.D2")

genes_tree_cut.prom <- cutree(hclust_genes.prom, numclust)
gene_pal.prom <- brewer.pal(length(unique(genes_tree_cut.prom)), "Paired")
cells_tree_cut.prom <- cutree(hclust_cells.prom,numclust)
cell_pal.prom <- brewer.pal(length(unique(cells_tree_cut.prom)), "Paired")

LSI_out.prom <-  t(t(sk_diag.prom %*% t(LLL.prom$dk)) %*% t(LLL.prom$tk))
LSI_out.prom <- t(scale(t(LSI_out.prom)))

scale_max <- 1.5
scale_min <- -1.5
LSI_out.prom[LSI_out.prom > scale_max] <- scale_max
LSI_out.prom[LSI_out.prom < scale_min] <- scale_min

jpeg(paste0(outheatmap,".prom.jpg"), width=4, height=5, units="in", res=600,type="cairo")
heatmap.2(LSI_out.prom,
          col=hmcols,
          ColSideColors=cell_pal.prom[as.factor(cells_tree_cut.prom)],
          RowSideColors=gene_pal.prom[as.factor(genes_tree_cut.prom)],
          Rowv = as.dendrogram(hclust_genes.prom), Colv = as.dendrogram(hclust_cells.prom),
          labRow=FALSE, labCol=FALSE, trace="none",  scale="none",
          useRaster=TRUE)
dev.off()

pdf(paste0(compplotsout,".prom.pdf"))
par(mfrow=c(3,3))
plot(LLL.prom$dk[,1],log10(colSums(ncounts.prom)),pch=20,col="dodgerblue2",xlab="PC1",ylab="Read Depth")
for(i in seq(1,dim(LLL.prom$dk)[2],2)){
    plot(LLL.prom$dk[,i],LLL.prom$dk[,i+1],pch=20,col="mediumseagreen",xlab=i,ylab=i+1)
}
dev.off()

set.seed(0) 
tsnetfidf = Rtsne(LLL.prom$dk,pca=F)

pdf(paste0(tsneout,".prom.pdf"))
plot(tsnetfidf$Y,pch=20,col=cell_pal.prom[as.factor(cells_tree_cut.prom)],main="Standard Components")

set.seed(0)
LLLtsne <- my_lsa(tf_idf_counts.prom,50)
sk_diagtsne <- matrix(0, nrow=length(LLLtsne$sk), ncol=length(LLLtsne$sk))
diag(sk_diagtsne) <- LLLtsne$sk
LLL_dtsne <- t(sk_diagtsne %*% t(LLLtsne$dk))
set.seed(0)
tsnetfidf = Rtsne(LLL_dtsne,pca=F)
plot(tsnetfidf$Y,pch=20,col=cell_pal.prom[as.factor(cells_tree_cut.prom)],main="50 Components")
dev.off()

sitebeds.prom = strsplit2(annot.ncounts.prom,"-|:")
site_clusters.prom = cbind(annot.ncounts.prom,genes_tree_cut.prom)
site_clusters.prom = cbind(sitebeds.prom,site_clusters.prom)
colnames(site_clusters.prom) = c("Chrom","Start","End","Site","Cluster")
cell_clusters.prom = cbind(colnames(ncounts.prom),cells_tree_cut.prom)
colnames(cell_clusters.prom) = c("Cell","Cluster")
write.table(site_clusters.prom,paste0(sitesout,".prom.txt"),row.names=F,col.names=F,sep="\t",quote=F)
write.table(cell_clusters.prom,paste0(cellsout,".prom.txt"),row.names=F,col.names=F,sep="\t",quote=F)


distal.bin = bigmat.bin[-promind,]
#
#distal.quant = sortmat.nonzero[-promind,]
#
label.distal = label[-promind]

num_cells_ncounted.distal = rowSums(distal.bin)
num_sites_ncounted.distal = colSums(distal.bin)

pdf(paste0(outhists,".distal.pdf"))
hist(log10(num_cells_ncounted.distal),main="Number of Cells Each Site is Observed In",breaks=50)
abline(v=log10(dim(distal.bin)[2]*siteperccutoff),lwd=2,col="indianred")
hist(log10(num_sites_ncounted.distal),main="Number of Sites Each Cell Uses",breaks=50)
abline(v=log10(quantile(num_sites_ncounted.distal,probs=cellcutoff)),lwd=2,col="indianred")
#dev.off()

annot.ncounts.distal = label.distal[num_cells_ncounted.distal >= dim(distal.bin)[2]*siteperccutoff]
ncounts.distal = distal.bin[num_cells_ncounted.distal >= dim(distal.bin)[2]*siteperccutoff,]
#
#nquants.distal = distal.quant[num_cells_ncounted.distal >= num_cells_ncounted.distal[order(num_cells_ncounted.distal,decreasing=T)[sitecutoff]],]
#
new_counts = colSums(ncounts.distal)
ncounts.distal = ncounts.distal[,new_counts >= quantile(new_counts,probs=cellcutoff)]
#
#nquants.distal = nquants.distal[,new_counts >= quantile(new_counts,probs=cellcutoff)]
#
hist(log10(new_counts),main="Number of Sites Each Cell Uses",breaks=50)
abline(v=log10(quantile(new_counts,probs=cellcutoff)),lwd=2,col="indianred")
dev.off()

annot.ncounts.distal = annot.ncounts.distal[rowSums(ncounts.distal) > 0]
#
#nquants.distal = nquants.distal[rowSums(nquants.distal) > 0,]
#
ncounts.distal = ncounts.distal[rowSums(ncounts.distal) > 0,]
ncounts.distal = Matrix(ncounts.distal,sparse=T)

nfreqs.distal <- t(t(ncounts.distal) / colSums(ncounts.distal))
tf_idf_counts.distal <- nfreqs.distal * log(1 + ncol(ncounts.distal) / rowSums(ncounts.distal))
#
#nnot.ncounts.promnfreqs.distal <- t(t(nquants.distal) / colSums(nquants.distal))
#tf_idf_counts.distal <- nfreqs.distal * log(1 + ncol(nquants.distal) / rowSums(nquants.distal))
#
set.seed(0)
LLL.distal <- my_lsa(tf_idf_counts.distal,numcomponents)

sk_diag.distal <- matrix(0, nrow=length(LLL.distal$sk), ncol=length(LLL.distal$sk))
diag(sk_diag.distal) <- LLL.distal$sk
sk_diag.distal[1,1] = 0
LLL_d.distal <- t(sk_diag.distal %*% t(LLL.distal$dk))
hclust_cells.distal <- hclust(proxy::dist(LLL_d.distal, method="cosine"), method="ward.D2")
#hclust_cells.distal <- hclust(proxy::dist(LLL_d.distal), method="ward.D2")
TTT_d.distal <- t(sk_diag.distal %*% t(LLL.distal$tk))
hclust_genes.distal <- hclust(proxy::dist(TTT_d.distal, method="cosine"), method="ward.D2")
#hclust_genes.distal <- hclust(proxy::dist(TTT_d.distal), method="ward.D2")

genes_tree_cut.distal <- cutree(hclust_genes.distal, numclust)
gene_pal.distal <- brewer.pal(length(unique(genes_tree_cut.distal)), "Paired")
cells_tree_cut.distal <- cutree(hclust_cells.distal, numclust)
cell_pal.distal <- brewer.pal(length(unique(cells_tree_cut.distal)), "Paired")

LSI_out.distal <-  t(t(sk_diag.distal %*% t(LLL.distal$dk)) %*% t(LLL.distal$tk))
LSI_out.distal <- t(scale(t(LSI_out.distal)))

LSI_out.distal[LSI_out.distal > scale_max] <- scale_max
LSI_out.distal[LSI_out.distal < scale_min] <- scale_min

jpeg(paste0(outheatmap,".distal.jpg"), width=4, height=5, units="in", res=600,type="cairo")
heatmap.2(LSI_out.distal,
          col=hmcols,
          ColSideColors=cell_pal.distal[as.factor(cells_tree_cut.distal)],
          RowSideColors=gene_pal.distal[as.factor(genes_tree_cut.distal)],
          Rowv = as.dendrogram(hclust_genes.distal), Colv = as.dendrogram(hclust_cells.distal),
          labRow=FALSE, labCol=FALSE, trace="none",  scale="none",
          useRaster=TRUE)
dev.off()

pdf(paste0(compplotsout,".distal.pdf"))
par(mfrow=c(3,3))
plot(LLL.distal$dk[,1],log10(colSums(ncounts.distal)),pch=20,col="dodgerblue2",xlab="PC1",ylab="Read Depth")
for(i in seq(1,dim(LLL.distal$dk)[2],2)){
    plot(LLL.distal$dk[,i],LLL.distal$dk[,i+1],pch=20,col="mediumseagreen",xlab=i,ylab=i+1)
}
dev.off()

set.seed(0)
tsnetfidf = Rtsne(LLL.distal$dk,pca=F)

pdf(paste0(tsneout,".distal.pdf"))
plot(tsnetfidf$Y,pch=20,col=cell_pal.distal[as.factor(cells_tree_cut.distal)],main="Standard Components")

set.seed(0)
LLLtsne <- my_lsa(tf_idf_counts.distal,50)
sk_diagtsne <- matrix(0, nrow=length(LLLtsne$sk), ncol=length(LLLtsne$sk))
diag(sk_diagtsne) <- LLLtsne$sk
LLL_dtsne <- t(sk_diagtsne %*% t(LLLtsne$dk))
set.seed(0)
tsnetfidf = Rtsne(LLL_dtsne,pca=F)
plot(tsnetfidf$Y,pch=20,col=cell_pal.distal[as.factor(cells_tree_cut.distal)],main="50 Components")
dev.off()

sitebeds.distal = strsplit2(annot.ncounts.distal,"-|:")
site_clusters.distal = cbind(annot.ncounts.distal,genes_tree_cut.distal)
site_clusters.distal = cbind(sitebeds.distal,site_clusters.distal)
colnames(site_clusters.distal) = c("Chrom","Start","End","Site","Cluster")
cell_clusters.distal = cbind(colnames(ncounts.distal),cells_tree_cut.distal)
colnames(cell_clusters.distal) = c("Cell","Cluster")
write.table(site_clusters.distal,paste0(sitesout,".distal.txt"),row.names=F,col.names=F,sep="\t",quote=F)
write.table(cell_clusters.distal,paste0(cellsout,".distal.txt"),row.names=F,col.names=F,sep="\t",quote=F)




