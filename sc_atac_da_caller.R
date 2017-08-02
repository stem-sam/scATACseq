#library(Matrix)
library(limma)
#library(stringr)
#library(proxy)
library(data.table)
library(L1Graph,lib.loc="/net/shendure/vol1/home/cusanovi/R/x86_64-unknown-linux-gnu-library/3.2/")
library(simplePPT,lib.loc="/net/shendure/vol1/home/cusanovi/R/x86_64-unknown-linux-gnu-library/3.2/")
library(DDRTree,lib.loc="/net/shendure/vol1/home/cusanovi/R/x86_64-unknown-linux-gnu-library/3.2/")
library(monocle,lib.loc="/net/shendure/vol1/home/cusanovi/R/x86_64-unknown-linux-gnu-library/3.2/")
source("/net/shendure/vol10/projects/mouse_development/scripts/monocle_patch.R")

gz = FALSE
GREAT = TRUE
inmat = "/net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/matrix/scatac_hEB.q10.sort.true.nodups_peaks.dhsmatrix.txt"
print(inmat)
#cellassignments = "/net/shendure/vol10/projects/mouse_atlas/nobackup/tissues/Cerebellum_62216/SC_atac_combined.Cerebellum_62216_P1.true.nodups.5kbwindows.cellclusters.5percnew.txt"
cellassignments = "/net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/tSNE/scatac.peaks.densitypeaks.cellclusters.rho50.delta2.5_merged_clusters_7.txt"
#readreport = "/net/shendure/vol10/projects/mouse_atlas/nobackup/tissues/Cerebellum_62216/SC_atac_combined.Cerebellum_62216_P1.true.report.txt"
readreport = "/net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/fastq/scatac_hEB.q10.sort.true.nodups.bam_H9_DHS.txt"
tissue = "Human_EB"
#outprefix = strsplit(inmat,"[.]dhsmatrix[.]txt")[[1]]
outprefix = "/net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/matrix/scatac_hEB.q10.sort.true.nodups_peaks_merged_clusters"

print("Loading data...")
if(gz == TRUE){
        bigmat = fread(paste0("zcat ",inmat),header=T)
} else {
        bigmat = fread(inmat,header=T)
}
bigmat = as.data.frame(bigmat)

mastersites = bigmat[,4]
sortmat = as.matrix(bigmat[,-c(1:4)])
nonzeroes = rowSums(sortmat) > 0

sortmat.nonzero = sortmat[nonzeroes,]
mastersites.nonzero = mastersites[nonzeroes]
#save(sortmat.nonzero, file=outobj)

cells = read.table(cellassignments)
goodcells = match(cells[,1],colnames(sortmat.nonzero))
sortmat.nonzero.assigned = sortmat.nonzero[,goodcells]

readcounts = read.table(readreport,header=T)
readcounts.assigned = readcounts[match(cells[,1],rownames(readcounts)),]

print("Reformatting data...")
currnonzs = unlist(apply(sortmat.nonzero.assigned,2,function(x) which(x > 0)))
namers = gsub('[[:digit:]]+', '', names(currnonzs))
namenums = c(1:length(unique(namers)))[unlist(as.factor(namers))]
mastermat = cbind(currnonzs,namenums)
namelist = unique(namers)

bigmat.bin = Matrix::sparseMatrix(i=mastermat[,1],j=mastermat[,2],x=rep(1,times=dim(mastermat)[1]))
rownames(bigmat.bin) = mastersites.nonzero
colnames(bigmat.bin) = namelist
sites = strsplit2(mastersites.nonzero,"_")
label = paste0(sites[,1],":",sites[,2],"-",sites[,3])

rm(bigmat)
rm(sortmat)
rm(sortmat.nonzero)
rm(sortmat.nonzero.assigned)
rm(sites)
rm(readcounts)
rm(nonzeroes)
gc()


num_cells_ncounted = rowSums(bigmat.bin)
#num_sites_ncounted = colSums(bigmat.bin)

bigmat.bin.final = bigmat.bin[num_cells_ncounted > 9,]
label.final = label[num_cells_ncounted > 9]

rm(bigmat.bin)
rm(num_cells_ncounted)
rm(label)
gc()

print("Building Monocle object...")
fda <- as.data.frame(as.character(label.final))
fda <- new("AnnotatedDataFrame", data = fda)

clusterfactor = rep(0,dim(cells)[1])
clusterind = grep("1",cells[,2])
clusterfactor[clusterind] = 1
clusterfactor = factor(clusterfactor,levels=c("0","1"))
pda = as.data.frame(cbind(as.character(cells[,1]),as.numeric(log10(readcounts.assigned[,2])),clusterfactor))
names(pda) <- c("CellID", "ReadDepth","CellCluster")
pda <- new("AnnotatedDataFrame", data = pda)

rownames(bigmat.bin.final) = NULL
colnames(bigmat.bin.final) = NULL
submat_cds <-  newCellDataSet(as(bigmat.bin.final, "sparseMatrix"),
        featureData = fda,
        phenoData = pda,
        expressionFamily=binomialff(),
        lowerDetectionLimit=1)

pData(submat_cds)$CellID = as.character(cells[,1])
pData(submat_cds)$ReadDepth = as.numeric(log10(readcounts.assigned[,2]))
pData(submat_cds)$CellCluster = factor(clusterfactor,levels=c("0","1"))
pData(submat_cds)$Size_Factor = 1

rm(readcounts.assigned)
rm(label.final)
rm(fda)
rm(pda)
rm(bigmat.bin.final)
gc()

print("Running models...")
for(i in 1:length(table(cells[,2]))){
    print(paste0("Analyzing cluster ",i,"..."))
    clusterfactor = rep(0,dim(cells)[1])
    clusterind = which(cells[,2] == i)
    clusterfactor[clusterind] = 1
    clusterfactor = factor(clusterfactor,levels=c("0","1"))
    print(table(clusterfactor))
    pData(submat_cds)$CellCluster = factor(clusterfactor,levels=c("0","1"))
    differtest <- differentialGeneTest(submat_cds, fullModelFormulaStr = "~CellCluster + ReadDepth", reducedModelFormulaStr = "~ReadDepth", cores=15)
#    core = 15 original, and sumat_cds
    print(head(differtest))
#    a = rowSums(bigmat.bin.final[,pData(submat_cds)[,3] == 1])
#    b = sum(pData(submat_cds)[,3] == 1) - a
#    c = rowSums(bigmat.bin.final[,pData(submat_cds)[,3] == 0])
#    d = sum(pData(submat_cds)[,3] == 0) - c
#    logits = log((a/c)/(b/d))
#    logitsalt = logits
#    maxer = max(abs(logits[logits < Inf & logits > -Inf]))
#    logitsalt[logits == Inf] = maxer + 1
#    logitsalt[logits == -Inf] = -maxer - 1
#    newind = as.numeric(rownames(differtest))
#    tots = a + c
#    differtest$logits = logits[newind]
#    differtest$logitsalt = logitsalt[newind]
#    differtest$tots = tots[newind]
    sigcol = rep("black",times = length(differtest$qval))
    sigcol[differtest$qval < 0.01] = "indianred"
    png(paste0(outprefix,".volcano.cluster",i,".png"), width=4, height=5, units="in", res=600,type="cairo")
    plot(differtest$beta,-log10(differtest$pval),col=sigcol,xlab="Cluster Coeff",ylab="-Log10(P-value)",pch=20)
    dev.off()
    pdf(paste0(outprefix,".volcano.cluster",i,".pdf"))
    plot(differtest$beta,-log10(differtest$pval),col=sigcol,xlab="Cluster Coeff",ylab="-Log10(P-value)",pch=20)
    dev.off()
    write.table(differtest,paste0(outprefix,".datests.cluster",i,".txt"),row.names=F,quote=F,sep="\t")
    sigup = strsplit2(differtest[which(differtest$qval < 0.01 & differtest$beta > 0),6],":|-")
    sigup = cbind(sigup,signif(differtest[which(differtest$qval < 0.01 & differtest$beta > 0),5],6))
    sigdown = strsplit2(differtest[which(differtest$qval < 0.01 & differtest$beta < 0),6],":|-")
    sigdown = cbind(sigdown,signif(differtest[which(differtest$qval < 0.01 & differtest$beta < 0),5],6))
    write.table(sigup[order(as.numeric(sigup[,4])),],paste0(outprefix,".sigopen.cluster",i,".txt"),row.names=F,col.names=F,sep="\t",quote=F)
    write.table(sigdown[order(as.numeric(sigdown[,4])),],paste0(outprefix,".sigclosed.cluster",i,".txt"),row.names=F,col.names=F,sep="\t",quote=F)
    if(GREAT == TRUE){
        write.table(sigup[,1:3],paste0("/net/shendure/vol10/www/content/members/sgregala/human_EB/great_bed/",tissue,".bed"),row.names=F,col.names=F,sep="\t",quote=F)
        system(paste0("wget -O ",outprefix,"sigopen.cluster",i,".great_enrichments.txt \"http://bejerano.stanford.edu/great/public/cgi-bin/greatStart.php?outputType=batch&requestSpecies=hg19&requestName=CurrentResults&requestSender=Client+A&requestURL=http%3A%2F%2Fkrishna.gs.washington.edu%2Fcontent%2Fmembers%2Fcusanovich%2Fmouse_atlas%2Fgreat_bed%2F",tissue,".bed\""))
    }
    rm(differtest)
    rm(sigup)
    rm(sigdown)
    rm(sigcol)
    gc()
}





