args = commandArgs(TRUE)
report = read.table(args[1],header=T)
pdf(sub("txt","pdf",args[1]))
plot((report[,2]-report[,3])/report[,2],log10(report[,3]),xlim=c(0,1),xlab="Fraction PCR Duplicate",ylab="Log10(Unique Reads)",pch=20)
dev.off()

