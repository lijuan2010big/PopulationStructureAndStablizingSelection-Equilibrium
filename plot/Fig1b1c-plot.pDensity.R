
library(RColorBrewer)


args=commandArgs(trailingOnly = TRUE)

Ns=c(0.1, 1, 5)

critNm=c(0.624, 0.089, 0.046)
Nm=c(2, 0,0.022)


colors=brewer.pal(n = 8, name = "Dark2")


pdf("undDomi_alleleFrqDist.pdf", width=3, height=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)

Ns=0.1
Nm=c(0.15, 0.53, 0.624, 1)
linety=c("solid", "solid","11","solid")
col1=c(colors[2], colors[1], "black", colors[3])

plot(NULL, xlim=c(0,1), ylim=c(0, 3), xlab=expression("Allele frequency, "*italic(p)), ylab="Probability density", 
     main=substitute(italic("Ns") * " = " * Ns, list(Ns=Ns)), cex.axis=0.8,  cex.lab=1.2)

legend(0.45, 3, legend=c(Nm[1],Nm[2], paste(Nm[3], "(crit)", sep=""), Nm[4]), lty=linety, 
       col=col1, seg.len = 1 , title = expression(italic("Nm")), cex=0.8, bty="n")

for(m in 1:length(Nm)){
  data=read.table(paste(args[1], "/pDens_Nu0.01_Nm", Nm[m], "_Ns", Ns, ".dat", sep=""))
  lines(data[,1], data[,2], col=col1[m], lty=linety[m])
}

Ns=1
Nm=c(0.05, 0.089, 0.15, 0.5)
linety=c("solid", "11", "solid","solid")
col1=c(colors[2], "black", colors[1], colors[3])

plot(NULL, xlim=c(0,1), ylim=c(0, 3), xlab=expression("Allele frequency, "*italic(p)), ylab="Probability density", 
     main=substitute(italic("Ns") * " = " * Ns, list(Ns=Ns) ), cex.axis=0.8,  cex.lab=1.2)


legend(0.45, 3, legend=c(Nm[1], paste(Nm[2], "(crit)", sep=""),Nm[3], Nm[4]), lty=linety, 
       col=col1, seg.len = 1 , title = expression(italic("Nm")), cex=0.8, bty="n")


for(m in 1:length(Nm)){
  data=read.table(paste(args[1], "/pDens_Nu0.01_Nm", Nm[m], "_Ns", Ns, ".dat", sep=""))
  lines(data[,1], data[,2], col=col1[m], lty=linety[m])
}


dev.off()
