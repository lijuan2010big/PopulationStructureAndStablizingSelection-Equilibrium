library(RColorBrewer)
library(fields)
library(plotrix)

args=commandArgs(trailingOnly = TRUE)

L=200
D=100
N=200
Nm=2

Nu=c(0.001, 0.01, 0.001)
Ns=c(0.1, 1, 1)

colors=brewer.pal(n=9,name = "Set1")[2:4]

pdf(paste(args[1],"/stabSel_propVg_Nm2.uBs.pdf", sep=""), width=3, height=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)

  plot(NULL, xlim=log10(c(0.2, 250)), ylim=c(0,1), 
       ylab="", , xlab="", xaxt="n", yaxt="n",
       main=substitute(italic("Nm")*"="*x, list(x=Nm)))
  posLab=c(0.2, 1, 5,  25, 250)
  axis(1, at=log10(posLab), labels = posLab, cex.axis=1)
  mtext(expression(italic(v)), side=1, line = 1.4, cex=1.2)
  posLab=seq(0,1,0.2)
  axis(2, at=posLab, labels = posLab)
  mtext(expression("Proportion of total "*italic(V[W]^"(gen)")), side=2, line=1.8, las=0, cex=1.2)
  
for(i in 1:3){
  data=read.table(paste(args[1], "/icdfVg_Nu",Nu[i],"_Nm", Nm, "_Ns",Ns[i],".eff.dat", sep=""))
  lines(log10(data[,1]), data[,2]/data[1,2], col=colors[i])
}



dev.off()