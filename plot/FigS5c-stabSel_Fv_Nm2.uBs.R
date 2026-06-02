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

pdf("stabSel_Fv_Nm2.uBs.pdf", width=3, height=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)

plot(NULL, xlim=log10(c(0.2,250)), ylim=c(0, 0.5), ylab="", xlab="", xaxt="n",yaxt="n",
       main=substitute(italic("Nm") * " = " * m, list(m=round(Nm,3))))
posLab=c(0.2, 1, 5,  25, 250)
axis(1, at=log10(posLab), labels = posLab)
posLab=seq(0,0.5,0.1)
axis(2, at=posLab, labels = posLab)
  
mtext(expression(italic(v)*" (in units of "*italic("v")["*"]*")"), side=1, line = 1.6, cex=1.2)
mtext(expression("F("*italic(v)*")"), side=2, line=2, las=0, cex=1.2)
  
legend("topright", legend= c(substitute(italic(mu *"/s = ")*y*italic(", N"*mu*" = ")*x, list(x=Nu[1], y=Nu[1]/Ns[1])), 
                             substitute(italic(mu *"/s = ")*y*italic(", N"*mu*" = ")*x, list(x=Nu[2], y=Nu[2]/Ns[2])),
                             substitute(italic(mu *"/s = ")*y*italic(", N"*mu*" = ")*x, list(x=Nu[3], y=Nu[3]/Ns[3]))), 
       lty=rep("solid",3), col=colors, seg.len = 1, cex=0.8, bty="n")

for(i in 1:3){
  expData=read.table(paste(args[1], "/Fy_Nu", Nu[i], "_Nm", Nm, "_Ns", Ns[i], ".eff.dat", sep=""))
  lines(log10(expData[,1]), expData[,2], col=colors[i])
}

dev.off()