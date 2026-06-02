library(RColorBrewer)
library(plotrix)

args=commandArgs(trailingOnly = TRUE)

colors=brewer.pal(n=9,name = "Set1")[2:4]

Ns=c(0.1, 1, 5)
Nu=0.01
numLoc=seq(100, 1000,100)

pdf(paste(args[1],"stabSel_NuL.pdf", sep=""), height=3, width=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)

plot(NULL, xlim=c(0,10), ylim=c(1, 1.5), xlab="", ylab="", main=expression(italic(N)*mu*" = "*italic(L)*mu*" = 0.01, "*italic("Ns")* " = 5"))
mtext(expression(italic("N")*mu*italic("L")), side=1, line = 1.5, cex=1.2)
mtext("Scaled genic variance", side=2, line=2.2, las=0, cex=1.2)

  
for(j in 3:3){
  
  dat=read.table(paste(args[2], "/HetDeme_D1_Nu0.01_uL0.01_Ns.dat", sep=""))
  tmp=dat[dat[,1]==Ns[j],]
  lines(c(0,11),rep(tmp[1,2]/2*Ns[j]/Nu,2), lty="11")
  datEff=read.table(paste(args[2], "/HetDeme_D1_Nu0.01_uL0.01_Ns.Eff.dat", sep=""))
  tmp=datEff[datEff[,1]==Ns[j],]
  lines(c(0,11), rep(tmp[1,2]/2*Ns[j]/Nu, 2), lty="solid")

  for(i in 1:10){
    L=numLoc[i]
    HT=vector()
    numPop=L
  for(r in seq(1, 20, 1)){
      data=read.table(paste(args[1], "/het_trait_Nu0.01_Ns", Ns[j], "_N", numPop, "_D1_L", L, "_r", r, sep=""), header = TRUE)
      l=length(data[,1])
      HT=c(HT, mean(data[seq((l-300),l,4), 2]))
  }
  
  VT=Ns[j]/Nu*HT/2
  points(Nu*L, mean(VT), cex=0.7, pch=19)
  lines(rep(Nu*L, 2),c(mean(VT)-1.96*std.error(VT),  mean(VT)+1.96*std.error(VT)))
}
}

dev.off()