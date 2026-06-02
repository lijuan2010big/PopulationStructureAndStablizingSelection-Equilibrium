
library(RColorBrewer)

args=commandArgs(trailingOnly = TRUE)

Nm=c(0.1, 0.2, 0.5)
pmean=c(0.04523373134770202, 0.008386037921965258, 0.002644693586148858)
pmeanEff=c(0.49999998312228405,0.012077551463072865, 0.003250840240476799)
HTexp=2*pmean*(1-pmean)
HTexpEff=2*pmeanEff*(1-pmeanEff)
intfrq=0.5

colours=c("grey70", "grey40", "grey20")

lwidth=c(3, 1.5, 0.8)

pdf(paste("stabSel_undDomi.diffInit.time.het.pdf", sep=""), height=3, width=6)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1, xpd=TRUE)

plot(NULL, xlim=c(0,1000000), ylim=c(0,0.5), xlab="", ylab=expression(italic(H[T])), main="Single locus, IC1")
mtext("Time in generations", side=1, line=1.5)

legend("topleft", legend=Nm, text.col=colours, title = expression(italic("Nm")), title.col = "black", cex=0.8, bty="n", horiz=TRUE)
legend("topleft", inset=c(0,0.2), title="Approximation", legend=c("LE", "LD"), lty=c("11", "solid"), seg.len = 1, cex=0.8, bty="n", horiz = TRUE)


for(i in 1:3){
  m=Nm[i]
  #simulation
    x=read.table(paste(args[1], "/init0_frq_Nm", m, "_Ns5_d100_N100", sep=""), header = FALSE)
    l=length(x[1,])
    HetMeta=apply(x[,2:l], 1, mean)
    HT=2*HetMeta*(1-HetMeta)
    lines(x[,1], HT, col=colours[i], lwd=1)
}

for(i in 1:3){
  m=Nm[i]
  lines(c(0,1030000), rep(HTexp[i], 2), lty="11", col=colours[i])
  #lines(c(0,1030000), rep(HTexpEff[i], 2), lty="solid", col=colours[i])
  }

plot(NULL, xlim=c(0,1000000), ylim=c(0, 0.5), xlab="", ylab=expression(italic(H[T])), main="Single locus, IC2")
mtext("Time in generations", side=1, line=1.5)

for(i in 1:3){
  m=Nm[i]
  
  x=read.table(paste(args[1], "/init0.5_frq_Nm", m, "_Ns5_d100_N100", sep=""), header = FALSE)
  l=length(x[1,])
  HetMeta=apply(x[,2:l], 1, mean)
  HT=2*HetMeta*(1-HetMeta)
  lines(x[,1], HT, col=colours[i], lwd=1)
  
}



for(i in 1:3){
  m=Nm[i]
  lines(c(0,1030000), rep(HTexp[i], 2), lty="11", col=colours[i])
  #lines(c(0,1030000), rep(HTexpEff[i], 2), lty="solid", col=colours[i])
  }


plot(NULL, xlim=c(0,1000000), ylim=c(0, 0.5), xlab="", ylab=expression(italic(H[T])), main="Polygenic, IC1")
mtext("Time in generations", side=1, line=1.5)


for(i in 1:3){
  m=Nm[i]
  #Simulation
  x=read.table(paste(args[1], "/het_trait_frqB0_Nu0.01_Nm", m, "_Ns5_N200_D100_L200", sep=""), header = TRUE)
  lines(x[,1], x[,2], col=colours[i], lwd=1)

}

for(i in 1:3){
  m=Nm[i]
  lines(c(0,1030000), rep(HTexp[i], 2), lty="11", col=colours[i])
  lines(c(0,1030000), rep(HTexpEff[i], 2), lty="solid", col=colours[i])
  }

plot(NULL, xlim=c(0,1000000), ylim=c(0, 0.5), xlab="", ylab=expression(italic(H[T])), main="Polygenic, IC2")
mtext("Time in generations", side=1, line=1.5)

for(i in 1:3){
  m=Nm[i]
  x=read.table(paste(args[1], "/het_trait_frqB1_Nu0.01_Nm", m, "_Ns5_N200_D100_L200", sep=""), header = TRUE)
  lines(x[,1], x[,2], col=colours[i], lwd=1)
}

for(i in 1:3){
  m=Nm[i]
  lines(c(0,1030000), rep(HTexp[i], 2), lty="11", col=colours[i])
  lines(c(0,1030000), rep(HTexpEff[i], 2), lty="solid", col=colours[i])
}



dev.off()
