library(RColorBrewer)
library(plotrix)
options(digits=10)

args=commandArgs(trailingOnly = TRUE)
colors=brewer.pal(n=9,name = "Set1")[2:4]
Ns=c(0.1, 1, 1)
Nu=c(0.001, 0.01, 0.001)

pdf(paste("stabSel.uBs.pdf", sep=""), height=3, width=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)

lNm=0.005

#### genic VT as a function of Nm.
plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(0, log10(27)), xlab="", ylab="", yaxt="n", xaxt="n",
     cex.axis=0.8, main="Full population")
### set up y axis
posLab=c(1, 2, 5, 10, 25)
axis(side=2, labels=posLab, at=log10(posLab))
mtext(expression("Scaled genic variance, "*italic(V[T]^"(gen)"/V["*"])), side=2, line=1.6, las=0, cex=0.8)
### set up x axis
axis(side=1, labels=c(0, 0.02, 0.1, 1, 10), at=log10(c(0.005, 0.02, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.03, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
mtext(expression(italic(Nm)), side=1, line = 1.6, cex=1.2)
### set up break x axis
axis.break(axis=1, breakpos=log10(0.012)) 

abline(h=0, lwd=0.7, col="grey50")

for(i in 1:2){
  #### Analytical results
  seff=read.table(paste(args[1], "/Nu", Nu[i],"_Ns", Ns[i], "_Nm_NmEff_NsEff_pmean_HW.uBs.dat", sep=""))
  VT=Ns[i]/Nu[i]*seff[,4]*(1-seff[,4])
  lines(log10(seff[seff[,1]>=0.015,1]), log10(VT[seff[,1]>=0.015]),col=colors[i],lwd=1 )    
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(VT[1:8]), col=colors[i], lwd=1 )
}


legend("topright", legend=expression(italic(mu*"/s")*" = 0.01"), cex=0.8, bty="n")

legend("topright", inset=c(0.05, 0.15), legend=Nu[1:2], text.col=c(colors), pch=c(NA, NA), title = expression( italic("N"*mu) ), cex=0.8, title.col="black", title.adj = 0.7, bty="n")

### genic VW as a function of Nm.
plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(log10(0.9),log10(10)), xlab="", ylab="", yaxt="n", xaxt="n",
     cex.axis=0.8, main="Within demes")
posLab=c(0.1, 0.5, 1, 2, 5, 10)
### set up y axis
axis(side=2, labels=posLab, at=log10(posLab))
mtext(expression("Scaled genic variance, "*italic(V[W]^"(gen)"/V["*"])), side=2, line=1.6, las=0, cex=0.8)
### set up x axis
axis(side=1, labels=c(0, 0.02, 0.1, 1, 10), at=log10(c(0.005, 0.02, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.03, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
mtext(expression(italic(Nm)), side=1, line = 1.6, cex=1.2)
### set up break x axis
axis.break(axis=1, breakpos=log10(0.012)) 

abline(h=0, lwd=0.7, col="grey50")

for(i in 1:2){
  seff=read.table(paste(args[1], "/Nu", Nu[i],"_Ns", Ns[i], "_Nm_NmEff_NsEff_pmean_HW.uBs.dat", sep=""))
  VW=Ns[i]/(Nu[i]*2)*seff[,5]
  lines(log10(seff[seff[,1]>=0.015,1]), log10(VW[seff[,1]>=0.015]),col=colors[i], lwd=1)    
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(VW[1:8]), col=colors[i], lwd=1.5 )
}

if(FALSE){
#### FST
plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(0,1), xlab="", ylab="", cex.axis=0.8, xaxt="n")
### set up y axis
mtext(expression(italic("F"["ST"])), side=2, line=1.8, las=0, cex=1.2)
### set up x axis
axis(side=1, labels=c(0, 0.02, 0.1, 1, 10), at=log10(c(0.005, 0.02, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.03, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
mtext(expression(italic(Nm)), side=1, line = 1.6, cex=1.2)
### set up break x axis
axis.break(axis=1, breakpos=log10(0.012)) 

for(i in 1:2){
  seff=read.table(paste(args[1], "/Nu", Nu[i],"_Ns", Ns[i], "_Nm_NmEff_NsEff_pmean_HW.uBs.dat", sep=""))
  HW=seff[,5]
  HT=2*seff[,4]*(1-seff[,4])
  FST=(HT-HW)/HT
  lines(log10(seff[seff[,1]>=0.015,1]), FST[seff[,1]>=0.015],  col=colors[i], lwd=1)    
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), FST[1:8], col=colors[i],, lwd=1.5)
}

Nu=0.01
Nm1=seq(0.015,15,0.01)
lines(log10(Nm1), 1/(4*Nm1+1+8*Nu),  col="black", lty="21")
Nm1=seq(0,0.007,0.001)
lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), 1/(4*Nm1+1+8*Nu),  col="black", lty="21")
}
dev.off()
