library(RColorBrewer)
library(plotrix)
options(digits=10)

args=commandArgs(trailingOnly = TRUE)

colors=rev(brewer.pal(n = 11, name = "RdBu")[c(2, 4,9,11)])
#colors1=col2rgb(colors)/255

Nm=c(0, 0.1, 0.5, 2)
Vs=1
Nu=0.01
N=200
L=200
D=100

tmp=matrix(,nrow=length(Nm), ncol=L)
VT=tmp
VW=tmp
effectSize=tmp
FST=tmp

for(i in 1:length(Nm)){
  HT=matrix(,ncol=L, nrow=10)
  HW=matrix(,ncol=L, nrow=10)
  tmpFST=matrix(,ncol=L, nrow=10)
  for(r in 1:10){
    data=read.table(paste(args[1], "/simu/Nu0.01_Vs1_aExpNuL_Nm", Nm[i], "_L200_N200_D100_r", r, ".freq", sep=""))
    effectSize[i,]=as.matrix(data[1,2:length(data[1,])])
    for(e in 1:L){
        Ns=N*effectSize[i, e]^2/(2*Vs) 
        tmpHT=vector()
        tmpHW=vector()
        for(gen in seq(max(data[,1])-1e4, max(data[,1]), 2e3)){
          freq=data[data[,1]==gen, e+1]
          pmean=mean(freq)
          tmpHT=c(tmpHT, 2*pmean*(1-pmean))
          tmpHW=c(tmpHW, mean(as.matrix(2*freq*(1-freq))))
        }
        HT[r,e]=mean(tmpHT)
        HW[r,e]=mean(tmpHW)
        tmpFST[r,e]=1-HW[r,e]/HT[r,e]
    }
  }
  
  for(e in 1:L){
    Ns=N*effectSize[i, e]^2/(2*Vs) 
    VT[i, e]=mean(HT[,e])*Ns/(2*Nu)
    VW[i, e]=mean(HW[,e])*Ns/(2*Nu) 
    FST[i, e]=mean(tmpFST[,e])
  }
}


pdf(paste(args[1], "/stabSel_mixEffSizeExp.eff.pdf", sep=""), height=3, width=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)


### Genic variance in Full population
plot(NULL, xlim=log10(c(0.01, 10)), ylim=c(log10(0.1),log10(500)), 
     , xlab="", ylab="", main="Full population",
     type="l", lty=1, xaxt="n", yaxt="n") 
### x axis
mtext(expression(italic(Ns)), side=1, line = 1.5, cex=1.2)
axis(side=1, labels=c(0.01, 0.1, 1, 10), at=log10(c(0.01, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.01, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
### y axis
axis(side=2, labels=c( 0.1, 1, 10, 100, 500), at=log10(c(0.1, 1, 10, 100, 500)))
axis(side=2, at=log10(c(seq(0.1, 1, 0.1), seq(1,10,1), seq(10,100,10), seq(100,500,100))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
mtext(expression("Scaled genic variance, "* italic(v[T]^"(gen)"/v["*"])), side=2, line=1.9, las=0, cex=0.9)

######
for(i in 1:length(Nm)){
  for(e in 1:L){
    Ns=N*effectSize[i,e]^2/(2*Vs)
    points(log10(Ns), log10(VT[i, e]), pch=4, cex=0.4, col=colors[i], lwd=0.4)
  }
}

for(i in 1:length(Nm)){
  data=read.table(paste(args[1], "/predict/Nu0.01_Nm", Nm[i], "_Ns_pmean_HW.dat", sep=""))
  lines(log10(data[,1]), log10(data[,1]/Nu*(1-data[,2])*data[,2]), ,lty="11", col=colors[i])
  
  data=read.table(paste(args[1], "/predict/Nu0.01_uL0.01_Nm", Nm[i], "_mixEffSize_Ns_NmEff_NsEff_pmean_HW.i5.dat", sep=""))
  lines(log10(data[,1]), log10(data[,1]/Nu*data[,4]*(1-data[,4])), col=colors[i], lwd=1)
}

lines(c(log10(0.004),log10(12)), rep(0,2), lwd=0.5, col="grey90")
lines(c(0,0), c(log10(0.004),log10(600)), lwd=0.5, col="grey90")

#### Isolated pop.

data=read.table(paste(args[1], "/predict/Nu0.01_uL0.01_NmInfinity_mixEffSize_Ns_het.dat", sep=""))
lines(log10(data[,1]), log10(data[,1]/Nu*data[,2]/2), col="black", lty="11",lwd=1.5)

legend("topleft", inset=c(0.0,0), legend=c("LE", "LD"), cex=0.7, lty=c("11", "solid"), seg.len = 1, title="Approximation", bty="n", title.col = "black", horiz = TRUE)
legend("topleft", inset=c(0.0,0.2), legend=c(Nm, expression(infinity)), cex=0.7, text.col=c(colors, "black"), title=expression(italic("Nm")), title.col = "black", title.adj = 0.6, bty="n")


### Genic variance within demes
plot(NULL, xlim=log10(c(0.01, 10)), ylim=c(log10(0.01),log10(10)), 
     , xlab="", ylab="", main="Within demes",
     type="l", lty=1, xaxt="n", yaxt="n")
### x axis
mtext(expression(italic(Ns)), side=1, line = 1.5, cex=1.2)
axis(side=1, labels=c(0.01, 0.1, 1, 10), at=log10(c(0.01, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.01, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
### y axis
axis(side=2, labels=c(0.01, 0.1, 1, 10), at=log10(c(0.01, 0.1, 1, 10)))
axis(side=2, at=log10(c(seq(0.01, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
mtext(expression("Scaled genic variance, "* italic(v[W]^"(gen)"/v["*"])), side=2, line=1.9, las=0, cex=0.9)

##### 
for(i in 1:length(Nm)){
  for(e in 1:L){
    Ns=N*effectSize[i,e]^2/(2*Vs)
    points(log10(Ns), log10(VW[i,e]), pch=4, cex=0.4, col=colors[i], lwd=0.4)
  }
}

for(i in 1:length(Nm)){
  data=read.table(paste(args[1], "/predict/Nu0.01_Nm", Nm[i], "_Ns_pmean_HW.dat", sep=""))
  lines(log10(data[,1]), log10(data[,1]/Nu*data[,3]/2), col=colors[i], lty="11")
  
  data=read.table(paste(args[1], "/predict/Nu0.01_uL0.01_Nm", Nm[i], "_mixEffSize_Ns_NmEff_NsEff_pmean_HW.i5.dat", sep=""))
  lines(log10(data[,1]), log10(data[,1]/Nu*data[,5]/2),  col=colors[i], lwd=1)
}

lines(c(log10(0.004),log10(12)), rep(0,2), lwd=0.5, col="grey90")
lines(c(0,0), c(log10(0.004),log10(12)), lwd=0.5, col="grey90")

#### Isolated pop.

data=read.table(paste(args[1], "/predict/Nu0.01_uL0.01_NmInfinity_mixEffSize_Ns_het.dat", sep=""))
lines(log10(data[,1]), log10(data[,1]/Nu*data[,2]/2), col="black", lty="11",lwd=1.5)

### FST
plot(NULL, xlim=log10(c(0.01, 10)), ylim=c(0,1), 
     , xlab="", ylab="",
     type="l", lty=1, xaxt="n")
### x axis
mtext(expression(italic(Ns)), side=1, line = 1.5, cex=1.2)
axis(side=1, labels=c(0.01, 0.1, 1, 10), at=log10(c(0.01, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.01, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
### y axis
mtext(expression(italic(F[ST])), side=2, line=1.9, las=0, cex=0.9)

for(i in 1:length(Nm)){
  for(e in 1:L){
    Ns=N*effectSize[i,e]^2/(2*Vs)
    points(log10(Ns), FST[i,e], pch=4, cex=0.4, col=colors[i], lwd=0.4)
  }
}

for(i in 1:length(Nm)){
  data=read.table(paste(args[1], "/predict/Nu0.01_Nm", Nm[i], "_Ns_pmean_HW.dat", sep=""))
  lines(log10(data[,1]), 1-data[,3]/(2*data[,2]*(1-data[,2])),lty="11", col=colors[i])

  data=read.table(paste(args[1], "/predict/Nu0.01_uL0.01_Nm", Nm[i], "_mixEffSize_Ns_NmEff_NsEff_pmean_HW.i5.dat", sep=""))
  lines(log10(data[,1]), 1-data[,5]/(2*data[,4]*(1-data[,4])),  col=colors[i], lwd=1)

}


dev.off()
