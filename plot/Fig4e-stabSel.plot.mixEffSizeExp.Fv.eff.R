library(RColorBrewer)
library(plotrix)
options(digits=10)

args=commandArgs(trailingOnly = TRUE)

colors=rev(brewer.pal(n = 11, name = "RdBu")[c(2, 3, 4,9,11)])

Nm=c(0, 0.1, 0.5, 2,10)
Vs=1
Nu=0.01
N=200
L=200
D=100
threshold=2

meanNs=2
delta=0.01

pdf(paste(args[1], "/stabSel_mixEffSizeExp.fv.eff.pdf", sep=""), height=3, width=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)


### Genic variance in Full population
plot(NULL, xlim=log10(c(0.01, 10)), ylim=c(1, 80), 
     , xlab="", ylab="", main="",
     type="l", lty=1, xaxt="n") 
### x axis
mtext(expression(italic(Ns)), side=1, line = 1.5, cex=1.2)
axis(side=1, labels=c(0.01, 0.1, 1, 10), at=log10(c(0.01, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.01, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
### y axis
mtext(expression("Number of demes with "*italic(tilde(v)>2)), side=2, line=1.9, las=0, cex=0.9)


legend("topleft", inset=c(0.0,0), legend=c("LE", "LD"), cex=0.7, lty=c("11", "solid"), seg.len = 1, title="Approximation", bty="n", title.col = "black", horiz = TRUE)
legend("topleft", inset=c(0.0,0.2), legend=c(Nm), cex=0.7, text.col=c(colors), title=expression(italic("Nm")), title.col = "black", title.adj = 0.6, bty="n")


for(m in 1:length(Nm)){
  dataExp=read.table(paste(args[1], "/predict/Fy_Nu0.01_Nm", Nm[m], "_Ns_Fv_Prob_mixExp2.v2.dat", sep=""))
  lines(log10(dataExp[,1]), dataExp[,6]*D, col=colors[m], lwd=1)
  lines(log10(dataExp[,1]), dataExp[,3]*D, col=colors[m], lty="11")
  
  freqDeme=data.frame()  
  
  for (r in 1:10) {
    
  data=read.table(paste(args[1], "/simu/Nu0.01_Vs1_aExpNuL_Nm", Nm[m], "_L200_N200_D100_r",r, ".freq", sep=""))
  
  effectSize=data[1, 2:length(data[1,])]

  Ns=as.numeric(N*effectSize^2/(2*Vs))
    
  for(gen in seq(max(data[,1])-1e4, max(data[,1]), 5e3)){
    freq=data[data[,1]==gen, 2:201]
    
    tmpVW=matrix(, ncol=200, nrow=D)
    
    for(d in 1:D){
      tmpFrq=as.numeric(freq[d,])
      tmpVW[d,]=tmpFrq*(1-tmpFrq)*Ns/Nu
    }
    
    tmpFreqDeme=apply(tmpVW, 2, function(x){length(which(x>2))})
    freqDeme=rbind(freqDeme,tmpFreqDeme)
  }
  }
    
    meanFreqDeme=apply(freqDeme, 2, mean)
    sdFreqDeme=apply(freqDeme, 2, sd)
    points(log10(Ns), meanFreqDeme, pch=4, cex=0.4, lwd=0.4, col=colors[m])
}

dev.off()
