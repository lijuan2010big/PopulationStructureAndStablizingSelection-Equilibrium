library(RColorBrewer)
library(plotrix)
options(digits=10)

args=commandArgs(trailingOnly = TRUE)

colors=rev(brewer.pal(n = 11, name = "RdBu")[c(2, 3, 4, 9 , 11)])

Nm=c(0, 0.1, 0.5, 2, 10)
Vs=1
Nu=0.01
N=200
L=200
D=100

tmp=matrix(,nrow=length(Nm), ncol=L)
VT=tmp
VW=tmp
effectSize=tmp


pdf(paste(args[1], "/stabSel_mixEffSizeExp.pairwise.eff.1.pdf", sep=""), height=3, width=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)

### Genic variance in Full population
plot(NULL, xlim=c(0.5, 5.5), ylim=c(0,0.6), 
     , xlab="", ylab="", main="",
     type="l", lty=1, xaxt="n", yaxt="n") 
### x axis
mtext(expression(italic(Nm)), side=1, line = 1.5, cex=1.2)
axis(side=1, at=1:5, labels = Nm)
### y axis
mtext(expression("Proportion of "* italic(V[W]^"(gen)") *" in deme B"), side=2, line=2.35, las=0, cex=0.8)
mtext(expression("explained by loci with " * italic(tilde(v)>2)* " in deme A"), side=2, line=1.5, las=0, cex=0.8)
axis(side=2, at=seq(0,1,0.1), labels = seq(0,1,0.1), cex.axis=0.8)

for(i in 1:length(Nm)){
  data=read.table(paste(args[1], "/predict/Nu0.01_Nm", Nm[i],"_Ns_Vg_Fv1_Fv2_mixExp2.v2.sim.dat", sep=""))

  rgbValue=col2rgb(colors[i])/255
  Vg=sum(data[,6])
  VgV2=sum(data[,6]*data[,7])

  Vg=sum(data[,6])
  VgV2=sum(data[,6]*data[,8])
  points(i, VgV2/Vg, pch=8, col=colors[i], cex=0.5, lwd=1)
  
}

for(i in 1:length(Nm)){
  proportion1=c()
  proportion2=c()
  for(r in 1:10){
    
    data=read.table(paste(args[1], "/simu/Nu0.01_Vs1_aExpNuL_Nm", Nm[i], "_L200_N200_D100_r1.freq", sep=""))
    
    effectSize=as.matrix(data[1,2:length(data[1,])])
    NsNu=as.vector(N*effectSize^2/(2*Vs)/Nu)

    for(gen in seq(max(data[,1])-1e4, max(data[,1]), 5e3)){
      freq=data[data[,1]==gen, 2:(L+1)]

      Vw=t(t(as.matrix(freq*(1-freq)))*NsNu)
      
      for(d in seq(1, (D-2), 1)){
        tmp=Vw[(d+1):D,]
        VW=apply(tmp, 1, sum)
        
        VW1=apply(tmp, 1, function(x){sum(x[Vw[d,]>0])})
        proportion1=c(proportion1, VW1/VW)
        
        VW2=apply(tmp, 1, function(x){sum(x[Vw[d,]>2])})
        proportion2=c(proportion2, VW2/VW)
      }
    }
  }
    qProp=quantile(proportion2, probs=c(0.05, 0.25, 0.5, 0.75, 0.95 ))
    boxplot(qProp, at=i, boxwex=0.4, add=TRUE, range=Inf, lty="solid", col=NA, border=colors[i], xaxt="n", yaxt="n", boxlwd=1.2, medlwd = 1.2, whisklwd = 1.2)
}


dev.off()
