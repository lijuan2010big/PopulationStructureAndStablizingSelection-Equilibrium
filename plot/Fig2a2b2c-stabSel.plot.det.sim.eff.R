library(RColorBrewer)
library(plotrix)
options(digits=10)

args=commandArgs(trailingOnly = TRUE)

Nm=c(0, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 2 , 5, 10)
Ns=c(0.1, 1, 5)
Vs=c(5, 0.5, 0.1)
colors=brewer.pal(n=9,name = "Set1")[2:4]
Nu=0.01

critNm=c(0.6236608707517198, 0.08876520761800123, 0.046083603164263603)

lwidth=c(0.5, 1, 2)

tmpInit=matrix(, nrow=length(Nm), ncol=3)
VT_mean=tmpInit
VT_sd=tmpInit
VW_mean=tmpInit
VW_sd=tmpInit
FST_mean=tmpInit
FST_sd=tmpInit

for(i in 1:length(Nm)){
  m=Nm[i]
  for(j in 1:3){
    s=Ns[j]
    HT=vector()
    HW=vector()
    FST=vector()
    for(r in seq(1, 10, 1)){
      #### simulation
      data=read.table(paste(args[1], "/simu/het_trait_frqB1_Nu0.01_Nm", m, "_Ns", s, "_N200_D100_L200_r", r, sep=""), header = TRUE)
      l=length(data[,1])
      HTtmp=data[seq((l-50),l, 25), 2]
      HT=c(HT,HTtmp)
      HWtmp=apply(data[seq((l-50),l, 25), seq(5, length(data[10,]), 3)], 1, mean)
      HW=c(HW, HWtmp)
      FST=c(FST, 1-HWtmp/HTtmp)
    }
    
    VT=s/(Nu*2)*HT
    VW=s/(Nu*2)*HW
    VT_mean[i,j]=mean(VT)
    VT_sd[i,j]=1.96*std.error(VT)
    VW_mean[i,j]=mean(VW)
    VW_sd[i,j]=1.96*std.error(VW)
    FST_mean[i,j]=mean(FST)
    FST_sd[i,j]=1.96*std.error(FST)
    print(c(sd(FST),std.error(FST)))
  }
}

pdf(paste("stabSel_det_sim.eff.pdf", sep=""), height=3, width=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)

lNm=0.005

#### genic VT as a function of Nm.
plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(0, log10(125)), xlab="", ylab="", yaxt="n", xaxt="n",
     cex.axis=0.8, main="Full population")
### set up y axis
posLab=c(1, 2, 5, 10, 20, 50, 125)
axis(side=2, labels=posLab, at=log10(posLab))
mtext(expression("Scaled genic variance, "*italic(V[T]^"(gen)"/V["*"])), side=2, line=1.6, las=0, cex=0.8)
### set up x axis
axis(side=1, labels=c(0, 0.02, 0.1, 1, 10), at=log10(c(0.005, 0.02, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.03, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
mtext(expression(italic(Nm)), side=1, line = 1.6, cex=1.2)
### set up break x axis
axis.break(axis=1, breakpos=log10(0.012)) 

for(j in 1:3){
  for(i in 1:length(Nm)){
    m=if(Nm[i]==0) lNm else Nm[i]
    up=VT_mean[i,j]+VT_sd[i,j]
    down=VT_mean[i,j]-VT_sd[i,j]
    
      points(log10(m), log10(VT_mean[i,j]), col=colors[j], pch=1, cex=0.7)
      if(down<0){
        lines(log10(c(m, m)), c(-4,log10(up)), col=colors[j])
      }
      else{
        lines(log10(c(m, m)), c(log10(down),log10(up)), col=colors[j])
      }

}
  #### Analytical results
  pmean=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_pmean_HW.dat", sep=""))
  VT=Ns[j]/Nu*pmean[,2]*(1-pmean[,2])
  lines(log10(pmean[pmean[,1]>=0.015,1]), log10(VT[pmean[,1]>=0.015]),col=colors[j], lty="11" )    
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(VT[1:8]), col=colors[j],lty="11")
  
  ## s_eff
  seff=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_NmEff_NsEff_pmean_HW.dat", sep=""))
  VT=Ns[j]/Nu*seff[,4]*(1-seff[,4])
  lines(log10(seff[seff[,1]>=0.015,1]), log10(VT[seff[,1]>=0.015]),col=colors[j],lwd=1 )    
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(VT[1:8]), col=colors[j], lwd=1 )
}


legend("topright", title="Approximation", legend=c("LE", "LD"), lty=c("11", "solid"), seg.len = 1, cex=0.8, bty="n", horiz = TRUE)
legend("topright", inset=c(0.1, 0.2), legend=c(1/(2*Vs)), text.col=c(colors), pch=c(NA, NA, NA), title = expression( italic("Ns") ), title.col="black", title.adj = 0.7,cex=0.8, bty="n")

### genic VW as a function of Nm.
plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(-1,log10(10)), xlab="", ylab="", yaxt="n", xaxt="n",
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

for(j in 1:3){
  s=Vs[j]
for(i in 1:length(Nm)){
  m=Nm[i]
    down=VW_mean[i,j]-VW_sd[i,j]
    up=VW_mean[i,j]+VW_sd[i,j]
    
    if(m==0){
      points(log10(lNm), log10(VW_mean[i,j]),col=colors[j], pch=1, cex=0.7)
      if(down<0){
        lines(log10(c(lNm, lNm)), c(-4,log10(up)),col=colors[j])
      }
      else{
        lines(log10(c(lNm, lNm)), c(log10(down),log10(up)), col=colors[j])
      }
    }
    else{
    points(log10(m), log10(VW_mean[i,j]),col=colors[j], pch=1, cex=0.7)
    if(down<0){
      lines(log10(c(m, m)), c(-4,log10(up)),col=colors[j])
    }
    else{
      lines(log10(c(m, m)), c(log10(down),log10(up)), col=colors[j])
    }
    }
  }
  #### Analytical results
  HW=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_pmean_HW.dat", sep=""))
  VW=Ns[j]/(Nu*2)*HW[,3]
  lines(log10(HW[HW[,1]>=0.015,1]), log10(VW[HW[,1]>=0.015]),col=colors[j], lty="11")    
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(VW[1:8]), col=colors[j],lty="11")
  
  seff=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_NmEff_NsEff_pmean_HW.dat", sep=""))
  VW=Ns[j]/(Nu*2)*seff[,5]
  lines(log10(seff[seff[,1]>=0.015,1]), log10(VW[seff[,1]>=0.015]),col=colors[j], lwd=1)    
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(VW[1:8]), col=colors[j], lwd=1 )
}

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

for(j in 1:3){
  s=Ns[j]
for(i in 1:length(Nm)){
  m=Nm[i]

  if(m==0){
    points(log10(lNm), FST_mean[i,j],  pch=1, col=colors[j], cex=0.7)
    lines(log10(c(lNm, lNm)), c(FST_mean[i,j]-FST_sd[i,j],FST_mean[i,j]+FST_sd[i,j]), col=colors[j])
  }else{
    points(log10(m), FST_mean[i,j],  pch=1, col=colors[j], cex=0.7)
    lines(log10(c(m, m)), c(FST_mean[i,j]-FST_sd[i,j],FST_mean[i,j]+FST_sd[i,j]), col=colors[j])
  }
}
    #### Analytical results
  pmeanHW=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_pmean_HW.dat", sep=""))
  HW=pmeanHW[,3]
  HT=2*pmeanHW[,2]*(1-pmeanHW[,2])
  FST=(HT-HW)/HT
  lines(log10(pmeanHW[pmeanHW[,1]>=0.015,1]), FST[pmean[,1]>=0.015], col=colors[j],lty="11")    
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), FST[1:8], col=colors[j], lty="11")
  
  seff=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_NmEff_NsEff_pmean_HW.dat", sep=""))
  HW=seff[,5]
  HT=2*seff[,4]*(1-seff[,4])
  FST=(HT-HW)/HT
  lines(log10(seff[seff[,1]>=0.015,1]), FST[seff[,1]>=0.015],  col=colors[j], lwd=1)    
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), FST[1:8], col=colors[j],, lwd=1)
  
  
}

Nu=0.01
Nm1=seq(0.015,15,0.01)
lines(log10(Nm1), 1/(4*Nm1+1+8*Nu),  col="black", lty="21")
Nm1=seq(0,0.007,0.001)
lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), 1/(4*Nm1+1+8*Nu),  col="black", lty="21")


dev.off()
