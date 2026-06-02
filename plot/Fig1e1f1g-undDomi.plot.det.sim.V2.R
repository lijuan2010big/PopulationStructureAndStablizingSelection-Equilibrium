library(RColorBrewer)
library(plotrix)
options(digits=10)

args=commandArgs(trailingOnly = TRUE)
colors=brewer.pal(n=9,name = "Set1")[2:4]

Nm=c(0, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 2 , 5, 10)
Ns=c(0.1, 1, 5)

critNm=c(0.624, 0.089, 0.046)

lwidth=c(0.5, 1, 2)

intfrq=0.5

numRow=length(Nm)
initTmp=matrix(rep(0,numRow), ncol=3, nrow=numRow)
HT_mean=initTmp
HT_sd=initTmp
HS_mean=initTmp
HS_sd=initTmp
FST_mean=initTmp
FST_sd=initTmp


for(i in 1:length(Nm)){
  m=Nm[i]
  for(j in 1:3){
    s=Ns[j]
    #data=data.frame()
    HT=vector()
    HS=vector()
    FST=vector()
    for(r in seq(1, 10, 1)){
      #### simulation
      x=read.table(paste(args[1], "/simu/frq_", intfrq,"_Nm", m, "_Ns", s, "_r", r, sep=""), header = FALSE)
      l=length(x[,1])
      d=length(x[1,])
      tmp=apply(x[(l-100):l,2:d], 1, mean)
      HTtmp=mean(2*tmp*(1-tmp))
      HT=c(HT, HTtmp)
      tmp=apply(2*x[(l-100):l,2:d]*(1-x[(l-100):l,2:d]), 1, mean)
      HStmp=mean(tmp)
      HS=c(HS, HStmp)
      FST=c(FST, 1-HStmp/HTtmp)
    }
    
    HT_mean[i,j]=mean(HT)
    HT_sd[i,j]=1.96*std.error(HT)
    HS_mean[i,j]=mean(HS)
    HS_sd[i,j]=1.96*std.error(HS)
    FST_mean[i,j]=mean(FST)
    FST_sd[i,j]=1.96*std.error(FST)
  }
}


pdf(paste("undDomi_det_sim_", intfrq, ".V2.pdf", sep=""), height=3, width=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)
lNm=0.005
#### HT as a functioni of Nm

plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(-3, log10(0.51)), xlab="", ylab="", yaxt="n", xaxt="n",
     cex.axis=0.8, main="Full population")

posLab=c(0.001, 0.01, 0.1, 0.5)
axis(side=2, labels=posLab, at=log10(posLab), cex.axis=0.8)

axis(side=1, labels=c(0, 0.02, 0.1, 1, 10), at=log10(c(0.005, 0.02, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.03, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
mtext(expression(italic(Nm)), side=1, line = 1.6, cex=1.2)
mtext(expression("Heterozygosity, "*italic(H[T])), side=2, line=2, las=0, cex=1.2)
axis.break(axis=1, breakpos=log10(0.012)) 

legend("bottomleft", legend=c(Ns, "Neutral"), col=c(colors, "black"), title = expression(italic(Ns)), lty=c(rep("solid", 3), "21"), seg.len = 1, cex=0.8, bty="n")



for(j in 1:3){
  s=Ns[j]
for(i in 1:length(Nm)){
  m=Nm[i]
    if(m==0){
      points(log10(0.005), log10(HT_mean[i,j]), col=colors[j], pch=1, cex=0.7)
      if(HT_mean[i,j]-HT_sd[i,j]<0){
        lines(log10(c(lNm, lNm)), c(-4,log10(HT_mean[i,j]+HT_sd)), col=colors[j])
      }
      else if(HT_mean[i,j]+HT_sd[i,j]>0.5){
        lines(log10(c(lNm, lNm)), c(log10(HT_mean[i,j]-HT_sd[i,j]),log10(0.5)), col=colors[j])
      }
      else{
        lines(log10(c(lNm, lNm)), c(log10(HT_mean[i,j]-HT_sd[i,j]),log10(HT_mean[i,j]+HT_sd[i,j])), col=colors[j])
      }
      
    }
    else{
    points(log10(m), log10(HT_mean[i,j]), col=colors[j], pch=1, cex=0.7)
    if(HT_mean[i,j]-HT_sd[i,j]<0){
      lines(log10(c(m, m)), c(-4,log10(HT_mean[i,j]+HT_sd)), col=colors[j])
    }
    else if(HT_mean[i,j]+HT_sd[i,j]>0.5){
      lines(log10(c(m, m)), c(log10(HT_mean[i,j]-HT_sd[i,j]),log10(0.5)), col=colors[j])
    }
    else{
      lines(log10(c(m, m)), c(log10(HT_mean[i,j]-HT_sd[i,j]),log10(HT_mean[i,j]+HT_sd[i,j])), col=colors[j])
    }
    }
  }
  #### Analytical results
  pmean=read.table(paste(args[2], "/Nu0.01_Ns", s, "_Nm_pmean.dat", sep=""))
  HT=2*pmean[,2]*(1-pmean[,2])
  lines(log10(pmean[pmean[,1]>=0.015,1]), log10(HT[pmean[,1]>=0.015]),col=colors[j])    
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(HT[1:8]), col=colors[j])
  
  lines(rep(log10(critNm[j]), 2), log10(c(0.00001, 0.5)), col=colors[j], lwd=0.7)
}

pmean=read.table(paste(args[2], "/Nu0.01_Ns0_Nm_pmean.dat", sep=""))
HT=2*pmean[,2]*(1-pmean[,2])
lines(log10(pmean[pmean[,1]>=0.015,1]), log10(HT[pmean[,1]>=0.015]),col="black", lty="21", lwd=1.5)    
lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(HT[1:8]), col="black", lty="21", lwd=1.5)


###HW as a function of Nm.
plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(-3, log10(0.51)), xlab="", ylab="", yaxt="n", xaxt="n",
     cex.axis=0.8, main="Within demes")
posLab=c(0.001, 0.01, 0.1, 0.5)
axis(side=2, labels=posLab, at=log10(posLab), cex.axis=0.8)
axis(side=1, labels=c(0, 0.02, 0.1, 1, 10), at=log10(c(0.005, 0.02, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.03, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
mtext(expression(italic(Nm)), side=1, line = 1.6, cex=1.2)
mtext(expression("Heterozygosity, "*italic(H[W])), side=2, line=2, las=0, cex=1.2)
axis.break(axis=1, breakpos=log10(0.012)) 


for(j in 1:3){
  s=Ns[j]
for(i in 1:length(Nm)){
  m=Nm[i]
    if(m==0){
      points(log10(0.005), log10(HS_mean[i,j]), col=colors[j], pch=1, cex=0.7)
      if(HS_mean[i,j]-HS_sd[i,j]<0){
        lines(log10(c(lNm, lNm)), c(-4,log10(HS_mean[i,j]+HS_sd[i,j])), col=colors[j])
      }
      else{
        lines(log10(c(lNm, lNm)), c(log10(HS_mean[i,j]-HS_sd[i,j]),log10(HS_mean[i,j]+HS_sd[i,j])), col=colors[j])
      }
    }
    else{
    points(log10(m), log10(HS_mean[i,j]), col=colors[j], pch=1, cex=0.7)
    if(HS_mean[i,j]-HS_sd[i,j]<0){
      lines(log10(c(m, m)), c(-4,log10(HS_mean[i,j]+HS_sd[i,j])), col=colors[j])
    }
    else{
      lines(log10(c(m, m)), c(log10(HS_mean[i,j]-HS_sd[i,j]),log10(HS_mean[i,j]+HS_sd[i,j])), col=colors[j])
    }
    }
    

  }
  #### Analytical results
  HW=read.table(paste(args[2], "/Nu0.01_Ns", s, "_Nm_HW.dat", sep=""))
  lines(log10(HW[HW[,1]>=0.015,1]), log10(HW[HW[,1]>=0.015,2]),col=colors[j])    
  options(digits=10)
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(HW[1:8, 2]), col=colors[j])
  
  lines(rep(log10(critNm[j]), 2), log10(c(0.00001, max(HW[,2]))), col=colors[j], lwd=0.7)
}

HW=read.table(paste(args[2], "/Nu0.01_Ns0_Nm_HW.dat", sep=""))
lines(log10(HW[HW[,1]>=0.015,1]), log10(HW[HW[,1]>=0.015,2]),col="black", lty="21", lwd=1.5)    
options(digits=10)
lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(HW[1:8, 2]), col="black", lty="21", lwd=1.5)

#### FST
plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(0,1), xlab="", ylab="", cex.axis=0.8, xaxt="n")
# x axis
axis(side=1, labels=c(0, 0.02, 0.1, 1, 10), at=log10(c(0.005, 0.02, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.03, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
mtext(expression(italic(Nm)), side=1, line = 1.6, cex=1.2)
axis.break(axis=1, breakpos=log10(0.012)) 
# y axis
mtext(expression(italic("F"["ST"])), side=2, line=1.8, las=0, cex=1.2)

for(j in 1:3){
  s=Ns[j]
for(i in 1:length(Nm)){
  m=Nm[i]
    up = FST_mean[i,j]+FST_sd[i,j]
    down=FST_mean[i,j]-FST_sd[i,j]
    
  
    if(m==0){
      points(log10(0.005), FST_mean[i,j], col=colors[j], pch=1,  cex=0.7)
      if( up > 1 ){
        lines(log10(c(lNm, lNm)), c(down,1), col=colors[j])
      }
      else if( down < 0){
        lines(log10(c(lNm, lNm)), c(0,up), col=colors[j])
      }
      else{
        lines(log10(c(lNm, lNm)), c(down, up), col=colors[j])
      }
    }
    else{
      points(log10(m), FST_mean[i,j], col=colors[j], pch=1,  cex=0.7)
    if( up > 1 ){
      lines(log10(c(m, m)), c(down,1), col=colors[j])
    }
    else if( down < 0){
      lines(log10(c(m, m)), c(0,up), col=colors[j])
    }
    else{
      lines(log10(c(m, m)), c(down, up), col=colors[j])
    }
    }
}
  #### Analytical results
  pmean=read.table(paste(args[2], "/Nu0.01_Ns", s, "_Nm_pmean.dat", sep=""))
  HW=read.table(paste(args[2], "/Nu0.01_Ns", s, "_Nm_HW.dat", sep=""))
  HT=2*pmean[,2]*(1-pmean[,2])
  FST=(HT-HW[,2])/HT
  lines(log10(pmean[pmean[,1]>=0.015,1]), FST[pmean[,1]>=0.015], lty=1, col=colors[j])    
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), FST[1:8], lty=1, col=colors[j])
  
  ### Add Nm_crit lines
  lines(rep(log10(critNm[j]), 2), c(-0.1, 1-max(HW[,2])/(2*0.5*0.5)), col=colors[j], lwd=0.7)
}

Nu=0.01
Nm1=seq(0.015,10,0.01)
lines(log10(Nm1), 1/(4*Nm1+1+8*Nu),  col="black", lty="21")
Nm1=seq(0,0.007,0.001)
lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), 1/(4*Nm1+1+8*Nu),  col="black", lty="21")


dev.off()
