library(RColorBrewer)
library(plotrix)
options(digits=10)

args=commandArgs(trailingOnly = TRUE)
colors=brewer.pal(n=9,name = "Set1")[2:4]

Nm=c(0, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 2 , 5, 10)
#Nm=c(0, 0.1, 0.5,  2)
Nu=0.01
L=200
a=0.1
N=200
Vs=c(10, 1, 0.2)
Ns=c(0.1, 1, 5)

Vstar=4*Nu*L/N*Vs

tmpInit=matrix(, nrow=length(Nm), ncol=3)
QST_mean=tmpInit
QST_sd=tmpInit
GT_mean=tmpInit
GT_sd=tmpInit
GW_mean=tmpInit
GW_sd=tmpInit

for(i in 1:length(Nm)){
  m=Nm[i]
  for(j in 1:3){
    s=Ns[j]
    GW=vector()
    GT=vector()
    QST=vector()
    for(r in seq(1, 4, 1)){
      #### simulation het_trait_frqB_Nu0.01_Nm1_Ns1_N100_D100_L100_r1
      data=read.table(paste(args[1], "/simu/het_trait_frqB1_Nu0.01_Nm", m, "_Ns", s, "_N200_D100_L200_r", r, sep=""), header = TRUE)
      l=length(data[,1])
      GWtmp=apply(data[seq((l-50),l,25), seq(7, length(data[1,]), 3)], 1, mean)
      GTtmp=mean(data[seq((l-50),l,25), 4])
      GBtmp=apply((data[seq((l-50),l,25), seq(6, length(data[1,]), 3)]-data[seq((l-50),l,25),3])^2, 1, mean)
      GW=c(GW, GWtmp)
      GT=c(GT, GTtmp)
      QSTtmp=GBtmp/(GBtmp+2*GWtmp)
      QST=c(QST, QSTtmp)
    }

    GW_mean[i,j]=mean(GW)
    GW_sd[i,j]=1.96*std.error(GW)
    GT_mean[i,j]=mean(GT)
    GT_sd[i,j]=1.96*std.error(GT)
    QST_mean[i,j]=mean(QST)
    QST_sd[i,j]=1.96*std.error(QST)
  }
}
    
geneticVW=function(Vg, Vs){
  return((Vg-Vs+sqrt(Vg^2+6*Vg*Vs+Vs^2))/4)
}

pdf(paste("stabSel_trait_qst.eff.pdf", sep=""), height=3, width=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)

lNm=0.005
### GT

if(FALSE){
plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(log10(0.007), log10(1)), xlab="", ylab="", yaxt="n", xaxt="n",
     cex.axis=0.8, main="Full population")
### set up y axis
posLab=c(0.005, 0.01, 0.05, 0.2 ,0.5,1)
axis(side=2, labels=posLab, at=log10(posLab), cex.axis=0.8)
mtext(expression("Genetic variance, "*italic("V"["T"])), side=2, line=1.8, las=0, cex=1.2)
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
      points(log10(lNm), log10(GT_mean[i,j]),  pch=1, col=colors[j], cex=0.7)
      lines(log10(c(lNm, lNm)), log10(c(GT_mean[i,j]-GT_sd[i,j],GT_mean[i,j]+GT_sd[i,j])), col=colors[j])
    }
    else{
      points(log10(m), log10(GT_mean[i,j]),  pch=1, col=colors[j], cex=0.7)
      lines(log10(c(m, m)), log10(c(GT_mean[i,j]-GT_sd[i,j],GT_mean[i,j]+GT_sd[i,j])), col=colors[j])
    }
    points(log10(m), log10(GT_mean[i,j]),  pch=1, col=colors[j], cex=0.7)
    lines(log10(c(m, m)), log10(c(GT_mean[i,j]-GT_sd[i,j],GT_mean[i,j]+GT_sd[i,j])), col=colors[j])
  }
    #### Analytical results
    HW=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_pmean_HW.dat", sep=""))
    VW=geneticVW(HW[,3]*L*a^2, Vs[j])
    VT=VW+VW/(2*HW[,1]+2*N*VW/Vs[j])
    lines(log10(HW[HW[,1]>0.015,1]), log10(VT[HW[,1]>0.015]), col=colors[j], lty="11")
    lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(VT[1:8]), col=colors[j], lty="11")
    
    seff=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_NmEff_NsEff_pmean_HW.dat", sep=""))
    VW=geneticVW(seff[,5]*L*a^2, Vs[j])
    VT=VW+VW/(2*seff[,2]+2*N*VW/Vs[j])
    lines(log10(seff[seff[,1]>0.015,1]), log10(VT[seff[,1]>0.015]), col=colors[j], lwd=1)
    lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(VT[1:8]), col=colors[j], lwd=1)
    
  
}


#### GW
plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(log10(0.007), log10(1)), xlab="", ylab="", yaxt="n", xaxt="n",cex.axis=0.8, main="Within demes")
### set up y axis
posLab=c(0.005, 0.01,  0.05,  0.2 ,0.5,1)
axis(side=2, labels=posLab, at=log10(posLab), cex.axis=0.8)
mtext(expression("Genetic variance, "*italic("V"["W"])), side=2, line=1.8, las=0, cex=1.2)
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
      points(log10(lNm), log10(GW_mean[i,j]),  pch=1, col=colors[j], cex=0.7)
      lines(log10(c(lNm, lNm)), log10(c(GW_mean[i,j]-GW_sd[i,j],GW_mean[i,j]+GW_sd[i,j])), col=colors[j])
    }
    else{
      points(log10(m), log10(GW_mean[i,j]),  pch=1, col=colors[j], cex=0.7)
      lines(log10(c(m, m)), log10(c(GW_mean[i,j]-GW_sd[i,j],GW_mean[i,j]+GW_sd[i,j])), col=colors[j])
    }
    
  }
    #### Analytical results
    HW=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_pmean_HW.dat", sep=""))
    VW=geneticVW(HW[,3]*L*a^2, Vs[j])
    lines(log10(HW[HW[,1]>0.015,1]),log10(VW[HW[,1]>0.015]), col=colors[j], lty="11")
    lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(VW[1:8]), col=colors[j],  lty="11")

    seff=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_NmEff_NsEff_pmean_HW.dat", sep=""))
    VW=geneticVW(seff[,5]*L*a^2, Vs[j])
    lines(log10(seff[seff[,1]>0.015,1]),log10(VW[seff[,1]>0.015]), col=colors[j], lwd=1)
    lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(VW[1:8]), col=colors[j],lwd=1)
    
}

}

#### QST
plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(0,0.2),, xlab="", ylab="", cex.axis=0.8, xaxt="n")
### set up y axis
mtext(expression(italic(Q[ST])), side=2, line=1.8, las=0, cex=1.2)
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
      points(log10(lNm), QST_mean[i,j],  pch=1, col=colors[j], cex=0.7)
      lines(log10(c(lNm, lNm)), c(QST_mean[i,j]-QST_sd[i,j],QST_mean[i,j]+QST_sd[i,j]), col=colors[j])
    }
    else{
      points(log10(m), QST_mean[i,j],  pch=1, col=colors[j], cex=0.7)
      lines(log10(c(m, m)), c(QST_mean[i,j]-QST_sd[i,j],QST_mean[i,j]+QST_sd[i,j]), col=colors[j])
    }
  }
  
  #### Analytical results
  HW=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_pmean_HW.dat", sep=""))
  VW=geneticVW(L*HW[,3]*a^2, Vs[j])
  QST=1/(1+4*HW[,1]+4*N*(VW/Vs[j]))
  lines(log10(HW[HW[,1]>0.015,1]), QST[HW[,1]>0.015], col=colors[j], lty="11")
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), QST[1:8], col=colors[j], lty="11")
  
  seff=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_NmEff_NsEff_pmean_HW.dat", sep=""))
  VW=geneticVW(L*seff[,5]*a^2, Vs[j])
  QST=1/(1+4*seff[,2]+4*N*(VW/Vs[j]))
  lines(log10(seff[seff[,1]>0.015,1]), QST[seff[,1]>0.015], col=colors[j],lwd=1 )
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), QST[1:8], col=colors[j],lwd=1)

}

Nu=0.01
Nm1=seq(0.015,15,0.01)
lines(log10(Nm1), 1/(4*Nm1+1+8*Nu),  col="black", lty="21")
Nm1=seq(0,0.007,0.001)
lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), 1/(4*Nm1+1+8*Nu),  col="black", lty="21")


### GT
plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(log10(0.15), log10(10)), xlab="", ylab="", yaxt="n", xaxt="n",
     cex.axis=0.8, main="Full population")
### set up y axis
posLab=c(0.1, 0.2, 0.4,  1,  2, 5, 10)
axis(side=2, labels=posLab, at=log10(posLab), cex.axis=0.8)
mtext(expression("Scaled genetic variance, "*italic("V"["T"]/V["*"])), side=2, line=1.8, las=0, cex=1)
### set up x axis
axis(side=1, labels=c(0, 0.02, 0.1, 1, 10), at=log10(c(0.005, 0.02, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.03, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
mtext(expression(italic(Nm)), side=1, line = 1.6, cex=1.2)
### set up break x axis
axis.break(axis=1, breakpos=log10(0.012)) 

legend("topright", legend=Ns, text.col=c(colors), pch=c(NA, NA, NA), title = expression( italic("Ns") ), title.col="black", title.adj = 0.7,cex=0.8, bty="n")


for(j in 1:3){
  s=Ns[j]
  for(i in 1:length(Nm)){
    m=Nm[i]
    if(m==0){
      points(log10(lNm), log10(GT_mean[i,j]/Vstar[j]),  pch=1, col=colors[j], cex=0.7)
      lines(log10(c(lNm, lNm)), log10(c(GT_mean[i,j]-GT_sd[i,j],GT_mean[i,j]+GT_sd[i,j])/Vstar[j]), col=colors[j])
    }
    else{
      points(log10(m), log10(GT_mean[i,j]/Vstar[j]),  pch=1, col=colors[j], cex=0.7)
      lines(log10(c(m, m)), log10(c(GT_mean[i,j]-GT_sd[i,j],GT_mean[i,j]+GT_sd[i,j])/Vstar[j]), col=colors[j])
    }
    points(log10(m), log10(GT_mean[i,j]/Vstar[j]),  pch=1, col=colors[j], cex=0.7)
    lines(log10(c(m, m)), log10(c(GT_mean[i,j]-GT_sd[i,j],GT_mean[i,j]+GT_sd[i,j])/Vstar[j]), col=colors[j])
  }
  #### Analytical results
  HW=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_pmean_HW.dat", sep=""))
  VW=geneticVW(HW[,3]*L*a^2, Vs[j])
  VT=VW+VW/(2*HW[,1]+2*N*VW/Vs[j])
  lines(log10(HW[HW[,1]>0.015,1]), log10(VT[HW[,1]>0.015]/Vstar[j]), col=colors[j], lty="11")
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(VT[1:8]/Vstar[j]), col=colors[j], lty="11")
  
  seff=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_NmEff_NsEff_pmean_HW.dat", sep=""))
  VW=geneticVW(seff[,5]*L*a^2, Vs[j])
  VT=VW+VW/(2*seff[,2]+2*N*VW/Vs[j])
  lines(log10(seff[seff[,1]>0.015,1]), log10(VT[seff[,1]>0.015]/Vstar[j]), col=colors[j], lwd=1)
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(VT[1:8]/Vstar[j]), col=colors[j], lwd=1)
  
}


#### GW
plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(log10(0.15), log10(10)), xlab="", ylab="", yaxt="n", xaxt="n",cex.axis=0.8, main="Within demes")
### set up y axis
posLab=c(0.1, 0.2, 0.4,  1,  2, 5, 10)
axis(side=2, labels=posLab, at=log10(posLab), cex.axis=0.8)
mtext(expression("Scaled genetic variance, "*italic("V"["W"]/V["*"])), side=2, line=1.8, las=0, cex=1)
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
    if(m==0) m=lNm else m=Nm[i]

    points(log10(m), log10(GW_mean[i,j]/Vstar[j]),  pch=1, col=colors[j], cex=0.7)
    lines(log10(c(m, m)), log10(c(GW_mean[i,j]-GW_sd[i,j],GW_mean[i,j]+GW_sd[i,j])/Vstar[j]), col=colors[j])
    
  }
  #### Analytical results
  HW=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_pmean_HW.dat", sep=""))
  VW=geneticVW(HW[,3]*L*a^2, Vs[j])
  lines(log10(HW[HW[,1]>0.015,1]),log10(VW[HW[,1]>0.015]/Vstar[j]), col=colors[j], lty="11")
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(VW[1:8]/Vstar[j]), col=colors[j],  lty="11")
  
  seff=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_NmEff_NsEff_pmean_HW.dat", sep=""))
  VW=geneticVW(seff[,5]*L*a^2, Vs[j])
  lines(log10(seff[seff[,1]>0.015,1]),log10(VW[seff[,1]>0.015]/Vstar[j]), col=colors[j], lwd=1)
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), log10(VW[1:8]/Vstar[j]), col=colors[j],lwd=1)
  
}

dev.off()