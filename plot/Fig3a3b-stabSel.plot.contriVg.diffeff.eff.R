library(RColorBrewer)
library(plotrix)
options(digits=10)

args=commandArgs(trailingOnly = TRUE)
colors=brewer.pal(n=9,name = "Set1")[2:4]

Nm=c(0, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 2 , 5, 10)
Vs=c(0.5, 1)
numLoci=c(100, 200)
popSize=c(100, 200)
symL=c(1, 1)
Nu=0.01
effectSize=c(0.05, 0.1, 0.2)
Ns=c(0.25, 1, 4)
propo=c(0.6, 0.3, 0.1)
lNm=0.005 #. the lower boundary of axis

### Row-Nm; Column-Ns&L
tmpMatrix=matrix(, nrow=length(Nm), ncol=length(effectSize)*length(numLoci))
proVT_mean=tmpMatrix
proVT_sd=tmpMatrix
proVW_mean=tmpMatrix
proVW_sd=tmpMatrix

for(l in c(length(numLoci))){ #number of loci
  for(i in 1:length(Nm)){ # migration
    m=Nm[i]
    data=data.frame()
    for(r in seq(1, 10, 1)){
      #### simulation
      x=read.table(paste(args[1], "/simu/het_trait_frqB1_Nm", m, "_Vs", Vs[l], "_N", popSize[l], "_D100_L", numLoci[l], "_r", r, sep=""), header = TRUE)
      generationSample=length(x[,1])
      tmp=apply(x[seq((generationSample-50),generationSample, 25), ], 2, mean)
      data=rbind(data, tmp)
    }
    
    # selection
    HT=data[,2:4]
    VT=t(propo*Ns*t(HT))
    totalVT=apply(VT, 1, sum)
    
    for(j in 1:3){
      proVT_mean[i, 3*(l-1)+j]=mean(VT[,j]/totalVT)
      proVT_sd[i,3*(l-1)+j]=1.96*std.error(VT[,j]/totalVT)
      
      proVWperDeme=numLoci[l]*propo[j]*effectSize[j]^2*data[, seq(j+6, length(data[1,]), 6)]/data[, seq(12, length(data[1,]), 6)]
      proVW=apply(proVWperDeme, 1, mean)
      proVW_mean[i,3*(l-1)+j]=mean(proVW)
      proVW_sd[i,3*(l-1)+j]=1.96*std.error(proVW)
    }
  }
}


pdf(paste("stabSel_contriVg_diffeff.eff.pdf", sep=""), height=3, width=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)

### proportion of genic VT 
plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(0, 0.8), xlab="", ylab="", xaxt="n",cex.axis=0.8, main="Full population")
### y axis
mtext(expression("Proportion of total "*italic(V[T]^"(gen)")), side=2, line=1.9, las=0, cex=1.2)
### x axis
axis(side=1, labels=c(0, 0.02, 0.1, 1, 10), at=log10(c(0.005, 0.02, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.03, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
mtext(expression(italic(Nm)), side=1, line = 1.6, cex=1.2)
axis.break(axis=1, breakpos=log10(0.012)) 

legend("topleft", inset=c(-0.05,0),legend=expression(italic(Ns)), cex=0.8, bty="n", horiz=TRUE)
legend("topleft", inset=c(0.1,0), legend=c(0.25, 1, 4), 
       text.col=colors, pch=c(NA, NA, NA), cex=0.8,  bty="n", horiz=TRUE)

# proportion of genic for each effect sizes
for(i in 1:3){
  lines(log10(c(0.005, 10)), c(propo[i], propo[i]), col="grey30", lty=2, lwd=0.5)
}
# simulations
for(j in 1:3){
  for(l in 1:2){
    for(i in 1:length(Nm)){
      m=if(Nm[i]==0) lNm else Nm[i]
      meanVal=proVT_mean[i,3*(l-1)+j]
      u=meanVal+proVT_sd[i,3*(l-1)+j]
      d=meanVal-proVT_sd[i,3*(l-1)+j]
      
      points(log10(m), meanVal, col=colors[j], cex=0.7, pch=symL[l])
      lines(log10(c(m, m)), c(d, u), col=colors[j])

    }
  }
}

# Analytical results
for(j in 1:3){
  pmean=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_pmean_HW.dat", sep=""))
  seff=read.table(paste(args[1], "/predict/Nu0.01_uL0.01_Nm_NmEff_NsEff_pmean_HW.dat", sep="")) 
  
  HT=2*pmean[,2]*(1-pmean[,2])
  HTeff=2*seff[,5+j]*(1-seff[,5+j])
  
  if(j==1){
    VT=effectSize[j]^2*propo[j]*HT
    Nm1=pmean[,1]
    VTeff=effectSize[j]^2*propo[j]*HTeff
    Nmeff=seff[,1]
  }
  else{
    VT=cbind(VT, effectSize[j]^2*propo[j]*HT)
    VTeff=cbind(VTeff, effectSize[j]^2*propo[j]*HTeff)
  }
}

VTtotal=apply(VT, 1, sum)
VTeffTotal=apply(VTeff, 1, sum)

for(j in 1:3){
  tmpPro=VT[,j]/VTtotal
  lines(log10(Nm1[Nm1>=0.015]), tmpPro[Nm1>=0.015],col=colors[j], lty="11")    
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), tmpPro[1:8], col=colors[j],lty="11")
  
  tmpPro=VTeff[,j]/VTeffTotal
  lines(log10(Nmeff[Nmeff>=0.015]), tmpPro[Nmeff>=0.015],col=colors[j], lwd=1)    
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), tmpPro[1:8], col=colors[j], lwd=1)
}

### proportion of genic VW
plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(0, 0.8), xlab="", ylab="", xaxt="n", cex.axis=0.8, main="Within demes")
# x axis
axis(side=1, labels=c(0, 0.02, 0.1, 1, 10), at=log10(c(0.005, 0.02, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.03, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
mtext(expression(italic(Nm)), side=1, line = 1.6, cex=1.2)
axis.break(axis=1, breakpos=log10(0.012)) 
# y axis
mtext(expression("Proportion of total "*italic(V[W]^"(gen)")), side=2, line=1.9, las=0, cex=1.2)

# prop. of mutational target size for each effect sizes
for(i in 1:3){
  lines(log10(c(0.005, 10)), c(propo[i], propo[i]), col="grey30", lty=2, lwd=0.5)
}

# simulation
for(j in 1:3){
  for(l in 1:2){
    for(i in 1:length(Nm)){
      m=if(Nm[i]==0) lNm else Nm[i]
      meanVal=proVW_mean[i,3*(l-1)+j]
      u=meanVal+proVW_sd[i,3*(l-1)+j]
      d=meanVal-proVW_sd[i,3*(l-1)+j]
      
      points(log10(m), meanVal, col=colors[j], cex=0.7, pch=symL[l])
      lines(log10(c(m, m)), c(d, u), col=colors[j])
    }
  }
}


# Analytical results
for(j in 1:3){
  HW=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_pmean_HW.dat", sep=""))
  seff=read.table(paste(args[1], "/predict/Nu0.01_uL0.01_Nm_NmEff_NsEff_pmean_HW.dat", sep="")) 
  if(j==1){
    VW=effectSize[j]^2*propo[j]*HW[, 3]
    Nm1=HW[,1]
    
    VWeff=effectSize[j]^2*propo[j]*seff[, 8+j]
    Nmeff=seff[,1]
  }
  else{
    VW=cbind(VW, effectSize[j]^2*propo[j]*HW[, 3])
    VWeff=cbind(VWeff, effectSize[j]^2*propo[j]*seff[, 8+j])
  }
}

VWtotal=apply(VW, 1, sum)
VWeffTotal=apply(VWeff, 1, sum)

for(j in 1:3){
  tmpPro=VW[,j]/VWtotal
  lines(log10(Nm1[Nm1>=0.015]), tmpPro[Nm1>=0.015],col=colors[j], lty="11")    
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), tmpPro[1:8], col=colors[j], lty="11")
  
  tmpPro=VWeff[,j]/VWeffTotal
  lines(log10(Nmeff[Nmeff>=0.015]), tmpPro[Nmeff>=0.015],col=colors[j], lwd=1)    
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), tmpPro[1:8], col=colors[j],  lwd=1 )
}

dev.off()
