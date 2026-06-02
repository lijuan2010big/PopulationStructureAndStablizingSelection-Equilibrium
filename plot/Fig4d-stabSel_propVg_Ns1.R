library(RColorBrewer)
library(fields)
library(plotrix)

args=commandArgs(trailingOnly = TRUE)

Nm=c(0, 0.25, 2, 10)
colors=rev(brewer.pal(n = 11, name = "RdBu")[c(2, 3, 9, 11)])
Ns=c(0.1, 1, 5)
Nu=0.01
L=200
N=200
D=100

critNm=c(0.62366087075, 0.088765207618, 0.046083603)

pdf("stabSel_propVg_Ns1.pdf", width=3, height=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)

# Analytical prediction
for(j in c(2)){
  s=Ns[j]
  #substitute(italic("Ns") * " = " * s *","*italic("Nm")[crit]*"="*mc, list(s=s, mc=round(critNm[j], digits = 3)))
  plot(NULL, xlim=c(0, Ns[j]/(4*Nu)), ylim=c(0,1), ylab="", main=substitute(italic("Ns") * " = " * s , list(s=s)), 
       xlab="", xaxt="n")
  axis(side=1, at=seq(0,25,5), label=seq(0,25,5), line=0, cex.axis=0.9)
  mtext(expression(italic(v)), side=1, line = 1.4, cex=1.2)
  mtext(expression("Proportion of total" * italic(V[W]^"(gen)")), side=2, line=2.1, las=0, cex=1.2)
  
  for(i in 1:4){
    m=Nm[i]
    data=read.table(paste(args[1], "/N200/icdfVg_Nu0.01_Nm", Nm[i], "_Ns",Ns[j],".eff.dat", sep=""))
    lines(data[,1], data[,2]/data[1,2], col=colors[i])
    
    data=read.table(paste(args[1], "/N200/icdfVg_Nu0.01_Nm", Nm[i], "_Ns",Ns[j],".dat", sep=""))
    lines(data[,1], data[,2]/data[1,2], col=colors[i], lty="11")
  }
  
  data=read.table(paste(args[1], "/N200/icdfVg_Nu0.01_NmCrt_Ns",Ns[j],".eff.dat", sep=""))
  lines(data[,1], data[,2]/data[1,2], col="grey50")
  
  data=read.table(paste(args[1], "/N200/icdfVg_Nu0.01_NmCrt_Ns",Ns[j],".dat", sep=""))
  lines(data[,1], data[,2]/data[1,2], lty="11", col="grey50") 
  

  
}

propQuantInv <- function(distr, value) {
  prop=sum(distr[distr>=value])/sum(distr)
  return(prop)
}

Nm=c(Nm,0.089)
colors=c(colors, "grey50")

for(i in 1:5){
for(j in c(2)){
  
  qValue=seq(0, 0.6, 0.1)
  
  low=Ns[j]/Nu * 1/400 * (1-1/400)
  up=Ns[j]/(4*Nu)
  qValue=c(seq(low, up-1, 2), up)
  
  
  qProp=matrix(, ncol=length(qValue), nrow=10)
  
  for(r in 1:10){
    # simulation
    freq=read.table(paste(args[1], "/simu/frqB1_Nu0.01_Nm", Nm[i],"_Ns", Ns[j],"_N200_D100_L200_r",r, ".freq", sep=""))
    
    qPropR=matrix(, ncol=length(qValue), nrow=3)
    g=1
    
    if((Nm[i]<0.5 && Nm[i]>0) || (Nm[i]==0 && r<5)) lowGen=490000 else lowGen=90000
    for(gen in seq(lowGen, lowGen+10000, 5000)){
      frq1=c(as.matrix(freq[freq[,1]==gen,2:(L+1)]))
      tmp=frq1*(1-frq1)*Ns[j]/Nu
      vgenic=tmp[tmp>0]
      for(x in 1:length(qValue)){
        qPropR[g,x]=propQuantInv(vgenic, qValue[x])
      }
      g=g+1
    }
    
    qProp[r,]=apply(qPropR, 2, mean)
  }
  
  qPropMean=apply(qProp, 2, mean)
  qPropSD=apply(qProp, 2, std.error)*1.96

  points(qValue, qPropMean,  col=colors[i], pch=1, cex=0.5)
  
  for(numQ in seq(1, length(qValue), 1)){
    lines(rep(qValue[numQ],2), c(qPropMean[numQ]+qPropSD[numQ], qPropMean[numQ]-qPropSD[numQ]),  col=colors[i])
  }
  
}
}

if(FALSE){
#### VT

plot(NULL, xlim=log10(c(0.049,125)), ylim=c(0,1), ylab="", , xlab="", xaxt="n", yaxt="n",
     main=substitute("In full population, "*italic("Nm")*"="*x, list(x=m)))
posLab=c(0.1, 1, 5,  25, 125)
axis(1, at=log10(posLab), labels = posLab, cex.axis=1)
mtext(expression(italic(v)), side=1, line = 1.4, cex=1.2)
posLab=seq(0,1,0.2)
axis(2, at=posLab, labels = posLab)
mtext(expression("Proportion of total "*italic(V[W]^"(gen)")), side=2, line=1.8, las=0, cex=1.2)

propQuantInv <- function(distr, value) {
  prop=sum(distr[distr>=value])/sum(distr)
  return(prop)
}

for(i in 1:5){
for(j in c(2)){
  
  qValue=seq(0, 0.6, 0.1)
  
  low=Ns[j]/Nu * 1/(D*N) * (1-1/(D*N))
  up=Ns[j]/(4*Nu)
  qValue=c(low-0.00001, 10^seq(log10(low)+0.1, log10(up)-0.1, 0.1), up)
  
  qProp=matrix(, ncol=length(qValue), nrow=10)
  
  for(r in 1:10){
    freq=read.table(paste(args[1], "/simu/frqB1_Nu0.01_Nm", Nm[i], "_Ns", Ns[j],"_N200_D100_L200_r",r, ".freq", sep=""))
    
    frq=c()
    if((Nm[i]<0.5 && Nm[i]>0) || (Nm[i]==0 && r<5))  lowGen=490000 else lowGen=90000
    for(gen in seq(lowGen, lowGen+10000, 5000)){
      frq1=as.matrix(freq[freq[,1]==gen,2:(L+1)])
      frq=c(frq, apply(frq1, 2, mean))
    }
    tmp=frq*(1-frq)*Ns[j]/Nu
    vgenic=tmp[tmp>0]
    
    for(x in 1:length(qValue)){
      qProp[r,x]=propQuantInv(vgenic, qValue[x])
    }
    
  }
  
  qPropMean=apply(qProp, 2, mean)
  qPropSD=apply(qProp, 2, std.error)*1.96
  
  points(log10(qValue), qPropMean,  col=colors[i], pch=1, cex=0.5)
  
  for(numQ in seq(1, length(qValue), 1)){
    lines(rep(log10(qValue[numQ]),2), c(qPropMean[numQ]+qPropSD[numQ], qPropMean[numQ]-qPropSD[numQ]),  col=colors[i])
  }
  
}
}
}

dev.off()