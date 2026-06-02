library(RColorBrewer)
library(fields)
library(plotrix)

args=commandArgs(trailingOnly = TRUE)


Ns=c(0.2, 1, 5)
Nm=c(0.1, 2)
L=200
Nu=0.01
D=100
N=200

colors=brewer.pal(n=9,name = "Set1")[2:4]

pdf(paste(args[1],"/stabSel_propVg_Nm2.pdf", sep=""), width=3, height=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)

for(i in c(2)){
  m=Nm[i]
  plot(NULL, xlim=log10(c(0.049,125)), ylim=c(0,1), 
       ylab="", , xlab="", xaxt="n", yaxt="n",
       main=substitute(italic("Nm")*"="*x, list(x=m)))
  posLab=c(0.1, 1, 5,  25, 125)
  axis(1, at=log10(posLab), labels = posLab, cex.axis=1)
  mtext(expression(italic(v)), side=1, line = 1.4, cex=1.2)
  posLab=seq(0,1,0.2)
  axis(2, at=posLab, labels = posLab)
  mtext(expression("Proportion of total "*italic(V[W]^"(gen)")), side=2, line=1.8, las=0, cex=1.2)
  
  # Analytical prediction
  for(j in 1:3){
    s=Ns[j]
    data=read.table(paste(args[1], "/N200/icdfVg_Nu0.01_Nm", Nm[i], "_Ns",Ns[j],".eff.dat", sep=""))
    lines(log10(data[,1]), data[,2]/data[1,2], col=colors[j])
    data=read.table(paste(args[1], "/N200/icdfVg_Nu0.01_Nm", Nm[i], "_Ns",Ns[j],".dat", sep=""))
    lines(log10(data[,1]), data[,2]/data[1,2], col=colors[j], lty="11")
  }
  
}

propQuantInv <- function(distr, value) {
  prop=sum(distr[distr>=value])/sum(distr)
  return(prop)
  }

  for(j in 1:3){

    qValue=seq(0, 0.6, 0.1)
    
    low=Ns[j]/Nu * 1/400 * (1-1/400)
    up=Ns[j]/(4*Nu)
    qValue=c(low-0.00001, 10^seq(log10(low)+0.2, log10(up)-0.1, 0.2), up)
    

    qProp=matrix(, ncol=length(qValue), nrow=10)
    
    for(r in 1:10){
    freq=read.table(paste(args[1], "/simu/frqB1_Nu0.01_Nm2_Ns", Ns[j],"_N200_D100_L200_r",r, ".freq", sep=""))
    
    qPropR=matrix(, ncol=length(qValue), nrow=3)
    g=1
    if(Nm[i]<0.5 && Nm[i]>0)  lowGen=490000 else lowGen=90000
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

    points(log10(qValue), qPropMean,  col=colors[j], pch=1, cex=0.5)
    
    for(numQ in seq(1, length(qValue), 1)){
      lines(rep(log10(qValue[numQ]),2), c(qPropMean[numQ]+qPropSD[numQ], qPropMean[numQ]-qPropSD[numQ]),  col=colors[j])
    }
    
  }


if(FALSE){
#### VT

m=2
plot(NULL, xlim=log10(c(0.049,125)), ylim=c(0,1), ylab="", , xlab="", xaxt="n", yaxt="n",
     main=substitute("Full population, "*italic("Nm")*"="*x, list(x=m)))
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

for(j in 1:3){
  
  qValue=seq(0, 0.6, 0.1)
  
  low=Ns[j]/Nu * 1/(D*N) * (1-1/(D*N))
  up=Ns[j]/(4*Nu)
  qValue=c(low-0.00001, 10^seq(log10(low)+0.1, log10(up)-0.1, 0.1), up)
  
  qProp=matrix(, ncol=length(qValue), nrow=4)
  
  for(r in 1:4){
    freq=read.table(paste(args[1], "/simu/frqB1_Nu0.01_Nm2_Ns", Ns[j],"_N200_D100_L200_r",r, ".freq", sep=""))
    
    frq=c()
    if(Nm[i]<0.5 && Nm[i]>0)  lowGen=490000 else lowGen=90000
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
  
  points(log10(qValue), qPropMean,  col=colors[j], pch=1, cex=0.5)
  
  for(numQ in seq(1, length(qValue), 1)){
    lines(rep(log10(qValue[numQ]),2), c(qPropMean[numQ]+qPropSD[numQ], qPropMean[numQ]-qPropSD[numQ]),  col=colors[j])
  }
  
}
}

dev.off()