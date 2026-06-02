library(RColorBrewer)
library(fields)
library(plotrix)

args=commandArgs(trailingOnly = TRUE)

Ns=c(0.1, 1, 5)

Nm=c(0, 0.25, 2, 10)
colors=rev(brewer.pal(n = 11, name = "RdBu")[c(2, 3, 9, 11)])

critNm=c(0.62366087075, 0.088765207618, 0.046083603)

pdf("stabSel_Fv_Ns1.pdf", width=3, height=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)

for(i in c(2)){
  s=Ns[i]
  xupper=s/0.04
  if(xupper>25){
    xupper=25
  }
  plot(NULL, xlim=c(0,xupper), ylim=c(0,1), ylab="", xlab="", xaxt="n",yaxt="n", 
       main=substitute(italic("Ns") * " = " * s , list(s=s)))
  posLab=seq(0,1, 0.2)
  axis(1, at=posLab*xupper, labels = posLab*xupper, cex.axis=0.8)
  posLab=seq(0,1,0.2)
  axis(2, at=posLab, labels = posLab)
  
  mtext(expression(italic(v)*" (in units of "*italic("v")["*"]*")"), side=1, line = 1.4, cex=1.2)
  mtext(expression("F("*italic(v)*")"), side=2, line=2, las=0, cex=1.2)
  
  # Analytical prediciton
  j=1
  for(m in Nm){
    data=read.table(paste(args[1], "/N200/Fy_Nu0.01_Nm", m, "_Ns", s, ".eff.dat", sep=""))
    lines(data[,1], data[,2], col=colors[j])
    data=read.table(paste(args[1], "/N200/Fy_Nu0.01_Nm", m, "_Ns", s, ".dat", sep=""))
    lines(data[,1], data[,2], col=colors[j], lty="11")
    j=j+1
  }
  
  if(s==1){
    legend("topright", inset=c(0.0,0.2), legend=c(0, "0.089(crit)", "1/4", 2, 10), cex=0.7, text.col=c(colors[1], "grey50", colors[2:4]), title=expression(italic("Nm")), title.col = "black", title.adj=0.4, bty="n")
    legend("topright", inset=c(0.0,0), legend=c("LE", "LD"), cex=0.7, lty=c("11", "solid"), seg.len = 1, title="Approximation", bty="n", horiz = TRUE)  }
  
  data=read.table(paste(args[1], "/N200/Fy_Nu0.01_Nmcrt_Ns", s, ".dat", sep=""))
  lines(data[,1], data[,2], col="grey50", lty="11", lwd=1.2)
  
  data=read.table(paste(args[1], "/N200/Fy_Nu0.01_Nmcrt_Ns", s, ".eff.dat", sep=""))
  lines(data[,1], data[,2], col="grey50", lwd=1.2)
  i=i+1
}

Nu=0.01
N=200
L=200
D=100
quantInv <- function(distr, value) ecdf(distr)(value)

Nm=c(Nm,0.089)
colors=c(colors, "grey50")

for(i in 1:5){

for(j in c(2)){

low=Ns[j]/Nu *(1/(N*2))*(1-1/(N*2))
up=Ns[j]/(4*Nu)

qValue=c(seq(low, up-1, 2), up)

qProb=matrix(, ncol=length(qValue), nrow=10)

for(r in 1:10){
  freq=read.table(paste(args[1], "/simu/frqB1_Nu0.01_Nm",Nm[i], "_Ns", Ns[j],"_N200_D100_L200_r",r, ".freq", sep=""))
  
  if((Nm[i]<0.5 && Nm[i]>0) || (Nm[i]==0 && r<5))  lowGen=490000 else lowGen=90000

  frq1=c(as.matrix(freq[freq[,1]==lowGen,2:(L+1)]))

  for(gen in seq(lowGen+5000, lowGen+10000, 5000)){
    frq1=c(frq1,as.matrix(freq[freq[,1]==gen,2:(L+1)]))
  }
  
  vgenic=Ns[j]/Nu*frq1*(1-frq1)
  
  qProb[r,]=quantInv(vgenic, qValue-0.000001)
}

qProbMean=apply(qProb, 2, mean)
qProbSD=apply(qProb, 2, std.error)*1.96

points(qValue, 1-qProbMean,  col=colors[i], pch=1, cex=0.5)

for(numQ in seq(1, length(qValue), 1)){
  lines(rep(qValue[numQ],2), c(1-qProbMean[numQ]+qProbSD[numQ], 1-qProbMean[numQ]-qProbSD[numQ]),  col=colors[i])
}


}
}

if(FALSE){
##########In the full population  
plot(NULL, xlim=log10(c(0.049,125)), ylim=c(0, 1), ylab="", xlab="", xaxt="n",yaxt="n",
     main=substitute("In full population, "*italic("Nm") * " = " * m, list(m=round(m,3))))
posLab=c(0.1, 1, 5,  25, 125)
axis(1, at=log10(posLab), labels = posLab)
posLab=seq(0,1,0.2)
axis(2, at=posLab, labels = posLab)

mtext(expression(italic(v)), side=1, line = 1.4, cex=1.2)
mtext(expression("F("*italic(v)*")"), side=2, line=2, las=0, cex=1.2)

for(i in 2:5){
  
  for(j in c(2)){
  #expData=read.table(paste(args[2], "/Fy_Nu0.01_Nm", m, "_Ns", Ns[j], ".dat", sep=""))
  #lines(log10(expData[,1]), expData[,2], col=colors[j])
  
  low=Ns[j]/Nu *(1/(D*N*2))*(1-1/(D*N*2))
  up=Ns[j]/(4*Nu)
  
  print(round(log10(low*2)))
  
  qValue=c(10^seq(log10(low), log10(up)-0.1, 0.1), up)
  
  qProb=matrix(, ncol=length(qValue), nrow=10)
  
  print("start")
  for(r in 1:10){
    freq=read.table(paste(args[1], "/simu/frqB1_Nu0.01_Nm", Nm[i],"_Ns", Ns[j],"_N200_D100_L200_r",r, ".freq", sep=""))
    
    pmean=c()
    
    if((Nm[i]<0.5 && Nm[i]>0) || (Nm[i]==0 && r<5)) lowGen=490000 else lowGen=90000

    for(gen in seq(lowGen, lowGen+10000, 5000)){
      frq1=as.matrix(freq[freq[,1]==gen,2:(L+1)])
      pmean=c(pmean, apply(frq1, 2, mean))
    }
    
    vgenic=Ns[j]/Nu*pmean*(1-pmean)

    qProb[r,]=quantInv(vgenic, qValue-0.000001)
  }
  
  qProbMean=apply(qProb, 2, mean)
  qProbSD=apply(qProb, 2, std.error)*1.96
  
  
  points(log10(qValue), 1-qProbMean,  col=colors[i], pch=1, cex=0.5)
  
  for(numQ in seq(1, length(qValue), 1)){
    lines(rep(log10(qValue[numQ]),2), c(1-qProbMean[numQ]+qProbSD[numQ], 1-qProbMean[numQ]-qProbSD[numQ]),  col=colors[i])
  }
  }
}
}



dev.off()
