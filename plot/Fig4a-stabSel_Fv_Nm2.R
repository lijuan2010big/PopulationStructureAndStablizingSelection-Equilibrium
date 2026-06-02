library(RColorBrewer)
library(fields)
library(plotrix)

args=commandArgs(trailingOnly = TRUE)


Ns=c(0.2, 1, 5)
Nm=c(0, 0.1, 2)
L=200
Nu=0.01
D=100
N=200

colors=brewer.pal(n=9,name = "Set1")[2:4]

quantInv <- function(distr, value) ecdf(distr)(value)

pdf("stabSel_Fv_Nm2.pdf", width=3, height=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)


m=2
  
  plot(NULL, xlim=log10(c(0.049,125)), ylim=c(0, 1), ylab="", xlab="", xaxt="n",yaxt="n",
       main=substitute(italic("Nm") * " = " * m, list(m=round(m,3))))
  posLab=c(0.1, 1, 5,  25, 125)
  axis(1, at=log10(posLab), labels = posLab)
  posLab=seq(0,1,0.2)
  axis(2, at=posLab, labels = posLab)
  
  mtext(expression(italic(v)*" (in units of "*italic("v")["*"]*")"), side=1, line = 1.4, cex=1.2)
  mtext(expression("F("*italic(v)*")"), side=2, line=2, las=0, cex=1.2)
  
  if(m>0){
    legend("topright", legend=Ns, text.col=c(colors), title = expression( italic("Ns")), title.col = "black", title.adj = 0.7, seg.len = 1, cex=0.8, bty="n")
  }
  

  for(j in 1:3){
    expData=read.table(paste(args[1], "/N200/Fy_Nu0.01_Nm", m, "_Ns", Ns[j], ".eff.dat", sep=""))
    lines(log10(expData[,1]), expData[,2], col=colors[j])
    
    expData=read.table(paste(args[1], "/N200/Fy_Nu0.01_Nm", m, "_Ns", Ns[j], ".dat", sep=""))
    lines(log10(expData[,1]), expData[,2], col=colors[j], lty="11")

    
    low=Ns[j]/Nu *(1/(N*2))*(1-1/(N*2))
    up=Ns[j]/(4*Nu)
    
    print(round(log10(low*2)))

    qValue=c(10^seq(log10(low), log10(up)-0.1, 0.2), up)

    qProb=matrix(, ncol=length(qValue), nrow=10)
    
    for(r in 1:10){
    freq=read.table(paste(args[1], "/simu/frqB1_Nu0.01_Nm2_Ns", Ns[j],"_N200_D100_L200_r",r, ".freq", sep=""))
    
    frq1=c(as.matrix(freq[freq[,1]==90000,2:(L+1)]))
    
    for(gen in seq(95000, 100000, 5000)){
      frq1=c(frq1,as.matrix(freq[freq[,1]==gen,2:(L+1)]))
    }
    
    vgenic=Ns[j]/Nu*frq1*(1-frq1)

    qProb[r,]=quantInv(vgenic, qValue-0.000001)
    }
    
    qProbMean=apply(qProb, 2, mean)
    qProbSD=apply(qProb, 2, std.error)*1.96

    
    points(log10(qValue), 1-qProbMean,  col=colors[j], pch=1, cex=0.5)
    
    for(numQ in seq(1, length(qValue), 1)){
      lines(rep(log10(qValue[numQ]),2), c(1-qProbMean[numQ]+qProbSD[numQ], 1-qProbMean[numQ]-qProbSD[numQ]),  col=colors[j])
    }
    
  
    
  }

if(FALSE){
  ##########In the full population  
  
  m=2
  
  plot(NULL, xlim=log10(c(0.049,125)), ylim=c(0, 1), ylab="", xlab="", xaxt="n",yaxt="n",
       main=substitute("Full population, "*italic("Nm") * " = " * m, list(m=round(m,3))))
  posLab=c(0.1, 1, 5,  25, 125)
  axis(1, at=log10(posLab), labels = posLab)
  posLab=seq(0,1,0.2)
  axis(2, at=posLab, labels = posLab)
  
  mtext(expression(italic(v)*" (in units of "*italic("v")["*"]*")"), side=1, line = 1.4, cex=1.2)
  mtext(expression("F("*italic(v)*")"), side=2, line=2, las=0, cex=1.2)

  
  if(m>0){
    legend("topright", legend=Ns, col=c(colors), lty=c(rep("solid",3)), title = expression( italic("Ns")), seg.len = 1, cex=0.8, bty="n")
  }
  
  
  for(j in 1:3){
    low=Ns[j]/Nu *(1/(D*N*2))*(1-1/(D*N*2))
    up=Ns[j]/(4*Nu)
    
    print(round(log10(low*2)))
    
    qValue=c(10^seq(log10(low), log10(up)-0.1, 0.1), up)
    
    qProb=matrix(, ncol=length(qValue), nrow=10)
    
    print("start")
    for(r in 1:10){
      freq=read.table(paste(args[1], "/Simu/frqB1_Nu0.01_Nm2_Ns", Ns[j],"_N200_D100_L200_r",r, ".freq", sep=""))
      
      pmean=c()
      if(m<0.5 && m>0)  lowGen=490000 else lowGen=90000
      for(gen in seq(lowGen, lowGen+10000, 5000)){
        frq1=as.matrix(freq[freq[,1]==gen,2:(L+1)])
        pmean=c(pmean, apply(frq1, 2, mean))
      }
      
      vgenic=Ns[j]/Nu*pmean*(1-pmean)
      
      qProb[r,]=quantInv(vgenic, qValue-0.000001)
    }
    
    qProbMean=apply(qProb, 2, mean)
    qProbSD=apply(qProb, 2, std.error)*1.96
    
    
    points(log10(qValue), 1-qProbMean,  col=colors[j], pch=1, cex=0.5)
    
    for(numQ in seq(1, length(qValue), 1)){
      lines(rep(log10(qValue[numQ]),2), c(1-qProbMean[numQ]+qProbSD[numQ], 1-qProbMean[numQ]-qProbSD[numQ]),  col=colors[j])
    }
  }
}


dev.off()