library(RColorBrewer)
library(fields)
library(plotrix)

args=commandArgs(trailingOnly = TRUE)


Ns=c(0.1, 0.2, 1, 5)
Nm=c(0, 0.1, 2)
L=200
Nu=0.01
D=100
N=200



colors=brewer.pal(n=9,name = "Set1")[2:4]

quantInv <- function(distr, value) ecdf(distr)(value)

pdf("stabSel_Fv.metapop.numDeme.eff.pdf", width=3, height=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)

D=c(50, 100, 200)
  m=2
  
  pbar=c(0.206794, 0.089373, 0.0156406, 0.00264005)
  
  for(j in c(2,3,4)){

    
  plot(NULL, xlim=log10(c(0.1,25)), ylim=c(0, 1), ylab="", xlab="", xaxt="n",yaxt="n",
       main=substitute("Full population, "*italic("Ns") * " = " * x, list(x=Ns[j])))
  posLab=c(0.1, 1, 5,  25)
  axis(1, at=log10(posLab), labels = posLab, cex.axis=1)
  posLab=seq(0,1,0.2)
  axis(2, at=posLab, labels = posLab)
  
  if(j==1){
    mtext(expression(italic(v)*" (in units of "*italic("v")["*"]*")"), side=1, line = 1.6, cex=1.2)
  }
  else{
  mtext(expression(italic(v)), side=1, line = 1.4, cex=1.2)}
  mtext(expression("F("*italic(v)*")"), side=2, line=2, las=0, cex=1.2)
  
  print(Ns[j]/Nu*pbar[j]*(1-pbar[j]))
  abline(v=log10(Ns[j]/Nu*pbar[j]*(1-pbar[j])), col="grey50", lwd=1)
  
  if(j==2){
    legend("topright", legend=D, col=c(colors), lty=c(rep("solid",3)), title = expression( italic("D")), seg.len = 1, cex=0.8, bty="n")
  }
  


  for(d in 1:3){
    low=Ns[j]/Nu *(1/(D[d]*N*2))*(1-1/(D[d]*N*2))
    up=Ns[j]/(4*Nu)
    
    print(round(log10(low*2)))
    
    qValue=c(10^seq(log10(low), log10(up)-0.1, 0.1), up)
    
    qProb=matrix(, ncol=length(qValue), nrow=4)
    
    print("start")
    for(r in 1:4){
      freq=read.table(paste(args[1], "/frqB_Nu0.01_Nm2_Ns", Ns[j],"_N200_D", D[d],"_L200_r",r, ".freq", sep=""))
      
      pmean=c()
      
      for(gen in seq(40000, 50000, 2000)){
        frq1=as.matrix(freq[freq[,1]==gen,2:(L+1)])
        pmean=c(pmean, apply(frq1, 2, mean))
      }
      
      vgenic=Ns[j]/Nu*pmean*(1-pmean)
      
      qProb[r,]=quantInv(vgenic, qValue-0.000001)
    }
    
    qProbMean=apply(qProb, 2, mean)
    qProbSD=apply(qProb, 2, std.error)*1.96
    
    
    points(log10(qValue), 1-qProbMean,  col=colors[d], pch=1, cex=0.5)
    
    for(numQ in seq(1, length(qValue), 1)){
      lines(rep(log10(qValue[numQ]),2), c(1-qProbMean[numQ]+qProbSD[numQ], 1-qProbMean[numQ]-qProbSD[numQ]),  col=colors[d])
    }
  }

}

dev.off()