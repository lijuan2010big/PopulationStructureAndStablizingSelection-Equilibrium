library(plotrix)
library(RColorBrewer)

args=commandArgs(trailingOnly = TRUE)

Ns=c(0.1, 1, 5)
colors=brewer.pal(n=9,name = "Set1")[2:4]

critNm=c(0.6236608707517198, 0.08892506364329157, 0.046083603164263603)

pdf("pmean.V1.pdf", width=3, height=3)

par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1, xpd=FALSE)


### Plot Pmean as a function of Nm

plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(0,1), xlab=expression(italic(Nm)), ylab=expression(italic(bar(p))), 
     type="l", lty=1, main="Full population", xaxt="n")
axis(side=1, labels=c(0, 0.02, 0.1, 1, 10), at=log10(c(0.005, 0.02, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.03, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
axis.break(axis=1, breakpos=log10(0.012)) 

legend("topleft", legend=Ns, text.col=c(colors), pch=c(NA, NA, NA), title = expression( italic("Ns") ), title.col="black", title.adj = 0.7,cex=0.8, bty="n")


for(i in 1:3){
  pmean=read.table(paste("Nu0.01_Ns", Ns[i], "_Nm_pmean.dat", sep=""))
  lines(log10(pmean[pmean[,1]>=0.011,1]), pmean[pmean[,1]>=0.011,2], lty=1, col=colors[i], lwd=2)    
  lines(log10(pmean[pmean[,1]>=0.011,1]), 1-pmean[pmean[,1]>=0.011,2], lty=1, col=colors[i], lwd=2)   
  points(log10(critNm[i]), 0.5, pch=8, col=colors[i], cex=0.7)
  options(digits=10)
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), pmean[1:8, 2], lty=1, col=colors[i], lwd=2)
}



#### Plot Pmean with estimation

plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(0,1), xlab=expression(italic(Nm)), ylab=expression(italic(bar(p))), 
     type="l", lty=1, main="Full population", xaxt="n")
axis(side=1, labels=c(0, 0.02, 0.1, 1, 10), at=log10(c(0.005, 0.02, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.03, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)

for(i in 1:3){
  pmean=read.table(paste("Nu0.01_Ns", Ns[i], "_Nm_pmean.dat", sep=""))
  lines(log10(pmean[pmean[,1]>=0.011,1]), pmean[pmean[,1]>=0.011,2], lty=1, col=colors[i], lwd=1)    
  lines(log10(pmean[pmean[,1]>=0.011,1]), 1-pmean[pmean[,1]>=0.011,2], lty=1, col=colors[i], lwd=1)   
  points(log10(critNm[i]), 0.5, pch=8, col=colors[i], cex=0.7)
  options(digits=10)
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), pmean[1:8, 2], lty=1, col=colors[i], lwd=1)
}

axis.break(axis=1, breakpos=log10(0.012)) 

Nu=0.01

pdif=function(Nu, Ns, Nm){
  return(sqrt(1-Nu/Ns*(1+2*Nm)*(1+4*Nm)/(2*Nm^2)))
}
lines(log10(seq(critNm[1],10, 0.01)), 1/2*(1-pdif(Nu, Ns[1], seq(critNm[1], 10, 0.01))), lty="11", col=colors[1], lwd=3)
lines(log10(seq(critNm[1],10, 0.01)), 1/2*(1+pdif(Nu, Ns[1], seq(critNm[1], 10, 0.01))), lty="11", col=colors[1], lwd=3)


lines(log10(seq(critNm[2], 10, 0.01)), 1/2*(1-pdif(Nu, Ns[2], seq(critNm[2], 10, 0.01))), lty="11", lwd=3, col=colors[2])
lines(log10(seq(critNm[2], 10, 0.01)), 1/2*(1+pdif(Nu, Ns[2], seq(critNm[2], 10, 0.01))), lty="11", lwd=3, col=colors[2])


pbar=function(Nu, Ns, Nm){
  return((Nu/(Ns-Nm))*((4*Ns)^(4*Nm)+4*(Ns-Nm)*gamma(4*Nm))/(((4*Ns)^(4*Nm)-4*Nm)*gamma(4*Nm)))
}

pbar(0.01, 5, 1)
lines(log10(seq(critNm[3], 10, 0.01)), pbar(Nu, Ns[3], seq(critNm[3], 10, 0.01)), lty="21", lwd=3, col=colors[3])
lines(log10(seq(critNm[3], 10, 0.01)), 1-pbar(Nu, Ns[3], seq(critNm[3], 10, 0.01)), lty="21", lwd=3, col=colors[3])


#### Plot Fst as a function of Nm

plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(0,1), xlab=expression(italic(Nm)), ylab=expression(italic(F[ST])), 
     type="l", lty=1, main="", xaxt="n")
axis(side=1, labels=c(0, 0.02, 0.1, 1, 10), at=log10(c(0.005, 0.02, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.03, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)

for(i in 1:3){
  pmean=read.table(paste("Nu0.01_Ns", Ns[i], "_Nm_pmean.dat", sep=""))
  HW=read.table(paste("Nu0.01_Ns", Ns[i], "_Nm_HW.dat", sep=""))
  HT=2*pmean[,2]*(1-pmean[,2])
  FST=(HT-HW[,2])/HT
  print(FST)
  lines(log10(pmean[pmean[,1]>=0.015,1]), FST[pmean[,1]>=0.015], lty=1, col=colors[i], lwd=1)    
  options(digits=10)
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), FST[1:8], lty=1, col=colors[i], lwd=1)
}

axis.break(axis=1, breakpos=log10(0.012)) 

Nu=0.01


weakFst1=function(Nu, Ns, Nm){
  return((1/(1+4*Nm))*(1+8*Nm*Ns/((1+4*Nm)*(3+4*Nm))-8*Nu/(1+4*Nm)))
}
weakFst2=function(Nu, Ns, Nm){
  return((1/(1+4*Nm))*(1+8*Nm*(1-4*Nm)*Ns/((1+2*Nm)*(1+4*Nm)*(3+4*Nm))+64*Nm*Nu/((1+4*Nm)*(3+4*Nm))))
}

Ns=0.1
weakFst1(Nu, Ns, 1)

lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), weakFst1(Nu, Ns, seq(0, 0.007, 0.001)), lty="11", col=colors[1], lwd=3)
lines(log10(seq(0.015, critNm[1], 0.01)), weakFst1(Nu, Ns, seq(0.015, critNm[1], 0.01)), lty="11", col=colors[1], lwd=3)
lines(log10(seq(critNm[1],10, 0.01)), weakFst2(Nu, Ns, seq(critNm[1],10, 0.01)), lty="11", col=colors[1], lwd=3)

Ns=1
lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), weakFst1(Nu, Ns, seq(0, 0.007, 0.001)), lty="11", col=colors[2], lwd=3)
lines(log10(seq(0.015, critNm[2], 0.01)), weakFst1(Nu, Ns, seq(0.015, critNm[2], 0.01)), lty="11", col=colors[2], lwd=3)
lines(log10(seq(critNm[2],10, 0.01)), weakFst2(Nu, Ns, seq(critNm[2],10, 0.01)), lty="11", col=colors[2], lwd=3)


strongFst1=function(Nu, Ns, Nm){
  return(1 + ((Nm + 2*Nu) *(1 + 2*Nm - 4*Ns + 4*Nu))/(2*Ns^2))
}
strongFst2=function(Nu, Ns, Nm){
  return(1/(4*Ns)*(((4*Ns)^(4*Nm)+4*(Nm+4*(Ns-Nm)^2)*gamma(4*Nm))/((4*Ns)^(4*Nm)+4*(Ns-Nm)*gamma(4*Nm))
  ))
}

Ns=5
lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), strongFst2(Nu, Ns, seq(0, 0.007, 0.001)),lty="21", col=colors[3], lwd=3)
lines(log10(seq(0.015, critNm[3], 0.01)), strongFst1(Nu, Ns, seq(0.015, critNm[3], 0.01)), lty="21", col=colors[3], lwd=3)
lines(log10(seq(critNm[3],10, 0.01)), strongFst2(Nu, Ns, seq(critNm[3],10, 0.01)), lty="21", col=colors[3], lwd=3)

legend("topright", legend=c("Weak sel.", "Strong sel.", "Diffusion"), title="Approx.", lwd=c(3,3,1), lty=c("11", "21", "solid"), seg.len = 2, cex=0.7, bty="n")

legend("topright", inset = c(0.1,0.3),  legend=NA, title = expression(italic(Ns)), pch=c(NA), bty="n", cex=0.8)
legend("topright", inset = c(0.1,0.35), legend=c(0.1, 1, 5),  text.col = colors, pch=c(NA, NA, NA), bty="n", cex=0.8)
 

dev.off()