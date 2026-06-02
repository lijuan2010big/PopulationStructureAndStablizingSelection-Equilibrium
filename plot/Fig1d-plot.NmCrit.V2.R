library(plotrix)
library(RColorBrewer)

colors=brewer.pal(n=9,name = "Set1")[2:4]

weakSelNmCrit=function(Nu, Ns){
  ratio=Nu/Ns
  return(((3*ratio+sqrt((ratio)^2+2*(ratio)))/(2*(1-4*ratio))))
}


strongSelNmCrit=function(Nu, Ns){
  return(sqrt(Nu/(log(4 *Ns) + 0.57721)))
}

pdf("critNm_Approx_Numerical.V2.pdf", width=3, height=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.4, 0.5, 0), las=1, lwd=1, cex=1, xpd=FALSE)

#######################################
#### Plot n u/s
colors=brewer.pal(n=9,name = "Set1")[2:4]
Ns=c(0.1, 1, 5)

plot(NULL, xlim=log10(c(0.0005,0.25)), ylim=c(log10(0.01),log10(10)), 
     , ylab="", xlab="", 
     type="l", lty=1, xaxt="n", yaxt="n")
### y axis
axis(side=2, labels=c(0.01, 0.1, 1, 10), at=log10(c(0.01, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=2, at=log10(c(seq(0.01, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
mtext(expression(italic(Nm)[crit]), side=2, las=0, line=1.5)
### y axis
xlabs=c(0.001, 0.01, 0.05, 0.25)
axis(side=1, labels=xlabs, at=log10(xlabs))
mtext(expression(italic(mu/s)),side=1, las=1, line=1.5)

data1=read.table("NmCritUSweak.fixedNs.dat")
data2=read.table("NmCritUSstrong.fixedNs.dat")

for(i in 1:3){
  s=Ns[i]
  us=c(seq(0.0001,0.001,0.0001), seq(0.001,0.01,0.001), seq(0.01,0.25,0.01))
  mutR=us*s
  
  if(s<=1){
    lines(log10(us), log10(weakSelNmCrit(mutR, s)), lty=1, type = "l", col="black", lwd=1.5)
  }
  else{
    lines(log10(us), log10(strongSelNmCrit(mutR, s)),lty="21", type="l",  lwd=1.5, col=colors[3])
  }
  points(log10(data1[,1]), log10(data1[,(i+1)]), cex=0.7,  pch=i-1, lwd=1, col=colors[i])
  points(log10(data2[,1]), log10(data2[,(i+1)]), cex=0.7,  pch=i-1, lwd=1, col=colors[i])
}

legend("topleft", inset = c(0,0),  title="Approximation", legend=c("Weak sel.", "Strong sel."), lty=c("solid", "21"), col=c("black", colors[3]), bty="n", cex=0.7,pt.cex=0.5,lwd=1.5, seg.len = 1)
legend("topleft",inset = c(0.05,0.25), title=expression(italic("Ns")), legend=Ns, pch=c(0,1,2), col=colors, cex=0.7,  bty="n")


dev.off()