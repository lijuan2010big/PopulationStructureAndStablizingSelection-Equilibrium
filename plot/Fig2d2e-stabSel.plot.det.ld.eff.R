library(RColorBrewer)
library(plotrix)
options(digits=10)

args=commandArgs(trailingOnly = TRUE)

colors=brewer.pal(n=9,name = "Set1")[2:4]

Nm=c(0, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 2 , 5, 10)
Ns=c(0.1, 1, 5)
Vs=c(10, 1, 0.2)
Nu=0.01
a=0.1
L=200
N=200

tmpInit=matrix(, nrow=length(Nm), ncol=3)
VTr_mean=tmpInit
VTr_sd=tmpInit
VWr_mean=tmpInit
VWr_sd=tmpInit

for(i in 1:length(Nm)){
  m=Nm[i]
  for(j in 1:3){
    s=Ns[j]
    VTr=vector()
    VWr=vector()
    for(r in seq(1, 4, 1)){
      #### simulation het_trait_frqB_Nu0.01_Nm1_Ns1_N100_D100_L100_r1
      data=read.table(paste(args[1], "/simu/het_trait_frqB1_Nu0.01_Nm", Nm[i], "_Ns", Ns[j], "_N200_D100_L200_r", r, sep=""), header = TRUE)
      for(r in 1:4){
        l=length(data[,1])
        Vg=data[seq((l-50),l,25), 2]*a^2*L
        traitVar=data[seq((l-50),l,25),4]
        VTr=c(VTr,traitVar/Vg)
        
        Vg=data[seq((l-50),l,25), seq(5, length(data[1,]), 3)]*a^2*L
        traitVarDeme=data[seq((l-50),l,25),seq(7, length(data[1,]), 3)]
        traitVarMean=mean(as.matrix(traitVarDeme))
        VgMean=mean(as.matrix(Vg))
        VWr=c(VWr, traitVarMean/VgMean)
      }
    }
    
    VTr_mean[i,j]=mean(VTr)
    VTr_sd[i,j]=1.96*std.error(VTr)
    VWr_mean[i,j]=mean(VWr)
    VWr_sd[i,j]=1.96*std.error(VWr)
  }
}


propV=function(Vg, Vs){
  x=Vs/Vg
  return(-(3+x-sqrt(1+6*x+x^2))/4)
}

geneticVW=function(Vg, Vs){
  return((Vg-Vs+sqrt(Vg^2+6*Vg*Vs+Vs^2))/4)
}

pdf("stabSel_LD_det.eff.pdf", height=3, width=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)
lNm=0.005
### metapopulation

plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(0,1), xlab="", ylab="", xaxt="n",
     cex.axis=0.8, main="Full population")
### set up y axis
mtext(expression(italic(V[T]/V[T]^{(gen)})), side=2, line=2, las=0, cex=1.2)
### set up x axis
axis(side=1, labels=c(0, 0.02, 0.1, 1, 10), at=log10(c(0.005, 0.02, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.03, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
mtext(expression(italic(Nm)), side=1, line = 1.6, cex=1.2)
### set up break x axis
axis.break(axis=1, breakpos=log10(0.012)) 

legend("topleft", title="Approximation", legend=c("QG", "QG(eff)"), lty=c("11", "solid"), seg.len = 1, cex=0.8, bty="n", horiz = FALSE)


# simulation
for(j in 1:3){
  for(i in 1:length(Nm)){
    m=if(Nm[i]==0) lNm else Nm[i]
    d=VTr_mean[i,j]-VTr_sd[i,j]
    u=VTr_mean[i,j]+VTr_sd[i,j]
    points(log10(m),VTr_mean[i,j], col=colors[j], pch=1, cex=0.7)
    lines(log10(rep(m, 2)), c(d, u), col=colors[j], pch=1, cex=0.7)
  }
}

# deterministic
for(j in 1:3){
  pmeanHW=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_pmean_HW.dat", sep=""))
  m=pmeanHW[,1]
  VW=geneticVW(pmeanHW[,3]*a^2*L,Vs[j])
  HT=2*pmeanHW[,2]*(1-pmeanHW[,2])
  VT=HT*a^2*L
  VB=VW/(2*m+2*N*(VW/Vs[j]))
  
  ratio=(VW+VB)/VT
  lines(log10(m[m>0.015]), ratio[m>0.015], col=colors[j],  lty="11")
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), ratio[1:8], col=colors[j], lty="11")
  
  seff=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_NmEff_NsEff_pmean_HW.dat", sep=""))
  m=seff[,2]
  VW=geneticVW(seff[,5]*a^2*L,Vs[j])
  HT=2*seff[,4]*(1-seff[,4])
  VT=HT*a^2*L
  VB=VW/(2*m+2*N*(VW/Vs[j]))
  
  ratio=(VW+VB)/VT
  lines(log10(m[m>0.015]), ratio[m>0.015], col=colors[j], lwd=1)
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), ratio[1:8], col=colors[j], lwd=1)
}


### deme
plot(NULL, xlim=log10(c(0.005, 10)), ylim=c(0.78, 1), xlab="", ylab="", xaxt="n",cex.axis=0.8, main="Within demes")
### set up y axis
mtext(expression(italic(V[W]/V[W]^{(gen)})), side=2, line=2, las=0, cex=1.2)
### set up x axis
axis(side=1, labels=c(0, 0.02, 0.1, 1, 10), at=log10(c(0.005, 0.02, 0.1, 1, 10)), cex.axis=1, tck=-0.04)
axis(side=1, at=log10(c(seq(0.03, 0.1, 0.01), seq(0.1, 1, 0.1), seq(1,10,1))), labels=FALSE,  lwd.ticks=0.2, tck=-0.03)
mtext(expression(italic(Nm)), side=1, line = 1.6, cex=1.2)
### set up break x axis
axis.break(axis=1, breakpos=log10(0.012)) 


for(j in 1:3){
  for(i in 1:length(Nm)){
    m=if(Nm[i]==0) lNm else Nm[i]
    d=VWr_mean[i,j]-VWr_sd[i,j]
    u=VWr_mean[i,j]+VWr_sd[i,j]
    points(log10(m), VWr_mean[i,j], col=colors[j], pch=1, cex=0.7)
    lines(log10(c(m, m)), c(d, u), col=colors[j], pch=1, cex=0.7)
  }
}



for(j in 1:3){
  HW=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_pmean_HW.dat", sep=""))
  VW=HW[,3]*a^2*L
  lines(log10(HW[HW[,1]>0.015,1]), (propV(VW, Vs[j])+1)[HW[,1]>0.015], col=colors[j], lty="11" )
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), (propV(VW, Vs[j])+1)[1:8], col=colors[j], lty="11")
  
  seff=read.table(paste(args[1], "/predict/Nu0.01_Ns", Ns[j], "_Nm_NmEff_NsEff_pmean_HW.dat", sep=""))
  VW=seff[,5]*a^2*L
  lines(log10(seff[seff[,1]>0.015,1]), (propV(VW, Vs[j])+1)[seff[,1]>0.015], col=colors[j], lwd=1)
  lines(seq(log10(0.005),log10(0.01), by=(log10(0.01)-log10(0.005))/7), (propV(VW, Vs[j])+1)[1:8], col=colors[j],lwd=1 )
}


dev.off()
