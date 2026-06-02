library(RColorBrewer)
library(plotrix)


colors=brewer.pal(n=9,name = "Set1")[2:3]

colors1=col2rgb(colors, alpha = FALSE)/255

pmean=c(0.5, 0.168)
het_meta=2*pmean*(1-pmean)
het_deme=c(0.33, 0.25)

Nm=c(0.5, 2)

args=commandArgs(trailingOnly = TRUE)

numDeme=c(25, 50, 100, 250, 500)
demeSize=c(2000, 1000, 500, 200, 100)
fold=c(5,1)
Ns=c(2,1,0.5,0.2,0.1)
Nu=Ns/10

pchtype=c(1,0)

het1=c(0.331098, 0.332908, 0.33322, 0.333304, 0.333329, 0.333332, 0.333329)
het5=c(0.294011, 0.32629,0.331098,0.332688,0.33322, 0.333304,0.33322)

intfrq=0.5

pdf("underDomi_numDeme_NuD.pdf", width=3, height=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)

for(m in c(1)){

plot(NULL, xlim=c(log(20),log(550)), ylim=c(log(0.1), log(0.5)), xlab="Number of demes", ylab="", xaxt="n", yaxt="n",
     main=expression("Fixed "* italic("N")*mu*italic("D")))
labnum=c(25, 50, 100, 250, 500)
axis(side=1, at=log(labnum), labels = labnum,  cex.axis=1)
ylabnum=c(0.5,0.4,0.3, 0.2, 0.1, 0.05, 0.02)
axis(side=2, at=log(ylabnum), labels = ylabnum,  cex.axis=1)
mtext("Heterozygosity", side=2, line=2, las=0)


text(log(20), log(0.1),adj=c(0,0), cex=0.7,labels= expression(mu*"="*10^-4*", "*italic(s)*"="*10^-3*", "*italic(Nm)*"=0.5"), bty="n")


  lines(c(0,log(600)), c(log(het_meta[m]), log(het_meta[m])), lty=1, col=colors[1])
  lines(log(c(10, numDeme, 1000)), log(het1),,lty="21", col=rgb(colors1[1,1], colors1[2,1],colors1[3,1],alpha=0.5))
  lines(log(c(10, numDeme, 1000)), log(het5),,lty="solid", col=rgb(colors1[1,1], colors1[2,1],colors1[3,1],alpha=0.5))
  
  points(log(c(10, numDeme, 1000)), log(het1),pch=19, cex=0.2,col=rgb(colors1[1,1], colors1[2,1],colors1[3,1],alpha=0.5))
  points(log(c(10, numDeme, 1000)), log(het5),pch=19, cex=0.2,col=rgb(colors1[1,1], colors1[2,1],colors1[3,1],alpha=0.5))
  
  
  for(f in 1:2){ 
  for(i in 1:5){
    d=numDeme[i]
    HT=vector()
    HS=vector()
    for(r in 1:10){
      fileData=read.table(paste(args[1], "/frq0.5_Nm", Nm[m], "_D", numDeme[i], "_N", demeSize[i]/fold[f], "_Ns", Ns[i]/fold[f], "_Nu", Nu[i]/fold[f], "_r", r, sep=""), header = FALSE)
      l=length(fileData[,1])
      x=fileData[(l-100):l, 2:(numDeme[i]+1)]
      HetMeta=apply(x, 1, mean)
      HT_r=2*HetMeta*(1-HetMeta)
      HT_r_mean=mean(HT_r)
      HT=c(HT, HT_r_mean)
      HS_r=apply(2*x*(1-x), 1, mean)
      HS_r_mean=mean(HS_r)
      HS=c(HS, HS_r_mean)
    }

  HT_mean=mean(HT)
  HT_sd=std.error(HT)*1.96
  HS_mean=mean(HS)
  HS_sd=std.error(HS)*1.96
  
  points(log(d)+0.03*(f-1.5), log(HT_mean), cex=0.8,pch=pchtype[f], col=colors[1])
  lines(log(c(d,d))+0.03*(f-1.5), c(log(HT_mean-HT_sd), log(HT_mean+HT_sd)), col=colors[1])
    

  
  points(log(d)+0.03*(f-1.5), log(HS_mean), pch=pchtype[f], cex=0.8,col=rgb(colors1[1,1], colors1[2,1],colors1[3,1],alpha=0.5))
  lines(log(c(d,d))+0.03*(f-1.5), c(log(HS_mean-HS_sd), log(HS_mean+HS_sd)), col=rgb(colors1[1,1], colors1[2,1],colors1[3,1],alpha=0.5))

  }
}
}


legend("bottomleft", inset = c(-0.07, 0.1), legend=expression(italic("N")*mu*italic("D")),pch=NA, cex=0.8, bty="n", horiz = TRUE)
legend("bottomleft", inset = c(0.25, 0.1), title="", legend=c(1, 5), pch=c(pchtype), 
       lty=c("21", "solid"), cex=0.8, seg.len = 2, bty="n", horiz = TRUE)

legend("bottomleft", inset = c(0, 0.2),legend=c(expression(italic(H[T])), expression(italic(H[W]))), pch=c(NA, NA), 
       text.col=c(colors[1], rgb(colors1[1,1], colors1[2,1],colors1[3,1],alpha=0.5), colors[2], rgb(colors1[1,2], colors1[2,2],colors1[3,2],alpha=0.5)), 
       cex=0.8, bty="n", horiz = TRUE)


dev.off()