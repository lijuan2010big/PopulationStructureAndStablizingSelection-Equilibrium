library(RColorBrewer)
library(plotrix)


colors=brewer.pal(n=9,name = "Set1")[2:3]

colors1=col2rgb(colors, alpha = FALSE)/255


het_meta=c(0.5, 0.0475866)
het_deme=c(0.3333039222620644, 0.033664)

args=commandArgs(trailingOnly = TRUE)

numDeme=c(25, 50, 100, 250, 500)

Ns=c(0.1, 1)

intfrq=0.5

lwidth=c(3, 1.5, 0.8)


pdf("underDomi_numDeme_Nu.pdf", width=3, height=3)
par(mfrow=c(1, 1), mai=c(0.5, 0.7, 0.3, 0.1), mgp=c(1.6, 0.5, 0), las=1, lwd=1.5, cex=1)

j=1
plot(NULL, xlim=c(log(20),log(500)), ylim=c(log(0.1), log(0.5)), xlab="Number of demes", ylab="", xaxt="n", yaxt="n" 
     , main=expression("Fixed"*italic(N)*mu))
labnum=c(25, 50, 100, 250, 500)
axis(side=1, at=log(labnum), labels = labnum,  cex.axis=1)
ylabnum=c(0.5,0.4,0.3, 0.2, 0.1, 0.05, 0.02)
axis(side=2, at=log(ylabnum), labels = ylabnum,  cex.axis=1)
mtext("Heterozygosity", side=2, line=2, las=0)

text(log(20), log(0.1),adj=c(0,0), cex=0.7,labels= expression(italic(N)*"=100, "*italic(N)*mu*"=0.01, "*italic(Ns)*"=0.1, "*italic(Nm)*"=0.5"), bty="n")

legend("bottomleft", inset = c(0, 0.1),legend=c(expression(italic(H[T])), expression(italic(H[W]))), pch=c(NA, NA), 
       text.col=c(colors[1], rgb(colors1[1,1], colors1[2,1],colors1[3,1],alpha=0.5), colors[2], rgb(colors1[1,2], colors1[2,2],colors1[3,2],alpha=0.5)), 
       cex=0.8, bty="n", horiz = TRUE)

for(i in c(1)){
  lines(c(0,log(600)), c(log(het_meta[i]), log(het_meta[i])), lty=1, col=colors[i])
  lines(c(0,log(600)), c(log(het_deme[i]), log(het_deme[i])), col=rgb(colors1[1,i], colors1[2,i],colors1[3,i],alpha=0.5), lty=1)
  for(d in numDeme){
    HT=vector()
    HS=vector()
    for(r in 1:10){
      #frq0.5_Nm0.5_D500_N100_Ns0.1_Nu0.01_r6
      fileData=read.table(paste(args[1], "/frq0.5_Nm0.5_D", d, "_N100_Ns0.1_Nu0.01_r", r, sep=""), header = FALSE)
      l=length(fileData[,1])
      x=fileData[(l-99):l, 2:(numDeme[i]+1)]
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
  print(c(HT_mean, HS_mean))
  points(log(d)-0.05, log(HT_mean), pch=19, col=colors[i])
  lines(log(c(d,d))-0.05, c(log(HT_mean-HT_sd), log(HT_mean+HT_sd)), col=colors[i])
  points(log(d)+0.05, log(HS_mean), pch=19, col=rgb(colors1[1,i], colors1[2,i],colors1[3,i],alpha=0.5))
  lines(log(c(d,d))+0.05, c(log(HS_mean-HS_sd), log(HS_mean+HS_sd)), col=rgb(colors1[1,i], colors1[2,i],colors1[3,i],alpha=0.5))
}
}



dev.off()