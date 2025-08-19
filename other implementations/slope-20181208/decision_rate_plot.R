### Plot decision rate

source("utilities.R")

uvnd_uv = c(1:46,48:65, 67:100) ## old  81
uvk = uv[uvnd_uv]
resMain = read.table("results/slope.tab", sep="\t", header=T)
drMain = decision_rate_meta(resMain[uvnd_uv,], uv=uvk)
resMainNeg = read.table("results/slope_p_val_det.tab", sep="\t", header=T)
drMainNeg = decision_rate_meta(resMainNeg[uvnd_uv,], uv=uvk)
resResit = read.table("results/resit_gp.tab", sep="\t", header=T)
drResit = decision_rate_meta_pos(resResit[uvnd_uv,], uv=uvk)
resANM = read.table("results/GPIpHISC.tab", sep="\t", header=T)
drANM = decision_rate_meta_pos(resGPIpHISC[uvnd_uv,], uv=uvk)
resIGCI = read.table("results/igci.tab", sep="\t", header=T)
drIGCI = decision_rate_meta(resIGCI[uvnd_uv,], uv=uvk)
resCure = read.table("results/cure.tab", sep="\t", header=T)
drCure = decision_rate_meta(resCure[uvnd_uv,], uv=uvk)

pplus = rep(0.0, length(uvnd_uv))
pminus = rep(0.0, length(uvnd_uv))
sum = 0
pos = rep(0,length(uvnd_uv))
library("Hmisc")
#add = 66 / 79
add = 34.4979 / 98
for(i in 1:length(uvnd_uv)){
    sum = sum + add
    r = binconf(x=sum/2,n=sum, method="exact")
    pplus[i] = r[2]
    pminus[i] = r[3]
}
pos = 1:(length(uvnd_uv)) / (length(uvnd_uv))
pos = c(0, pos)
pplus = c(0,pplus)
pminus = c(1,pminus)

pdf(width = 12, height=6, file="results/decision_rate_restricted.pdf")
colors = c("darkblue","darkgreen", "darkred", "orange", "black", "green")
plot(c(0,1), c(0,1), type="n", xlab="decision rate", ylab="% of correct decisions", bty="n", cex.lab=1.3, cex.axis=1.3)
polygon(c(pos, rev(pos)), c(pplus, rev(pminus)), col = "lightgrey", border = NA)
lines(drMain$S, drMain$D, col=colors[1], pch=20, lty=1, lwd=3)
lines(drMainNeg$S, drMainNeg$D, col=colors[2], pch=20, lty=1, lwd=2)
lines(drResit$S, drResit$D, col=colors[3], pch=20, lty=1, lwd=2)
lines(drANM$S, drANM$D, col=colors[4], pch=20, lty=1, lwd=2)
lines(drIGCI$S, drIGCI$D, col=colors[5], pch=20, lty=1, lwd=2)
lines(drCure$S, drCure$D, col=colors[5], pch=20, lty=1, lwd=2)
abline(h=0.9,col="darkgrey",lty=2)
abline(h=0.8,col="darkgrey",lty=2)
abline(h=0.7,col="darkgrey",lty=2)
abline(h=0.6,col="darkgrey",lty=2)

legend(0.8, 0.425, c("Slope", "Slope-", "Resit", "ANM", "IGCI", "Cure"), cex=1.3, col=colors, lty=1, lwd=c(3,2,2), bty="n", y.intersp=1.2)
dev.off()

df = data.frame(P=pos, PP=pplus, PM=pminus, SlopePos=drMain$S, Slope=drMain$D, CurePos=drCure$S, Cure=drCure$D, ResitPos=drResit$S, Resit=drResit$D, IGCIPos=drIGCI$S, IGCI=drIGCI$D, GPIPos=drANM$S, GPI=drANM$D, SlopeDPos=drMainNeg$S, SlopeD=drMainNeg$D)
write.table(df, file="decision_rate_weighted_.dat", row.names=F, quote=F, col.names = T)


uvnd_uv = c(1:46,48:65, 67:100) ## old  81
uvk = uv[uvnd_uv]
resMain = read.table("results/slope_p_val.tab", sep="\t", header=T)
drMain = decision_rate_w(resMain[uvnd_uv,])
resMainNeg = read.table("results/slope_p_val.tab", sep="\t", header=T)
resMainNeg$Eps = -log10(resMainNeg$PV)
resMainNeg$Eps[resMainNeg$Eps > 300] = 300
drMainNeg = decision_rate_w(resMainNeg[uvnd_uv,])
resResit = read.table("results/resit_gp.tab", sep="\t", header=T)
drResit = decision_rate_w_pos(resResit[uvnd_uv,])
resANM = read.table("results/GPIpHISC.tab", sep="\t", header=T)
drANM = decision_rate_w_pos(resGPIpHISC[uvnd_uv,])
resIGCI = read.table("results/igci.tab", sep="\t", header=T)
drIGCI = decision_rate_w(resIGCI[uvnd_uv,])
resCure = read.table("results/cure.tab", sep="\t", header=T)
drCure = decision_rate_w(resCure[uvnd_uv,])

pplus = rep(0.0, length(uvnd_uv))
pminus = rep(0.0, length(uvnd_uv))
sum = 0
pos = rep(0,length(uvnd_uv))
library("Hmisc")
add = 1
for(i in 1:length(uvnd_uv)){
    sum = sum + add
    r = binconf(x=sum/2,n=sum, method="exact")
    pplus[i] = r[2]
    pminus[i] = r[3]
}
pos = 1:(length(uvnd_uv)) / (length(uvnd_uv))
pos = c(0, pos)
pplus = c(0,pplus)
pminus = c(1,pminus)

pdf(width = 12, height=6, file="results/decision_rate_no_weights.pdf")
colors = c("darkblue","darkgreen", "darkred", "orange", "black", "green")
plot(c(0,1), c(0,1), type="n", xlab="decision rate", ylab="% of correct decisions", bty="n", cex.lab=1.3, cex.axis=1.3)
polygon(c(pos, rev(pos)), c(pplus, rev(pminus)), col = "lightgrey", border = NA)
lines(pos, drMain, col=colors[1], pch=20, lty=1, lwd=3)
lines(pos, drMainNeg, col=colors[2], pch=20, lty=1, lwd=2)
lines(pos, drResit, col=colors[3], pch=20, lty=1, lwd=2)
lines(pos, drANM, col=colors[4], pch=20, lty=1, lwd=2)
lines(pos, drIGCI, col=colors[5], pch=20, lty=1, lwd=2)
lines(pos, drCure, col=colors[5], pch=20, lty=1, lwd=2)
abline(h=0.9,col="darkgrey",lty=2)
abline(h=0.8,col="darkgrey",lty=2)
abline(h=0.7,col="darkgrey",lty=2)
abline(h=0.6,col="darkgrey",lty=2)

legend(0.8, 0.425, c("Slope", "SlopeP", "Resit", "ANM", "IGCI", "Cure"), cex=1.3, col=colors, lty=1, lwd=c(3,2,2), bty="n", y.intersp=1.2)
dev.off()

df = data.frame(P=pos, PP=pplus, PM=pminus, Slope=drMain, Cure=drCure, Resit=drResit, IGCI=drIGCI, GPI=drANM, SlopeP=drMainNeg)
write.table(df, file="decision_rate_minus.dat", row.names=F, quote=F, col.names = T)
