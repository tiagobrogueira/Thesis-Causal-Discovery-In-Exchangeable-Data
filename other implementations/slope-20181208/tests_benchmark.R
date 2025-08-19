#!/usr/bin/Rscript

source("utilities.R")
source("Slope.R")

### Test Slope
callF = Slope
pair = rep("", length(uv))
eps = rep(0, length(uv))
cds = rep("--", length(uv))
time = rep(0, length(uv))
dXY = rep(0, length(uv))
dYX = rep(0, length(uv))
pval = rep(0, length(uv))
for(i in 1:1){
    print(uv[i])
    t1 = Sys.time()
    t = readI(uv[i])[,1:2]
    res = callF(t, mixedFunctions=T, alpha=1.01)#, mixedFunctions=T, nof=8)     ## don't use p-value
    t2 = Sys.time()
    elapsed = as.numeric(difftime(t2,t1), units="secs")
    pair[i] = uv[i]
    eps[i] = res$eps
    cds[i] = res$cd
    time[i] = elapsed
    dXY[i] = res$sc[1]
    dYX[i] = res$sc[2]
    pval[i] = res$p.value
    print(res$cd)
}
corr = rep(0,length(uv))
corr[ref.uv$V6 == cds] = 1
print(sum(corr))
print(sum(cds == "--"))
resSlope = data.frame(Correct=corr, Eps=eps, Cds=cds, T=time, DXY=dXY, DYX=dYX, PV=pval)
write.table(resSlope, file="results/slope_p_val.tab", row.names=F, quote=F, sep="\t")
print("done")

### test resit
source("utilities.R")
source("resit/code/startups/startupICML.R", chdir = TRUE)
pair = rep("", length(uv))
eps = rep(0, length(uv))
cds = rep("--", length(uv))
p_table = data.frame(pl=rep(1, length(uv)), pr=rep(1, length(uv)))
for(i in 5:5){
    print(uv[i])
    t = readI(uv[i])[,1:2]
    if(uv[i] != 70){
        res = tryCatch({
            r = ICML(t, model = train_gp, indtest = indtestHsic, output = FALSE)
        }, error = function(e) {
            print(e)
            return(list(Eps = 1.0, Cd = "--", P_vals=c(1,1)))
        })
        eps[i] = res$Eps
        cds[i] = res$Cd
        if(length(res$P_vals) > 1){
            p_table$pl[i] = res$P_vals[1]
            p_table$pr[i] = res$P_vals[2]
        }
    }else{
        eps[i] = 1.0
        cds[i] = "--"
    }
    pair[i] = uv[i]
}
corr = rep(0,length(uv))
corr[ref.uv$V6 == cds] = 1
print(sum(corr))
print(sum(cds == "--"))
resResit = data.frame(Correct=corr, Eps=eps, Cds=cds)
write.table(resResit, file="results/resit_gp.tab", row.names=F, quote=F, sep="\t")
write.table(p_table, file="results/resit_gp_p_values.tab", row.names=F, quote=F, sep="\t")
print(sum(p_table$pl < 0.05 & p_table$pr < 0.05))

### test lingram
source("resit/code/startups/startupLINGAM.R", chdir = TRUE)
lingamWrapper = function(t){
    res = tryCatch({
        r = lingamWrap(t)
    }, error = function(e) {
        print(e)
        return(list(B=matrix(c(0,0,0,0), nrow=2, ncol=2), Adj=matrix(c(F,F,F,F), nrow=2, ncol=2)))
    })
    C = 1 * res$Adj
    if(sum(C) != 1){
        causd = "--"
        p_val = 0
    }else{
        if(C[2,1] == 1){
            causd = "<-"
            p_val = res$B[1,2]
        }else if(C[1,2] == 1){
            causd = "->"
            p_val = res$B[2,1]
        }else{
            causd = "--"
            p_val = 0
        }
    }
    print(res)
    return(list(cd=causd, epsilon=p_val))
}
pair = rep("", length(uv))
eps = rep(0, length(uv))
cds = rep("--", length(uv))
for(i in 1:length(uv)){
    print(uv[i])
    t = readI(uv[i])[,1:2]
    res = lingamWrapper(t)
    eps[i] = res$Eps
    cds[i] = res$Cd
    pair[i] = uv[i]
}
corr = rep(0,length(uv))
corr[ref.uv$V6 == cds] = 1
print(sum(corr))
print(sum(cds == "--"))
resLingam = data.frame(Correct=corr, Eps=eps, Cds=cds)
write.table(resResit, file="results/lingam_gp.tab", row.names=F, quote=F, sep="\t")

### test igci
pair = rep("", length(uv))
eps = rep(0, length(uv))
cds = rep("--", length(uv))
for(i in 1:length(uv)){
    print(uv[i])
    t = readI(uv[i])[,1:2]
    res = IGCI(t)
    eps[i] = res$eps
    cds[i] = res$cd
    pair[i] = uv[i]
}
corr = rep(0,length(uv))
corr[ref.uv$V6 == cds] = 1
print(sum(corr))
print(sum(cds == "--"))
resIGCI = data.frame(Correct=corr, Eps=eps, Cds=cds)
write.table(resIGCI, file="results/igci.tab", row.names=F, quote=F, sep="\t")

#### Check simulated data from Mooij et al. JMLR 2016
# requires df with Eps, GT (0,1 or differen)
# (optional) W (weights)
getTableTF = function(res, pos=1, rev=F){
    ll = dim(res)[1]
    if(!("W" %in% colnames(res))){
        res = data.frame(res, W=rep(1,ll))
    }
    res = res[with(res, order(Eps)),]
    if(rev){
        res = res[with(res, order(-Eps)),]
    }
    tabTF = data.frame(TP=rep(0,ll+1), FP=rep(0,ll+1), FN=rep(0,ll+1), TN=rep(0,ll+1))
    tp = 0
    fp = 0
    fn = sum(res$W[res$GT == pos])
    tn = sum(res$W) - fn
    tabTF$TP[1] = tp
    tabTF$FP[1] = fp
    tabTF$FN[1] = fn
    tabTF$TN[1] = tn
    for(i in 1:ll){
        if(res$GT[i] == pos){
            tp = tp + res$W[i]
            fn = fn - res$W[i]
        }else{
            fp = fp + res$W[i]
            tn = tn - res$W[i]
        }
        tabTF$TP[i+1] = tp
        tabTF$FP[i+1] = fp
        tabTF$FN[i+1] = fn
        tabTF$TN[i+1] = tn
    }
    return(tabTF)
}
getROCAUC = function(tabTF){
    ll = dim(tabTF)[1]
    tpr = rep(0, ll)
    fpr = rep(0, ll)
    for(i in 1:ll){
        tpr[i] = tabTF$TP[i] / (tabTF$TP[i] + tabTF$FN[i])
        fpr[i] = tabTF$FP[i] / (tabTF$FP[i] + tabTF$TN[i])
    }
    roc = integrate(splinefun(fpr, tpr), 0, 1)
    return(roc$value)
}
getPRAUC = function(tabTF){
    ll = dim(tabTF)[1]
    pre = rep(0, ll)
    rec = rep(0, ll)
    for(i in 1:ll){
        pre[i] = tabTF$TP[i] / (tabTF$TP[i] + tabTF$FP[i])
        rec[i] = tabTF$TP[i] / (tabTF$TP[i] + tabTF$FN[i])
    }
    prac = integrate(splinefun(rec, pre), 0, 1)
    return(prac$value)
}
getPRACC = function(decR, xx=NULL){
    if(is.null(xx)){
        xx = 1:length(decR) / length(decR)
    }
    prac = integrate(splinefun(xx, decR), 0, 1)
    return(prac$value)
}

dfRP = data.frame(SimS=rep(0,5), SimLnS=rep(0,5), SimGS=rep(0,5), SimS2=rep(0,5), SimLnS2=rep(0,5), SimGS2=rep(0,5))
pos = 1:100 / 100
dfDr = data.frame(Pos=pos, SimS=rep(0,100), SimLnS=rep(0,100), SimGS=rep(0,100), SimS2=rep(0,100), SimLnS2=rep(0,100), SimGS2=rep(0,100))
dird = "data/Benchmark_simulated/SIM-G"
pmeta = read.table(paste(dird, "/pairmeta.txt", sep=""))
readS = function(i){
    f = paste(c(dird, "/pair0", pmeta$V1[i], ".txt"), collapse="")
    if(i < 10){
        f = paste(c(dird, "/pair000", pmeta$V1[i], ".txt"), collapse="")
    }else if(i < 100){
        f = paste(c(dird, "/pair00", pmeta$V1[i], ".txt"), collapse="")
    }
    t = read.table(f, sep="", header=F, stringsAsFactors = F)
    return(t)
}
runSimTest = function(callF, nof=6, mixed=F){
    pmeta = read.table(paste(dird, "/pairmeta.txt", sep=""))
    pair = rep("", 100)
    eps = rep(0, 100)
    cds = rep("--", 100)
    pv = rep(0, 100)
    for(i in 1:100){
        print(i)
        t = readS(i)[,1:2]
        res = callF(t, alpha=2, nof=nof, mixedFunctions=mixed)
        pair[i] = uv[i]
        pv[i] = res$p.value
        eps[i] = res$eps
        cds[i] = res$cd
    }
    corr = rep(0,100)
    corr[pmeta$V2 == 1 & cds == "->"] = 1
    corr[pmeta$V2 == 2 & cds == "<-"] = 1
    print(sum(corr))
    print(sum(cds == "--"))
    simSlope = data.frame(Correct=corr, Eps=eps, Cds=cds, P=pv)
    drsimS = decision_rate_w(simSlope)
    # AUC
    tXtY = getTableTF(data.frame(simSlope, GT=pmeta$V2), pos=1, rev=F)
    rocrX = getROCAUC(tXtY)
    praucX = getPRAUC(tXtY)
    tYtX = getTableTF(data.frame(simSlope, GT=pmeta$V2), pos=2, rev=T)
    rocrY = getROCAUC(tYtX)
    praucY = getPRAUC(tYtX)
    # AUC ACC
    accAUC = getPRACC(drsimS)
    return(list(dr=drsimS[2:101], roc=c(rocrX,rocrY), pr=c(praucX,praucY), drauc=accAUC))
}
getROCPRBenchm = function(res, uvnd_uv){
    gt = rep(0,length(meta.uv$V6))
    gt[ref.uv$V6 == "->"] = 1
    df = data.frame(res[uvnd_uv,], GT=gt[uvnd_uv], W=meta.uv$V6[uvnd_uv])
    tXtY = getTableTF(df, pos=1, rev=F)
    rocrX = getROCAUC(tXtY)
    praucX = getPRAUC(tXtY)
    tYtX = getTableTF(df, pos=0, rev=T)
    rocrY = getROCAUC(tYtX)
    praucY = getPRAUC(tYtX)
    return(list(roc=c(rocrX,rocrY), pr=c(praucX,praucY)))
}
getROCPRBenchmNW = function(res, uvnd_uv){
    gt = rep(0,length(meta.uv$V6))
    gt[ref.uv$V6 == "->"] = 1
    df = data.frame(res[uvnd_uv,], GT=gt[uvnd_uv])
    tXtY = getTableTF(df, pos=1, rev=F)
    rocrX = getROCAUC(tXtY)
    praucX = getPRAUC(tXtY)
    tYtX = getTableTF(df, pos=0, rev=T)
    rocrY = getROCAUC(tYtX)
    praucY = getPRAUC(tYtX)
    return(list(roc=c(rocrX,rocrY), pr=c(praucX,praucY)))
}
dird = "data/Benchmark_simulated/SIM"
res = runSimTest(Slope)
dfRP[1:2,1] = res$roc
dfRP[3:4,1] = res$pr
dfRP[5,1] = res$drauc
dfDr[,2] = res$dr

dird = "data/Benchmark_simulated/SIM-ln"
res = runSimTest(Slope)
dfRP[1:2,2] = res$roc
dfRP[3:4,2] = res$pr
dfRP[5,2] = res$drauc
dfDr[,3] = res$dr

dird = "data/Benchmark_simulated/SIM-G"
res = runSimTest(Slope)
dfRP[1:2,3] = res$roc
dfRP[3:4,3] = res$pr
dfRP[5,3] = res$drauc
dfDr[,4] = res$dr
#################################
dird = "data/Benchmark_simulated/SIM"
res = runSimTest(Slope, nof=8, mixed=T)
dfRP[1:2,4] = res$roc
dfRP[3:4,4] = res$pr
dfRP[5,4] = res$drauc
dfDr[,5] = res$dr

dird = "data/Benchmark_simulated/SIM-ln"
res = runSimTest(Slope, nof=8, mixed=T)
dfRP[1:2,5] = res$roc
dfRP[3:4,5] = res$pr
dfRP[5,5] = res$drauc
dfDr[,6] = res$dr

dird = "data/Benchmark_simulated/SIM-G"
res = runSimTest(Slope, nof=8, mixed=T)
dfRP[1:2,6] = res$roc
dfRP[3:4,6] = res$pr
dfRP[5,6] = res$drauc
dfDr[,7] = res$dr

pdf(width = 12, height=6, file="results/decision_rate.pdf")
colors = c("black","red","darkgreen")
plot(c(0,1), c(0,1), type="n", xlab="decision rate", ylab="% of correct decisions", bty="n", cex.lab=1.3, cex.axis=1.3)
lines(pos, dfDr[,2], col=colors[1], pch=20, lty=1, lwd=2)
lines(pos, dfDr[,3], col=colors[2], pch=20, lty=1, lwd=2)
lines(pos, dfDr[,4], col=colors[3], pch=20, lty=1, lwd=2)
lines(pos, dfDr[,5], col=colors[1], pch=20, lty=2, lwd=2)
lines(pos, dfDr[,6], col=colors[2], pch=20, lty=2, lwd=2)
lines(pos, dfDr[,7], col=colors[3], pch=20, lty=2, lwd=2)
abline(h=0.95,col="darkgrey",lty=2)
abline(h=0.9,col="darkgrey",lty=2)
abline(h=0.85,col="darkgrey",lty=2)
abline(h=0.8,col="darkgrey",lty=2)
legend(0.8, 0.425, c("sim", "simln", "simg", "sim2", "simln2", "simg2"), cex=1.3, col=colors, lty=1, lwd=c(3,2,2), bty="n", y.intersp=1.2)
dev.off()

write.table(dfDr, file="results/decision_rates_sim.tab", sep="\t", quote=F, row.names=F)
write.table(dfRP, file="results/sim_roc_pr_mean.tab", sep="\t", quote=F, row.names=F)

dfb98RP = data.frame(Slope=rep(0,4), SlopeM=rep(0,4), Cure=rep(0,4), Resit=rep(0,4), IGCI=rep(0,4), ANM=rep(0,4))
resSlope = read.table("results/slope.tab", sep="\t", header=T)
resSlopeM = read.table("results/slope_mixed.tab", sep="\t", header=T)
resCure = read.table("results/cure.tab", sep="\t", header=T)
resResit = read.table("results/resit_gp.tab", sep="\t", header=T)
resIGCI = read.table("results/igci.tab", sep="\t", header=T)
resANM = read.table("results/GPIpHISC.tab", sep="\t", header=T)
uvnd_uv = c(1:46,48:65, 67:100)
r = getROCPRBenchm(resSlope, uvnd_uv)
dfb98RP[,1] = c(r$roc, r$pr)
r = getROCPRBenchm(resSlopeM, uvnd_uv)
dfb98RP[,2] = c(r$roc, r$pr)
r = getROCPRBenchm(resCure, uvnd_uv)
dfb98RP[,3] = c(r$roc, r$pr)
r = getROCPRBenchm(resResit, uvnd_uv)
dfb98RP[,4] = c(r$roc, r$pr)
r = getROCPRBenchm(resIGCI, uvnd_uv)
dfb98RP[,5] = c(r$roc, r$pr)
r = getROCPRBenchm(resANM, uvnd_uv)
dfb98RP[,6] = c(r$roc, r$pr)

dfb79RP = data.frame(Slope=rep(0,4), SlopeM=rep(0,4), Cure=rep(0,4), Resit=rep(0,4), IGCI=rep(0,4), ANM=rep(0,4))
uvnd_uv = c(1:46,48:65,67:81)
r = getROCPRBenchmNW(resSlope[1:79,], uvnd_uv)
dfb79RP[,1] = c(r$roc, r$pr)
r = getROCPRBenchmNW(resSlopeM[1:79,], uvnd_uv)
dfb79RP[,2] = c(r$roc, r$pr)
r = getROCPRBenchmNW(resCure[1:79,], uvnd_uv)
dfb79RP[,3] = c(r$roc, r$pr)
r = getROCPRBenchmNW(resResit[1:79,], uvnd_uv)
dfb79RP[,4] = c(r$roc, r$pr)
r = getROCPRBenchmNW(resIGCI[1:79,], uvnd_uv)
dfb79RP[,5] = c(r$roc, r$pr)
r = getROCPRBenchmNW(resANM[1:79,], uvnd_uv)
dfb79RP[,6] = c(r$roc, r$pr)

write.table(dfb98RP, file="results/tueb_98_rocpr.tab", sep="\t", quote=F, row.names=F)
write.table(dfb79RP, file="results/tueb_79_rocpr.tab", sep="\t", quote=F, row.names=F)
