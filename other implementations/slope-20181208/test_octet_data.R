#!/usr/bin/Rscript

source("utilities.R")
source("Slope.R")
source("resit/code/startups/startupICML.R", chdir = TRUE)

testOctetData = function(testP=NULL){
    pairs = c("pair0001", "pair0002", "pair0003", "pair0004", "pair0005", "pair0006", "pair0007", "pair0008", "pair0009", "pair0010")
    gt=rep("<-", 10)
    w=c(0.5,0.5,0.34,0.33,0.34,0.33,1,0.33,0.33,1)
    decision=rep("--", 10)
    confidence=rep(0,10)
    correct=rep(0,10)
    for(i in 1:10){
        fileI = paste(c("data/octet/", pairs[i], ".txt"), collapse="")
        t = read.table(file=fileI, header=F, sep="\t")
        res = testP(t)
        decision[i] = res$cd
        confidence[i] = res$eps
        if(decision[i] == gt[i]){
            correct[i] = 1
        }
    }
    df = data.frame(Pairs=pairs, GT=gt, W=w, Decision=decision, Confidence=confidence, Correct=correct)
    return(df)
}

icmlWrapper = function(t){
    res = tryCatch({
        r = ICML(t, model = train_gp, indtest = indtestHsic, output = FALSE)
    }, error = function(e) {
        print(e)
        return(list(Eps = 1.0, Cd = "--", P_vals=c(1,1)))
    })
    return(list(eps=res$Eps, cd=res$Cd, pvalues=res$P_vals))
}

resSlope = testOctetData(testP=Slope)
write.table(data.frame(Correct=resSlope$Correct, Eps=resSlope$Confidence, Cds=resSlope$Decision), file="results/octet_slope.tab", row.names=F, quote=F, sep="\t")
resResit = testOctetData(testP=icmlWrapper)
write.table(data.frame(Correct=resResit$Correct, Eps=resResit$Confidence, Cds=resResit$Decision), file="results/octet_resit.tab", row.names=F, quote=F, sep="\t")
resIGCI = testOctetData(testP=IGCI)
write.table(data.frame(Correct=resIGCI$Correct, Eps=resIGCI$Confidence, Cds=resIGCI$Decision), file="results/octet_igci.tab", row.names=F, quote=F, sep="\t")
