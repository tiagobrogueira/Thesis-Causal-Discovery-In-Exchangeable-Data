#!/usr/bin/Rscript

source("utilities.R")
source("Slope.R")
source("resit/code/startups/startupICML.R", chdir = TRUE)

source("test_synthetic_data.R")

icmlWrapper = function(t){
    res = tryCatch({
        r = ICML(t, model = train_gp, indtest = indtestHsic, output = FALSE)
    }, error = function(e) {
        print(e)
        return(list(Eps = 1.0, Cd = "--", P_vals=c(1,1)))
    })
    return(list(eps=res$Eps, cd=res$Cd, pvalues=res$P_vals))
}

crackWrapper = function(t, name="test"){
    write.table(t, file="temp/current.txt", sep=" ", quote=F, row.names=F, col.names=F)
    system(paste(c("./crack -s 3 -o temp/crack_ -x 1 -c -d ' ' -t 'i' -i temp/current.txt -a ", name), collapse=""))
    result_t = read.table("temp/crack_results.tab", header=F, stringsAsFactors=F, sep="\t")
    row = dim(result_t)[1]
    return(list(cd=result_t[row,3], eps=result_t[row,2]))
}

### Synthetic for linear an cubic and reciprocal

resSlope = data.frame(test(testP=Slope))
resSlope = rbind(resSlope,test(testP=Slope, cubic=2))
resSlope = rbind(resSlope,test(testP=Slope, neg=1))
write.table(resSlope, "results/synthetic_data_slope.tab", row.names=F, quote=F)

resIGCI = data.frame(test(testP=IGCI))
resIGCI = rbind(resIGCI,test(testP=IGCI, cubic=2))
resIGCI = rbind(resIGCI,test(testP=IGCI, neg=1))
write.table(resIGCI, "results/synthetic_data_igci.tab", row.names=F, quote=F)

resResit = data.frame(test(testP=icmlWrapper))
resResit = rbind(resResit,test(testP=icmlWrapper, cubic=2))
resResit = rbind(resResit,test(testP=icmlWrapper, neg=1))
write.table(resResit, "results/synthetic_data_resit.tab", row.names=F, quote=F)


### Synthetic -- test number of produced bins
df = testBins()
perc = df$B / df$P
mean = 1000 / df$V
resBins = data.frame(df, perc, mean)
write.table(resBins, "results/synthetic_bins.dat", sep=" ", col.names=F, row.names=F, quote=F)

### Pvalues and Confidence
testBinsAndPvalues()

### Consitency test with WMW
testConsistency()

## eval confidene
conf = read.table("results/confidence_test.tab", header=T, sep="\t")
conf.ANMs = -log(-conf$ANM_C)
conf.ANMs[conf.ANMs == Inf] = 1000
box.df.anm = data.frame(S100=conf.ANMs[conf$Samples == 100], S250=conf.ANMs[conf$Samples == 250], S500=conf.ANMs[conf$Samples == 500], S1000=conf.ANMs[conf$Samples == 1000])
box.df.slope = data.frame(S100=conf$Slope_C[conf$Samples == 100], S250=conf$Slope_C[conf$Samples == 250], S500=conf$Slope_C[conf$Samples == 500], S1000=conf$Slope_C[conf$Samples == 1000])

anm_bounds = rbind(boxplot.stats(box.df.anm[,1])$stats, boxplot.stats(box.df.anm[,2])$stats, boxplot.stats(box.df.anm[,3])$stats, boxplot.stats(box.df.anm[,4])$stats)
slope_bounds = rbind(boxplot.stats(box.df.slope[,1])$stats, boxplot.stats(box.df.slope[,2])$stats, boxplot.stats(box.df.slope[,3])$stats, boxplot.stats(box.df.slope[,4])$stats)

write.table(anm_bounds, "results/confidence_resit.dat", sep=" ", row.names=F, col.names=F, quote=F)
write.table(slope_bounds, "results/confidence_slope.dat", sep=" ", row.names=F, col.names=F, quote=F)

## colors
boxplot(box.df.anm)
boxplot(box.df.slope)

## eval pvalues
pval = read.table("results/pvalues_test.tab", header=T, sep="\t")
alpha = -0.05
p100OL = sum(pval$ANM_LP1[pval$Samples == 100] > alpha | pval$ANM_LP2[pval$Samples == 100] > alpha) / 100
p100BL = sum(pval$ANM_LP1[pval$Samples == 100] > alpha & pval$ANM_LP2[pval$Samples == 100] > alpha) / 100
p250OL = sum(pval$ANM_LP1[pval$Samples == 250] > alpha | pval$ANM_LP2[pval$Samples == 250] > alpha) / 100
p250BL = sum(pval$ANM_LP1[pval$Samples == 250] > alpha & pval$ANM_LP2[pval$Samples == 250] > alpha) / 100
p500OL = sum(pval$ANM_LP1[pval$Samples == 500] > alpha | pval$ANM_LP2[pval$Samples == 500] > alpha) / 100
p500BL = sum(pval$ANM_LP1[pval$Samples == 500] > alpha & pval$ANM_LP2[pval$Samples == 500] > alpha) / 100
p1000OL = sum(pval$ANM_LP1[pval$Samples == 1000] > alpha | pval$ANM_LP2[pval$Samples == 1000] > alpha) / 100
p1000BL = sum(pval$ANM_LP1[pval$Samples == 1000] > alpha & pval$ANM_LP2[pval$Samples == 1000] > alpha) / 100

p100OC = sum(pval$ANM_CP1[pval$Samples == 100] > alpha | pval$ANM_CP2[pval$Samples == 100] > alpha) / 100
p100BC = sum(pval$ANM_CP1[pval$Samples == 100] > alpha & pval$ANM_CP2[pval$Samples == 100] > alpha) / 100
p250OC = sum(pval$ANM_CP1[pval$Samples == 250] > alpha | pval$ANM_CP2[pval$Samples == 250] > alpha) / 100
p250BC = sum(pval$ANM_CP1[pval$Samples == 250] > alpha & pval$ANM_CP2[pval$Samples == 250] > alpha) / 100
p500OC = sum(pval$ANM_CP1[pval$Samples == 500] > alpha | pval$ANM_CP2[pval$Samples == 500] > alpha) / 100
p500BC = sum(pval$ANM_CP1[pval$Samples == 500] > alpha & pval$ANM_CP2[pval$Samples == 500] > alpha) / 100
p1000OC = sum(pval$ANM_CP1[pval$Samples == 1000] > alpha | pval$ANM_CP2[pval$Samples == 1000] > alpha) / 100
p1000BC = sum(pval$ANM_CP1[pval$Samples == 1000] > alpha & pval$ANM_CP2[pval$Samples == 1000] > alpha) / 100

p.values.distrib = data.frame(Samples=c(100,250,500,1000), OneL=c(p100OL, p250OL, p500OL, p1000OL), BothL=c(p100BL, p250BL, p500BL, p1000BL), OneC=c(p100OC, p250OC, p500OC, p1000OC), BothC=c(p100BC, p250BC, p500BC, p1000BC))
write.table(p.values.distrib, "results/pvalue_distribution.tab", sep=" ", row.names=F, col.names=F, quote=F)


