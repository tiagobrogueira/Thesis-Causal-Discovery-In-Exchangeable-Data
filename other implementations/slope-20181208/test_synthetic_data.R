#!/usr/bin/Rscript

logg = function(x){
    if(x == 0){
        return(0)
    }else{
        return(log2(x))
    }
}

getYWithGauss = function(x, Sd=1, q=1, b=0, c=0,d=1,e=0,f=0, roundIt=F){
    xe = x
    xe[x == 0] = 0.01
    noise = rnorm(length(x), mean=0, sd=Sd)
    if(q != 1){
        sq0 = noise < 0
        noise = abs(noise)^q
        noise[sq0] = -noise[sq0]
    }
    if(roundIt){
        noise = round(noise)
    }
    y = b * (x^3) + c * (x^2) + e * exp(x) + d * x + f * xe^(-1) + noise
    df = data.frame(X=x, Y=y)
    return(df)
}

getYWithUnif = function(x, t=0.5, b=0, d=1,e=0,f=0, roundIt=F){
    xe = x
    xe[x == 0] = 0.01
    noise = runif(length(x), min=-t, max=t)
    if(roundIt){
        noise = round(noise)
    }
    y = b * (x^3) + e * exp(x) + d * x + f * xe^(-1) + noise
    df = data.frame(X=x, Y=y)
    return(df)
}

getYWithLocal = function(x, b=0, d=1,e=0,f=0, roundIt=F){
    xd = round(x * 10000)
    xs = unique(sort(xd))
    y = x
    for(s in xs){
        pos = (xd == s)
        val = s / 10000
        l = sum(pos)
        xhead = 1:sum(pos) - sum(pos) / 2
        par = runif(1, min=0.1, max=0.5)
        yn = xhead * d + val * d + val^2 * b
        y[pos] = y[pos] + yn
    }
    df = data.frame(X=x, Y=y)
    return(df)
}

getYWithNonAdditive = function(x, v=0.25, sigma=1, b=0, c=0, d=1,e=0,f=0){
    xe = x
    xe[x == 0] = 0.01
    y = b * (x^3) + e * exp(x) + c * (x^2) + d * x + f * xe^(-1)
    for(i in 1:length(x)){
        X = x[i]
        E = sigma * rnorm(1) * abs(sin(2*pi*v*X)) + 0.25 * sigma * rnorm(1) * abs(sin(20 * pi * v * X))
        y[i] = y[i] + E
    }
    df = data.frame(X=x, Y=y)
    return(df)
}

### TODO
## --> load function that should be tested
## OR source this file and apply the test function

### Fix "samples" and "iterations"
test = function(testP=NULL, samples = 1000, iterations=100, cubic=0, expo=0, neg=0, rou=F){
    results = list(uu=0,ug=0,un=0,gu=0,gg=0,gn=0,bu=0,bg=0,bn=0,pu=0,pg=0,pn=0)
    bb = cubic
    
    ######### uniform
    
    print("X~unif, N~unif")
    set.seed(1)
    sum_corr = 0
    for(i in 1:iterations){
        tX = runif(1, min=1, max=10)
        tN = runif(1, min=0.01, max=5)
        x = runif(samples, min=-tX, max=tX)
        df = getYWithUnif(x, t=tN, b=bb, e=expo, f=neg, roundIt=rou)
        cd = testP(df)$cd
        if(cd == "->"){
            sum_corr = sum_corr + 1
        }
    }
    results$uu = sum_corr
    
    print("X~unif, N~gauss")
    set.seed(1)
    sum_corr = 0
    for(i in 1:iterations){
        tX = runif(1, min=1, max=10)
        sdN = runif(1, min=0, max=5)
        x = runif(samples, min=-tX, max=tX)
        df = getYWithGauss(x, Sd=sdN, b=bb, e=expo, f=neg, roundIt=rou)
        cd = testP(df)$cd
        if(cd == "->"){
            sum_corr = sum_corr + 1
        }
    }
    sum_corr
    results$ug = sum_corr
    
    print("X~unif, N~non-additive")
    set.seed(1)
    sum_corr = 0
    for(i in 1:iterations){
        tX = runif(1, min=1, max=10)
        sdN = runif(1, min=0.25, max=1.1)
        x = runif(samples, min=-tX, max=tX)
        df = getYWithNonAdditive(x, v=sdN, b=bb, e=expo, f=neg)
        cd = testP(df)$cd
        if(cd == "->"){
            sum_corr = sum_corr + 1
        }
    }
    results$un = sum_corr
    
    ######### gaussian
    
    print("X~gauss, N~unif")
    set.seed(1)
    sum_corr = 0
    for(i in 1:iterations){
        tX = runif(1, min=1, max=10)
        x = rnorm(samples, sd=tX)
        sq0 = x < 0
        x = abs(x)^(0.7)
        x[sq0] = -x[sq0]
        tN = runif(1, min=0, max=max(x)/2)
        df = getYWithUnif(x, t=tN, b=bb, f=neg)
        cd = testP(df)$cd
        if(cd == "->"){
            sum_corr = sum_corr + 1
        }
    }
    results$gu = sum_corr
    
    print("X~gauss, N~gauss")
    set.seed(1)
    sum_corr = 0
    for(i in 1:iterations){
        tX = runif(1, min=1, max=10)
        x = rnorm(samples, sd=tX)
        sq0 = x < 0
        x = abs(x)^(0.7)
        x[sq0] = -x[sq0]
        sdN = runif(1, min=0, max=max(x)/2)
        df = getYWithGauss(x, Sd=sdN, b=bb, f=neg)
        cd = testP(df)$cd
        if(cd == "->"){
            sum_corr = sum_corr + 1
        }
    }
    results$gg = sum_corr
    
    print("X~gauss, N~non-additive")
    set.seed(1)
    sum_corr = 0
    for(i in 1:iterations){
        sdN = runif(1, min=0.25, max=1.1)
        tX = runif(1, min=1, max=10)
        x = rnorm(samples, sd=tX)
        sq0 = x < 0
        x = abs(x)^(0.7)
        x[sq0] = -x[sq0]
        df = getYWithNonAdditive(x, v=sdN, b=bb, f=neg)
        cd = testP(df)$cd
        if(cd == "->"){
            sum_corr = sum_corr + 1
        }
    }
    results$gn = sum_corr
    
    ######### binomial
    
    print("X~binomial, N~unif")
    set.seed(1)
    sum_corr = 0
    for(i in 1:iterations){
        tX = runif(1, min=1, max=10)
        tP = runif(1, min=0.1, max=0.9)
        x = rbinom(samples, size=ceiling(tX), prob=tP)
        tN = runif(1, min=0, max=max(x))
        df = getYWithUnif(x, t=tN, b=bb, e=expo, f=neg, roundIt=rou)
        cd = testP(df)$cd
        if(cd == "->"){
            sum_corr = sum_corr + 1
        }
    }
    results$bu = sum_corr
    
    print("X~binomial, N~gauss")
    set.seed(1)
    sum_corr = 0
    for(i in 1:iterations){
        tX = runif(1, min=1, max=10)
        tP = runif(1, min=0.1, max=0.9)
        x = rbinom(samples, size=ceiling(tX), prob=tP)
        sdN = runif(1, min=0, max=max(x)/2)
        df = getYWithGauss(x, Sd=sdN, b=bb, e=expo, f=neg, roundIt=rou)
        cd = testP(df)$cd
        if(cd == "->"){
            sum_corr = sum_corr + 1
        }
    }
    results$bg = sum_corr
    
    print("X~binomial, N~non-additive")
    set.seed(1)
    sum_corr = 0
    for(i in 1:iterations){
        tX = runif(1, min=1, max=10)
        tP = runif(1, min=0.1, max=0.9)
        x = rbinom(samples, size=ceiling(tX), prob=tP)
        sdN = runif(1, min=0.25, max=1.1)
        df = getYWithNonAdditive(x, v=sdN, b=bb, e=expo, f=neg)
        cd = testP(df)$cd
        if(cd == "->"){
            sum_corr = sum_corr + 1
        }
    }
    results$bn = sum_corr
    
    ######### poisson
    
    print("X~poisson, N~unif")
    set.seed(1)
    sum_corr = 0
    for(i in 1:iterations){
        tX = runif(1, min=1, max=10)
        x = rpois(samples, lambda=round(tX))
        tN = runif(1, min=0, max=max(x))
        df = getYWithUnif(x, t=tN, b=bb, e=expo, f=neg, roundIt=rou)
        cd = testP(df)$cd
        if(cd == "->"){
            sum_corr = sum_corr + 1
        }
    }
    results$pu = sum_corr
    
    print("X~poisson, N~gauss")
    set.seed(1)
    sum_corr = 0
    for(i in 1:iterations){
        tX = runif(1, min=1, max=10)
        x = rpois(samples, lambda=round(tX))
        sdN = runif(1, min=0, max=max(x)/2)
        df = getYWithGauss(x, Sd=sdN, b=bb, e=expo, f=neg, roundIt=rou)
        cd = testP(df)$cd
        if(cd == "->"){
            sum_corr = sum_corr + 1
        }
    }
    results$pg = sum_corr
    
    print("X~poisson, N~non-additive")
    set.seed(1)
    sum_corr = 0
    for(i in 1:iterations){
        tX = runif(1, min=1, max=10)
        x = rpois(samples, lambda=round(tX))
        sdN = runif(1, min=0.25, max=1.1)
        df = getYWithNonAdditive(x, v=sdN, b=bb, e=expo, f=neg)
        cd = testP(df)$cd
        if(cd == "->"){
            sum_corr = sum_corr + 1
        }
    }
    results$pn = sum_corr
    return(results)
}

#### Test bins
testBins = function(iterations=100, cubic=0){
    samples= 1000
    different_values = c(200, 150, 100,75, 50, 40, 30, 20, 10)
    lv = length(different_values)
    potentials = rep(0, lv)
    used_bins = rep(0, lv)
    for(k in 1:lv){
        num_vals = different_values[k]
        print(num_vals)
        bb = 0
        pp = 0
        for(i in 1:iterations){
            sdN = runif(1, min=1, max=10)
            x = round(runif(samples, min=1, max=num_vals))
            df = getYWithGauss(x, Sd=sdN, b=cubic)
            bins = Slope(df)$bins
            bb = bb + (bins[1] - 1)
            pp = pp + (bins[3] - 1)
        }
        potentials[k] = pp
        used_bins[k] = bb
    }
    df = data.frame(V=different_values, B=used_bins, P=potentials)
    return(df)
}

#### Test confidence
testBinsAndPvalues = function(iterations=100){
    df.conf = data.frame(Samples=rep(0,iterations*4), ANM_L=rep(0,iterations*4), ANM_C=rep(0,iterations*4), Slope_L=rep(0,iterations*4), Slope_C=rep(0,iterations*4))
    df.pvalues = data.frame(Samples=rep(0,iterations*4), ANM_LP1=rep(0,iterations*4), ANM_LP2=rep(0,iterations*4), ANM_CP1=rep(0,iterations*4), ANM_CP2=rep(0,iterations*4))
    
    row = 1
    testInstances = c(100,250,500,1000)
    for(samples in testInstances){
        for(i in 1:iterations){
            tX = runif(1, min=1, max=10)
            x = rnorm(samples, sd=tX)
            tN = runif(1, min=0, max=max(x)/2)
            df_L = getYWithUnif(x, t=tN)
            df_C = getYWithUnif(x, t=tN, b=2)
            
            # note down samples
            df.conf$Samples[row] = samples
            df.pvalues$Samples[row] = samples
            
            # apply anm
            anmL = icmlWrapper(df_L)
            df.conf$ANM_L[row] = anmL$eps
            df.pvalues$ANM_LP1[row] = anmL$pvalues[1]
            df.pvalues$ANM_LP2[row] = anmL$pvalues[2]
            anmC = icmlWrapper(df_C)
            df.conf$ANM_C[row] = anmC$eps
            df.pvalues$ANM_CP1[row] = anmC$pvalues[1]
            df.pvalues$ANM_CP2[row] = anmC$pvalues[2]
            
            # apply Slope
            slopeL = Slope(df_L)
            df.conf$Slope_L[row] = slopeL$eps
            slopeC = Slope(df_C)
            df.conf$Slope_C[row] = slopeC$eps
            
            # next row
            row = row + 1
        }
    }
    write.table(df.conf, "results/confidence_test.tab", row.names=F, quote=F, sep="\t")
    write.table(df.pvalues, "results/pvalues_test.tab", row.names=F, quote=F, sep="\t")
}

testConsistency = function(){
    iterations=100
    samples = 1000
    print("X~unif, N~unif")
    set.seed(1)
    uu = rep(-1,100)
    bb=0
    for(i in 1:iterations){
        print(i)
        if(bb == 0){
            bb = 2
        }else{
            bb = 0
        }
        tX = runif(1, min=1, max=10)
        tN = runif(1, min=1, max=10)
        x = runif(samples, min=-tX, max=tX)
        df = getYWithUnif(x, t=tN, b=bb)
        x = normX(df[,1],1)
        y = normX(df[,2],1)
        uu[i] = fitTestWMW(y,x)
    }
    
    print("X~unif, N~gauss")
    set.seed(1)
    ug = rep(-1,100)
    bb=0
    for(i in 1:iterations){
        print(i)
        if(bb == 0){
            bb = 2
        }else{
            bb = 0
        }
        tX = runif(1, min=1, max=10)
        sdN = runif(1, min=0, max=5)
        x = runif(samples, min=-tX, max=tX)
        df = getYWithGauss(x, Sd=sdN, b=bb)
        x = normX(df[,1],1)
        y = normX(df[,2],1)
        ug[i] = fitTestWMW(y,x)
    }
    
    print("X~unif, N~non-additive")
    set.seed(1)
    un = rep(-1,100)
    bb=0
    for(i in 1:iterations){
        print(i)
        if(bb == 0){
            bb = 2
        }else{
            bb = 0
        }
        tX = runif(1, min=1, max=10)
        sdN = runif(1, min=0.25, max=1.1)
        x = runif(samples, min=-tX, max=tX)
        df = getYWithNonAdditive(x, v=sdN, b=bb)
        x = normX(df[,1],1)
        y = normX(df[,2],1)
        un[i] = fitTestWMW(y,x)
    }
    results = data.frame(uu, ug, un)
    write.table(results, "results/consistency_test.tab", row.names=F, quote=F, sep="\t")
}
