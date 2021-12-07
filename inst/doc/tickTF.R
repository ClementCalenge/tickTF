## ----setup, include=FALSE, cache=FALSE--------------------
# set global chunk options
library('knitr')
opts_chunk$set(fig.path="caperpy-",
               fig.align="center",
               fig.show="hold",
               echo=TRUE,
               results="markup",
               fig.width=9,
               fig.height=9, out.width='0.8\\linewidth',
               out.height='0.8\\linewidth',
               cache=FALSE,
               dev='png',
               concordance=TRUE,
               error=FALSE)
opts_knit$set(aliases = c(h = 'fig.height',
              w = 'fig.width',
              wo='out.width',
              ho='out.height'))
options(replace.assign=TRUE,width=60)
load("extdata.rda")
set.seed(9567)


## ----eval=FALSE-------------------------------------------
## ## If devtools is not yet installed, type
## install.packages("devtools")
## 
## ## Install the package caperpyogm
## devtools::install_github("ClementCalenge/tickTF", ref="main")


## ----package-loading--------------------------------------
library(tickTF)

## Our dataset
head(fticks)


## ----relation-age-bsa-------------------------------------
plot(fticks$age, fticks$bsa, xlab="Age", ylab="Body surface area")
ag <- seq(0,8,length=100)
lines(ag, predict(loess(bsa~age, data=fticks),
                  newdata=data.frame(age=ag)), col="red", lwd=2)


## ----distribution-Yj-nimble-------------------------------
library(nimble)

## Function giving the probability of Y_j (here passed as argument x)
## given the value of Omega (here passed as argument lambda) and phi

dticks <- nimbleFunction(
    run=function(x=integer(0), lambda=double(0, default=1),
                 phi=double(0, default=1),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    proba <- phi/(phi+lambda)
    p1 <- pnbinom(9, size=phi, prob=proba,
                  lower.tail = TRUE, log.p = TRUE)
    p3 <- pnbinom(20, size=phi, prob=proba,
                  lower.tail = FALSE, log.p = TRUE)
    p2 <- log(0.0000001+
              pnbinom(20, size=phi, prob=proba,
                      lower.tail = TRUE, log.p = FALSE)-
              pnbinom(9, size=phi, prob=proba,
                      lower.tail = TRUE, log.p = FALSE))
    lpro <- (c(p1,p2,p3)[x])
    if (log) return(lpro)
    else return(exp(lpro))    
})


## Function needed to simulate random values from the probability distribution
## of Y_j, as a function of Omega (here passed as lambda) and phi.
rticks <- nimbleFunction(
    run = function(n = integer(0), lambda=double(0, default=1),
                   phi=double(0, default=1)) {
        returnType(double(0))
        if(n != 1) print("rtiques only allows n = 1; using n = 1.")
        u <- rnbinom(1, size=phi, prob=phi/(phi+lambda))
        if (u<10) {
            return(1)
        }
        if(u<21) {
            return(2)
        }
        return(3)
    })


## ----nimble-code-for-fit----------------------------------
codePaperModel <- nimbleCode({
    
    ## Priors
    intercept ~ dnorm(0,sd=20)
    gammah ~ dnorm(0,sd=20)
    gammav ~ dnorm(0,sd=20)
    beta0p ~ dnorm(0,sd=20)
    beta1p ~ dnorm(0,sd=20)
    beta2p ~ dnorm(0,sd=20)
    for (i in 1:3) {
        sigma_u[i]~dunif(0,100)
    }
    sigma_s~dunif(0,100)
    phi~dunif(0,100);

    
    ## Year random effects
    for (i in 1:nyears) {
        delta_u[i] ~ dnorm(intercept, sd=sigma_u[pery[i]]);
    }

    ## Because we work on centred data for body surface area
    ## (for better mixing)
    beta0 <- beta0p - beta1p*xbar + beta2p*pow(xbar,2.0)
    beta1 <- beta1p - 2*beta2p*xbar
    beta2 <- beta2p

    
    ## Integration of body surface area
    for (i in 1:nobs) {
        
        ## Expected bsa
        yhat[i] <- beta0p + beta1p*AGEc[i] + beta2p*pow(AGEc[i],2.0)
        SURFACE[i] ~ dnorm(yhat[i], sd=sigma_s)
        epsilon[i] <- SURFACE[i] - yhat[i]

        ## Integration
        Ej[i] <- step(5-AGE[i])*(beta0*AGE[i] + (beta1/2.0)*pow(AGE[i],2.0) +
                                 (beta2/3.0)*pow(AGE[i],3.0) +
                                 epsilon[i]*AGE[i]) +
            (1-step(5-AGE[i]))*((5.0*beta0 + 12.5*beta1 + (125.0/3.0)*beta2 + 5.0*epsilon[i]) +
                                exp(gammav)*(beta0*AGE[i] + (beta1/2.0)*pow(AGE[i],2.0) +
                                             (beta2/3.0)*pow(AGE[i],3.0) +
                                             epsilon[i]*AGE[i]) -
                                exp(gammav)*(5.0*beta0 + 12.5*beta1 +
                                             (125.0/3.0)*beta2 + 5.0*epsilon[i]))

        ## Expectation
        Lambdaj[i] <- Ej[i] * exp(delta_u[YEAR[i]]+gammah*HUMIDITY[i])

        ## Likelihood
        TICK_CLASS[i]~dticks(Lambdaj[i], phi)
    }

})


## Year as an integer variable
int_year <- as.integer(fticks$year-(min(fticks$year)-1))

## period for each year
yl <- unique(fticks[,c("year","lothar")])
pery <- yl$lothar[order(yl$year)]


## Constants
TicksConst <- list(nobs=nrow(fticks), nyears=max(int_year), YEAR=int_year, pery=pery)

## Data
TicksData <- list(TICK_CLASS=as.integer(fticks$ticks),
                  AGE=fticks$age,
                  xbar=mean(fticks$age),
                  AGEc=fticks$age-mean(fticks$age),
                  SURFACE=fticks$bsa,
                  HUMIDITY=fticks$humidity-mean(fticks$humidity))

## Starting values
inits <- list(beta0p = 0.17, beta1p = 0.01, beta2p = -0.0007,
              gammav = -0.8, gammah=0, intercept = 2.8, phi = 0.79,
              sigma_u = c(0.64, 0.64, 0.64), sigma_s = 0.02,
              delta_u = c(2.07, 2.91, 2.79, 
                          2.82, 2.57, 3.47, 2.85, 2.77, 2.59, 2.4, 1.45, 3.01, 2.5, 
                          2.68, 2.76, 2.76, 2.71, 3.41, 2.43, 2.4, 2.59, 2.19, 2.63, 
                          1.72, 2.12, 2.59, 2.19))

## We use the same for the 4 chains
linits <- list(inits,inits,inits,inits)


## ----main-fit-model-nimble, eval=FALSE--------------------
## set.seed(123)
## resultsModel <-
##     nimbleMCMC(code = codePaperModel, constants = TicksConst,
##                data = TicksData, inits = linits,
##                nchains = 4, niter = 51000, nburnin=1000, thin=20,
##                summary = TRUE, samplesAsCodaMCMC=TRUE,
##                monitors = c("delta_u","gammav", "gammah", "sigma_u",
##                             "intercept",
##                             "phi","sigma_s","beta0p","beta1p","beta2p"))


## ----graphical-display-chains,eval=FALSE------------------
## library(bayesplot)
## mcmc_trace(resultsModel$samples,
##            regex_pars=c("intercept", "phi", "sigma_s","sigma_u","gammav","gammah",
##                         "beta0p","beta1p","beta2p"))


## ----gelman-main-results----------------------------------
library(coda)
gelman.diag(resultsModel$samples)


## ----simulate-data, eval=FALSE----------------------------
## simulatedData <- simulatePaperModel(resultsModel$samples, TicksData, TicksConst)


## ----number-in-each-class---------------------------------
## Observed:
sapply(1:3,function(i) mean(TicksData$TICK_CLASS==i))

## Simulated credible intervals
sapply(1:3, function(i) {quantile(rowMeans(simulatedData==i),c(0.025, 0.975))})


## ----proportion-in-each-class-----------------------------
## Class 1: are all the observed proportion within the limits of the 95% CI?
observed <- (tapply(TicksData$TICK_CLASS==1, TicksData$AGE, mean))
sim <- apply(simulatedData,1,function(x) tapply(x==1, TicksData$AGE, mean))
qs <- apply(sim,1,function(x) quantile(x, c(0.025,0.975)))
(observed>=qs[1,]&observed<=qs[2,])

## Class 2: are all the observed proportion within the limits of the 95% CI?
observed <- (tapply(TicksData$TICK_CLASS==2, TicksData$AGE, mean))
sim <- apply(simulatedData,1,function(x) tapply(x==2, TicksData$AGE, mean))
qs <- apply(sim,1,function(x) quantile(x, c(0.025,0.975)))
(observed>=qs[1,]&observed<=qs[2,])

## Class 3: are all the observed proportion within the limits of the 95% CI?
observed <- (tapply(TicksData$TICK_CLASS==3, TicksData$AGE, mean))
sim <- apply(simulatedData,1,function(x) tapply(x==3, TicksData$AGE, mean))
qs <- apply(sim,1,function(x) quantile(x, c(0.025,0.975)))
(observed>=qs[1,]&observed<=qs[2,])


## ---------------------------------------------------------
## Class 1
ry <- apply(simulatedData, 1, function(x) tapply(x, TicksConst$YEAR, function(y) mean(y==1)))
rg <- apply(ry,1,quantile,c(0.025,0.975))
ta <- tapply(TicksData$TICK_CLASS, TicksConst$YEAR, function(y) mean(y==1))
mean(ta>=rg[1,]&ta<=rg[2,])

## Class 2
ry <- apply(simulatedData, 1, function(x) tapply(x, TicksConst$YEAR, function(y) mean(y==2)))
rg <- apply(ry,1,quantile,c(0.025,0.975))
ta <- tapply(TicksData$TICK_CLASS, TicksConst$YEAR, function(y) mean(y==2))
mean(ta>=rg[1,]&ta<=rg[2,])

## Class 3
ry <- apply(simulatedData, 1, function(x) tapply(x, TicksConst$YEAR, function(y) mean(y==3)))
rg <- apply(ry,1,quantile,c(0.025,0.975))
ta <- tapply(TicksData$TICK_CLASS, TicksConst$YEAR, function(y) mean(y==3))
mean(ta>=rg[1,]&ta<=rg[2,])


## ----summary-coefs----------------------------------------
ro <- round(resultsModel$summary$all.chains,3)
ro[-grep("delta",rownames(ro)),]


## ---------------------------------------------------------
codePaperModel0 <- nimbleCode({
    
    ## Priors
    intercept ~ dnorm(0,sd=20)
    gammav ~ dnorm(0,sd=20)
    beta0p ~ dnorm(0,sd=20)
    beta1p ~ dnorm(0,sd=20)
    beta2p ~ dnorm(0,sd=20)
    for (i in 1:3) {
        sigma_u[i]~dunif(0,100)
    }
    sigma_s~dunif(0,100)
    phi~dunif(0,100);

    
    ## Year random effects
    for (i in 1:nyears) {
        delta_u[i] ~ dnorm(intercept, sd=sigma_u[pery[i]]);
    }

    ## Because we work on centred data for body surface area
    ## (for better mixing)
    beta0 <- beta0p - beta1p*xbar + beta2p*pow(xbar,2.0)
    beta1 <- beta1p - 2*beta2p*xbar
    beta2 <- beta2p

    
    ## Integration of body surface area
    for (i in 1:nobs) {
        
        ## Expected bsa
        yhat[i] <- beta0p + beta1p*AGEc[i] + beta2p*pow(AGEc[i],2.0)
        SURFACE[i] ~ dnorm(yhat[i], sd=sigma_s)
        epsilon[i] <- SURFACE[i] - yhat[i]

        ## Integration
        Ej[i] <- step(5-AGE[i])*(beta0*AGE[i] + (beta1/2.0)*pow(AGE[i],2.0) +
                                 (beta2/3.0)*pow(AGE[i],3.0) +
                                 epsilon[i]*AGE[i]) +
            (1-step(5-AGE[i]))*((5.0*beta0 + 12.5*beta1 + (125.0/3.0)*beta2 + 5.0*epsilon[i]) +
                                exp(gammav)*(beta0*AGE[i] + (beta1/2.0)*pow(AGE[i],2.0) +
                                             (beta2/3.0)*pow(AGE[i],3.0) +
                                             epsilon[i]*AGE[i]) -
                                exp(gammav)*(5.0*beta0 + 12.5*beta1 +
                                             (125.0/3.0)*beta2 + 5.0*epsilon[i]))

        ## Expectation
        Lambdaj[i] <- Ej[i] * exp(delta_u[YEAR[i]])

        ## Likelihood
        TICK_CLASS[i]~dticks(Lambdaj[i], phi)
    }

})


TicksData0 <- TicksData
TicksData0$HUMIDITY <- NULL

## Starting values
inits <- list(beta0p = 0.17, beta1p = 0.01, beta2p = -0.0007,
              gammav = -0.8, intercept = 2.8, phi = 0.79,
              sigma_u = c(0.64, 0.64, 0.64), sigma_s = 0.02,
              delta_u = c(2.07, 2.91, 2.79, 
                          2.82, 2.57, 3.47, 2.85, 2.77, 2.59, 2.4, 1.45, 3.01, 2.5, 
                          2.68, 2.76, 2.76, 2.71, 3.41, 2.43, 2.4, 2.59, 2.19, 2.63, 
                          1.72, 2.12, 2.59, 2.19))

## We use the same for the 4 chains
linits <- list(inits,inits,inits,inits)



## ----main-fit-model-nimble-0, eval=FALSE------------------
## set.seed(123)
## resultsModel0 <-
##     nimbleMCMC(code = codePaperModel0, constants = TicksConst,
##                data = TicksData0, inits = linits,
##                nchains = 4, niter = 51000, nburnin=1000, thin=20,
##                summary = TRUE, samplesAsCodaMCMC=TRUE,
##                monitors = c("delta_u","gammav",  "sigma_u",
##                             "intercept",
##                             "phi","sigma_s","beta0p","beta1p","beta2p"))
## 


## ----waic-humidity, eval=FALSE----------------------------
## withHumidity <- waicTicks(resultsModel$samples, TicksData, TicksConst)
## withoutHumidity <- waicTicks(resultsModel0$samples, TicksData0, TicksConst)


## ---------------------------------------------------------
withHumidity$waic
withoutHumidity$waic


## ---------------------------------------------------------
plot(resultsModel$summary$all.chains[-31,2],
     resultsModel0$summary$all.chains[,2],
     xlim=c(0,5), ylim=c(0,5), asp=1,
     xlab="Parameters of the model with humidity",
     ylab="Parameters of the model without humidity")
abline(0,1)


## ---------------------------------------------------------
hist(apply(simulatedData,1,function(x) cor(x, TicksData$TICK_CLASS)),
     main="Spearman Correlation",
     xlab="Correlation coefficient")


## ---------------------------------------------------------
rand <- t(sapply(1:999, function(i) {
    x <- sample(TicksData$TICK_CLASS)
    c(mean(x[TicksData$TICK_CLASS==1]==1),
      mean(x[TicksData$TICK_CLASS>2]>2),
      mean(x[TicksData$TICK_CLASS==3]==3))
}))

obe <- t(apply(simulatedData,1,function(x) {
    c(mean(x[TicksData$TICK_CLASS==1]==1),
      mean(x[TicksData$TICK_CLASS>2]>2),
      mean(x[TicksData$TICK_CLASS==3]==3))      
}))

sran <- apply(round(rand,2),2,
              function(x)
    paste0(round(mean(x),2), " (SE = ",
           round(sd(x),2), ")"))

sobe <- apply(round(obe,2),2,
              function(x)
    paste0(round(mean(x),2), " (SE = ",
           round(sd(x),2), ")"))
data.frame(Class=1:3,
           observed=sobe,
           UnderRandomAssumption=sran)


## ----eval=FALSE-------------------------------------------
## en <- predictEN(resultsModel$samples, TicksData, TicksConst)


## ---------------------------------------------------------
lam <- en$EN
lr <- range(unlist(lam))
age <- en$age
plot(age, lam[1,], ty="n", ylim=lr, xlab="Age (days)", ylab="Mean number of ticks")
tmp <- lapply(1:nrow(lam), function(i) lines(age, lam[i,], col=rgb(0.1,0.1,0.1, 0.02)))
mo <- apply(lam,2,mean)
lines(age, mo, lwd=5)
lines(age, mo, lwd=3, col="yellow")


## ---------------------------------------------------------
library(ggplot2)
sam <- do.call(rbind,resultsModel$samples)
delt <- apply(sam[,grep("delta",colnames(sam))],2,function(x) x-sam[,"intercept"])
dfdel <- data.frame(x=factor(rep(1992:(1992+26), each=nrow(delt))), y=as.vector(delt))

plo <- ggplot(dfdel, ggplot2::aes(y = y, 
                                           x = x)) +
    geom_violin(draw_quantiles = c(0.1,0.9),trim=TRUE) +
    geom_boxplot(width = 0.2, fill = "grey", 
                 outlier.shape = NA)+ggplot2::ylim(c(-3,3))+
    ylab("Random effects") + 
    xlab("Year")+
    geom_vline(xintercept=c(8.5,18.5),col="red", lwd=1.5)+
    theme_bw()
    
suppressWarnings(print(plo))


## ----code-one-more-variable-------------------------------
codePaperModelP1 <- nimbleCode({
    
    ## Priors
    intercept ~ dnorm(0,sd=20)
    gammah ~ dnorm(0,sd=20)
    gammav ~ dnorm(0,sd=20)
    gamman ~ dnorm(0,sd=20)
    beta0p ~ dnorm(0,sd=20)
    beta1p ~ dnorm(0,sd=20)
    beta2p ~ dnorm(0,sd=20)
    for (i in 1:3) {
        sigma_u[i]~dunif(0,100)
    }
    sigma_s~dunif(0,100)
    phi~dunif(0,100);

    
    ## Year random effects
    for (i in 1:nyears) {
        delta_u[i] ~ dnorm(intercept, sd=sigma_u[pery[i]]);
    }

    ## Because we work on centred data for body surface area
    ## (for better mixing)
    beta0 <- beta0p - beta1p*xbar + beta2p*pow(xbar,2.0)
    beta1 <- beta1p - 2*beta2p*xbar
    beta2 <- beta2p

    
    ## Integration of body surface area
    for (i in 1:nobs) {
        
        ## Expected bsa
        yhat[i] <- beta0p + beta1p*AGEc[i] + beta2p*pow(AGEc[i],2.0)
        SURFACE[i] ~ dnorm(yhat[i], sd=sigma_s)
        epsilon[i] <- SURFACE[i] - yhat[i]

        ## Integration
        Ej[i] <- step(5-AGE[i])*(beta0*AGE[i] + (beta1/2.0)*pow(AGE[i],2.0) +
                                 (beta2/3.0)*pow(AGE[i],3.0) +
                                 epsilon[i]*AGE[i]) +
            (1-step(5-AGE[i]))*((5.0*beta0 + 12.5*beta1 + (125.0/3.0)*beta2 + 5.0*epsilon[i]) +
                                exp(gammav)*(beta0*AGE[i] + (beta1/2.0)*pow(AGE[i],2.0) +
                                             (beta2/3.0)*pow(AGE[i],3.0) +
                                             epsilon[i]*AGE[i]) -
                                exp(gammav)*(5.0*beta0 + 12.5*beta1 +
                                             (125.0/3.0)*beta2 + 5.0*epsilon[i]))

        ## Expectation
        Lambdaj[i] <- Ej[i] * exp(delta_u[YEAR[i]]+gammah*HUMIDITY[i]+gamman*VARIABLES[i])


        ## Likelihood
        TICK_CLASS[i]~dticks(Lambdaj[i], phi)
    }

})


## Data
TicksDataP1 <- TicksData
TicksDataP1$VARIABLES <- fticks$temperature-mean(fticks$temperature)

## Starting values
inits$gammah <- 0
inits$gamman <- 0

## We use the same for the 4 chains
linits <- list(inits,inits,inits,inits)



## ----eval=FALSE-------------------------------------------
## set.seed(123)
## resultsModelTemperature <-
##     nimbleMCMC(code = codePaperModelP1, constants = TicksConst,
##                data = TicksDataP1, inits = linits,
##                nchains = 4, niter = 51000, nburnin=1000, thin=20,
##                summary = TRUE, samplesAsCodaMCMC=TRUE,
##                monitors = c("delta_u","gammav", "sigma_u","gamman","gammah",
##                             "intercept","phi","sigma_s","beta0p",
##                             "beta1p","beta2p"))
## 


## ---------------------------------------------------------
round(resultsModelTemperature$summary$all.chains["gamman",],3)


## ----echo=FALSE-------------------------------------------
abun <- tapply(fticks$density, fticks$year, mean)
plot(as.numeric(names(abun)),
     abun,
     xlab="Year",ylab="Roe Deer Abundance", ty="b")


## ----eval=FALSE-------------------------------------------
## TicksDataP2 <- TicksData
## TicksDataP2$VARIABLES <- fticks$density-mean(fticks$density)
## 
## set.seed(123)
## resultsModelDensity <-
##     nimbleMCMC(code = codePaperModelP1, constants = TicksConst,
##                data = TicksDataP2, inits = linits,
##                nchains = 4, niter = 51000, nburnin=1000, thin=20,
##                summary = TRUE, samplesAsCodaMCMC=TRUE,
##                monitors = c("delta_u","gammav", "sigma_u","gamman",
##                             "intercept","phi","sigma_s","beta0p",
##                             "beta1p","beta2p"))
## 


## ---------------------------------------------------------
resultsModelDensity$summary$all.chains["gamman",]


## ---------------------------------------------------------
codePaperModelHab <- nimbleCode({
    
    ## Priors
    intercept ~ dnorm(0,sd=20)
    gammah ~ dnorm(0,sd=20)
    gammav ~ dnorm(0,sd=20)
    beta0p ~ dnorm(0,sd=20)
    beta1p ~ dnorm(0,sd=20)
    beta2p ~ dnorm(0,sd=20)
    for (i in 1:3) {
        sigma_u[i]~dunif(0,100)
    }
    for (h in 1:5) {
        gammahab[h]~dnorm(0,sigma_hab)
    }
    sigma_s~dunif(0,100)
    sigma_hab~dunif(0,100)
    phi~dunif(0,100);

    
    ## Year random effects
    for (i in 1:nyears) {
        delta_u[i] ~ dnorm(intercept, sd=sigma_u[pery[i]]);
    }

    ## Because we work on centred data for body surface area
    ## (for better mixing)
    beta0 <- beta0p - beta1p*xbar + beta2p*pow(xbar,2.0)
    beta1 <- beta1p - 2*beta2p*xbar
    beta2 <- beta2p

    
    ## Integration of body surface area
    for (i in 1:nobs) {
        
        ## Expected bsa
        yhat[i] <- beta0p + beta1p*AGEc[i] + beta2p*pow(AGEc[i],2.0)
        SURFACE[i] ~ dnorm(yhat[i], sd=sigma_s)
        epsilon[i] <- SURFACE[i] - yhat[i]

        ## Integration
        Ej[i] <- step(5-AGE[i])*(beta0*AGE[i] + (beta1/2.0)*pow(AGE[i],2.0) +
                                 (beta2/3.0)*pow(AGE[i],3.0) +
                                 epsilon[i]*AGE[i]) +
            (1-step(5-AGE[i]))*((5.0*beta0 + 12.5*beta1 + (125.0/3.0)*beta2 + 5.0*epsilon[i]) +
                                exp(gammav)*(beta0*AGE[i] + (beta1/2.0)*pow(AGE[i],2.0) +
                                             (beta2/3.0)*pow(AGE[i],3.0) +
                                             epsilon[i]*AGE[i]) -
                                exp(gammav)*(5.0*beta0 + 12.5*beta1 +
                                             (125.0/3.0)*beta2 + 5.0*epsilon[i]))

        ## Expectation
        Lambdaj[i] <- Ej[i] * exp(delta_u[YEAR[i]]+gammah*HUMIDITY[i]+gammahab[HABITAT[i]])


        ## Likelihood
        TICK_CLASS[i]~dticks(Lambdaj[i], phi)
    }

})


fticksh <- fticks[!is.na(fticks$habitat),]

## Year as an integer variable
int_yearh <- as.integer(fticksh$year-(min(fticksh$year)-1))

## period for each year
ylh <- unique(fticksh[,c("year","lothar")])
peryh <- ylh$lothar[order(ylh$year)]


## Constants
TicksConsth <- list(nobs=nrow(fticksh), nyears=max(int_yearh), YEAR=int_yearh, pery=peryh,
                    HABITAT=as.integer(fticksh$habitat))

## Data
TicksDatah <- list(TICK_CLASS=as.integer(fticksh$ticks),
                   AGE=fticksh$age,
                   xbar=mean(fticksh$age),
                   AGEc=fticksh$age-mean(fticksh$age),
                   SURFACE=fticksh$bsa,
                   HUMIDITY=fticksh$humidity-mean(fticksh$humidity))

## Starting values
inits <- list(beta0p = 0.17, beta1p = 0.01, beta2p = -0.0007, gammahab=rep(0,5),
              gammav = -0.8, gammah=0, intercept = 2.8, phi = 0.79,
              sigma_u = c(0.64, 0.64, 0.64), sigma_s = 0.02,
              delta_u = c(2.07, 2.91, 2.79, 
                          2.82, 2.57, 3.47, 2.85, 2.77, 2.59, 2.4, 1.45, 3.01, 2.5, 
                          2.68, 2.76, 2.76, 2.71, 3.41, 2.43, 2.4, 2.59, 2.19, 2.63, 
                          1.72, 2.12, 2.59, 2.19)[1:23])

## We use the same for the 4 chains
linits <- list(inits,inits,inits,inits)



## ----eval=FALSE-------------------------------------------
## set.seed(123)
## resultsModelHabitat <-
##     nimbleMCMC(code = codePaperModelHab, constants = TicksConsth,
##                data = TicksDatah, inits = linits,
##                nchains = 4, niter = 51000, nburnin=1000, thin=20,
##                summary = TRUE, samplesAsCodaMCMC=TRUE,
##                monitors = c("delta_u","gammav", "sigma_u","gammahab",
##                             "intercept","phi","sigma_s","beta0p",
##                             "beta1p","beta2p"))
## 


## ---------------------------------------------------------
re <- resultsModelHabitat$summary$all.chains
round(re[grep("gammahab",rownames(re)),],2)

