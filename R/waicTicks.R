waicTicks <-
function(mcmcList, TicksData, TicksConst){
    if (!inherits(mcmcList, "mcmc.list"))
        stop("mcmcList should be of class \"mcmc.list\"")
    if (any(colnames(mcmcList[[1]])=="gammahab[1]"))
        stop("waicTicks cannot be used on the model with habitat")

    samf <- do.call(rbind, mcmcList)

    res <- sapply(1:nrow(samf), function(i) {

        x <- samf[i,]

        beta0p <- x["beta0p"]
        beta1p <- x["beta1p"]
        beta2p <- x["beta2p"]
        gammav <- x["gammav"]
        delta_u <- x[grep("delta_u", colnames(samf))]
        addit <- rep(0,length(TicksData$TICK_CLASS))

        if (any(names(x)=="gammah")) {
            gammah <- x["gammah"]
            addit <- addit+gammah*TicksData$HUMIDITY
        }

        if (any(names(x)=="gamman")) {
            gamman <- x["gamman"]
            addit <- addit+gamman*TicksData$VARIABLES
        }

        intercept <- x["intercept"]
        phi <- x["phi"]
        sigma_u <- c(x["sigma_u[1]"],x["sigma_u[2]"],x["sigma_u[3]"])
        sigma_s <- x["sigma_s"]
        xbar <- TicksData$xbar
        beta0 <- beta0p - beta1p*xbar + beta2p*(xbar^2)
        beta1 <- beta1p - 2*beta2p*xbar
        beta2 <- beta2p

        epsilon <- TicksData$SURFACE-(beta0p+beta1p*TicksData$AGEc+
                                      beta2p*(TicksData$AGEc^2))

        AGE <- TicksData$AGE

        coefintegr <- as.numeric((5-AGE)>=0)*(beta0*AGE + (beta1/2.0)*(AGE^2) +
                                              (beta2/3.0)*(AGE^3) + epsilon*AGE) +
            ((5-AGE)<0)*((5.0*beta0 + 12.5*beta1 + (125.0/3.0)*beta2 + 5.0*epsilon) +
                         exp(gammav)*(beta0*AGE + (beta1/2.0)*(AGE^2) +
                                      (beta2/3.0)*(AGE^3) +
                                      epsilon*AGE) - exp(gammav)*
                         (5.0*beta0 + 12.5*beta1 + (125.0/3.0)*beta2 + 5.0*epsilon))

        lambda <- coefintegr * exp(delta_u[TicksConst$YEAR]+addit)

        pr <- dticksModel(TicksData$TICK_CLASS, lambda, phi)

        return(pr)
    })

    lppdiln <- log(colMeans(t(res)))
    pwaic2ln <- apply(log(t(res)),2,var)

    return(list(lppd=lppdiln, pwaic=pwaic2ln,
                waic=c(WAIC=-2*(sum(lppdiln)-sum(pwaic2ln)),
                       SE = sd(-2*(lppdiln-pwaic2ln)))))
}
