predictEN <-
function(mcmcList, TicksData, TicksConst){

    if (!inherits(mcmcList, "mcmc.list"))
        stop("mcmcList should be of class \"mcmc.list\"")
    samf <- do.call(rbind, mcmcList)

    res <- sapply(1:nrow(samf), function(i) {

        x <- samf[i,]
        beta0p <- x["beta0p"]
        beta1p <- x["beta1p"]
        beta2p <- x["beta2p"]
        gammav <- x["gammav"]
        intercept <- x["intercept"]
        agb <- seq(0,8,length=100)
        xbar <- TicksData$xbar
        beta0 <- beta0p - beta1p*xbar + beta2p*(xbar^2)
        beta1 <- beta1p - 2*beta2p*xbar
        beta2 <- beta2p
        epsilon <- 0

        uu <- sapply(agb, function(AGE) {
            coefintegr <- as.numeric((5-AGE)>=0)*(beta0*AGE + (beta1/2.0)*(AGE^2) +
                                                  (beta2/3.0)*(AGE^3) + epsilon*AGE) +
                ((5-AGE)<0)*((5.0*beta0 + 12.5*beta1 + (125.0/3.0)*beta2 + 5.0*epsilon) +
                             exp(gammav)*(beta0*AGE + (beta1/2.0)*(AGE^2) +
                                               (beta2/3.0)*(AGE^3) +
                                          epsilon*AGE) - exp(gammav)*
                             (5.0*beta0 + 12.5*beta1 + (125.0/3.0)*beta2 + 5.0*epsilon))

            LAMBDA <- coefintegr * exp(intercept)
            return(LAMBDA)
        })

        return(uu)
    })
    return(list(age=seq(0,8,length=100), EN=t(res)))
}
