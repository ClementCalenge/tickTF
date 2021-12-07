dticksModel <-
    function(x, lambda, phi) {
        if (!all(x%in%c(1:3)))
            stop("x should be an integer value comprised between 1 and 3")
        if (length(x)!=length(lambda))
            stop("lambda should have the same length as x")
        if (length(phi)>1)
            stop("phi should be of length 1")

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
        if (length(x)==1) {
            lpro <- (c(p1,p2,p3)[x])
        } else {
            cb <- cbind(p1,p2,p3)
            lpro <- sapply(1:length(x), function(i) cb[i,x[i]])
        }
        return(exp(lpro))

}
