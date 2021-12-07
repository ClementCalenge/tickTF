\name{predictEN}
\alias{predictEN}
\title{
  Expectation of the Number of Ticks on a Roe Deer as a Function of Age
}
\description{
  This function calculates the expected number of ticks on an
  average-sized roe deer fawn of a given age in average environmental
  conditions, at various ages, as 
  predicted by an infection model.
}
\usage{
predictEN(mcmcList, TicksData, TicksConst)
}
\arguments{
\item{mcmcList}{
    An object of class \code{mcmc.list} containing the results of the
    fit of a model.
}
\item{TicksData}{
  The list containing the data passed to \code{nimbleMCMC}. WARNING: the
  function does not check that this list actually corresponds to the
  list passed to \code{nimbleMCMC} and trusts the user.
}
  \item{TicksConst}{
    The list containing the constants passed to \code{nimbleMCMC}. WARNING: the
  function does not check that this list actually corresponds to the
  list passed to \code{nimbleMCMC} and trusts the user.
}
}
\value{
  A list containing the following elements:
  \item{age}{a numeric vector giving the ages, from 1 to 8 days, at
    which the expected number is calculated}
  \item{EN}{a matrix with R rows (corresponding to the R MCMC
    iterations) and N columns (corresponding to the N values of age for
    which the expected number is calculated) giving the expected number
    of ticks.}
}
\references{
    Bariod L., Said S., Calenge C., Chabot S., Badeaud V. \& Bourgoin
  G. in prep. Infection of Roe Deer Fawns by Ticks from 1992 to 2018 in the
  Trois-Fontaines Forest (France).  

}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\examples{
\dontrun{
    en <- predictEN(resultsModel$samples, TicksData, TicksConst)

    lam <- en$EN
    lr <- range(unlist(lam))
    age <- en$age
    plot(age, lam[1,], ty="n", ylim=lr, xlab="Age (days)", ylab="Mean number of ticks")
    tmp <- lapply(1:nrow(lam), function(i) lines(age, lam[i,], col=rgb(0.1,0.1,0.1, 0.02)))
    mo <- apply(lam,2,mean)
    lines(age, mo, lwd=5)
    lines(age, mo, lwd=3, col="yellow")

}

}
\keyword{model }
