\name{simulatePaperModel}
\alias{simulatePaperModel}
\title{
  Simulate Tick-Infection Datasets Using Fitted Models
}
\description{
  This function can be used to simulate one dataset per MCMC iteration
  using the final model describing the tick infection of roe deer fawns
  in the Trois-Fontaines forest.
}
\usage{
simulatePaperModel(mcmcList, TicksData, TicksConst)
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
  A matrix with N rows (the N fawns of the dataset) and R columns (the R
  MCMC iterations) containing, for each MCMC iteration, a simulation of
  the tick-burden class (from 1 to three, see \code{?fticks}) according
  to the fitted model.
}
\references{
  Bariod L., Said S., Calenge C., Chabot S., Badeaud V. \& Bourgoin
  G. in prep. Infection of Roe Deer Fawns by Ticks from 1992 to 2018 in the
  Trois-Fontaines Forest (France).  
}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\seealso{
  \code{\link{fticks}} for a description of the dataset used, and
  \code{\link{resultsModel}} for a description of the final models (and
  alternative models).  See also the vignette of the package for an
  example of use.
}
\examples{

## see vignette.

}
\keyword{ model }
