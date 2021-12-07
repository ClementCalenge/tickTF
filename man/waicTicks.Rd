\name{waicTicks}
\alias{waicTicks}
\title{
  WAIC of Models Fitted to Describe the Infection of Roe Deer Fawns by Ticks
}
\description{
  This function is used to calculate the WAIC of the models fitted in
  Bariod et al. (in prep).
}
\usage{
waicTicks(mcmcList, TicksData, TicksConst)
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
  \item{lppd}{the log predictive pointwise density for each fawn
    in the dataset}
  \item{pwaic }{the contribution of each fawn to the WAIC in the
    dataset}
  \item{waic}{the value of the WAIC as well as its standard error}
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
    ## WAIC for the model with humidity
    withHumidity <- waicTicks(resultsModel$samples, TicksData, TicksConst)

    ## WAIC for the model without humidity
    withoutHumidity <- waicTicks(resultsModel0$samples, TicksData0, TicksConst)
}
}
\keyword{model }
