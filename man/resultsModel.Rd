\name{resultsModel}
\alias{resultsModel}
\alias{resultsModel0}
\alias{resultsModelTemperature}
\alias{resultsModelDensity}
\alias{resultsModelHabitat}
\docType{data}
\title{
  Datasets containing Models Fitted by Bariod et al. (in prep.)
}
\description{
  These datasets contain the MCMC iterations corresponding to the
  various models fitted to describe the infection process of roe deer
  fawns by ticks in the Trois-Fontaines forest (France).
}
\usage{

data("resultsModel")

data("resultsModel0")

data("resultsModelTemperature")

data("resultsModelDensity")

data("resultsModelHabitat")

}
\format{
  All these objects were returned by the function \code{nimbleMCMC} from
  the package nimble (see the vignette of the package).
}
\details{
  Each dataset contains the MCMC iterations obtained by fitting one of
  the models studied by Bariod et al. to the dataset \code{fticks}. The
  vignette of the package describes how these datasets were built. More
  precisely:

  \code{resultsModel} contains the  results of the fit of the final
  model (presented in the paper), including the humidity as a variable
  affecting the infectious detection rate.

  \code{resultsModel0} contains the results of the fit of the final model
  (presented in the paper), but without the humidity as a variable
  affecting the infectious detection rate.

  \code{resultsModel0} contains the results of the fit of the same model
  as \code{resultsModel}, but with the mean temperature as an additional
  variable affecting the infectious detection rate.

  \code{resultsModelDensity} contains the results of the fit of the same
  model as \code{resultsModel}, but with the roe deer density of the
  year of capture as an additional variable affecting the infectious
  detection rate.

  \code{resultsModelHabitat} contains the results of the fit of the same
  model as \code{resultsModel}, but with the habitat type of the site of
  capture as an additional random effect affecting the infectious
  detection rate.

}
\source{
  The original Dataset has been kindly provided by the Office Francais
  de la Biodiversite:
  
  Office francais de la biodiversite, Direction de la Recherche et de
  l'Appui Scientifique, Office Francais de la Biodiversite, Birieux,
  France

}
\references{
  Bariod L., Said S., Calenge C., Chabot S., Badeaud V. \& Bourgoin
  G. in prep. Infection of Roe Deer Fawns by Ticks from 1992 to 2018 in the
  Trois-Fontaines Forest (France).
}
\examples{
## Structure of the dataset
str(resultsModel)

## See the vignette for a description of how these results where
## obtained
}
\keyword{datasets}
