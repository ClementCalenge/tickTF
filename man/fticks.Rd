\name{fticks}
\alias{fticks}
\docType{data}
\title{
  Tick Burden in Roe Deer Fawns in the Trois-Fontaines Forest (France)
}
\description{
  This data.frame contains the data used in Bariod et al. (in prep.) to
  model the infection process of roe deer fawns by ticks in the
  Trois-Fontaines forest (Northeastern France).
}
\usage{data("fticks")}
\format{
  A data frame with 965 observations on the following 9 variables.
  \describe{
    \item{\code{ticks}}{an integer vector containing the tick-burden
      class of each deer (1 = less than 10 ticks; 2 = between 10 and 20
      ticks; 3 = more than 10 ticks).}
    \item{\code{lothar}}{an integer vector containing the period to
      which the year of capture belong 
      (defined on the basis of the forest structure following the
      hurricane Lothar). Three periods are defined: 1 corresponds to the
      year prior to the Lothar hurricane in 1999; 2 corresponds to the 10
      years following the Lothar hurricane; 3 corresponds to the period
      2010--2018;}
    \item{\code{year}}{a numeric vector containing the year of capture
      of the fawns}
    \item{\code{age}}{a numeric vector containing the age of the
      captured fawns}
    \item{\code{bsa}}{a numeric vector containing the body surface area
      of the captured fawns.}
    \item{\code{temperature}}{a numeric vector containing the mean
      temperature over the 6 days prior to the capture.}
    \item{\code{humidity}}{a numeric vector containing the mean
      humidity over the 6 days prior to the capture.}
    \item{\code{habitat}}{an integer vector containing the habitat where
      the fawns were captured (1 = low scrub; 2 = high scrub; 3 = clear
      forest; 4 = dense and young forest; 5 = dense and older forest
      with continuous canopy) }
    \item{\code{density}}{a numeric vector containing the estimated
      roe deer density for the year of capture.}
  }
}
\source{
  Office francais de la biodiversite, Direction de la Recherche et de
  l'Appui Scientifique, Office Francais de la Biodiversite, Birieux,
  France 
}
\references{
  Bariod L., Said S., Calenge C., Chabot S., Badeaud V. \& Bourgoin
  G. in prep. Infection of Roe Deer Fawns by Ticks from 1992 to 2018 in the
  Trois-Fontaines Forest (France)
}
\examples{
data(fticks)
head(fticks)
}
\keyword{datasets}
