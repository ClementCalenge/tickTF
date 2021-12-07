\name{dticksModel}
\alias{dticksModel}
\title{
  Probability That a Fawn Belongs to a Given Tick Burden Class
}
\description{
  This function gives the probability that a fawn belongs to one of the
  three tick-burden classes (1 = less that 10 ticks; 2 = between 10 and
  20 ticks; 3 = more than 20 ticks), given that the number of ticks
  follows a negative binomial distribution.
}
\usage{
dticksModel(x, lambda, phi)
}
\arguments{
  \item{x}{
    a vector of integer values comprised between 1 and 3, giving the
    tick burden classes for the roe deer fawn.
  }
  \item{lambda}{
    a vector of the same length as \code{x} containing the mean of the
    negative binomial distribution for each value of \code{x}.
}
  \item{phi}{
    a real number. The dispersion parameter of the negative binomial
    distribution.
  }
}
\value{
  A real value corresponding to the probability of \code{x}
}
\references{
  Bariod L., Said S., Calenge C., Chabot S., Badeaud V. \& Bourgoin
  G. in prep. Infection of Roe Deer Fawns by Ticks from 1992 to 2018 in the
  Trois-Fontaines Forest (France)
}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\examples{


dticksModel(c(1,2,1), lambda=c(8,14,3), phi=0.9)

}
\keyword{ distribution }
