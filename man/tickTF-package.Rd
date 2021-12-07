\name{tickTF-package}
\alias{tickTF-package}
\alias{tickTF}
\docType{package}
\title{
\packageTitle{tickTF}
}
\description{
\packageDescription{tickTF}
}
\details{

  This package contains the data and functions used to describe the
  infection process of the roe deer fawns by the ticks described by
  Bariod et al. (in prep.).  The dataset \code{fticks} contains the
  original data, and the vignette of the package describes how the
  functions of the package were used to fit the model. This package
  relies strongly on the package nimble for the fit.
  
}
\author{
\packageAuthor{tickTF}

Maintainer: \packageMaintainer{tickTF}
}
\references{
  Bariod L., Said S., Calenge C., Chabot S., Badeaud V. \& Bourgoin
  G. in prep. Infection of Roe Deer Fawns by Ticks from 1992 to 2018 in the
  Trois-Fontaines Forest (France).  
}
\keyword{ package }
\seealso{
  \code{\link{fticks}} for the dataset used; \code{vignette("tickTF")}
  for a description of the calculations carried out in Bariod et al. (in
  prep.). 
}
\examples{
## dataset
head(fticks)
}
