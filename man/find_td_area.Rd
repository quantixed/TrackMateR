% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateTrackDensity.R
\name{find_td_area}
\alias{find_td_area}
\title{Track Density - Find search area}
\usage{
find_td_area(r, xy, a, b)
}
\arguments{
\item{r}{radius of search circle}

\item{xy}{numeric vector (length 2)}

\item{a}{numeric vector (length 2)}

\item{b}{numeric vector (length 2)}
}
\value{
numeric variable
}
\description{
Find the area of the intersection of the circle centered at xy with radius r
and the radius with vertical sides at a and horizontal sides at b.
xy, a, and b must be vectors of length 2, and xy must lie within the rectangle.
}
\examples{
find_td_area(r=2, xy=c(4, 4), a=c(0, 8), b=c(0, 5))
}
\keyword{internal}
