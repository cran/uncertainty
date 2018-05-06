\name{uncertainty}
\alias{uncertainty}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates an uncertainty object
}
\description{
Builds an uncertainty estimation object using a measurand model and an uncertainty budget object
}
\usage{
uncertainty(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  an uncertainty budget object
}
  \item{\dots}{
  additional parameters
}
}
\details{
	Creates an uncertainty estimation object. Uses an uncertainty budget object to estimate the expected value and uncertainty of a measurand by applying a selected estimation method.
}
\value{
	An uncertainty estimation object
}
\references{
JCGM 100:2008. \emph{Guide to the expression of uncertainty of measurement}

JCGM 100:2005. \emph{Supplement 1 Propagation of distributions usign a Monte Carlo method}

EURACHEM/CITAC Guide CG 4. \emph{Quantifying Uncertainty in Analytical Measurement}

Becker, R.A., Chambers, J.M. and Wilks, A.R. (1988) \emph{The New S Language}. Wadsworth & Brooks/Cole.
}
\author{
H. Gasca-Aragon

Maintainer: H. Gasca-Aragon <hugo_gasca_aragon@hotmail.com>
}
\note{
	none
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{uncertainty.default}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ uncertainty }
% __ONLY ONE__ keyword per line
