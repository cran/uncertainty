\name{uncertaintyBudget}
\alias{uncertaintyBudget}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generic function for uncertainty budget object
}
\description{
Generic function for creating an uncertainty budget object
}
\usage{
uncertaintyBudget(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a list with the vector entries name, label, mean, u(uncertainty), distribution and dof, one for each quantity.
}
  \item{\dots}{
  additional parameters
}
}
\details{
	\code{uncertaintyBudget} is a generic function (under S3 protocol) for searching the default method.
}
\value{
	An uncertainty budget object with attributes:

	\code{name} the name of each input quantity

	\code{mean} the mean value of each input quantity

	\code{u} the uncertainty of each input quantity

	\code{dof} the degrees of freedom of each input quantity

	\code{label} the label of each input quantity

	\code{distribution} the distribution of each input quantity, valid values are (normal, unif, t, chisq, f, triangle, binomial, bernoulli, beta, gamma)

	\code{cor} the correlation matrix among the input quantities
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
	\code{\link{uncertaintyBudget.default}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ uncertainty }
\keyword{ budget }
% __ONLY ONE__ keyword per line
