\name{uncertaintyBudget.default}
\alias{uncertaintyBudget.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generic function for calling an uncertainty budget object
}
\description{
Creates an uncertainty budget.
}
\usage{
\method{uncertaintyBudget}{default}(x, y, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	a list with the vector entries name, label, mean, u(uncertainty), distribution and dof, one for each input quantity.
}
  \item{y}{
	a correlation matrix of the input quantities, interpreted in the same order of input quantities as the vector name
}
  \item{\dots}{
     additional parameters
}
}
\details{
	Creates an uncertainty budget object
}
\value{
	An uncertainty budget object with attributes:

	\code{name} the name of each input quantity

	\code{mean} the mean value of each input quantity

	\code{u} the uncertainty of each input quantity

	\code{dof}the degrees of freedom of each input quantity

	\code{label}the label of each input quantity

	\code{distribution} the distribution of each input quantity, valid values are (bernoulli, beta, binomial, cuachy, chisq, exp, f, gamma, lognormal, poission, normal, unif, t, traingular, weibull, arcsine)

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
	\code{\link{uncertaintyBudget}}, \code{\link{uncertainty}}, \code{\link{print.uncertaintyBudget}}
}
\examples{
require(mvtnorm)

cor.mat<- matrix(c(1,-0.7,-0.7,1),2,2)

u.budget<- uncertaintyBudget(x=list(name=c("x0","x1"), 
	mean=c(10,20), u=c(1,5), dof=c(10,10),
	label=c("x[0]", "x[1]"), distribution=c("normal","normal")), y=cor.mat)
u.budget
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ uncertainty }
\keyword{ budget }
\keyword{ distribution }
% __ONLY ONE__ keyword per line
