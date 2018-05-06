\name{uncertainty.default}
\alias{uncertainty.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generic function for calling an uncertainty object
}
\description{
Creates an uncertainty estimation object using a measurand model and an uncertainty budget object
}
\usage{
\method{uncertainty}{default}(x, y, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	an uncertainty budget object
}
  \item{y}{
	a list with the measurand description and selected estimation method, the measurand description includes: measurand_name, measurand_model, measurand_label, alpha (significance level), method and method parameters.
	
	the valid methods are: GFO, GSO, MC, Kragten.

	currently the only method parameter implemented is the number of simulated samples (B) for the method MC.
}
  \item{\dots}{
     additional parameters
}
}
\details{
	Creates an uncertainty estimation object. Uses an uncertainty budget object to estimate the expected value and uncertainty of a measurand by applying a selected estimation method.
}
\value{
  An uncertainty estimation object with the structure:
  \code{method} selected estimating method, 
  \code{call} current call invocation,
  \code{uncertaintyBudget} an uncertainty budget object,
  \code{measurand} name, label, model describing the measurand,
  \code{mean} the estimated mean,
  \code{sd} the estimated standard deviation,
  \code{u} the estimated standard uncertainty,
  \code{alpha} the significante level used in the estimation,
  \code{dof} the estimated degrees of freedom,
  \code{U} the estimated expanded uncertainty,
  \code{lcl} the lower confidence interval,
  \code{ucl} the upper confidence interval,
  \code{variables} a vector with the input quantities,
  \code{contribution} a vector with the uncertainty contributions,
  \code{cor.contribution} the uncertainty contribution due to overall correlation,
  \code{partial} a vector of the partial derivatives of the measurand.model with respect to each input quantity,
  \code{coeff} a vector of the sensibility coefficients for each input quantity.
}
\references{
JCGM 100:2008. \emph{Guide to the expression of uncertainty of measurement}

JCGM 100:2005. \emph{Supplement 1 Propagation of distributions usign a Monte Carlo method}

EURACHEM/CITAC Guide CG 4. \emph{Quantifying Uncertainty in Analytical Measurement}
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
	\code{\link{uncertainty}}, \code{\link{uncertaintyBudget.default}}, \code{\link{print.uncertainty}}, \code{\link{plot.uncertainty}}, \code{\link{summary.uncertainty}}
}
\examples{
# create an uncertainty budget
cor.mat<- matrix(c(1,-0.7,-0.7,1),2,2)

u.budget<- uncertaintyBudget(x=list(name=c("x0","x1"), 
	mean=c(10,20), u=c(1,5), dof=c(10,10),
	label=c("x[0]", "x[1]"), distribution=c("normal","normal")), 
	y=cor.mat)
u.budget

# estimate the measurand uncertainty using an uncertainty budget,
# a measurand definition and a selected estimating method.
GFO.res<- uncertainty(x=u.budget, 
y=list(measurand_name="ratio.GFO", measurand_label="ratio[GFO]", 
measurand_model="x0/x1", method="GFO", alpha=0.05))

GFO.res
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ uncertainty }
\keyword{ budget }
\keyword{ estimating method }
\keyword{ measurand modeling }
% __ONLY ONE__ keyword per line
