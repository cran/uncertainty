\name{summary.uncertainty}
\alias{summary.uncertainty}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates an uncertainty summary object
}
\description{
Performs an uncertainty contribution estimation for the uncertainty object. 
The metric used to measure the contribution is percent of variance.
If correlation is present an additional entry is shown with the whole contribution due to correlated input quantities.
}
\usage{
\method{summary}{uncertainty}(object, ndigits = 3, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
	an uncerainty object
}
  \item{ndigits}{
	numeric, the number of digits for displaying.
}
  \item{\dots}{
	additional parameters
}
}
\details{
	none
}
\value{
  An uncertainty summary object: 
  \item{call}{the call invocation}
  \item{measurand.name}{name of the measurand}
  \item{measurand.label}{label of the measurand for displaying purposes}
  \item{budget}{a list with the name, mean, label, u(uncertainty), dof and uncertainty contribution for each input quantity plus a correlation entry if any}
}
\references{
JCGM 100:2008. \emph{Guide to the expression of uncertainty of measurement}

JCGM 100:2005. \emph{Supplement 1 Propagation of distributions usign a Monte Carlo method}

EURACHEM/CITAC Guide CG 4. \emph{Quantifying Uncertainty in Analytical Measurement}

Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}. Fourth edition. Springer.
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
	\code{\link{uncertainty.default}}, \code{\link{print.summary.uncertainty}}, \code{\link{summary}}
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

# create an uncertainty summary object
GFO.sum<- summary(GFO.res)

# implicit call to the print method
GFO.sum

# same as
print(GFO.sum)

# uncertainty summary structure
attributes(GFO.sum)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ uncertainty }
\keyword{ summary }
% __ONLY ONE__ keyword per line
