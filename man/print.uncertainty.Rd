\name{print.uncertainty}
\alias{print.uncertainty}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Displays the detailed content of a measurand model including its uncertainty estimate.
}
\description{
Displays the estimated value of the measurand, its standard deviation, its standard uncertainty, the degrees of freedom and
the significance level and an CI with that significance level.
}
\usage{
\method{print}{uncertainty}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	an uncertainty object
}
  \item{\dots}{
	additional parameters
}
}
\details{
	none
}
\value{
	None (invisible NULL)
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
	\code{\link{uncertainty.default}}, \code{\link{print}}
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

# implicit call to print method
GFO.res

# same as
print(GFO.res)

# structure of an uncertainty estimation object
attributes(GFO.res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }
\keyword{ uncertainty }
% __ONLY ONE__ keyword per line
