\name{print.summary.uncertainty}
\alias{print.summary.uncertainty}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Displays a list with the uncertainty contribution from each input quantity
}
\description{
For each input quantity (source of uncertainty) it shows the uncertainty contribution, 
measured in percent of variance of the measurand model.
}
\usage{
\method{print}{summary.uncertainty}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
     An uncertainty summary object
}
  \item{\dots}{
     Additional parameters
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
	\code{\link{summary.uncertainty}}, \code{\link{print}}
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
\keyword{ summary }
\keyword{ uncertainty }
\keyword{ print }
% __ONLY ONE__ keyword per line
