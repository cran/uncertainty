\name{print.uncertaintyBudget}
\alias{print.uncertaintyBudget}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prints an uncertainty budget object
}
\description{
Print the description of each uncertainty source
}
\usage{
\method{print}{uncertaintyBudget}(x, \dots)
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
	\code{\link{uncertaintyBudget.default}}, \code{\link{print}}
}
\examples{
cor.mat<- matrix(c(1,-0.7,-0.7,1),2,2)

u.budget<- uncertaintyBudget(x=list(name=c("x0","x1"), 
	mean=c(10,20), u=c(1,5), dof=c(10,10),
	label=c("x[0]", "x[1]"), distribution=c("normal","normal")), 
	y=cor.mat)

# implicitly calls the print method
u.budget

# same as
print(u.budget)

# uncertainty budget structure
attributes(u.budget)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ uncertainty }
\keyword{ budget }
\keyword{ print }
% __ONLY ONE__ keyword per line
