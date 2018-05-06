\name{plot.summary.uncertainty}
\alias{plot.summary.uncertainty}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots an uncertainty summary. It shows the uncertainty contribution from each involved quantity
}
\description{
Builds a barplot with a bar for each source of uncertainty.
If correlation is present then an additional entry is added.
The current metric used to display is % of total variance.
When correlation is present its contribution may be negative.
}
\usage{
\method{plot}{summary.uncertainty}(x, y = NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an uncertainty summary object
}
  \item{y}{
not used.
}
  \item{\dots}{
additional parameters to customize the plot
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

JCGM 100:2005. \emph{Supplement 1 Propagation of distributions usign a Monte Carlo method.}

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
	\code{\link{summary.uncertainty}}, \code{\link{plot}}
}
\examples{
# create an uncertainty budget
cor.mat<- matrix(c(1,-0.7,-0.7,1),2,2)

u.budget<- uncertaintyBudget(x=list(name=c("x0","x1"), 
	mean=c(10,20), u=c(1,5), dof=c(10,10),
	label=c("x[0]", "x[1]"), distribution=c("normal","normal")), 
	y=cor.mat)

# estimate the measurand uncertainty using an uncertainty budget,
# a measurand definition and a selected estimating method.
GFO.res<- uncertainty(x=u.budget, 
y=list(measurand_name="ratio.GFO", measurand_label="ratio[GFO]", 
measurand_model="x0/x1", method="GFO", alpha=0.05))

# create an uncertainty summary object
GFO.sum<- summary(GFO.res)

# display the chart
\dontrun{plot(GFO.sum)}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ barplot }
\keyword{ summary }
\keyword{ uncertainty }
\keyword{ plot }
% __ONLY ONE__ keyword per line
