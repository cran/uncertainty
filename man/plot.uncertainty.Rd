\name{plot.uncertainty}
\alias{plot.uncertainty}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots a probability density function related to the measurand model
}
\description{
Plot a probability density function attributed to the measurand, depending on the selected method to estimate the uncertainty.
}
\usage{
\method{plot}{uncertainty}(x, y = NULL, xlab = parse(text = x$measurand_label), 
main = "", ylab = "Probability density", from = x$mean - 4 * x$u, to = x$mean + 4 * x$u, 
lwd = 2, add = FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
     An uncertainty object
}
  \item{y}{
	not used, exists only for compatibility with the S3 generic function.
}
  \item{xlab}{
	string or expression, label for the x-axis.
}
  \item{main}{
	string or expression, label for the plot.
}
  \item{ylab}{
	string or expression, label for the y-axis.
}
  \item{from}{
	numeric, lower value of the x-axis to display.
}
  \item{to}{
	numeric, upper valur of the x-axis to display.
}
  \item{lwd}{
	numeric, line width.
}
  \item{add}{
	logic, decides to add the curve into an existing plot or to create a new plot.
}
  \item{\dots}{
	additional parameters.
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
	\code{\link{uncertainty.default}}, \code{\link{plot}}
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

# plot the estimated pdf of the measurand
\dontrun{plot(GFO.res)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ uncertainty }
\keyword{ plot }
% __ONLY ONE__ keyword per line
