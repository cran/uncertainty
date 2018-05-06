\name{uncertainty-package}
\alias{uncertainty-package}
\docType{package}
\title{
Uncertainty Estimation and Contribution Analysis
}
\description{
Uncertainty estimation and contribution analysis implemented by 4 methods:
the Gaussian method of first,
the Gaussian method of second order,
the Kragten numerical method and 
the Monte Carlo simulation method
}
\details{
\tabular{ll}{
Package: \tab uncertainty\cr
Type: \tab Package\cr
Version: \tab 0.1.1\cr
Date: \tab 2014-06-12\cr
License: \tab GPL (>=2)\cr
}
Define an "uncertainty budget" object, including all the involved variables.
Then estimate the "uncertainty" object by defining a measurand model, using the "uncertainty budget" and applying an estimation method.
Print or plot the measurand estimates or create a "summary uncertainty" object to print or plot the uncertainty contributions to the measurand model.
}
\author{
H. Gasca-Aragon

Maintainer: H. Gasca-Aragon <hugo_gasca_aragon@hotmail.com>
}
\references{
JCGM 100:2008. \emph{Guide to the expression of uncertainty of measurement}

JCGM 100:2005. \emph{Supplement 1 Propagation of distributions usign a Monte Carlo method}

EURACHEM/CITAC Guide CG 4. \emph{Quantifying Uncertainty in Analytical Measurement}

Becker, R.A., Chambers, J.M. and Wilks, A.R. (1988) \emph{The New S Language}. Wadsworth & Brooks/Cole.
}
\keyword{ measurand modeling }
\keyword{ uncertainty }
\keyword{ budget }
\keyword{ estimating method }
\seealso{
	\code{\link{uncertaintyBudget}}, \code{\link{print.uncertaintyBudget}}, \code{\link{uncertainty}}, \code{\link{print.uncertainty}}, \code{\link{plot.uncertainty}}, \code{\link{summary.uncertainty}}, \code{\link{print.summary.uncertainty}}, \code{\link{plot.summary.uncertainty}}
%#~~ Optional links to other man pages, e.g. ~~
%#~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{

require(mvtnorm)

cor.mat<- matrix(c(1,-0.7,-0.7,1),2,2)

u.budget<- uncertaintyBudget(x=list(name=c("x0","x1"), 
	mean=c(10,20), u=c(1,5), dof=c(10,10),
	label=c("x[0]", "x[1]"), distribution=c("normal","normal")), y=cor.mat)
u.budget


## Gaussian first order estimates
GFO.res<- uncertainty(x=u.budget,
y=list(measurand_name="ratio.GFO", 
measurand_label=expression(ratio[GFO]), 
measurand_model="x0/x1", 
method="GFO", alpha=0.05))

contr.GFO<- summary.uncertainty(GFO.res)

## Monte Carlo estimates
MC.res<- uncertainty(x=u.budget, 
y=list(measurand_name="ratio.MC", 
measurand_label=expression(ratio[MC]),
measurand_model="x0/x1", 
method="MC", alpha=0.05, B=1e5))

contr.MC<- summary.uncertainty(MC.res)

## print the estimates
MC.res
GFO.res

## print the uncertainty summary
contr.MC
contr.GFO

## Displaying both estimated distributions
\dontrun{
plot(MC.res, col=4, xlab=MC.res$measurand_model)
plot(GFO.res, lty=2, col=2, add=T)
legend(0.7, 2.5, legend=c("Monte Carlo", "Gaussian First Order"), 
lty=c(1,2), col=c(4,2), lwd=2, bg="white")
}

## Display both uncertainty summaries

\dontrun{
barplot(cbind(contr.GFO$budget$contrib, contr.MC$budget$contrib), 
beside=TRUE, horiz=TRUE, main="Uncertainty contribution by method",
xlab="percent Variance",
names.arg=c(GFO.res$measurand_label, MC.res$measurand_label))
}

##########################
## Example H.1 from GUM ##
##########################

# define the uncertainty budget

u.budget<- uncertaintyBudget(
  x=list(
    name=c("lambda.s", "alpha.s", "theta.bar", "Delta", "delta.alpha",
    "delta.theta", "d.bar", "d.cr",
	"d.cnr"),
    label=c("lambda[s]", "alpha[s]", "bar(theta)", "Delta", "delta[alpha]",
    "delta[theta]", "bar(d)", "d[cr]", "d[cnr]"),
    mean=c(50.000623,11.5e-6,-1e-1, 0, 0, 0, 2.15e-4, 0, 0),
    units=c("mm", "oC^-1","oC","oC", "oC^-1", "oC", "mm", "mm", "mm"),
    u=c(25e-6, 1.2e-6, 0.2, 0.35, 0.58e-6, 0.029, 5.8e-6, 3.9e-6, 6.7e-6),
    distribution=c("t","unif","unif","arcsine","unif","unif","t","t","t"),
    dof=c(18, 1, 1, 1, 50, 2, 24, 5, 8)
  ),
  y=diag(1, 9)
)

# define the measurand
measurand_name<- "lambda"
measurand_label<- "lambda"
measurand_model<- paste("(lambda.s*(1+alpha.s*(theta.bar+Delta+delta.theta))",
"+d.bar+d.cr+d.cnr)/(1+(alpha.s+delta.alpha)*(theta.bar+Delta))", sep="")

# estimate the measurand using the Gaussian First Order method (GUM)

u.GFO<- uncertainty(
	x=u.budget, 
	y=list(measurand_name=measurand_name, 
		measurand_label=measurand_label, 
		measurand_model=measurand_model,
		alpha=0.01,
		method="GFO"
	)
)

u.GFO
# same result as reported in Table H.1

# estimate the measurand using the Gaussian Second Order method

u.GSO<- uncertainty(
	x=u.budget, 
	y=list(measurand_name=measurand_name, 
		measurand_label=measurand_label, 
		measurand_model=measurand_model,
		alpha=0.01,
		method="GSO"
	)
) 

u.GSO
# same results as reported in section H.1.6, U(99) = 93 nm,
# the difference is due to rounding error.
# u = 34 nm, but dof are updated to 21 instead of keeping 16.

# estimate the measurand using the Monte Carlo method (GUM supplement 1)

u.MC<- uncertainty(
	x=u.budget, 
	y=list(measurand_name=measurand_name, 
		measurand_label=measurand_label, 
		measurand_model=measurand_model,
		alpha=0.01,
		method="MC", B=1e6
	)
) 

u.MC
# this result is not reported in the GUM

# estimate the measurand using the Kragten method

u.Kragten<- uncertainty(
	x=u.budget, 
	y=list(measurand_name=measurand_name, 
		measurand_label=measurand_label, 
		measurand_model=measurand_model,
		alpha=0.01,
		method="Kragten"
	)
) 

u.Kragten
# same as GFO results

}
