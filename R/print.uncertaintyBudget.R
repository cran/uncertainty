print.uncertaintyBudget <-
function(x, ...) {
	dat <- data.frame(quantity = x$name, mean = x$mean, u = x$u,
	                  distribution = x$distribution, dof = x$dof)
	cat("Uncertainty budget:\n")
	print(dat)
	cat("Correlation matrix:\n")
	print(x$cor)
}
