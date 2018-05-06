uncertainty.default <- function(x, y, ...) {
## this is the uncertainty estimator constructor
	if (class(x) != "uncertaintyBudget") {
	  stop("x must be an uncertainty budget object")
	}
	if (y$method == "MC") {
		res <- .internal_mc(name = y$measurand_name, label = y$measurand_label,
		                    model = y$measurand_model, ub = x, alpha = y$alpha,
		                    B = y$B)
	} else if (y$method == "GFO") {
		res <- .internal_gfo(name = y$measurand_name, label = y$measurand_label,
		                     model = y$measurand_model, ub = x, alpha = y$alpha)
	} else if (y$method == "GSO") {
		res <- .internal_gso(name = y$measurand_name, label = y$measurand_label,
		                     model = y$measurand_model, ub = x, alpha = y$alpha)
	} else if (y$method == "Kragten") {
		res <- .internal_kragten(name = y$measurand_name, label = y$measurand_label,
		                         model = y$measurand_model, ub = x, alpha = y$alpha)
	} else {
		stop("estimation methods available are (Kragten, GFO, GSO, MC)")
	}
	class(res) <- "uncertainty"
	return( res )
}
