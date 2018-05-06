summary.uncertainty <-
function(object, ndigits = 3, ...) {
	if (class(object) != "uncertainty") {
	  stop("object must be an uncertainty object.\n")
	}
	if (all(object$method != c("Kragten","GFO","GSO","MC"))) {
	  stop("object was created with an unknown method.\n")
	}
	res <- list(call = object$call, measurand_name = object$measurand_name,
	           measurand_label = object$measurand_label)
	if (object$cor_contribution != 0) {
		res$budget <- list(name = c(object$variables$name, "Correlations"),
			mean = c(object$variables$mean, 0),
			label = c(object$variables$label, "Correlations"),
			u = c(object$variables$u, 0),
			dof = c(object$variables$dof, 0),
			contrib = signif(c(object$contribution ^ 2, object$cor_contribution) /
			                   object$sd ^ 2 * 100, ndigits))
	} else {
		res$budget <- list(name = object$variables$name,
			mean = object$variables$mean,
			label = object$variables$label,
			u = object$variables$u,
			dof = object$variables$dof,
			contrib = signif(object$contribution ^ 2 / object$sd ^ 2 * 100, ndigits))
	}
	class(res)<- "summary.uncertainty"
	res
}
