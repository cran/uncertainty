uncertaintyBudget.default <-
function(x, y, ...) {
	colnames(y) <- rownames(y) <- x$name
	res <- list(name = x$name, mean = x$mean, u = x$u, dof = x$dof,
	            label = x$label, distribution = x$distribution, cor = y)
	class(res) <- "uncertaintyBudget"
	return(res)
}
