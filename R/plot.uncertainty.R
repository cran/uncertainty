plot.uncertainty <-
function(x, y = NULL,
	xlab = parse(text = x$measurand_label), main = "",
	ylab = "Probability density", from = x$mean - 4 * x$u, to = x$mean + 4 * x$u,
	lwd = 2, add = FALSE, ...) {
	mean <- x$mean
	sd <- x$u
	dof <- x$dof
	if (x$method == "MC") {
		if (add) {
			lines(density(x$simulated_values, from = from, to = to), ...)
		} else {
			plot(density(x$simulated_values, from = from, to = to), ylab = ylab,
			     lwd = lwd, main = main, xlab = xlab, ...)
		}
	} else {
		curve(dt( (x - mean) / sd, dof) / sd, xlab = xlab, ylab = ylab, main = main,
		      from = from, to = to, lwd = lwd, add = add, ...)
	}
}
