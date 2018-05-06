plot.summary.uncertainty <-
function(x, y = NULL, ...) {
	if (x$budget$contrib[length(x$budget$label)] > 0) {
		barplot(x$budget$contrib, horiz = TRUE,
		        names.arg = c(parse(text = x$budget$label)),
			      ylab = "Source", xlab = "% Variance", col = "dark magenta",
			      main = expression("Uncertainty Contribution"),
			      xlim = c(0, 100))
	} else {
		barplot(x$budget$contrib, horiz=TRUE,
		        names.arg = c(parse(text = x$budget$label)),
			      ylab = "Source", xlab = "% Variance", col = "dark magenta",
			      main  ="Uncertainty Contribution",
			      xlim = c(min(c(-100, x$budget$contrib * 1.1)),
			                                     max(c(200, x$budget$contrib * 1.1))))
	}
	lims <- range(x$budget$contrib * 1.1)
	if (lims[1] > 0) lims[1] <- 0
	if (lims[2] < 100) lims[2] <- 100
	signif((lims[2] - lims[1]) / 4, 1)
	abline(v = seq(from = signif(lims[1],1), to = signif(lims[2],1),
	               by = (signif( (lims[2] - lims[1]) / 4,1))), lty = 2)
	abline(v = 0, lty = 2)

	box()
}
