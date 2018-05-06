print.uncertainty <-
function(x, ...) {
	cat(paste(c("measurand:", x$measurand_name)), "\nmodel:", x$measurand_model)

	cat(paste(c("\nmean =", x$mean), sep = ""), "\n")
	cat(paste(c("sd    =", x$sd), sep = ""), "\n")
	cat(paste(c("u    =", x$u), sep = ""), "\n")
	cat(paste(c("dof  =", signif(x$dof, 3)), sep = ""), "\n")
	cat(paste(c("U(", signif(100 * (1 - x$alpha), 3),"%)=", x$U), sep = ""), "\n")
	cat(paste(100 * (1 - x$alpha), "% CI = (", x$lcl, ",", x$ucl, ")", sep = ""),
	    "\n")
}
