print.summary.uncertainty <-
function(x, ...) {
	cat("\nUncertainty summary contribution for ", x$measurand_name, "\n")
	amatrix <- matrix(0, length(x$budget$name), 1)
	amatrix[, 1] <- x$budget$contrib
	colnames(amatrix) <- c("%variance")
	rownames(amatrix) <- c(x$budget$name) #, "Correlations")

	my_line <- ""
	for (i in 1:max(nchar(x$budget$name) + 10)) {
		my_line <- paste(my_line, "-", sep = "")
	}

	cat(my_line, "\n")
	print(amatrix)
	cat(my_line, "\n")
}
