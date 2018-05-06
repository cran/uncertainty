.internal_build_envir <-
function(name, value) {
  n <- length(name)
  res <- paste(name[1], "=", value[1], sep = "")
  if (n > 1) {
    for (i in 2:n) {
      res <- paste(res, ", ", paste(name[i], "=", value[i], sep = ""), sep = "")
    }
  }
  res <- paste("list(", res, ")", sep = "")
  return(res)
}
.internal_dof_eff <-
function(c_u, dof, cor = diag(1, length(c_u))) {
	ct <- sqrt( t(c_u) %*% cor %*% c_u )
	res <- 0
	rr_1 <- sum(c_u ^ 4 / dof)
	rr_2 <- 0
	rr_3 <- 0
	n <- length(c_u)
	# this implements the Castrup.2010 algorithm
	for (i in 1:(n - 1)) {
		for (j in (i + 1):n) {
			rr_2 <- rr_2 + cor[i, j] ^ 2 * (c_u[i] ^ 2 / dof[i]) *
			  (c_u[j] ^ 2 / dof[j]) * (dof[i] + dof[j] + 1 / 2)
			rr_3 <- rr_3 + cor[i, j] * (c_u[i] ^ 3 / dof[i] * c_u[j] +
			                              c_u[i] * c_u[j] ^ 3 / dof[j])
		}
	}
	res <- ct ^ 4 / (rr_1 + rr_2 + 2 * rr_3)
	return (res)
}


.internal_gfo <-
function(name, label = name, model, ub, alpha = 0.05) {
	envir_str <- .internal_build_envir(ub$name, ub$mean)
	n <- length(ub$name)
	
	if (any(is.na(ub$cor))) {
	  ub$cor <- diag(1, n)
	}

	data_cov <- matrix(ub$u, n, n) * ub$cor * t(matrix(ub$u, n, n))

	# get the partial derivatives of the expression wrt each of the variables
	# compute the contribution to the uncertainty by evaluating the list of
	# variables as parameter of the environment.
	# this table contains the variables, the gradient of the function wrt the
	# variable and the uncertainty contribution for that variable.
	res <- data.frame(var = rep(NA, n), partial = rep(NA, n),
	                  coefficient = rep(NA, n), contribution = rep(NA, n))

	res$var <- ub$name
	for (i in 1:n) {
		res$partial[i] <- as.expression(D(parse(text = model), ub$name[i]))
		res$coefficient[i] <- eval(res$partial[i],
		                           envir = eval(parse(text = envir_str)))
	}

	res$contribution <- res$coefficient * ub$u

	# correlation contributions
	JJ <- res$coefficient %*% t(res$coefficient)

	GUM_mu <- eval(parse(text = model), envir = eval(parse(text = envir_str)))
	GUM_u <- sqrt(sum(JJ * data_cov))

	cor_contribution <- GUM_u ^ 2 - sum(res$contribution ^ 2)

	GUM_dof <- .internal_dof_eff(c_u = res$contribution, dof = ub$dof,
	                             cor = ub$cor)

#	print(GUM_dof)
	
	GUM_U <- qt(1 - alpha / 2, GUM_dof) * GUM_u

	su <- NA
	if (GUM_dof > 2) {
		su <- GUM_u * sqrt(GUM_dof / (GUM_dof - 2))
	}

	result_GUM <- list(method = "GFO",
		call = match.call(),
		measurand_name = name,
		measurand_model = model,
		measurand_label = label,
		mean = GUM_mu,
		sd = GUM_u,
		u = GUM_u,
		alpha = alpha,
		dof = GUM_dof,
		U = GUM_U,
		lcl = GUM_mu - GUM_U,
		ucl = GUM_mu + GUM_U,

		variables = list(name = ub$name, label = ub$label, mean = ub$mean, u = ub$u,
		                 dof = ub$dof, coeff = res$coefficient, cor = ub$cor),
		contribution = res$contribution,
		cor_contribution = cor_contribution,

		# additional ouput
		partial = res$partial,
		coeff = res$coefficient,
		su = su)

	return(result_GUM)
}
.internal_gso <-
function(name, label = name, model, ub, alpha = 0.05) {
	envir_str <- .internal_build_envir(ub$name, ub$mean)
	n <- length(ub$name)

	if (any(is.na(ub$cor))) ub$cor <- diag(1, n)

	data_cov <- matrix(ub$u, n, n) * ub$cor * t(matrix(ub$u, n, n))

	# get the partial derivatives of the expression wrt each of the variables.
	# compute the contribution to the uncertainty by evaluating the list of
	# variables as parameter of the environment.
	# this table contains the variables, the gradient of the function wrt the
	# variable and the uncertainty contribution for that variable.
	res <- data.frame(var = rep(NA, n), partial = rep(NA, n),
	                  coefficient = rep(NA, n), contribution = rep(NA, n))

	res$var <- ub$name
	for (i in 1:n) {
		res$partial[i] <- as.expression(D(parse(text = model), ub$name[i]))
		res$coefficient[i] <- eval(res$partial[i],
		                           envir = eval(parse(text = envir_str)))
	}
	res$contribution <- res$coefficient * ub$u

	JJ <- res$coefficient %*%t (res$coefficient)
	GFO_u <- sqrt(sum(JJ * data_cov))

	J <- rep(NA, n)
	H <- matrix(NA, n, n)
	for (ii in 1:n) {
		J[ii] <- eval(res$partial[ii], envir = eval(parse(text = envir_str)))
		for (jj in 1:n) {
			H[ii, jj] <- eval(D(res$partial[ii], ub$name[jj]),
			                 envir = eval(parse(text = envir_str)))
		}
	}

	# second order expected value
	mu_so <- eval(parse(text = model), envir = eval(parse(text = envir_str))) +
	  sum(diag(H %*% data_cov)) / 2
	# second order variance
	u_so <- sqrt(GFO_u ^ 2 + sum(diag(H %*% data_cov %*% H %*% data_cov)) / 2)

	dof_so <- .internal_dof_eff(c_u = res$contribution, dof = ub$dof,
	                            cor = ub$cor)

	cor_contribution <- u_so ^ 2 - sum(res$contribution ^ 2)

	U_so <- qt(1 - alpha / 2, dof_so) * u_so

	res <- list(method = "GSO",
		call = match.call(),
		measurand_name = name,
		measurand_model = model,
		measurand_label = label,
		mean = mu_so,
		sd = u_so,
		u = u_so,
		alpha = alpha,
		dof = dof_so,
		U = U_so,
		lcl = mu_so - U_so,
		ucl = mu_so + U_so,

		variables = list(name = ub$name, label = ub$label, mean = ub$mean, u = ub$u,
		                 dof = ub$dof, coeff = res$coefficient, cor = ub$cor),
		contribution = res$contribution,
		cor_contribution = cor_contribution,

		# additional ouput
		partial = res$partial,
		coeff = res$coefficient,
		su = u_so * sqrt(dof_so / (dof_so - 2)))

	return(res)
}
.internal_kragten <-
function(name, label = name, model, ub, alpha = 0.05) {
	envir_str <- .internal_build_envir(ub$name, ub$mean)
	n <- length(ub$name)

	f <- rep(NA, n + 1)

	f[n + 1] <- eval(parse(text = model), envir = eval(parse(text = envir_str)))

	for (i in 1:n) {
		# build a new point at (x[1],...,x[i]+u[i],...,x[n]) to evaluate the model
		new_mean <- ub$mean
		new_mean[i] <- new_mean[i] + ub$u[i]
		new_envir_str <- .internal_build_envir(ub$name, new_mean)
		f[i] <- eval(parse(text = model), envir = eval(parse(text = new_envir_str)))
	}

	Kragten_mean <- f[n + 1]
	Kragten_contribution <- (f[ - (n + 1)] - Kragten_mean)
	Kragten_u <- sqrt(sum(Kragten_contribution %*% ub$cor %*%
	                        Kragten_contribution))
	Kragten_cor_contribution <- Kragten_u ^ 2 - sum(Kragten_contribution ^ 2)
	Kragten_dof <- .internal_dof_eff(c_u = Kragten_contribution, dof = ub$dof,
	                                 cor = ub$cor)
	Kragten_U <- qt(1 - alpha / 2, Kragten_dof) * Kragten_u
	Kragten_coeff <- Kragten_contribution
	Kragten_coeff[ub$u > 0] <- Kragten_contribution[ub$u > 0] / ub$u[ub$u > 0]
	Kragten_coeff[ub$u == 0] <- NA

	result_Kragten <- list(method = "Kragten",
		call = match.call(),
		measurand_name = name,
		measurand_model = model,
		measurand_label = label,
		mean = Kragten_mean,
		sd = Kragten_u,
		u = Kragten_u,
		alpha = alpha,
		dof = Kragten_dof,
		U = Kragten_U,
		lcl = Kragten_mean - Kragten_U,
		ucl = Kragten_mean + Kragten_U,

		variables = list(name = ub$name, label = ub$label, mean = ub$mean, u = ub$u,
		                 dof = ub$dof, coeff = Kragten_coeff, cor = ub$cor),
		contribution = Kragten_contribution,
		cor_contribution = Kragten_cor_contribution,

		su = Kragten_u * sqrt(Kragten_dof / (Kragten_dof - 2))

		# no additional output
		)

	return(result_Kragten)
}
.internal_mc <-
function(name, label = name, model, ub, alpha = 0.05, B = 2000) {

	n <- length(ub$name)

	rvalues <- matrix(NA, B, n)

	var_env <- new.env()

	if (sum(abs(ub$cor - diag(1,n))) < sqrt(.Machine$double.eps)) {
	  # variables are uncorrelated
		# get the simulated variables
		for (i in 1:n) {
			if (ub$distribution[i] == "bernoulli") {
				rvalues[, i]<- rbinom(B, size = 1, ub$mean[i])
			} else if (ub$distribution[i] == "beta") {
				a <- ub$mean[i] * (ub$mean[i] * (1 - ub$mean[i]) - ub$u[i] ^ 2) /
				  ub$u[i] ^ 2
				b <- (1 - ub$mean[i]) * (ub$mean[i] * (1 - ub$mean[i]) - ub$u[i] ^ 2) /
				  ub$u[i] ^ 2
				rvalues[,i] <- rbeta(B, a, b)
			} else if (ub$distribution[i] == "binomial") {
				rvalues[, i] <- rbinom(B, size = ub$mean[i], ub$u[i])
			} else if (ub$distribution[i] == "cauchy") {
				rvalues[, i] <- rcauchy(B, location = ub$mean[i], scale = ub$u[i])
			} else if (ub$distribution[i] == "chisq") {
				rvalues[, i] <- rchisq(B, ub$dof[i])
			} else if (ub$distribution[i] == "exp") {
				rvalues[, i] <- rexp(B, rate = 1 / ub$mean[i])
			} else if (ub$distribution[i] == "f") {
				# this uses two dof as parameters
				rvalues[, i] <- rf(B, ub$mean[i], ub$u[i])
			} else if (ub$distribution[i] == "gamma") {
				rvalues[, i] <- rgamma(B, shape = (ub$mean[i] / ub$u[i]) ^ 2,
				                       scale = ub$u[i] ^ 2 / ub$mean[i])
			} else if (ub$distribution[i] == "lognormal") {
				rvalues[, i] <- exp(rnorm(B, log(ub$mean[i]) - 1 / 2 *
				                            log(1 + ub$u[i] ^ 2 / ub$mean[i] ^ 2),
				                          sqrt(log(1 + ub$u[i] ^ 2 / ub$mean[i] ^ 2))))
			} else if (ub$distribution[i] == "poisson") {
				rvalues[, i] <- rpois(B, ub$mean[i])
			} else if (ub$distribution[i] == "normal") {
				rvalues[, i] <- rnorm(B, ub$mean[i], ub$u[i])
			} else if (ub$distribution[i] == "unif") {
				rvalues[, i] <- runif(B, ub$mean[i] - sqrt(3) * ub$u[i],
				                      ub$mean[i] + sqrt(3) * ub$u[i])
			} else if (ub$distribution[i] == "t") {
				rvalues[, i] <- ub$mean[i] + rt(B, ub$dof[i]) * ub$u[i] /
				  sqrt(ub$dof[i] / (ub$dof[i] - 2))
			} else if (ub$distribution[i] == "triangular") {
				rvalues[, i] <- rtriangle(B, ub$mean[i] - sqrt(6) * ub$u[i],
				                          ub$mean[i] + sqrt(6) * ub$u[i], ub$mean[i])
			} else if (ub$distribution[i] == "weibull") {
				rvalues[, i] <- rweibull(B, shape = ub$mean[i], scale = ub$u[i])
			} else if (ub$distribution[i] == "arcsine") {
				rvalues[, i] <- ub$mean[i] + ub$u[i] * sin(2 * pi * runif(B))
			} else {
				rvalues[, i] <- rnorm(B, ub$mean[i], ub$u[i])
			}
			eval(parse(text = paste(ub$name[i]," <- rvalues[,", i, "]", sep = "")),
			     envir = var_env)
		}
	}  else {
	  # variables are correlated, generate approximate random variables
		z <- rmvnorm(B, mean = rep(0, n), sigma = ub$cor)
		z <- pnorm(z)

		for (i in 1:n) {
			if (ub$distribution[i] == "bernoulli") {
				rvalues[, i] <- qbinom(z[, i], size = 1, ub$mean[i])
			} else if (ub$distribution[i] == "beta") {
				a <- ub$mean[i] * (ub$mean[i] * (1 - ub$mean[i]) - ub$u[i] ^ 2) /
				  ub$u[i] ^ 2
				b <- (1 - ub$mean[i]) * (ub$mean[i] * (1 - ub$mean[i]) - ub$u[i] ^ 2) /
				  ub$u[i] ^ 2
				rvalues[, i] <- qbeta(z[, i], a, b)
			} else if (ub$distribution[i] == "binomial") {
				rvalues[, i] <- qbinom(z[, i], size = ub$mean[i], ub$u[i])
			} else if (ub$distribution[i] == "cauchy") {
				rvalues[, i] <- qcauchy(z[, i], location = ub$mean[i], scale = ub$u[i])
			} else if (ub$distribution[i] == "chisq") {
				rvalues[, i] <- qchisq(z[, i], ub$dof[i])
			} else if (ub$distribution[i] == "exp") {
				rvalues[, i] <- qexp(z[, i], rate = 1 / ub$mean[i])
			} else if (ub$distribution[i] == "f") {
				rvalues[, i] <- qf(z[, i], ub$mean[i], ub$u[i])
			} else if (ub$distribution[i] == "gamma") {
				rvalues[, i] <- qgamma(z[, i], shape = (ub$mean[i] / ub$u[i]) ^ 2,
				                       scale = ub$u[i] ^ 2 / ub$mean[i])
			} else if (ub$distribution[i] == "lognormal") {
				rvalues[, i] <- exp(qnorm(z[, i], log(ub$mean[i]) - 1 / 2 *
				                            log(1 + ub$u[i] ^ 2 / ub$mean[i] ^ 2),
				                          sqrt(log(1 + ub$u[i] ^ 2 / ub$mean[i] ^ 2))))
			} else if (ub$distribution[i] == "poisson") {
				rvalues[, i] <- qpois(z[, i], ub$mean[i])
			} else if (ub$distribution[i] == "normal") {
				rvalues[, i] <- qnorm(z[, i], ub$mean[i], ub$u[i])
			} else if (ub$distribution[i] == "unif") {
				rvalues[, i] <- qunif(z[, i], ub$mean[i] - sqrt(3) * ub$u[i],
				                      ub$mean[i] + sqrt(3) * ub$u[i])
			} else if (ub$distribution[i] == "t") {
				rvalues[, i] <- ub$mean[i]+qt(z[, i], ub$dof[i]) * ub$u[i] /
				  sqrt(ub$dof[i] / (ub$dof[i] - 2))
			} else if (ub$distribution[i] == "triangular") {
				rvalues[, i] <- qtriangle(z[, i], ub$mean[i] - sqrt(6) * ub$u[i],
				                          ub$mean[i] + sqrt(6) * ub$u[i], ub$mean[i])
			} else if (ub$distribution[i] == "weibull") {
				rvalues[, i] <- qweibull(z[, i], shape = ub$mean[i], scale = ub$u[i])
			} else if (ub$distribution[i] == "arcsine") {
				rvalues[, i] <- ub$mean[i] + ub$u[i] * sin(2 * pi * z[,i])
			} else {
				rvalues[, i] <- qnorm(z[, i], ub$mean[i], ub$u[i])
			}
			eval(parse(text = paste(ub$name[i], " <- rvalues[,", i, "]", sep = "")),
			     envir = var_env)
		}
	}

	f_mc <- eval(parse(text = model), envir = var_env)

	# infinite values are handled as NA
	f_mc[f_mc == Inf] <- NA
	f_mc[f_mc == -Inf] <- NA

	res_mean <- mean(f_mc, na.rm = TRUE)
	res_u <- sd(f_mc, na.rm = TRUE)

	######################################
	## uncertainty contribution estimation
	######################################
	res_contribution <- rep(1, n)
	# backup the simulated inputs with extension ".rv"
	for (i in 1:n) {
		eval(parse(text = paste(ub$name[i], ".rv", "<- ", ub$name[i], sep = "")),
		     envir = var_env)
	}

	# for each variable compute its contribution
	for (i in 1:n) {
		# restore the scalar values to the environment variables
		for (j in 1:n) {
			eval(parse(text = paste(ub$name[j], "<- ", ub$mean[j])), envir = var_env)
		}
		# now replace the ith variable with the simulated values into the environment
		eval(parse(text = paste(ub$name[i], "<- ", ub$name[i], ".rv", sep = "")),
		     envir = var_env)
		# compute the function while changing the ith variable
		f_i <- eval(parse(text = model), envir = var_env)
		# compute the sd of the values as the contribution for that variable
		res_contribution[i] <- sd(f_i - res_mean, na.rm = TRUE)
		if (is.na(res_contribution[i])) res_contribution[i] <- 0
	}

	res_cor_contribution <- res_u^2 - sum(res_contribution^2)
	res_coeff <- res_contribution
	res_coeff[ub$u > 0] <- res_contribution[ub$u > 0] / ub$u[ub$u > 0]
	res_coeff[ub$u == 0] <- NA

	## estimate the dof for the estimated uncertainty
	## update this code with the new criteria
	bb <- floor(sqrt(B))
	vi <- rep(NA, bb)
	for (j in 1:bb) {
		vi[j] <- var(f_mc[ ( (j - 1) * bb + 1):(j * bb)], na.rm = TRUE)
	}
	res_dof_ws <- 2 * res_u^4 / var(vi)
	res_U <- (f_mc[order(f_mc)[round(B * (1 - alpha / 2))]] -
	            f_mc[order(f_mc)[round(B * (alpha / 2))]]) / 2

	result_MC <- list(method = "MC",
		call = match.call(),
		measurand_name = name,
		measurand_model = model,
		measurand_label = label,
		mean = res_mean,
		sd = res_u,
		u = res_u,
		alpha = alpha,
		dof = res_dof_ws,
		U = res_U,
		lcl = f_mc[order(f_mc)[round(B * (alpha / 2))]],
		ucl = f_mc[order(f_mc)[round(B * (1 - alpha / 2))]],

		variables = list(name = ub$name, label = ub$label, mean = ub$mean, u = ub$u,
		                 dof = ub$dof, coeff = res_coeff, cor = ub$cor),
		contribution = res_contribution,
		cor_contribution = res_cor_contribution,

		# additional ouput
		ncycles = B,
		na_count = sum(is.na(f_mc)),
		simulated_values = f_mc,
		simulated_data = rvalues)

	return(result_MC)
}
