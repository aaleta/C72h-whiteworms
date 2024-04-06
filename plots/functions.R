R0 = function(beta, epsilon, gamma = 1, mu = 1) {
	epp = epsilon / gamma
	return((beta / mu) * (epp / (epp + 1)))
}

protectedF = function(x, r0) {
	return(abs(x - 1.0 + exp(-r0 * x)))
}


computeProtected = function(data, gamma = 1, mu = 1) {
	
	betas = seq(0, 10, by = 0.15)
	result = data.frame()
		
	epsilons = unique(data$epsilon)
	for(epsilon in epsilons) {
		for(beta in betas) {
			result = rbind(result,
				data.frame(beta = beta,
						epsilon = epsilon,
						protected = optimize(protectedF, interval = c(0, 1), r0 = R0(beta, epsilon), tol = 0.0001)$minimum
					)
			)
		}
	}

	return(result)
}

computeBotnetThreshold = function(beta_W_values, epsilon, gamma, mu) {
	return(0.5 * (-(epsilon + mu + gamma) + sqrt((epsilon + gamma - mu)**2 + 4 * beta_W_values * epsilon)))
}

computeBotnetTime = function(data, beta_W, gamma = 1, mu = 1) {
	
	result = data.frame()
		
	epsilons = unique(data$forced) * gamma
	epsilons = epsilons[seq(1, length(epsilons), by = 2)]
	for(epsilon in epsilons) {
		result = rbind(result,
			data.frame(forced = epsilon / gamma,
					protected = optimize(protectedF, interval = c(0, 1), r0 = R0(beta_W, epsilon, gamma), tol = 0.0001)$minimum
				)
		)
	}

	return(result)
}
