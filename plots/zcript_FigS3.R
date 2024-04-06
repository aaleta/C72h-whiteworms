source("functions.R")
library(tidyverse)
library(hrbrthemes)

plot_protected = function(generate_pdf = FALSE) {
	# Data computed in the notebook homogeneous_approximation.ipynb
	data = read.csv("data/homogeneous_protection_stochastic.csv")
	data = data %>% filter(epsilon %in% c(0.1, 0.15, 0.2, 0.3, 0.6, 1, 2))
	analytical = computeProtected(data)

	p = ggplot() +
		geom_line(data = analytical, aes(x = beta, y = protected, color = format(epsilon, digits = 2)), size = 1) +	
		geom_point(data = data, aes(x = beta, y = protected, color = format(epsilon, digits = 2)), size = 3) +
		geom_errorbar(data = data, aes(x = beta, ymin = low, ymax = high, color = format(epsilon, digits = 2)), width = .2) +
		geom_line(data = analytical, aes(x = beta, y = protected, color = format(epsilon, digits = 2)), size = 1) +
		theme_ipsum(base_family = "Helvetica",
				axis_title_size = 18, 
				axis_text_size = 18,
				axis_title_just = "rt") +
		scale_color_ipsum() +
		guides(color=guide_legend(title=expression(paste("\u03F5", "/", gamma)))) +
		theme(legend.text = element_text(size = 14),
			legend.key.width = unit(1.5, 'cm'),
			legend.title = element_text(size = 18),
			legend.position = c(0.12, 0.705),
			legend.box.background = element_rect(colour = "black")) +
		xlab(expression(beta[W])) +
		ylab("Fraction of protected devices")# +
		#scale_x_continuous(breaks = c(0, 7, 14, 21, 28))

	if(generate_pdf) {
		cairo_pdf("FigS3.pdf", width = 10, height = 6)
		plot(p)
		dev.off()
	}
	
	return(p)
}

plot_protected(TRUE)
