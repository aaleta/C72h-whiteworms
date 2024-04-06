source("functions.R")
library(tidyverse)
library(hrbrthemes)

plot_who_protects = function(generate_pdf = FALSE) {
	# Data computed in the notebook homogeneous_approximation.ipynb
	data = read.csv("data/homogeneous_protection_self_forced.csv")

	p = ggplot() +
		geom_line(data = data, aes(x = epsilon/gamma, y = self, color = format(beta, digits = 2), linetype = "self-protected"), linewidth = 1) +
		geom_line(data = data, aes(x = epsilon/gamma, y = protected, color = format(beta, digits = 2), linetype = "total protection"), linewidth = 1) +
		theme_ipsum(base_family = "Helvetica",
				axis_title_size = 18, 
				axis_text_size = 18,
				axis_title_just = "rt") +
		scale_color_ipsum() +
		scale_linetype_manual(values = c(2, 1)) +
		guides(color=guide_legend(title=expression(beta[W])),
			linetype=guide_legend(title="")) +
		theme(legend.text = element_text(size = 14),
			legend.key.width = unit(1.5, 'cm'),
			legend.title = element_text(size = 18),
			legend.position = "inside",
			legend.position.inside = c(0.14, 0.7),
			legend.box.background = element_rect(colour = "black"),
			legend.margin = margin(0.1,0.1,0.1,0.1, unit="cm")) +
		xlab(expression(paste("\u03F5", "/", gamma))) +
		ylab("Fraction of protected devices") +
		scale_x_log10(limits = c(0.05, 10))
		#scale_x_continuous(breaks = c(0, 7, 14, 21, 28))

	if(generate_pdf) {
		cairo_pdf("Fig2B.pdf", width = 10, height = 6)
		plot(p)
		dev.off()
	}
	
	return(p)
}
