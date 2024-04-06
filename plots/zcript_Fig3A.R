source("functions.R")
library(tidyverse)
library(hrbrthemes)

plot_competition = function(generate_pdf = FALSE) {
	# Data computed in the notebook homogeneous_approximation.ipynb
	data = read.csv("data/ER_botnet_size.csv")

	beta_W = seq(0, 10, by = 0.01)
	analytical = data.frame(beta_W = beta_W, beta_B = computeBotnetThreshold(beta_W, 1, 1, 1))
	analytical = analytical[analytical$beta_B > 0, ]

	p = ggplot() +
		geom_tile(data = data, aes(x = beta_w, y = beta_b, fill = botnet, color = botnet)) +
		geom_line(data = analytical, aes(x = beta_W, y = beta_B), color = "white", linewidth = 1, linetype = 2) +
		theme_ipsum(base_family = "Helvetica",
				axis_title_size = 18, 
				axis_text_size = 18,
				axis_title_just = "rt") +
		scale_color_gradientn(name = "Botnet size", colors = ipsum_pal()(3), limits = c(0, 1)) +
		scale_fill_gradientn(name = "Botnet size", colors = ipsum_pal()(3), limits = c(0, 1)) +
		guides(color = guide_colourbar(ticks.colour = "white",
						ticks.linewidth = 1,
						frame.colour = "white",
						frame.linewidth = 1,
						barheight = unit(.6,"npc"),
						barwidth = unit(.045,"npc")))+
		theme(legend.text = element_text(size = 14),
			legend.key.width = unit(1.5, 'cm'),
			legend.title = element_text(size = 18),
			legend.position = "right",
			plot.title = element_text(size = 24, hjust = 0.5, margin=margin(-20,0,0,0))) +
		ggtitle("ER network") +
		xlab(expression(beta[W])) +
		ylab(expression(beta[B])) +
		scale_x_continuous(limits = c(-0.1, 10.1)) +
		scale_y_continuous(limits = c(-0.1, 10.1))

	if(generate_pdf) {
		cairo_pdf("Fig3A.pdf", width = 10, height = 6)
		plot(p)
		dev.off()
	}

	return(p)
}

plot_competition(TRUE)
