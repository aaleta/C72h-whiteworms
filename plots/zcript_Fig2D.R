source("functions.R")
library(tidyverse)
library(hrbrthemes)

plot_time = function(generate_pdf = FALSE) {
	# Data computed in the notebook homogeneous_approximation.ipynb
	data = read.csv("data/homogeneous_botnet_time.csv")
	analytical = computeBotnetTime(data, beta_W = 2, gamma = 1, mu = 1)

	data$time = data$time * data$total_time
	data = data[data$time > 1, ]

	p = ggplot() +
		geom_tile(data = data, aes(x = forced_rate, y = threshold, fill = time, color = time)) +
		geom_line(data = analytical, aes(x = forced, y = 1 - protected), color = "black", linewidth = 1, linetype = 2) +
		theme_ipsum(base_family = "Helvetica",
				axis_title_size = 18, 
				axis_text_size = 18,
				axis_title_just = "rt") +
		scale_color_gradientn(name = "Timesteps above threshold", colors = ipsum_pal()(3), limits = c(1, 100), breaks = c(1, 25, 50, 75, 100)) +
		scale_fill_gradientn(name = "Timesteps above threshold", colors = ipsum_pal()(3), limits = c(1, 100), breaks = c(1, 25, 50, 75, 100)) +
		guides(color = guide_colourbar(ticks.colour = "white",
						ticks.linewidth = 1,
						frame.colour = "white",
						frame.linewidth = 1,
						barheight = unit(.6,"npc"),
						barwidth = unit(.045,"npc"),
						title.position = "right"))+
		theme(legend.text = element_text(size = 14),
			legend.key.width = unit(1.5, 'cm'),
			legend.title = element_text(angle = 90, size = 18, hjust = 0.5),
			legend.position = "right") +
		xlab(expression(paste("\u03F5", "/", gamma))) +
		ylab("Threshold botnet size") +
		scale_x_continuous(limits = c(-0.1, 10.1)) +
		scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, by = 0.25))

	if(generate_pdf) {
		cairo_pdf("Fig2D.pdf", width = 10, height = 6)
		plot(p)
		dev.off()
	}

	return(p)
}
