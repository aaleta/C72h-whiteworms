source("functions.R")
library(tidyverse)
library(hrbrthemes)

data = NULL
beta_W = seq(0, 20, by = 0.01)
for(epsilon in seq(0, 10, by = 0.1)) {
	data = rbind(data, data.frame(beta_W = beta_W, 
					epsilon = epsilon, 
					beta_B = computeBotnetThreshold(beta_W, epsilon, 1, 1)
				)
			)
}

colors = ipsum_pal()(2)

p = ggplot() +
	geom_line(data = data, aes(x = beta_W, y = beta_B, color = epsilon, group = epsilon)) +
	geom_line(data = data, aes(x = beta_W, y = beta_W - 1), color = "black", linetype = 2) +
	annotate("text", x = 14, y = 15, angle = 35, label = "paste(\"\u03F5\", \"\u2192\", infinity)", parse = TRUE, size = 8) +
	theme_ipsum(base_family = "Helvetica",
			axis_title_size = 18, 
			axis_text_size = 18,
			axis_title_just = "rt") +
	scale_color_gradient(name = expression("\u03F5"), low = colors[1], high = colors[2]) +
	guides(color = guide_colourbar(ticks.colour = "white",
					ticks.linewidth = 1,
					frame.colour = "white",
					frame.linewidth = 1,
					barheight = unit(.6,"npc"),
					barwidth = unit(.045,"npc")))+
	theme(legend.text = element_text(size = 14),
		legend.key.width = unit(1.5, 'cm'),
		legend.title = element_text(size = 18),
		legend.position = "right") +
	xlab(expression(beta[W])) +
	ylab(expression(beta[B])) +
	scale_y_continuous(limits = c(0, 20))

cairo_pdf("FigS1.pdf", width = 10, height = 6)
plot(p)
dev.off()
