library(cowplot)

source("zcript_Fig3A.R")
pA = plot_competition()

source("zcript_Fig3B.R")
pB = plot_competition()

cairo_pdf("Fig3.pdf", width = 20, height = 6)
plot_grid(pA, pB,
		axis = "b",
		align = "h", nrow = 1, ncol = 2,
		labels = c("a", "b"), label_size = 26,
		label_fontfamility = "Helvetica")
dev.off()
