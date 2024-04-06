source("zcript_Fig2A.R")
source("zcript_Fig2B.R")
source("zcript_Fig2C.R")
source("zcript_Fig2D.R")
library(cowplot)

pA = plot_protected()
#pA$layout$clip = "off"

pB = plot_who_protects()
#pB$layout$clip = "off"

pC = plot_competition()
#pC$layout$clip = "off"

pD = plot_time()
#pD$layout$clip = "off"

cairo_pdf("Fig2.pdf", width = 20, height = 12)
plot_grid(pA, pB, pC, pD,
		align = "h", nrow = 2, ncol = 2,
		labels = c("a", "b", "c", "d"), label_size = 26,
		label_fontfamility = "Helvetica")
dev.off()
