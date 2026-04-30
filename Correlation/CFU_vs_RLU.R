# Correlation and regression analysis
library(ggplot2)
library(RColorBrewer)

# Read the data "CFU_vs_RLU.txt", which is available at the GitHub page.
dat <- read.delim(file.choose(), header=T)

#calculation of regression and correlation coefficients
dat_lm <- summary(lm(CFU ~ RLU, data = dat))
r <- round(cor(dat[, "CFU"], dat[, "RLU"]), 2)
coef = round(dat_lm$coefficients["RLU", "Estimate"], 2)
intercept = round(dat_lm$coefficients["(Intercept)", "Estimate"], 2)
regression = paste("y =", coef, "x +", intercept, ",", "R =", r)

##genotype x treatment
Fig1D <- ggplot(dat, aes(x = RLU, y = CFU, color = Strain)) +
      scale_colour_brewer(palette = "Set1") +
      geom_point(alpha = 0.5, size = 3) + theme_classic() + ylim(c(3, 7.5)) +
      scale_x_continuous(limits = c(1.5, 6), breaks = c(1.5, 2.5, 3.5, 4.5, 5.5)) +
      xlab(expression(paste(log(RLU/cm^{2})))) + ylab(expression(paste(log(CFU/cm^{2})))) +
      stat_smooth(method = "lm", se = FALSE, colour = "black", size = 0.5) +
      theme(legend.position="bottom") + theme(legend.title=element_blank()) +
      labs(title = regression) + theme(plot.title = element_text(size = 12))

ggsave(file = "Fig1D.png", plot = Fig1D,  dpi = 300, width = 5, height = 5)
