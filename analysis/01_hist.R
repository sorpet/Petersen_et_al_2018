library(tidyverse)
library(cowplot)
library(here)
theme_set(theme_classic(30))

flow_cyt <- read_csv(here('data','flowcytometry.csv'))

cpalette = c("#f7fcfd","#e5f5f9","#ccece6","#99d8c9","#66c2a4","#41ae76","#238b45","#006d2c","#00441b","#00441b")

splot <- function(x) {
    plot <- flow_cyt %>% filter(sample_id == x) %>%
    ggplot(aes(gfp_norm)) +
    #scale_y_continuous(limits = c(0,0.1), breaks = c(0, 0.05, 0.1)) +
    scale_x_continuous(limits = c(10,300), breaks = c(50, 100, 150, 200,250,300)) +
    geom_density(aes(y=..scaled.., fill=factor(sample_id))) +
    scale_fill_manual(values = c(cpalette[[x]])) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks = element_line(size = 1.5),
          axis.line = element_line(size = 1.5),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          legend.position="none")
  return(plot)
}

splot1 <- function(x) {
  plot <- flow_cyt %>% filter(sample_id == x) %>%
    ggplot(aes(gfp_norm)) +
    #scale_y_continuous(limits = c(0,0.1), breaks = c(0, 0.05, 0.1)) +
    scale_x_continuous(limits = c(10,300), breaks = c(50, 100, 150, 200,250,300)) +
    geom_density(aes(y=..scaled.., fill=factor(sample_id))) +
    scale_fill_manual(values = c(cpalette[[x]])) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          #axis.ticks.x=element_blank())
          axis.ticks = element_line(size = 1.5),
          axis.line = element_line(size = 1.5),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          legend.position="none")
  return(plot)
}
plot_grid(splot(1),splot(2),splot(3),splot(4),splot(5),splot(6),splot(7),splot1(8), ncol=1)
ggsave(here('figures','Figure 1,3A_hist.pdf'))
