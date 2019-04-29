library(tidyverse)
library(stringr)
library(gtools)
library(cowplot)
library(here)
theme_set(theme_classic(30))

reorder_by_count <- function(x) {
  factor(x, levels = names(sort(table(x), decreasing = TRUE)))
  }


ngs <- read_csv(here('data','ngs.csv'))

tis_permutations <- as.tibble(permutations(4,6,c("A","T","C","G"),repeats=TRUE)) %>%
  mutate(tis = str_c(V1,V2,V3,V4,V5,V6)) %>%
  .[["tis"]]

ngs$tis <- as.factor(ngs$tis)
ngs$tis <- factor(ngs$tis, levels = tis_permutations)
ngs$sample <- as.factor(ngs$sample)

head(ngs)

plot1 <- function(x) {
  ngs %>%
    filter(sample == x) %>%
    ggplot() +
    geom_bar(aes(reorder_by_count(tis))) +
    scale_y_continuous(breaks = c(0, 1500, 3000)) +
    #scale_y_continuous(limits = c(0,5000), breaks = c(0, 2500, 5000)) +
    scale_x_discrete(drop=FALSE) +
    coord_cartesian(ylim=c(0, 3000)) +
    labs(y = "count") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  }

plot_grid(plot1(3), plot1(2), plot1(1), labels = c("A", "B", "C"), align = "v", ncol = 1)
ggsave(here('figures','Figure S2_counts.pdf'))

# if trusted tis's require >100 reads
# TEF1 4037-3193 = 844, 21 %
# RPL18b 4093 - 1920 = 2073, 51%
# REV1 4037 - 2317 = 1720, 42%

nuc <- ngs %>%
  filter(sample == "1") %>%
  extract(tis, into = c("-6","-5","-4","-3","-2","-1"), '(.)(.)(.)(.)(.)(.)') %>%
  gather(nucleotide_position,nucleotide,"-6","-5","-4","-3","-2","-1")

p1 <- nuc %>%
  ggplot() +
  geom_bar(mapping = aes(x = sample, fill = nucleotide), position = "fill") +
  labs(fill = "nuc.") +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p2 <- p1 + facet_wrap(~nucleotide_position)
plot_grid(p1, p2, labels = c("A", "B"), ncol = 2)
ggsave(here('figures','Figure S3_frequencies.pdf'), width = 10, height = 5, units = "in")
