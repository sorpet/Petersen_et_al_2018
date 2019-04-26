library(tidyverse)
library(broom)
library(cowplot)
library(here)
theme_set(theme_classic(30))
# data for paper --------------------------------------------------------


spect <- read_csv(here('data','spect.csv'),guess_max = 5000) %>%
  filter(meas_date=="2017-12-07" | meas_date=="2017-12-13" | meas_date=="2018-01-12" | meas_date=="2018-01-16") %>%
  mutate(tis = tis_no)

# Wrangle data ------------------------------------------------------------
spect_wt <- spect %>%
  filter(meas_date == "2017-12-13" & description=="CEN.PK113-7D") %>%
  select(main_gfp_median, main_rfp_median) %>%
  map(median, na.rm = TRUE)

spect <- spect %>%
  filter(!is.na(reporter), tis != 9, tis != 10) %>%
  select(meas_date,promoter, reporter, tis,chassis,carbon_source, main_od, main_gfp_median, main_rfp_median) %>%
  mutate(main_fluo_median = ifelse(reporter == "yegfp" | reporter == "ymukg1", main_gfp_median-spect_wt[[1]],
                                   ifelse(reporter == "mkate2", main_rfp_median-spect_wt[[2]],
                                          ifelse(reporter == "mgfp" | reporter == "zsgreen", main_gfp_median,NA)))) %>%
  group_by(meas_date, reporter) %>%
  mutate(main_fluo_median_norm = main_fluo_median/mean(main_fluo_median, na.rm = TRUE)*100) %>%
  ungroup() %>%
  select(-main_od, -main_gfp_median, -main_rfp_median, -main_fluo_median)
# Data types --------------------------------------------------------------
# Factor levels -----------------------------------------------------------
spect$promoter <- factor(spect$promoter, levels = c("tef1","adh2","ef1a"))
spect$reporter <- factor(spect$reporter, levels = c("ymukg1","mkate2","yegfp","mgfp","zsgreen"))
spect$chassis <- factor(spect$chassis, levels = c("yeast","cho"))
spect$carbon_source <- factor(spect$carbon_source, levels = c("glucose","ethanol"))

predictions <- read_csv(here('data','predictions.csv'))
predictions$promoter <- factor(predictions$promoter, levels = c("tef1","adh2","ef1a"))
predictions$reporter <- factor(predictions$reporter, levels = c("ymukg1","mkate2","yegfp","mgfp","zsgreen"))

df <- spect %>%
  filter((meas_date=="2017-12-13" | (meas_date=="2018-01-12" & promoter=="adh2"))) %>%
  group_by(promoter,reporter, tis) %>%
  summarise(median = median(main_fluo_median_norm, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(predictions, by = c("promoter","reporter", "tis"))

df %>%
  #group_by(promoter,reporter) %>%
  do(glance(lm(median ~ protein_abundance, data = .)))

plot2 <- function(data,xx,yy,zz,ww) {
  ggplot(data, mapping = aes_string(x=xx, y=yy)) +
    geom_point(size = 3) +
    stat_smooth(method=lm, se=FALSE, color="black", linetype = "dashed") +
    scale_x_continuous(limits = c(0,200), breaks = c(50, 100, 150, 200)) +
    scale_y_continuous(limits = c(35,105), breaks = c(40, 60, 80, 100)) +
    labs(subtitle = zz, x = "measured", y="predicted") +
    annotate("text", label=ww, parse=TRUE, size=10, x=100, y=102)
    #theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank())
}
plot21 <- function(data,xx,yy,zz,ww) {
  ggplot(data, mapping = aes_string(x=xx, y=yy)) +
    geom_point(size = 3) +
    
    stat_smooth(method=lm, se=FALSE, color="black", linetype = "dashed") +
    scale_x_continuous(limits = c(0,200), breaks = c(50, 100, 150, 200)) +
    scale_y_continuous(limits = c(1.5,5.1), breaks = c(2, 3, 4,5)) +
    labs(subtitle = zz, x = "measured", y="predicted") +
    annotate("text", label=ww, parse=TRUE, size=10, x=100, y=5)
    #theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank())
}

p16 <- df %>% filter(promoter == "tef1", reporter == "ymukg1") %>% plot2("median","efficiency","tef1 ymukg1","paste(R ^ 2, \" = 0.68\")")
p17 <- df %>% filter(promoter == "tef1", reporter == "yegfp") %>% plot2("median","efficiency","tef1 yegfp ","paste(R ^ 2, \" = 0.57\")")
p18 <- df %>% filter(promoter == "tef1", reporter == "mkate2") %>% plot2("median","efficiency","tef1 mkate2","paste(R ^ 2, \" = 0.86\")")
p19 <- df %>% filter(promoter == "adh2", reporter == "yegfp") %>% plot2("median","efficiency","adh2 yegfp ","paste(R ^ 2, \" = 0.54\")")
plot_grid(p16,p17,p18,p19, ncol = 2)
ggsave(here('figures','Figure S6 Noderer.pdf'), width = 10, height = 10, units = "in")

p20 <- df %>% filter(promoter == "tef1", reporter == "ymukg1") %>% plot21("median","protein_abundance","tef1 ymukg1","paste(R ^ 2, \" = 0.47\")")
p21 <- df %>% filter(promoter == "tef1", reporter == "yegfp") %>% plot21("median","protein_abundance","tef1 yegfp ","paste(R ^ 2, \" = 0.44\")")
p22 <- df %>% filter(promoter == "tef1", reporter == "mkate2") %>% plot21("median","protein_abundance","tef1 mkate2","paste(R ^ 2, \" = 0.77\")")
p23 <- df %>% filter(promoter == "adh2", reporter == "yegfp") %>% plot21("median","protein_abundance","adh2 yegfp ","paste(R ^ 2, \" = 0.44\")")
plot_grid(p20,p21,p22,p23, ncol = 2)
ggsave(here('figures','Figure S7 Decoene.pdf'), width = 10, height = 10, units = "in")

theme_set(theme_classic(20))

df %>% 
  #lm(MFE ~ median, data=.) %>% summary()
  ggplot(mapping = aes_string(x="MFE", y="median")) +
  geom_point(size = 3) +
  #labs(subtitle = "MFE") + #title = "A", 
  stat_smooth(method=lm, se=FALSE, color="black", linetype = "dashed") +
  scale_y_continuous(limits = c(0,200), breaks = c(50, 100, 150, 200)) +
  #scale_x_continuous(limits = c(-16, -8), breaks = c(-8, -10,-12,-14, -16)) +
  labs(x = "MFE, kcal/mol", y="Normalized mean fluorescence") +
  annotate("text", label="paste(R ^ 2, \" = 3.0e-3\")", parse=TRUE, size=10, y=182, x=-16, vjust = 0, hjust = 0) +
  annotate("text", label="paste('p', \" = 0.77\")", parse=TRUE, size=10, y=152, x=-16, vjust = 0, hjust = 0)
ggsave(here('figures','Figure S8 secondary structure.pdf'), width = 5, height = 5, units = "in")

