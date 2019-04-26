library(tidyverse)
library(ggthemes)
library(googlesheets)
library(GGally)
library(cowplot)
library(data.table)
library(stringr)
library(here)
theme_set(theme_classic(30))

# data for paper --------------------------------------------------------
spect <- read_csv(here('data','spect.csv'), guess_max = 5000) %>%
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
# EDA ---------------------------------------------------------------------

cpalette <- c("#252525","#ca0020", "#1a9641", "#0571b0")
plot1 <- function(data,arg2,arg3) {
  ggplot(data) +
    geom_point(aes_string(x=arg2, y="main_fluo_median_norm", color=arg3, shape=arg3), size = 3) +
    scale_y_continuous(limits = c(0,200), breaks = c(50, 100, 150, 200)) +
    labs(x = "TIS no.", y = "Normalized fluorescence") + #title = "B",
    scale_colour_manual(values= cpalette) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          legend.position = c(0.39,.90))
}

spect1 <- spect
spect1$tis <- factor(spect1$tis, levels = c("1","2","3","4","5","6","7","8"))

#PROMOTER
p1 <- spect1 %>% 
  filter(meas_date=="2018-01-12") %>%
  plot1("tis","promoter")
  # annotate("text",x=1,y=.1,label="a", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=2,y=.1,label="ab", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=3,y=.1,label=" bc", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=4,y=.1,label=" bc", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=5,y=.1,label="   d", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=6,y=.1,label="   def", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=7,y=.1,label="     ef", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=8,y=.1,label="      f", angle = 90, size = 6, vjust = 0.4, hjust = 0)
  #aov(main_fluo_median_norm ~ promoter * tis, data = .) %>%
  #summary()
  # aov(main_fluo_median_norm ~ tis, data = .) %>%
  # #summary() 
  # TukeyHSD()
  
# REPORTER
p2 <- spect1 %>%
  filter(meas_date=="2017-12-13") %>%
  plot1("tis","reporter")
  # annotate("text",x=1,y=.1,label="a", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=2,y=.1,label="ab", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=3,y=.1,label="  bc", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=4,y=.1,label="    cd", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=5,y=.1,label="      de", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=6,y=.1,label="        ef", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=7,y=.1,label="        ef", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=8,y=.1,label="          f", angle = 90, size = 6, vjust = 0.4, hjust = 0)
  # aov(main_fluo_median_norm ~ tis, data = .) %>%
  # #summary()
  # TukeyHSD()

# CARBON_SOURCE
p3 <- spect1 %>%
  filter(meas_date== "2017-12-13" & reporter == "yegfp" | meas_date=="2018-01-12" & promoter == "tef1") %>%
  plot1("tis","carbon_source")
  # annotate("text",x=1,y=.1,label="a", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=2,y=.1,label="ab", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=3,y=.1,label=" bc", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=4,y=.1,label=" bc", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=5,y=.1,label="   d", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=6,y=.1,label="   def", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=7,y=.1,label="     ef", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=8,y=.1,label="      f", angle = 90, size = 6, vjust = 0.4, hjust = 0)
  # aov(main_fluo_median_norm ~ tis, data = .) %>%
  # #summary()
  # TukeyHSD()
plot_grid(p1, p2, p3, ncol = 3)
ggsave(here('figures','Figure 3CDE.pdf'), width = 20.1, height = 5, units = "in")

# df1 = spect %>% filter(meas_date=="2018-01-12")
# 
# 
# df1 %>%
#   aov(main_fluo_median_norm ~ promoter * tis, data = .) %>%
#   summary()
  
# summary(res.aov)
# TukeyHSD(res.aov)

# df1 = spect %>% filter(meas_date=="2018-01-12" | meas_date=="2017-12-13")
# df1$tis <- factor(df1$tis, levels = c("1","2","3","4","5","6","7","8"))
# res.aov <- aov(main_fluo_median_norm ~ promoter * reporter * tis * carbon_source, data = df1)
# summary(res.aov)
# TukeyHSD(res.aov)

plot22 <- function(data,xx,yy,min1,max1,min2,max2,zz,ww,qq) {
  ggplot(data, mapping = aes_string(x=xx, y=yy)) +
    geom_point(size = 3) +
    geom_errorbar(aes_string(ymin=min1, ymax=max1),width = 0.2) +
    geom_errorbarh(aes_string(xmin=min2, xmax=max2),height = 0.2) +
    labs(subtitle = zz) + #title = "A", 
    stat_smooth(method=lm, se=FALSE, color="black", linetype = "dashed") +
    scale_y_continuous(limits = c(0,200), breaks = c(50, 100, 150, 200)) +
    scale_x_continuous(limits = c(0,200), breaks = c(50, 100, 150, 200)) +
    annotate("text", label=ww, parse=TRUE, size=8, x=5, y=180,hjust = 0) +
    annotate("text", label=qq, parse=TRUE, size=8, x=5, y=140,hjust = 0)
}

p4 <- spect %>%
  filter(meas_date=="2018-01-12" & reporter == "yegfp") %>%
  group_by(promoter, tis) %>%
  summarise(fluo_mean = mean(main_fluo_median_norm, na.rm = TRUE),
            fluo_se  = sd(main_fluo_median_norm, na.rm = TRUE) / sqrt(n())) %>%
  select(tis,promoter,fluo_mean, fluo_se) %>%
  setDT() %>% 
  dcast(tis ~ promoter, value.var = c("fluo_mean", "fluo_se")) %>%
  mutate(ymin_tef1 =fluo_mean_tef1 - fluo_se_tef1,
         ymax_tef1 =fluo_mean_tef1 + fluo_se_tef1,
         ymin_adh2 =fluo_mean_adh2 - fluo_se_adh2,
         ymax_adh2 =fluo_mean_adh2 + fluo_se_adh2) %>%
  plot22("fluo_mean_adh2","fluo_mean_tef1","ymin_tef1","ymax_tef1","ymin_adh2","ymax_adh2","promoter","paste(R ^ 2, \" = 0.90\")","paste(p, \" = 2.8e-4\")")
 
df <- spect %>%
  filter(meas_date=="2017-12-13") %>%
  group_by(reporter, tis) %>%
  summarise(fluo_mean = mean(main_fluo_median_norm, na.rm = TRUE),
            fluo_se  = sd(main_fluo_median_norm, na.rm = TRUE) / sqrt(n())) %>%
  select(tis,reporter,fluo_mean, fluo_se) %>%
  setDT() %>% 
  dcast(tis ~ reporter, value.var = c("fluo_mean", "fluo_se")) %>%
  #ggpairs()
  mutate(ymin_ymukg1 =fluo_mean_ymukg1 - fluo_se_ymukg1,
         ymax_ymukg1 =fluo_mean_ymukg1 + fluo_se_ymukg1,
         ymin_mkate2 =fluo_mean_mkate2 - fluo_se_mkate2,
         ymax_mkate2 =fluo_mean_mkate2 + fluo_se_mkate2,
         ymin_yegfp =fluo_mean_yegfp - fluo_se_yegfp,
         ymax_yegfp =fluo_mean_yegfp + fluo_se_yegfp)
p5 <- plot22(df,"fluo_mean_ymukg1","fluo_mean_yegfp","ymin_yegfp","ymax_yegfp","ymin_ymukg1","ymax_ymukg1","reporter1","paste(R ^ 2, \" = 0.92\")","paste(p, \" = 1.87e-4\")")
p6 <- plot22(df,"fluo_mean_mkate2","fluo_mean_ymukg1","ymin_ymukg1","ymax_ymukg1","ymin_mkate2","ymax_mkate2","reporter2","paste(R ^ 2, \" = 0.86\")","paste(p, \" = 8.26e-4\")")
p7 <- plot22(df,"fluo_mean_mkate2","fluo_mean_yegfp","ymin_yegfp","ymax_yegfp","ymin_mkate2","ymax_mkate2","reporter3","paste(R ^ 2, \" = 0.75\")","paste(p, \" = 5.65e-3\")")
p8 <- spect %>%
  filter((meas_date== "2017-12-13" & reporter == "yegfp") | (meas_date=="2018-01-12" & promoter == "adh2")) %>%
  group_by(carbon_source,tis) %>%
  summarise(fluo_mean = mean(main_fluo_median_norm, na.rm = TRUE),
            fluo_se  = sd(main_fluo_median_norm, na.rm = TRUE) / sqrt(n())) %>%
  select(tis,carbon_source,fluo_mean, fluo_se) %>%
  setDT() %>% 
  dcast(tis ~ carbon_source, value.var = c("fluo_mean", "fluo_se")) %>%
  mutate(ymin_glucose =fluo_mean_glucose - fluo_se_glucose,
         ymax_glucose =fluo_mean_glucose + fluo_se_glucose,
         ymin_ethanol =fluo_mean_ethanol - fluo_se_ethanol,
         ymax_ethanol =fluo_mean_ethanol + fluo_se_ethanol) %>%
  plot22("fluo_mean_ethanol","fluo_mean_glucose","ymin_glucose","ymax_glucose","ymin_ethanol","ymax_ethanol","carbon_source","paste(R ^ 2, \" = 0.98\")","paste(p, \" = 9.97e-7\")")
plot_grid(p4, p8, NULL, p5, p6, p7, ncol = 3)
ggsave(here('figures','Figure 4ABCDE.pdf'), width = 20.1, height = 10, units = "in")

plot2 <- function(data,xx,yy,zz,ww, qq) {
  ggplot(data, mapping = aes_string(x=xx, y=yy)) +
    geom_point(size = 3) +
    labs(subtitle = zz) + #title = "A", 
    stat_smooth(method=lm, se=FALSE, color="black", linetype = "dashed") +
    scale_y_continuous(limits = c(0,200), breaks = c(50, 100, 150, 200)) +
    #annotate("text", label=ww, parse=TRUE, size=10, x=80, y=170) +
    annotate("text", label=ww, parse=TRUE, size=8, x=5, y=180,hjust = 0) +
    annotate("text", label=qq, parse=TRUE, size=8, x=5, y=140,hjust = 0)
}

spect2 <- spect
spect2$tis <- factor(spect2$tis, levels = c("1","5","8","11"))

# CHASSIS
p9 <- spect %>%
  filter((meas_date=="2017-12-07" | meas_date=="2018-01-16" | meas_date=="2017-12-13") & (reporter == "yegfp" | reporter == "mgfp")) %>%
  filter(tis == 1 | tis == 5 | tis == 8| tis == 11) %>%
  plot1("tis","chassis") +
  # # annotate("text",x=1,y=.1,label="a", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=5,y=.1,label="a", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # # annotate("text",x=8,y=.1,label=" bc", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=11,y=.1,label="a", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  scale_x_continuous(limits = c(1,11), breaks = c(1, 5, 8,11))
spect2 %>%
  filter((meas_date=="2017-12-07" | meas_date=="2018-01-16" | meas_date=="2017-12-13") & (reporter == "yegfp" | reporter == "mgfp")) %>%
  filter(tis == 1 | tis == 5 | tis == 8| tis == 11) %>%
  aov(main_fluo_median_norm ~ tis, data = .) %>%
  TukeyHSD()

# REPORTER_CHO
p10 <- spect %>%
  filter((meas_date=="2017-12-07"| meas_date=="2018-01-16") & (tis == 1 | tis == 5 | tis == 8)) %>%
  plot1("tis","reporter") +
  # annotate("text",x=1,y=.1,label="a", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=5,y=.1,label="ab", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  # annotate("text",x=8,y=.1,label="  b", angle = 90, size = 6, vjust = 0.4, hjust = 0) +
  scale_x_continuous(limits = c(1,8), breaks = c(1, 5, 8))
spect2 %>%
  filter((meas_date=="2017-12-07"| meas_date=="2018-01-16") & (tis == 1 | tis == 5 | tis == 8)) %>%
  aov(main_fluo_median_norm ~ tis, data = .) %>%
  TukeyHSD()


p11 <- spect %>%
  filter((meas_date=="2017-12-07" | meas_date=="2018-01-16" | meas_date=="2017-12-13") & (reporter == "yegfp" | reporter == "mgfp")) %>%
  filter(tis == 1 | tis == 5 | tis == 8) %>%
  group_by(chassis,tis) %>%
  summarise(mean = mean(main_fluo_median_norm, na.rm = TRUE)) %>%
  spread(chassis, mean) %>%
  plot2("cho","yeast","chassis","paste(R ^ 2, \" = 0.98\")","paste(p, \" = 0.040\")")
p12 <- spect %>%
  filter(meas_date=="2017-12-07"| meas_date=="2018-01-16") %>%
  filter(tis == 1 | tis == 5 | tis == 8) %>%
  group_by(meas_date,reporter,tis) %>%
  #summarise(median = median(main_fluo_median_norm, na.rm = TRUE)) %>%
  spread(reporter, main_fluo_median_norm) %>% 
  #lm(zsgreen ~ mgfp, .) %>% summary()
  plot2("mgfp","zsgreen","reporter_cho","paste(R ^ 2, \" = 0.91\")","paste(p, \" = 4.02e-3\")")

plot_grid(p9,p10,p11,p12, ncol = 2)
ggsave(here('figures','Figure 5CDEF.pdf'), width = 13.4, height = 10, units = "in")

# correlations ------------------------------------------------------------
#PROMOTER
spect %>%
  filter(meas_date=="2018-01-12" & reporter == "yegfp") %>%
  group_by(promoter, tis) %>%
  summarise(mean = mean(main_fluo_median_norm, na.rm = TRUE)) %>%
  spread(promoter, mean) %>%
  # ggpairs()
  lm(tef1 ~ adh2, .) %>% summary()
#REPORTER
spect %>%
  filter(meas_date=="2017-12-13") %>%
  group_by(reporter,tis) %>%
  summarise(mean = mean(main_fluo_median_norm, na.rm = TRUE)) %>%
  spread(reporter, mean) %>%
  # lm(ymukg1 ~ mkate2, .) %>% summary()
  # lm(ymukg1 ~ yegfp, .) %>% summary()
  lm(mkate2 ~ yegfp, .) %>% summary()
  # ggpairs()
#CARBON_SOURCE
spect %>%
  filter(meas_date== "2017-12-13" & reporter == "yegfp" | meas_date=="2018-01-12" & promoter == "tef1") %>%
  group_by(carbon_source,tis) %>%
  summarise(mean = mean(main_fluo_median_norm, na.rm = TRUE)) %>%
  spread(carbon_source, mean) %>%
  lm(glucose ~ ethanol, .) %>% summary()
  # ggpairs()
#CHASSIS
spect %>%
  filter((meas_date=="2017-12-07" | meas_date=="2018-01-16" | meas_date=="2017-12-13") & (reporter == "yegfp" | reporter == "mgfp")) %>%
  filter(tis == 1 | tis == 5 | tis == 8) %>%
  group_by(chassis,tis) %>%
  summarise(mean = mean(main_fluo_median_norm, na.rm = TRUE)) %>%
  spread(chassis, mean) %>%
  lm(yeast ~ cho, .) %>% summary()
  #ggpairs()
#REPORTER CHO
spect %>%
  filter(meas_date=="2017-12-07"| meas_date=="2018-01-16") %>%
  filter(tis == 1 | tis == 5 | tis == 8) %>%
  group_by(meas_date,reporter,tis) %>%
  #summarise(median = median(main_fluo_median_norm, na.rm = TRUE)) %>%
  spread(reporter, main_fluo_median_norm) %>%
  lm(mgfp ~ zsgreen, .) %>% summary()
  #ggpairs()
    

