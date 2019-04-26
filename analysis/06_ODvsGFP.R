library(tidyverse)
library(here)
theme_set(theme_classic(25))

spect <- read_csv(here('data', 'spect.csv'),guess_max = 5000) %>%
  filter(!is.na(reporter), meas_date == "2017-05-26") %>%
  .[["main_od"]] %>%
  min()

read_csv(here('data', 'spect.csv'),guess_max = 5000) %>%
  filter(!is.na(reporter), meas_date == "2017-05-26") %>%
  select(main_od, main_gfp_median) %>%
  mutate(main_od = main_od-min(main_od)) %>%
  mutate(main_gfp_median_norm = main_gfp_median/mean(main_gfp_median, na.rm = TRUE)*100) %>%
  ggplot() +
  geom_point(mapping = aes(x=main_od, y=main_gfp_median), size = 3) +
  labs(x = "OD600", y = "Median fluorescence")
ggsave(here('figures', 'Figure S5_od_gfp.pdf'))
