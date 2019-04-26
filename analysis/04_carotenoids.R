library(tidyverse)
library(stringr)
library(here)
theme_set(theme_classic(30))

strain_description <- read_csv(here('data','strain_description.csv')) %>%
  filter(plate == "yp15", description != "empty", description != "media") %>%
  mutate(description = ifelse(description == "CEN.PK113-7D","0_0_0",description),
         well = str_c(row,column),
         erg9 = map_chr(description, ~ strsplit(.x, "_")[[1]][1]),
         crte = map_chr(description, ~ strsplit(.x, "_")[[1]][3]),
         erg9 = plyr::mapvalues(erg9, c("TGATAT","TAGGTT","TCGGTC"), c(1,5,8)),
         crte = plyr::mapvalues(crte, c("TGATAT","TAGGTT","TCGGTC"), c(1,5,8)),
         description = str_c(erg9,crte,sep="_")) %>%
  left_join(read_csv(here('data','platereader.csv')), by = c("hour","plate","row","column","well"))

by_well <- strain_description %>%
  left_join(read_csv(here('data','HPLC.csv')) %>% filter(date == "2018-05-16"), by = c("hour","plate","row","column", "well")) %>%
  mutate(od = od*10,
         amount = (amount * 5)/od)

strain_description$description <- factor(strain_description$description, levels = c("0_0","1_1","1_5","1_8","5_1","5_5","5_8","8_1","8_5","8_8"))
by_well %>%
  group_by(description) %>%
  summarise(amount_mean = mean(amount),
            amount_se  = sd(amount) / sqrt(n())) %>%
  ggplot(mapping = aes(x=description, y=amount_mean)) +
  geom_bar(stat="identity", fill = "white", colour = "black") +
  geom_errorbar(aes(ymin=amount_mean - amount_se, ymax=amount_mean + amount_se), width = 0.2) +
  labs(x = "TIS combination, (erg9 crte)", y = "amount, (µg/ml)/od") +
  coord_flip() +
  annotate("text",x=1,y=5, hjust = 0,label="a") +
  annotate("text",x=2,y=5, hjust = 0,label="b, c") +
  annotate("text",x=3,y=5, hjust = 0,label="b") +
  annotate("text",x=4,y=5, hjust = 0,label="a") +
  annotate("text",x=5,y=5, hjust = 0,label="b, c, d") +
  annotate("text",x=6,y=5, hjust = 0,label="b, c, d") +
  annotate("text",x=7,y=5, hjust = 0,label="a") +
  annotate("text",x=8,y=5, hjust = 0,label="d") +
  annotate("text",x=9,y=5, hjust = 0,label="d") +
  #annotate("text",x=10,y=5, hjust = 0,label="") +
  theme()
ggsave(here('figures','Figure 6D.pdf'), width = 7.5, height = 10, units = "in")

res.aov <- aov(amount ~ description, data = by_well)
summary(res.aov)
TukeyHSD(res.aov)

