library(tidyverse)
library(googlesheets)
library(stringr)
library(growthrates)
library(lattice)
library(here)


strain_description <- read_csv(here('data','strain_description.csv')) %>%
  filter(plate == "yp15", description != "empty", description != "media") %>%
  mutate(description = ifelse(description == "CEN.PK113-7D","0_0_0",description),
         well = str_c(row,column),
         erg9 = map_chr(description, ~ strsplit(.x, "_")[[1]][1]),
         crte = map_chr(description, ~ strsplit(.x, "_")[[1]][3]),
         erg9 = plyr::mapvalues(erg9, c("TGATAT","TAGGTT","TCGGTC"), c(1,5,8)),
         crte = plyr::mapvalues(crte, c("TGATAT","TAGGTT","TCGGTC"), c(1,5,8)),
         description = str_c(erg9,crte,sep="_"))

strains <- strain_description %>%
  left_join(read_csv(here('data','kinetic.csv'), guess_max = 30000), by = c("plate","row","column", "well")) %>%
  filter(od > 0) #%>%
  # filter(column == 10)

xyplot(od ~ time|description, data = strains, pch = 16, cex = 0.5)
many_easylinear <- all_easylinear(od ~ time|well, h= 10, data = strains, pch = 16, cex = 0.5)

# pdf("../figures/Figure SX growth fit ex.pdf", width = 9, height = 6)
# par(mfrow = c(2,3), mar = c(2,3,2,1), cex=1)
# plot(many_easylinear, log="y")
# par(mfrow = c(1,1), mar = c(2.5,4,2,1), cex=1)
# dev.off()
# 

pdf(here('figures','Figure S11 growth fit.pdf'), width = 6, height = 9)
par(mfrow = c(9,6), mar = c(2,3,2,1), cex=0.5)
plot(many_easylinear, log="y")
par(mfrow = c(1,1), mar = c(2.5,4,2,1), cex=1)
dev.off()

as.tibble(rsquared(many_easylinear)) %>% ggplot() + geom_histogram(aes(x=value))
library(GGally)
ggpairs(as.tibble(coef(many_easylinear)))
summary(as.tibble(coef(many_easylinear)))
#rownames_to_column(as.data.frame(coef(many_easylinear)), "well")[,c("well","mumax")]
#summarise(rownames_to_column(as.data.frame(coef(many_easylinear)), "well")[,c("well","mumax")], sd=sd(mumax))

strain_description$description <- factor(strain_description$description, levels = c("0_0","1_1","1_5","1_8","5_1","5_5","5_8","8_1","8_5","8_8"))
#strain_description$description <- factor(strain_description$description, levels = rev(c("1_8","8_5","8_1","8_8","5_8","1_1","5_5","5_1","1_5")))
# strain_description$description <- factor(strain_description$description, levels = rev(c("1_1","1_5","1_8","5_5","5_8","5_1","8_5","8_8","8_1","0_0")))
theme_set(theme_classic(30))
strain_description %>%
  #filter(!is.na(description)) %>%
  left_join(rownames_to_column(as.data.frame(coef(many_easylinear)), "well"), by = c("well")) %>%
  mutate(mumax = mumax * 60) %>%
  group_by(description) %>%
  summarise(u_mean = mean(mumax),
            u_se  = sd(mumax) / sqrt(n())) %>%
  ggplot(mapping = aes(x=description, y=u_mean)) +
  geom_bar(stat="identity", fill = "white", colour = "black") +
  geom_errorbar(aes(ymin=u_mean - u_se, ymax=u_mean + u_se), width = 0.2) +
  #scale_y_continuous(limits = c(0, 0.20)) +
  labs(x = "TIS combination, (erg9 crte)", y = "growth rate, h-1") +
  coord_flip() +
  #coord_flip(ylim=c(0.10, 0.25)) +
  #annotate("segment", x=c(1,1,2),xend=c(1,2,2), y= c(0.3,0.32,0.32), yend=c(0.32,0.32,0.3)) +
  #annotate("segment", x=c(2,2,3),xend=c(2,3,3), y= c(0.33,0.35,0.35), yend=c(0.35,0.35,0.33)) +
  #annotate("text",x=1.5,y=.4,label="p=0.012") +
  annotate("text",x=1,y=.1,label="a") +
  annotate("text",x=2,y=.1,label="a") +
  annotate("text",x=3,y=.1,label="a") +
  annotate("text",x=4,y=.1,label="") +
  annotate("text",x=5,y=.1,label="") +
  annotate("text",x=6,y=.1,label="") +
  annotate("text",x=7,y=.1,label="b") +
  annotate("text",x=8,y=.1,label="b") +
  annotate("text",x=9,y=.1,label="") +
  #annotate("text",x=10,y=.1,label="") +
  theme()
ggsave(here('figures','Figure 6C.pdf'), width = 7.5, height = 10, units = "in")
#ggsave("../figures/figure6C_u_ordered_test.pdf", width = 7.5, height = 10, units = "in")

summary(rsquared(many_easylinear))

# tests -------------------------------------------------------------------
df <- strain_description %>%
  left_join(rownames_to_column(as.data.frame(coef(many_easylinear)), "well"), by = c("well")) %>%
  mutate(mumax = mumax * 60) %>%
  select(description, mumax)
res.aov <- aov(mumax ~ description, data = df)
summary(res.aov)
TukeyHSD(res.aov)

  
  
