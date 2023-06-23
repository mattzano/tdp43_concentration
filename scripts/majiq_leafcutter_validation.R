#### validation of majiq and leafcutter hits

table_validation <- read_xlsx("~/Downloads/manual_annotations_shsy5y_matteo.xlsx")

# plot yes no unsure
table_validation[c(1:218),] %>%
  ggplot(aes(x = `manual annotation`)) +
  geom_bar() +
  theme_classic()

# same but with increasing experience
table_validation[c(1:218),] %>%
  ggplot(aes(x = `manual annotation`, y = juncID)) +
  geom_point() +
  coord_flip() +
  theme_classic()

# yes / no / unsure by junc cat
table_validation[c(1:218),] %>%
  ggplot(aes(x = `manual annotation`, y = junc_cat)) +
  #geom_point(position = position_dodge(width = 0.3)) +
  geom_count() +
  theme_classic()

# yes / no / unsure by psi_ctrl < 0.05 - similar to the one before

# yes / no / unsure by p value (is the software confidence similar to validation one) 
table_validation[c(1:218),] %>%
  ggplot(aes(x = `manual annotation`, y = -log10(pvalue))) +
  geom_point() +
  theme_classic()


### with psi
table_validation[c(1:218),] %>%
  ggplot(aes(x = `manual annotation`, y = dPSI)) +
  geom_point() +
  theme_classic()



