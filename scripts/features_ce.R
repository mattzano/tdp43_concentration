#manval <- manual_validation_sy
delta_filtereds_early_late <- big_delta %>%
  select(c(1:11,13,18:19,22:23,35:36)) %>%
  pivot_wider(names_from = .id, values_from = c(mean_dpsi_per_lsv_junction, probability_changing)) %>%
  filter(mean_dpsi_per_lsv_junction_dox0075 > 0.1 &
           probability_changing_dox0075 > 0.9) %>%
  mutate_if(is.numeric, replace_na, 0) %>%
  group_by(lsv_junc) %>%
  mutate(category = ifelse(mean_dpsi_per_lsv_junction_dox00125 > 0.2, "Early",
                           ifelse(mean_dpsi_per_lsv_junction_dox00187 > 0.2, "Eerlier",
                                  ifelse(mean_dpsi_per_lsv_junction_dox0021 > 0.2 |
                                           mean_dpsi_per_lsv_junction_dox0025 > 0.2, "Intermediate", "Late")))) %>%
  mutate(is_cryptic = ifelse(base_mean_psi < 0.05, "yes", "no")) %>%
  mutate(max_dpsi_change = mean_dpsi_per_lsv_junction_dox0075)


delta_filtereds_early_late <- big_delta %>% ##74067
  select(c(1:11,13,18:19,22:23,35:36)) %>%
  group_by(lsv_junc) %>%
  filter(n() == 5) %>% ##5925
  pivot_wider(names_from = .id, values_from = c(mean_dpsi_per_lsv_junction, probability_changing)) %>% ##1185
  filter(mean_dpsi_per_lsv_junction_dox0075 > 0.1 &
           probability_changing_dox0075 > 0.9) %>% #128
  group_by(lsv_junc) %>%
  mutate(category = ifelse(mean_dpsi_per_lsv_junction_dox00125 > 0.2, "Early",
                           ifelse(mean_dpsi_per_lsv_junction_dox00187 > 0.2, "Eerlier",
                                  ifelse(mean_dpsi_per_lsv_junction_dox0021 > 0.2 |
                                         mean_dpsi_per_lsv_junction_dox0025 > 0.2, "Intermediate", "Late")))) %>%
  mutate(is_cryptic = ifelse(base_mean_psi < 0.05, "yes", "no")) %>%
  mutate(max_dpsi_change = mean_dpsi_per_lsv_junction_dox0075)
table(delta_filtereds_early_late$category, delta_filtereds_early_late$is_cryptic)
length(unique(delta_filtereds_early_late$lsv_junc))

delta_filtereds_early_late_long <- delta_filtereds_early_late %>%
  pivot_longer(cols = c(16:25),
    names_to = ".id", values_to = "mean_dpsi_per_lsv_junction_probability_changing") %>%
  mutate(value = stringr::str_match(.id, '(.*)_(.*)')[,2],
         .id = stringr::str_match(.id, '(.*)_(.*)')[,3]) %>%
  pivot_wider(names_from = value, values_from = mean_dpsi_per_lsv_junction_probability_changing)
  

#plot max dpsi distribution by category
delta_filtereds_early_late_long %>%
  ggplot(aes(x = category, y = max_dpsi_change)) +
  geom_boxplot() +
  geom_point(aes(color = category)) +
  annotate("text",
           x = 1:length(table(delta_filtereds_early_late_long$category)),
           y = aggregate(max_dpsi_change ~ category, delta_filtereds_early_late_long, median)[ , 2],
           label = table(delta_filtereds_early_late_long$category),
           vjust = -12, hjust = -0.1) +
  scale_color_brewer(palette = "Set1") +
  theme_classic()

#plot max dpsi distribution by category and ce
delta_filtereds_early_late_long %>%
  ggplot(aes(x = category, y = max_dpsi_change)) +
  geom_boxplot(aes(fill = category)) +
  #geom_text(aes(label = count)) +
  geom_point(aes(color = category)) +
  facet_wrap(facets = vars(is_cryptic)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Pastel1") +
  theme_classic()

delta_filtereds_early_late_long %>%
  ggplot(aes(x = is_cryptic, y = max_dpsi_change)) +
  geom_boxplot(aes(fill = category)) +
  #geom_point(aes(color = category)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Pastel1") +
  #facet_wrap(facets = vars(is_cryptic)) +
  theme_classic()

#spaghettini
delta_filtereds_early_late_long %>%
  ggplot(aes(x = .id, y = mean_dpsi_per_lsv_junction, group = lsv_junc, color = category)) +
  geom_line() +
  facet_wrap(ncol = 4, facets = vars(is_cryptic, category)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))


#correlation between mean psi and dpsi
delta_filtereds_early_late_long %>%
  ggplot(aes(x = base_mean_psi, y = mean_dpsi_per_lsv_junction)) +
  geom_point(aes(color = category)) +
  geom_smooth(method = "lm", se = F, color = "black") +
  stat_cor() +
  facet_wrap(ncol = 1, facets = vars(category)) +
  theme_classic()

#nmd sensitivity
nmd_exp <- read.table("~/Desktop/nmd_chx.csv", sep = ",", header = T) %>%
  filter(.id == "Control_Control")

delta_filtereds_early_late_long_nmd <- delta_filtereds_early_late_long %>%
  left_join(nmd_exp[,c(12,19)])

table(delta_filtereds_early_late_long_nmd$color_gene_name)

delta_filtereds_early_late_long_nmd %>%
  ggplot(aes(x = category, y = max_dpsi_change)) +
  geom_boxplot(aes(fill = color_gene_name)) +
  #geom_point(aes(color = category)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Pastel1") +
  #facet_wrap(facets = vars(is_cryptic)) +
  theme_classic()

delta_filtereds_early_late_long_nmd %>%
  ggplot(aes(x = .id, y = mean_dpsi_per_lsv_junction, group = lsv_junc, color = category)) +
  geom_line() +
  facet_wrap(ncol = 4, facets = vars(color_gene_name, category)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))




