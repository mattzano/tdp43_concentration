#import chx data
###add info on gene expression first - i can try and use gene expression from curves
###label them as nmd or not (using the slope plots)
annotation_nmd_chx <- read.csv("~/Desktop/annotation_nmd.csv")

big_data_filtereds <- big_data[,c(1,2,5,9,11,17)] %>% #1701012/6 = 283502
  pivot_wider(names_from = .id, values_from = mean_psi_per_lsv_junction, values_fn = max) %>% ##232472
  left_join(annotation_nmd_chx) %>%
  filter(!is.na(Control_TDP43KD) & !is.na(Cycloheximide_TDP43KD)) %>%
  filter(no_dox < 0.05) %>%  #114616  ####to include all splicing events comment this line
  filter(dox_075-no_dox > 0.1) %>% #1669
  mutate(delta075025 = dox_075 - dox_025) %>%
  mutate(delta025021 = dox_025 - dox_021) %>%
  mutate(delta021018 = dox_021 - dox_0187) %>%
  mutate(delta018012 = dox_0187 - dox_0125) %>%
  mutate(delta012000 = dox_0125 - no_dox) %>%
  filter((delta075025   > - 0.1) &
           (delta025021 > - 0.1) &
           (delta021018 > - 0.1) &
           (delta018012 > - 0.1) &
           (delta012000 > - 0.1)) #1346
big_data_filtereds_dedup <- big_data_filtereds[!duplicated(big_data_filtereds$paste_into_igv_junction),] #1325

big_data_filtereds_normal <- big_data_filtereds_dedup %>%
  mutate(dox_075normal = (dox_075 - no_dox)  / (dox_075 - no_dox)) %>%
  mutate(dox_025normal = (dox_025 - no_dox)  / (dox_075 - no_dox)) %>%
  mutate(dox_021normal = (dox_021 - no_dox)  / (dox_075 - no_dox)) %>%
  mutate(dox_0187normal = (dox_0187 - no_dox)  / (dox_075 - no_dox)) %>%
  mutate(dox_0125normal = (dox_0125 - no_dox)  / (dox_075 - no_dox)) %>%
  mutate(no_doxnormal = (no_dox - no_dox)  / (dox_075 - no_dox)) %>%
  mutate(delta075025normal = dox_075normal - dox_025normal) %>%
  mutate(delta025021normal = dox_025normal - dox_021normal) %>%
  mutate(delta021018normal = dox_021normal - dox_0187normal) %>%
  mutate(delta018012normal = dox_0187normal - dox_0125normal) %>%
  mutate(delta012000normal = dox_0125normal - no_doxnormal) %>%
  filter((delta075025normal   > - 0.1) &
           (delta025021normal > - 0.1) &
           (delta021018normal > - 0.1) &
           (delta018012normal > - 0.1) &
           (delta012000normal > - 0.1))

plot_early_late <- big_data_filtereds %>%
  pivot_longer(cols = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075")) %>%
  mutate(name = factor(name, levels = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075"))) %>%
  #pivot_longer(cols = c("no_doxnormal","dox_0125normal", "dox_0187normal", "dox_021normal", "dox_025normal", "dox_075normal")) %>%
  #mutate(name = factor(name, levels = c("no_doxnormal","dox_0125normal", "dox_0187normal", "dox_021normal", "dox_025normal", "dox_075normal"))) %>%
  ggplot(mapping = aes(x = name, y = value, group = paste_into_igv_junction)) +
  geom_line(aes(color = color_gene_name, alpha = color_gene_name), show.legend = T) +
  #geom_text_repel(aes(label = label_junction), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  scale_color_manual(values = c("#3288BD","#D53E4F","gray")) +
  scale_alpha_manual(values = c(1,1,0.1)) +
  xlab("TDP-43 knockdown level") +
  ylab("Normalised PSI") +
  scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++")) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  
  #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_classic()
print(plot_early_late)

#count from chx and curves datasets and intersection (venn)



#verify expression and behaviour in curves






####nmd prediction early-late

###import darkcryptics

##




###eventually check expression - also to confirm dark cryptics downregulates genes


avg_read_counts <- featureCounts %>% 
  mutate(dox_0125 = rowMeans(dplyr::select(featureCounts, contains("DOX_0.0125")), na.rm = TRUE)) %>%
  mutate(dox_0187 = rowMeans(dplyr::select(featureCounts, contains("DOX_0.0187")), na.rm = TRUE)) %>%
  mutate(dox_021 = rowMeans(dplyr::select(featureCounts, contains("DOX_0.021")), na.rm = TRUE)) %>%
  mutate(dox_025 = rowMeans(dplyr::select(featureCounts, contains("DOX_0.025")), na.rm = TRUE)) %>%
  mutate(dox_075 = rowMeans(dplyr::select(featureCounts, contains("DOX_0.075")), na.rm = TRUE)) %>%
  mutate(no_dox = rowMeans(dplyr::select(featureCounts, contains("NT_0")), na.rm = TRUE)) %>%
  dplyr::select(c(20, 26, 21:25)) %>%
  pivot_longer(cols = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075"), 
               names_to = ".id", values_to = "avg_count") %>%
  mutate(.id = factor(.id, levels = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075"))) %>%
  group_by(gene_name, .id) %>%
  summarise(avg_count = max(avg_count))

avg_read_counts %>%
  mutate(.id2 = c(96.000000, 85.333333, 34.666667, 25.666667, 14.333333, 4.333333)) %>%
  dplyr::filter(gene_name == "CYFIP2") %>%
  ggplot(aes(x = .id2, y = avg_count)) + 
  geom_point() +
  #stat_summary(geom = "line", fun = mean)  +
  geom_smooth(method = "lm", se = F) + 
  scale_x_reverse() +
  xlab("TARDBP read counts") +
  ylab("CYFIP2 read counts") +
  theme_classic()
