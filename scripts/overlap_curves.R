##sy5y
input_list <- c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075")
datalist <- list()
for (i in input_list) {
  ####change path here
  input_splicing <- fread(paste0("/Users/matteozanovello/Documents/GitHub/tdp43_concentration/data/majiq_single/", i, "_parsed.csv"))
  datalist[[i]] <- input_splicing
}
big_data <- rbindlist(datalist, idcol = TRUE)
big_data$.id <- factor(big_data$.id, levels = input_list)


###filter
big_data_filtereds <- big_data[,c(1,2,5,9,11,12,17)] %>% #1701012/6 = 283502
  pivot_wider(names_from = .id, values_from = mean_psi_per_lsv_junction) %>% ##232472
  #filter(no_dox < 0.05) %>%  #114616  ####to include all splicing events comment this line
  filter(dox_075-no_dox > 0.1) %>% #1669
  filter((dox_075 - dox_025     > - 0.1) &
           (dox_025 - dox_021   > - 0.1) &
           (dox_021 - dox_0187  > - 0.1) &
           (dox_0187 - dox_0125 > - 0.1) &
           (dox_0125 - no_dox   > - 0.1)) %>%
  mutate(cryptic = ifelse(no_dox < 0.05, "cryptic", "no_cryptic"))#1346
big_data_filtereds_dedup <- big_data_filtereds[!duplicated(big_data_filtereds$paste_into_igv_junction),] #1325

big_data_filtereds_early_late_new <- big_data_filtereds_dedup %>%
  mutate(color_gene_name = as.factor(as.character(ifelse((dox_0125-no_dox > 0.1 | dox_0187-no_dox > 0.1) & 
                                                           dox_075-no_dox > 0.2, "early", 
                                                         ifelse(dox_025-no_dox < 0.1 & 
                                                                  dox_075-no_dox > 0.2, "late",
                                                                ifelse((dox_021-no_dox > 0.1 | dox_025-no_dox > 0.1) & 
                                                                         dox_075-no_dox > 0.2, "intermediate",
                                                                       "none"))))))
table(big_data_filtereds_early_late_new$color_gene_name, big_data_filtereds_early_late_new$cryptic)
 
plot_early_late <- big_data_filtereds_early_late_new %>%
  pivot_longer(cols = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075")) %>%
  mutate(name = factor(name, levels = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075"))) %>%
  ggplot(mapping = aes(x = name, y = value, group = paste_into_igv_junction)) +
  geom_line(aes(color = color_gene_name, alpha = color_gene_name), show.legend = F) +
  #geom_text_repel(aes(label = label_junction), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  scale_color_manual(values = c("#FAAC5E","#D53E4F","#3288BD","gray")) +
  scale_alpha_manual(values = c(1,1,1,0.1)) +
  #geom_hline(yintercept = 0.1, linetype = "dotted") +
  xlab("TDP-43 knockdown level") +
  ylab(expression(paste(Psi))) +
  scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++")) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  facet_wrap(cryptic ~ color_gene_name, nrow = 2) +
  #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_classic() +
  theme(text = element_text(size = 22))
print(plot_early_late)
#write.csv(big_data_filtereds_early_late, "~/Desktop/early_late.csv", row.names = F)


###dz
input_list_dz <- c("DZ_curves_0","DZ_curves_002", "DZ_curves_003", "DZ_curves_004", "DZ_curves_005", "DZ_curves_0075", "DZ_curves_01", "DZ_curves_1")
datalistas <- list()

for (i in input_list_dz) {
  input <- paste(i, "parsed", sep = "_")
  input_csv <- paste(input, "csv", sep = ".")
  input_splicing <- fread(file.path(here::here(), "data", "majiq_dz", input_csv))
  datalistas[[i]] <- input_splicing
}
big_data_dz <- rbindlist(datalistas, idcol = TRUE)
big_data_dz$.id <- factor(big_data_dz$.id, levels = input_list_dz)

#big_data_dz <- rbind(big_data_dz, big_data_dz[47996,])
big_data_dz_filtereds <- big_data_dz[,c(1,2,5,9,11,12,17)] %>% #1901543/6 = 316923.8
  pivot_wider(names_from = .id, values_from = mean_psi_per_lsv_junction) %>% ##195298
  #filter(DZ_curves_0 < 0.05) %>%  #91339 ##remove this to include all splicing changes
  filter(DZ_curves_1 - DZ_curves_0 > 0.1) %>% #1376
  filter((DZ_curves_1 - DZ_curves_01 > - 0.1) &
         (DZ_curves_01 - DZ_curves_0075 > - 0.1) &
         (DZ_curves_0075 - DZ_curves_005 > - 0.1) &
         (DZ_curves_005 - DZ_curves_004 > - 0.1) &
         (DZ_curves_004 - DZ_curves_003 > - 0.1) &
         (DZ_curves_003 - DZ_curves_002 > - 0.1) &
         (DZ_curves_002 - DZ_curves_0 > - 0.1)) %>%
  mutate(cryptic = ifelse(DZ_curves_0 < 0.05, "cryptic", "no_cryptic")) #915
big_data_dz_filtereds_dedup <- big_data_dz_filtereds[!duplicated(big_data_dz_filtereds$paste_into_igv_junction),] #901

big_data_dz_filtereds_early_late <- big_data_dz_filtereds_dedup %>%
  mutate(color_gene_name = as.factor(as.character(ifelse((DZ_curves_002-DZ_curves_0 > 0.1 | DZ_curves_003-DZ_curves_0 > 0.1) & 
                                                           DZ_curves_1-DZ_curves_0 > 0.2, "early", 
                                                         ifelse(DZ_curves_0075-DZ_curves_0 < 0.1 & 
                                                                  (DZ_curves_01-DZ_curves_0 > 0.1 | DZ_curves_1-DZ_curves_0 > 0.2), "late",
                                                                ifelse((DZ_curves_004-DZ_curves_0 > 0.1 | DZ_curves_005-DZ_curves_0 > 0.1 | DZ_curves_0075-DZ_curves_0 > 0.1) &
                                                                         DZ_curves_1-DZ_curves_0 > 0.2, "intermediate",
                                                                       "none")))))) 
table(big_data_dz_filtereds_early_late$color_gene_name, big_data_dz_filtereds_early_late$cryptic)

plot_dz_early_late <- big_data_dz_filtereds_early_late %>%
  pivot_longer(cols = c("DZ_curves_0","DZ_curves_002", "DZ_curves_003", "DZ_curves_004", "DZ_curves_005", "DZ_curves_0075", "DZ_curves_01", "DZ_curves_1")) %>%
  mutate(name = factor(name, levels = c("DZ_curves_0","DZ_curves_002", "DZ_curves_003", "DZ_curves_004", "DZ_curves_005", "DZ_curves_0075", "DZ_curves_01", "DZ_curves_1"))) %>%
  ggplot(mapping = aes(x = name, y = value, group = paste_into_igv_junction)) +
  geom_line(aes(color = color_gene_name, alpha = color_gene_name), show.legend = F) +
  #geom_text_repel(aes(label = label_junction), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  scale_color_manual(values = c("#FAAC5E","#D53E4F","#3288BD","gray")) +
  scale_alpha_manual(values = c(1,1,1,0.1)) +
  xlab("TDP-43 knockdown level") +
  ylab(expression(paste(Psi))) +
  scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++", "++++++", "+++++++")) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  facet_wrap(cryptic ~ color_gene_name, nrow = 2) +
  #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_classic()
print(plot_dz_early_late)
#write.csv(big_data_filtereds_early_late, "~/Desktop/early_late_dz.csv", row.names = F)

merged_curves <- full_join(big_data_filtereds_early_late_new, big_data_dz_filtereds_early_late, 
                           by = c("paste_into_igv_junction", "de_novo_junctions", "strand", "gene_name")) %>%
  mutate(overlap = ifelse(`color_gene_name.x` == `color_gene_name.y`, "yes", "no"))
overlap <- as.data.frame(table(merged_curves$cryptic.x, merged_curves$cryptic.y, merged_curves$color_gene_name.x, merged_curves$color_gene_name.y))
names(overlap) <- c("CE", "SY5Y", "DZ", "Freq")
#write.table(overlap, file = "~/Desktop/overlap.csv", row.names = F, sep = ",")

###import manual validation table
table_validation <- read_xlsx("~/Downloads/manual_annotations_shsy5y_matteo.xlsx") %>%
  filter(!is.na(`manual annotation`)) %>%
  mutate(paste_into_igv_junction = paste0(chr, ":", start, "-", end))

merged_curves_annotation <- left_join(merged_curves, table_validation, by = "paste_into_igv_junction")

