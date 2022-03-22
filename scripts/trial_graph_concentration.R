#####import dpsi table and filter for relevant events
base <- "noDox"
contrast_input <- c("dox00125", 
                    "dox00187", 
                    "dox0021", 
                    "dox0025", 
                    "dox0075")
datalista <- list()

for (i in contrast_input) {
  contrast <- paste(i, "annotated", "junctions", sep = "_")
  input_bc <- paste(base, contrast, sep = "-")
  input_bc_csv <- paste(input_bc, "csv", sep = ".")
  input_bc_splicing <- fread(file.path(here::here(), "data", "majiq_delta", input_bc_csv))
  names(input_bc_splicing)[13] <- "contrast_mean_psi"
  datalista[[i]] <- input_bc_splicing
}

big_delta <- rbindlist(datalista, idcol = TRUE)

big_delta_filter <- big_delta %>%
  filter(no_dox_mean_psi < 0.05 & contrast_mean_psi >= 0.1)
  #filter(mean_dpsi_per_lsv_junction > 0 & probability_changing > 0.9)

######import parsed data, merge them, normalize data from 0 (ctrl) to 1 (0.075) !!!
input_list <- c("no_dox",
           "dox_0125", 
           "dox_0187", 
           "dox_021", 
           "dox_025", 
           "dox_075")
datalist <- list()

for (i in input_list) {
  input <- paste(i, "parsed", sep = "_")
  input_csv <- paste(input, "csv", sep = ".")
  input_splicing <- fread(file.path(here::here(), "data", "majiq_single", input_csv))
  datalist[[i]] <- input_splicing
}

big_data <- rbindlist(datalist, idcol = TRUE)
big_data$.id <- factor(big_data$.id, levels = c("no_dox",
                                                "dox_0125", 
                                                "dox_0187", 
                                                "dox_021", 
                                                "dox_025", 
                                                "dox_075"))


big_data_filtereds <- big_data %>%
  filter(lsv_junc %in% big_delta$lsv_junc) %>% #try without cryptics
  group_by(lsv_junc) %>%
  filter(n() == 6) %>%
  filter(mean_psi_per_lsv_junction[.id == "no_dox"] < 0.05) %>%
  filter(mean_psi_per_lsv_junction[.id == "dox_025"] > 0.01) %>%
  filter((mean_psi_per_lsv_junction[.id == "dox_075"] -  mean_psi_per_lsv_junction[.id == "dox_025"]  > -0.1) &
         (mean_psi_per_lsv_junction[.id == "dox_025"] -  mean_psi_per_lsv_junction[.id == "dox_021"]  > -0.1) &
         (mean_psi_per_lsv_junction[.id == "dox_021"] -  mean_psi_per_lsv_junction[.id == "dox_0187"] > -0.1) &
         (mean_psi_per_lsv_junction[.id == "dox_0187"] - mean_psi_per_lsv_junction[.id == "dox_0125"] > -0.1) &
         (mean_psi_per_lsv_junction[.id == "dox_0125"] - mean_psi_per_lsv_junction[.id == "no_dox"]   > -0.1))

big_data_filtereds_a <- big_data_filtereds %>% filter(.id == "dox_075")
big_data_filtereds <- right_join(big_data_filtereds, big_data_filtereds_a[, c("lsv_junc","mean_psi_per_lsv_junction")], by = "lsv_junc")
colnames(big_data_filtereds)[5] <- "mean_psi_per_lsv_junction"
colnames(big_data_filtereds)[19] <- "mean_psi_per_lsv_junction_075"

big_data_filtereds <- big_data_filtereds %>%
  group_by(lsv_junc) %>%
  #filter(n() == 6) %>%
  mutate(mean_psi_per_lsv_junction_normal = (mean_psi_per_lsv_junction - mean_psi_per_lsv_junction[.id == "no_dox"]) 
                                                    / (mean_psi_per_lsv_junction_075 - mean_psi_per_lsv_junction[.id == "no_dox"])) %>%
  filter((mean_psi_per_lsv_junction_normal[.id == "dox_075"] - mean_psi_per_lsv_junction_normal[.id == "dox_025"] > -0.1) &
      (mean_psi_per_lsv_junction_normal[.id == "dox_025"] - mean_psi_per_lsv_junction_normal[.id == "dox_021"] > -0.1) &
      (mean_psi_per_lsv_junction_normal[.id == "dox_021"] - mean_psi_per_lsv_junction_normal[.id == "dox_0187"] > -0.1) &
      (mean_psi_per_lsv_junction_normal[.id == "dox_0187"] - mean_psi_per_lsv_junction_normal[.id == "dox_0125"] > -0.1) &
      (mean_psi_per_lsv_junction_normal[.id == "dox_0125"] - mean_psi_per_lsv_junction_normal[.id == "no_dox"] > -0.1))

big_data_filtereds_a <- big_data_filtereds %>% filter(.id == "dox_025")
big_data_filtereds <- right_join(big_data_filtereds, big_data_filtereds_a[, c("lsv_junc","mean_psi_per_lsv_junction_normal")], by = "lsv_junc")
colnames(big_data_filtereds)[20] <- "mean_psi_per_lsv_junction_normal"
colnames(big_data_filtereds)[21] <- "mean_psi_per_lsv_junction_normal_025"

big_data_filtereds_a <- big_data_filtereds %>% filter(.id == "dox_0187")
big_data_filtereds <- right_join(big_data_filtereds, big_data_filtereds_a[, c("lsv_junc","mean_psi_per_lsv_junction_normal")], by = "lsv_junc")
colnames(big_data_filtereds)[20] <- "mean_psi_per_lsv_junction_normal"
colnames(big_data_filtereds)[22] <- "mean_psi_per_lsv_junction_normal_0187"

  #mutate(color_gene_name = as.factor(as.character(ifelse(paste_into_igv_junction %in% c("chr19:17641556-17642414","chr8:79611214-79616822"), 1, 0))))
big_data_filtereds_b <- big_data_filtereds %>%
  group_by(lsv_junc) %>% 
  mutate(color_gene_name = as.factor(as.character(ifelse(paste_into_igv_junction %in% c("chr19:17641556-17642414","chr8:79611214-79616822"), 1, 0))))

big_data_filtereds <- big_data_filtereds %>%
  group_by(lsv_junc) %>% 
  dplyr::mutate(color_gene_name = as.factor(as.character(ifelse(paste_into_igv_junction == "chr19:17641556-17642414", 1, 
                                                                ifelse(paste_into_igv_junction == "chr19:17653320-17655274",2,
                                                                ifelse(paste_into_igv_junction == "chr19:17623216-17623542",3,4))))))
                                            
  #mutate(label_junction = case_when(((.id =="dox_0187") ~ gene_name, T ~ "")))

plot <- big_data_filtereds %>%
  dplyr::filter(paste_into_igv_junction == "chr8:79611214-79616822") %>%
  ggplot(mapping = aes(x = .id, y = mean_psi_per_lsv_junction, group = lsv_junc)) +
  #geom_line(color = "#3288BD", show.legend = F) +
  geom_bar(fill ="#3288BD", stat = "identity", width = 0.5) +
  geom_errorbar(aes(x = .id, ymin = mean_psi_per_lsv_junction - (stdev_psi_per_lsv_junction/sqrt(2)),
                             ymax = mean_psi_per_lsv_junction + (stdev_psi_per_lsv_junction/sqrt(2))), width = 0.3) +
  #geom_hline(aes(yintercept = 0.1)) +
  #geom_text_repel(aes(x = .id, y = mean_psi_per_lsv_junction, label = gene_name), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  #scale_color_manual(values = c("gray",
  #                              "#3288BD",
  #                              "#D53E4F",
  #                              "#FDAE61"
  #                              )) +
  #scale_alpha_manual(values = c(1)) +
  xlab("TDP-43 protein levels") +
  ylab("STMN2 cryptic exon PSI") +
  scale_x_discrete(labels = c("100%", "77%", "25%", "8%", "4%", "0%")) +
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1),
    labels = c("0%","20%","40%","60%","80%","100%")) +
  theme_classic()
print(plot)
ggsave("barplot_stmn.pdf", plot, width = 8, height = 5)

plot <- big_data_filtereds_c %>%
  #filter(gene_name == "STMN2") %>%
  ggplot(mapping = aes(x = .id, y = mean_psi_per_lsv_junction_normal, group = lsv_junc)) +
  geom_line(aes(color = color_gene_name, alpha = color_gene_name), show.legend = F) +
  #geom_text_repel(aes(x = .id, y = mean_psi_per_lsv_junction, label = label_junction), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  scale_color_manual(values = c("#3288BD",
                                "#D53E4F",
                                
                                
                                "#FDAE61",
                                "gray")) +
  scale_alpha_manual(values = c(1,1,1,0.2)) +
  xlab("TDP-43 knockdown level") +
  ylab("Normalised PSI") +
  scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++")) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_classic()
print(plot)

ggsave("spaghetti.png", plot, width = 8, height = 5)





big_data_matrix = big_delta %>%
  filter(junc_cat != "annotated") %>%
  mutate(is_a_cryptic = no_dox_mean_psi < 0.05 & contrast_mean_psi >= 0.1) %>%
           #abs(mean_dpsi_per_lsv_junction) > 0.15 & probability_changing > 0.9) %>% # 
  mutate(name = glue::glue("{gene_name}|{junc_cat}|{paste_into_igv_junction}")) %>%
  dplyr::select(name,is_a_cryptic,.id) %>%
  unique() %>%
  filter(is_a_cryptic == T)%>% 
  pivot_wider(values_from = "is_a_cryptic",
              id_cols = "name",
              names_from = ".id",
              values_fill = FALSE, values_fn = sum) %>%
  column_to_rownames('name') %>%
  mutate_all(as.logical)
  
upset(big_data_matrix * 1,
      order.by = "freq", keep.order = TRUE,
      mainbar.y.label = "Junction Intersections", 
      sets.x.label = "Junction Per Level of TDP-43 KD")




avg_read_counts <- featureCounts %>%
  mutate(dox_0125 = rowMeans(dplyr::select(featureCounts, contains("DOX_0.0125")), na.rm = TRUE)) %>%
  mutate(dox_0187 = rowMeans(dplyr::select(featureCounts, contains("DOX_0.0187")), na.rm = TRUE)) %>%
  mutate(dox_021 = rowMeans(dplyr::select(featureCounts, contains("DOX_0.021")), na.rm = TRUE)) %>%
  mutate(dox_025 = rowMeans(dplyr::select(featureCounts, contains("DOX_0.025")), na.rm = TRUE)) %>%
  mutate(dox_075 = rowMeans(dplyr::select(featureCounts, contains("DOX_0.075")), na.rm = TRUE)) %>%
  mutate(no_dox = rowMeans(dplyr::select(featureCounts, contains("NT_0")), na.rm = TRUE)) %>%
  dplyr::select(c(20, 26, 21:25)) %>%
  pivot_longer(cols = c("no_dox",
                        "dox_0125", 
                        "dox_0187", 
                        "dox_021", 
                        "dox_025", 
                        "dox_075"), 
               names_to = ".id", 
               values_to = "avg_count") %>%
  group_by(gene_name, .id) %>%
  summarise(avg_count = max(avg_count))

avg_read_counts$.id <- factor(avg_read_counts$.id, levels = c("no_dox",
                                                              "dox_0125", 
                                                              "dox_0187", 
                                                              "dox_021", 
                                                              "dox_025", 
                                                              "dox_075"))



avg_read_counts %>%
  filter(gene_name %in% c("STMN2","UNC13A","RANBP1","MYO9B", "STX2", "RELL1")) %>%
  #filter(avg_count[.id == "no_dox"] < 10000) %>%
  #group_by(gene_name) %>%
  #filter(n() == 6) %>%
  #filter(avg_count[.id == "no_dox"] > avg_count[.id == "dox_0187"] &
  #       avg_count[.id == "dox_0187"] > avg_count[.id == "dox_021"] &
  #       avg_count[.id == "dox_021"] > avg_count[.id == "dox_025"] &
  #       avg_count[.id == "dox_025"] > avg_count[.id == "dox_075"]) %>%
  #filter(avg_count[.id == "no_dox"] - avg_count[.id == "dox_075"] > 1000) %>%
  mutate(label_junction = case_when(.id ==" dox_0187"~ gene_name, T ~ "")) %>%
  ggplot(aes(x = .id, y = avg_count, group = gene_name)) +
  geom_point(aes(color = gene_name), show.legend = F) +
  geom_line(aes(color = gene_name), show.legend = T) +
  #geom_text_repel(aes(label = label_junction), max.overlaps = Inf, size=4, show.legend = F) +
  #scale_y_log10() +
  theme_classic()


plot_not_normal <- big_data_filtereds_b %>%
  filter(mean_psi_per_lsv_junction[.id == "dox_025"] < 0.05 &
         mean_psi_per_lsv_junction[.id == "dox_075"] > 0.1 &
         mean_psi_per_lsv_junction[.id == "no_dox"] < 0.05) %>%
  ggplot(mapping = aes(x = .id, y = mean_psi_per_lsv_junction, group = lsv_junc)) +
  geom_line(aes(color = gene_name), show.legend = T) +
  #geom_text_repel(aes(x = .id, y = mean_psi_per_lsv_junction, label = label_junction), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  #scale_color_manual(values = c("#3288BD",
  #                              "#D53E4F",
  #                              "#FDAE61",
  #                              "gray")) +
  scale_alpha_manual(values = c(1,1,1,0.2)) +
  xlab("TDP-43 knockdown level") +
  ylab("Normalised PSI") +
  scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++")) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_classic()
print(plot_not_normal)


unc_cryptic  <- big_data %>%
  dplyr::filter(gene_name == "UNC13A" & (.id == "no_dox" | .id == "dox_025")) %>%
  group_by(paste_into_igv_junction, exon_type) %>%
  filter(n() == 2) %>%
  #summarise(mean_psi_per_lsv_junction)
  dplyr::filter(mean_psi_per_lsv_junction[.id == "no_dox"] < 0.05 &
                  mean_psi_per_lsv_junction[.id == "dox_025"] > 0.1) %>%
  mutate(delta_psi = mean_psi_per_lsv_junction[.id == "dox_025"] -mean_psi_per_lsv_junction[.id == "no_dox"]) %>%
  dplyr::filter(delta_psi > 0.05 & .id == "dox_025")

unc_cryptic_df <- unc_cryptic %>%
  dplyr::select(c(2,5,10,11,13,14,19)) %>%
  mutate(type = "cryptic_exon")
names(unc_cryptic_df)[c(5,7,8)] <- c("seqnames", "start", "end")

write.table(unc_cryptic_df, "unc_conc.csv", quote = F, row.names = F, sep = ",")
