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
  filter(lsv_junc %in% big_delta_filter$lsv_junc) %>% #try without cryptics
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

big_data_filtereds_c <- big_data_filtereds %>%
  group_by(lsv_junc) %>% 
  mutate(color_gene_name = as.factor(as.character(ifelse(paste_into_igv_junction %in% c("chr19:17641556-17642414","chr8:79611214-79616822"), 1,
                                                     ifelse(mean_psi_per_lsv_junction_normal_025 < 0.1, 2, 
                                                         ifelse(mean_psi_per_lsv_junction_normal_0187 > 0.7, 3, 4))))))
                                                         #ifelse(mean_psi_per_lsv_junction_normal[.id == "dox_025"] < 0.1, 3, 4)))))
    #mean_psi_per_lsv_junction_normal[.id == "dox_0187"] > 0.6, "red", "orange"))))
                                                         #ifelse(mean_psi_per_lsv_junction_normal[.id == "dox_025"] < 0.1, "green", "orange")))))
#paste_into_igv_junction %in% c("chr19:17641556-17642414",
#                    "chr8:79611214-79616822"), 1, 
                                                                              #"chr20:63439708-63444659",
                                                                              #"chr10:3099819-3101365",
                                                                              #"chr2:241668985-241670726",
                                                                              #"chr9:128956174-128956350"), 1, 0)))) %>%

  #mutate(label_junction = case_when((paste_into_igv_junction %in% c("chr19:17641556-17642414",
  #                                                                  "chr8:79611214-79616822",
  #                                                                  "chr20:63439708-63444659",
  #                                                                  "chr10:3099819-3101365",
  #                                                                  "chr2:241668985-241670726",
  #                                                                  "chr9:128956174-128956350") &
  #                                     .id =="dox_0187") |
  #                                  (mean_psi_per_lsv_junction_normal[.id == "dox_025"] < 0.1 & .id == "dox_025") 
  #                                  ~ gene_name,
  #                                  T ~ ""))

plot <- big_data_filtereds_b %>%
  filter(paste_into_igv_junction %in% c("chr19:17641556-17642414","chr8:79611214-79616822")) %>%
  ggplot(mapping = aes(x = .id, y = mean_psi_per_lsv_junction_normal, group = lsv_junc)) +
  geom_line(aes(color = color_gene_name, alpha = color_gene_name), show.legend = F) +
  #geom_text_repel(aes(x = .id, y = mean_psi_per_lsv_junction, label = label_junction), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  scale_color_manual(values = c(#"gray",
                                "#3288BD")) +
  scale_alpha_manual(values = c(1)) +
  xlab("TDP-43 knockdown level") +
  ylab("Normalised PSI") +
  scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++")) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_classic()
print(plot)
ggsave("spaghetti_bw_unc_stmn.png", plot, width = 8, height = 5)

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

#shared_cryptic = tbl %>%
#  mutate(name = glue::glue("{gene_name}|{junc_cat}|{paste_into_igv_junction}")) %>%
#  dplyr::select(name,is_a_cryptic,comparison) %>%
#  unique() %>%
#  filter(is_a_cryptic == T) %>% 
#  mutate(event_present = is_a_cryptic) %>%
#  pivot_wider(values_from = "event_present",
#              id_cols = "name",
#              names_from = "comparison",
#              values_fill = FALSE,values_fn = sum) %>%
#  column_to_rownames('name') %>%
#  mutate_all(as.logical)

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

#big_data_matrix <- big_data %>%
#  filter(lsv_junc %in% big_delta_filter$lsv_junc) %>%
#  mutate(no_dox = ifelse(.id %in% "no_dox", 1, 0),
#         dox_0125 = ifelse(.id %in% "dox_0125", 1, 0),
#         dox_0187 = ifelse(.id %in% "dox_0187", 1, 0),
#         dox_021 = ifelse(.id %in% "dox_021", 1, 0),
#         dox_025 = ifelse(.id %in% "dox_025", 1, 0),
#         dox_075 = ifelse(.id %in% "dox_075", 1, 0))

#big_data_matrix_2 <- big_data_matrix[,c(2,12,19:24)] %>%
#  group_by(gene_name,lsv_junc) %>%
#  summarise_at(vars(contains('dox')), ~ max(as.numeric(.), na.rm = TRUE))

#big_data_matrix_2[,3:8] <- lapply(big_data_matrix_2[,3:8], as.integer)

#big_data_matrix_3 <- as.data.frame(big_data_matrix_2)

#upset(big_data_matrix_3, nsets = 6, sets = c("no_dox",
#                                             "dox_0125", 
#                                             "dox_0187", 
#                                             "dox_021", 
#                                             "dox_025", 
#                                             "dox_075"), order.by = "freq", keep.order = TRUE,
#      mainbar.y.label = "Junction Intersections", sets.x.label = "Junction Per Level of TDP-43 KD")
