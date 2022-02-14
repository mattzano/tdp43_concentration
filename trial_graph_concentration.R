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
#upset(big_data, nsets = 5, number.angles = 30, point.size = 3.5, line.size = 2)



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
  filter(mean_dpsi_per_lsv_junction > 0 &
           probability_changing > 0.9)

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

big_data_filtered <- big_data %>%
  filter(lsv_junc %in% big_delta_filter$lsv_junc) %>% #try without cryptics
  group_by(lsv_junc) %>%
  filter(n() == 6) %>%
  mutate(mean_psi_per_lsv_junction_normal = (mean_psi_per_lsv_junction - mean_psi_per_lsv_junction[.id == "no_dox"])
         / (mean_psi_per_lsv_junction[.id == "dox_075"] - mean_psi_per_lsv_junction[.id == "no_dox"]))

big_data_filtere <- big_data_filtered %>%
  group_by(lsv_junc) %>%
  filter(n() == 6) %>%
  filter(mean_psi_per_lsv_junction_normal[.id == "dox_075"] >= mean_psi_per_lsv_junction_normal[.id == "dox_025"] &
           mean_psi_per_lsv_junction_normal[.id == "dox_025"] >= mean_psi_per_lsv_junction_normal[.id == "dox_021"] &
           mean_psi_per_lsv_junction_normal[.id == "dox_021"] >= mean_psi_per_lsv_junction_normal[.id == "dox_0187"] &
           mean_psi_per_lsv_junction_normal[.id == "dox_0187"] >= mean_psi_per_lsv_junction_normal[.id == "dox_0125"]) %>%
  filter(mean_psi_per_lsv_junction_normal >= 0) 

big_data_filter <- big_data_filtere %>%
  group_by(lsv_junc) %>%
  filter(n() == 6) %>%
  filter(mean_psi_per_lsv_junction_normal[.id == "dox_025"] < 0.3 &
           mean_psi_per_lsv_junction_normal[.id == "dox_025"] > 0)

big_data_filtere %>%
  #filter(gene_name == "INSR") %>% 
  ggplot(mapping = aes(x = .id, y = mean_psi_per_lsv_junction_normal, group = lsv_junc, color = gene_name)) +
    geom_line() +
    theme_classic()





big_data_matrix <- big_data %>%
  filter(lsv_junc %in% big_delta_filter$lsv_junc) %>%
  mutate(no_dox = ifelse(.id %in% "no_dox", 1, 0),
         dox_0125 = ifelse(.id %in% "dox_0125", 1, 0),
         dox_0187 = ifelse(.id %in% "dox_0187", 1, 0),
         dox_021 = ifelse(.id %in% "dox_021", 1, 0),
         dox_025 = ifelse(.id %in% "dox_025", 1, 0),
         dox_075 = ifelse(.id %in% "dox_075", 1, 0))

big_data_matrix_2 <- big_data_matrix[,c(2,12,19:24)] %>%
  group_by(gene_name,lsv_junc) %>%
  summarise_at(vars(contains('dox')), ~ max(as.numeric(.), na.rm = TRUE))

big_data_matrix_2[,3:8] <- lapply(big_data_matrix_2[,3:8], as.integer)

big_data_matrix_3 <- as.data.frame(big_data_matrix_2)

upset(big_data_matrix_3, nsets = 6, sets = c("no_dox",
                                             "dox_0125", 
                                             "dox_0187", 
                                             "dox_021", 
                                             "dox_025", 
                                             "dox_075"), order.by = "freq", keep.order = TRUE,
      mainbar.y.label = "Junction Intersections", sets.x.label = "Junction Per Level of TDP-43 KD")


