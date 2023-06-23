#splicing results changing dox00187_2 with dox0021_2
#then make the 2 spaghetti plots

input_list <- c("no_dox", 
                #"dox_0125", 
                "dox_0187", 
                "dox_021")#, "dox_025", "dox_075")
datalistas <- list()

for (i in input_list) {
  input <- paste(i, "parsed", sep = "_")
  input_csv <- paste(input, "csv", sep = ".")
  input_splicing <- fread(file.path(here::here(), "data", "majiq_single", input_csv))
  datalistas[[i]] <- input_splicing
}
big_data <- rbindlist(datalistas, idcol = TRUE)
big_data$.id <- factor(big_data$.id, levels = input_list)

big_delta %>%
  dplyr::filter(probability_changing > 0.9 & no_dox_mean_psi < 0.05 & contrast_mean_psi > 0.1) %>%
  dplyr::filter(junc_cat %in% c("novel_acceptor", "novel_donor", "novel_combo", "novel_exon_skip")) %>%
  ggplot() +
  #geom_point(aes(x = .id, y = mean_dpsi_per_lsv_junction, color = junc_cat), position = "jitter") +
  geom_bar(aes(x = .id, fill = junc_cat, group = junc_cat), position = position_dodge()) +
  theme_classic()