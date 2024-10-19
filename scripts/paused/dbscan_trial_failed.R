library(dbscan)
library(umap)
library(tidyverse)

#big_delta_filter <- big_delta %>%
#  dplyr::filter(mean_dpsi_per_lsv_junction > 0.1 & probability_changing > 0.9)

manual_validation_sy5y <- readxl::read_xlsx("~/Documents/phd/research_lines/tdp-43 curves/manual_validation_curves/validation_curves.xlsx") %>%
  mutate(paste_into_igv_junction = paste0(chr, ":", start, "-", end)) %>%
  filter(validation == "yes" & cryptic == "yes")

big_data_filtereds <- big_data %>%
  filter(paste_into_igv_junction %in% manual_validation_sy5y$paste_into_igv_junction) %>%
  #filter(lsv_junc %in% big_delta_filter$lsv_junc) %>% #try without cryptics
  group_by(lsv_junc) %>% 
  filter(n() == 6) %>%
  filter(mean_psi_per_lsv_junction[.id == "no_dox"] < 0.05) %>% 
  filter(mean_psi_per_lsv_junction[.id == "dox_075"] > 0.15) %>%
  filter((mean_psi_per_lsv_junction[.id == "dox_075"]  - mean_psi_per_lsv_junction[.id == "dox_025"]  > - 0.05) &
           (mean_psi_per_lsv_junction[.id == "dox_025"]  - mean_psi_per_lsv_junction[.id == "dox_021"]  > - 0.05) &
           (mean_psi_per_lsv_junction[.id == "dox_021"]  - mean_psi_per_lsv_junction[.id == "dox_0187"] > - 0.05) &
           (mean_psi_per_lsv_junction[.id == "dox_0187"] - mean_psi_per_lsv_junction[.id == "dox_0125"] > - 0.05) &
           (mean_psi_per_lsv_junction[.id == "dox_0125"] - mean_psi_per_lsv_junction[.id == "no_dox"]   > - 0.05)) %>%
  dplyr::mutate(mean_psi_per_lsv_junction_normal = (mean_psi_per_lsv_junction - mean_psi_per_lsv_junction[.id == "no_dox"]) 
                / (mean_psi_per_lsv_junction[.id == "dox_075"] - mean_psi_per_lsv_junction[.id == "no_dox"]))

big_data_dbscan <- big_data_filtereds %>%
  ungroup() %>%
  dplyr::select(1,2,5,17) %>%
  #unique(c(.id,exon_type))
  #dplyr::mutate(id = paste0(.id,";",exon_type)) #%>%
  tidyr::pivot_wider(names_from = .id, values_from = mean_psi_per_lsv_junction, values_fn = max) %>%
  mutate(category = as.factor(as.character(ifelse((dox_0125-no_dox > 0.1 | dox_0187-no_dox > 0.1) & 
                                                           dox_075-no_dox > 0.15, "early", 
                                                         ifelse(dox_025-no_dox < 0.1 & 
                                                                  dox_075-no_dox > 0.15, "late",
                                                                ifelse((dox_021-no_dox > 0.1 | dox_025-no_dox > 0.1) & 
                                                                         dox_075-no_dox > 0.15, "intermediate",
                                                                       "none"))))))
  #dplyr::filter(no_dox < 0.05 & dox_075 > 0.1)
table(big_data_dbscan$category)
names(big_data_dbscan)[3] <- "dox_0"
big_data_dbscan %>%
  pivot_longer(cols = c(3:8), names_to = "level", values_to = "value") %>%
  #mutate(level = factor(level), levels = c("no_dox", "dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075")) %>%
  ggplot(aes(x = level, y = value, group = paste_into_igv_junction)) +
  geom_line(aes(color = category)) +
  theme_minimal()
  

#write.table(big_data_dbscan, "~/Desktop/data_for_ettore_normal_label.csv", quote = F, row.names = F, sep = ",")

penguins <- big_data_dbscan %>% 
  drop_na() %>%
  #dplyr::filter(color_gene_name != "none") %>%
  dplyr::mutate(ID=dplyr::row_number())

penguins_meta <- penguins %>%
  dplyr::select(ID, gene_name, category, no_dox, dox_075, paste_into_igv_junction) %>%
  mutate(delta_psi = dox_075-no_dox)

umap_fitt <- penguins %>%
  dplyr::select(where(is.numeric)) %>%
  tibble::column_to_rownames("ID") #%>%
  #scale()
write.table(umap_fitt, "~/Desktop/umap_fitt.csv", sep = ",", row.names = F, quote = F)

set.seed(420)
umap_fitt
umap_fit <- umap(umap_fitt[,c(1:6)], n_neighbors = 10, min_dist = 0.01)
#umap_fit <- dbscan(umap_fitt[,c(2:5)], eps = 0.05, minPts = 5)
#umap_fit$layout

umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  dplyr::rename(UMAP1="V1",
                UMAP2="V2") %>%
  dplyr::mutate(ID= dplyr::row_number())%>%
  dplyr::inner_join(penguins_meta, by="ID") #%>%
  #dplyr::filter(color_gene_name != "none")

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
             color = category,#)) +
             size = delta_psi)) +
  geom_point() +
  #geom_text_repel(aes(label = gene_name),
  #                max.overlaps = Inf,
  #                show.legend = F) +
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot") +
  scale_color_brewer(palette = "Set1") +
  theme_bw()





manual_validation_dz <- read_xlsx("~/Documents/phd/research_lines/tdp-43 curves/manual_validation_curves/dz_to_validate.xlsx") %>%
  mutate(paste_into_igv_junction = paste0(chr, ":", start, "-", end)) %>%
  filter(validation == "yes" & cryptic == "yes")

big_data_dz_filtereds <- big_data_dz %>% #1901543/6 = 316923.8
  dplyr::filter(paste_into_igv_junction %in% manual_validation_dz$paste_into_igv_junction) %>%
  group_by(lsv_junc) %>% 
  filter(n() == 8) %>%
  filter(mean_psi_per_lsv_junction[.id == "DZ_curves_0"] < 0.05) %>% 
  filter(mean_psi_per_lsv_junction[.id == "DZ_curves_1"] > 0.15) %>%
  filter((mean_psi_per_lsv_junction[.id == "DZ_curves_1"]  - mean_psi_per_lsv_junction[.id == "DZ_curves_01"]  > - 0.05) &
           (mean_psi_per_lsv_junction[.id == "DZ_curves_01"]  - mean_psi_per_lsv_junction[.id == "DZ_curves_0075"]  > - 0.05) &
           (mean_psi_per_lsv_junction[.id == "DZ_curves_0075"]  - mean_psi_per_lsv_junction[.id == "DZ_curves_005"] > - 0.05) &
           (mean_psi_per_lsv_junction[.id == "DZ_curves_005"] - mean_psi_per_lsv_junction[.id == "DZ_curves_004"] > - 0.05) &
           (mean_psi_per_lsv_junction[.id == "DZ_curves_004"] - mean_psi_per_lsv_junction[.id == "DZ_curves_003"]   > - 0.05) &
           (mean_psi_per_lsv_junction[.id == "DZ_curves_003"] - mean_psi_per_lsv_junction[.id == "DZ_curves_002"]   > - 0.05) &
           (mean_psi_per_lsv_junction[.id == "DZ_curves_002"] - mean_psi_per_lsv_junction[.id == "DZ_curves_0"]   > - 0.05)) %>%
  dplyr::mutate(mean_psi_per_lsv_junction_normal = (mean_psi_per_lsv_junction - mean_psi_per_lsv_junction[.id == "DZ_curves_0"]) 
                / (mean_psi_per_lsv_junction[.id == "DZ_curves_1"] - mean_psi_per_lsv_junction[.id == "DZ_curves_0"]))

big_data_dz_filtereds_early_late <- big_data_dz_filtereds %>%
  ungroup() %>%
  dplyr::select(1,2,5,17) %>%
  tidyr::pivot_wider(names_from = .id, values_from = mean_psi_per_lsv_junction, values_fn = max) %>%
  mutate(category = as.factor(as.character(ifelse((DZ_curves_002-DZ_curves_0 > 0.1 | DZ_curves_003-DZ_curves_0 > 0.1) & 
                                                           DZ_curves_1-DZ_curves_0 > 0.15, "early", 
                                                         ifelse(DZ_curves_004-DZ_curves_0 < 0.1 & 
                                                                  (DZ_curves_0075-DZ_curves_0 > 0.1 | 
                                                                     DZ_curves_01-DZ_curves_0 > 0.1 | DZ_curves_1-DZ_curves_0 > 0.15), "late",
                                                                ifelse((DZ_curves_005-DZ_curves_0 > 0.1 | DZ_curves_004-DZ_curves_0 > 0.1) &
                                                                         DZ_curves_1-DZ_curves_0 > 0.15, "intermediate",
                                                                       "none")))))) #%>%
  #pivot_longer(cols = c("DZ_curves_0","DZ_curves_002", "DZ_curves_003", "DZ_curves_004", "DZ_curves_005", "DZ_curves_0075", "DZ_curves_01", "DZ_curves_1")) %>%
  #mutate(name = factor(name, levels = c("DZ_curves_0","DZ_curves_002", "DZ_curves_003", "DZ_curves_004", "DZ_curves_005", "DZ_curves_0075", "DZ_curves_01", "DZ_curves_1")))
table(big_data_dz_filtereds_early_late$category)

penguins <- big_data_dz_filtereds_early_late %>% 
  drop_na() %>%
  #dplyr::filter(color_gene_name != "none") %>%
  dplyr::mutate(ID=dplyr::row_number())

penguins_meta <- penguins %>%
  dplyr::select(ID, gene_name, category, DZ_curves_0, DZ_curves_1, paste_into_igv_junction) %>%
  mutate(delta_psi = DZ_curves_1-DZ_curves_0)

umap_fitt <- penguins %>%
  dplyr::select(where(is.numeric)) %>%
  tibble::column_to_rownames("ID") #%>%
#scale()

umap_fitt
umap_fit <- umap(umap_fitt[,c(1:8)], n_neighbors = 10, min_dist = 0.01)
#umap_fit <- dbscan(umap_fitt[,c(2:5)], eps = 0.05, minPts = 5)
#umap_fit$layout

umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  dplyr::rename(UMAP1="V1",
                UMAP2="V2") %>%
  dplyr::mutate(ID= dplyr::row_number())%>%
  dplyr::inner_join(penguins_meta, by="ID") #%>%
#dplyr::filter(color_gene_name != "none")

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
             color = category,#)) +
             size = delta_psi)) +
  geom_point() +
  geom_text_repel(aes(label = gene_name),
                  max.overlaps = Inf,
                  show.legend = F) +
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot") +
  scale_color_brewer(palette = "Set1") +
  theme_bw()


####################################################################














#PCAtools::pca(t(umap_fitt)) %>% PCAtools::biplot()


#big_data_dbscan_clean <- big_data_dbscan %>%
#  tidyr::pivot_wider(names_from = .id, values_from = mean_psi_per_lsv_junction, values_fn = max)
#dplyr::distinct() %>%
#  dplyr::group_by(paste_into_igv_junction) %>%
#  dplyr::mutate(row = dplyr::row_number()) %>%
#  tidyr::pivot_wider(names_from = paste_into_igv_junction, values_from = mean_psi_per_lsv_junction) %>%
#  dplyr::select(-row)
#pivot_wider(., names_from = .id, values_from = mean_psi_per_lsv_junction) #, values_fn = mean(mean_psi_per_lsv_junction))

#which(duplicated(big_data_dbscan))

#duplicates_big <- big_data %>%
#  dplyr::group_by(paste_into_igv_junction, .id) %>%
#  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
#  dplyr::filter(n == 1)

#big_data_dbscan_clean <- big_data_dbscan %>%
#  tibble::column_to_rownames("paste_into_igv_junction") 
#  drop_na() %>%

#matrix_dbscan <- as.matrix(big_data_dbscan_clean) 
#plot(c(matrix_dbscan), pch='.')


big_data_dbscan_clean[is.na(big_data_dbscan_clean),]
c1 <- big_data_dbscan_clean[,c(2:5)] %>%
  dbscan(eps = 0.3, minPts = 5)
c1
#table(big_data_dbscan_long$names, c1$cluster)
plot(big_data_dbscan_clean, col=c1$cluster+1, pch='.')
#colors <- mapply(function(col, i) adjustcolor(col, alpha.f = cl$membership_prob[i]), 
#                 palette()[cl$cluster+1], seq_along(cl$cluster))
#points(big_data_dbscan_clean, col=colors, pch=20)

ggplot(big_data_dbscan_long, aes(x = names, y = values, group = paste_into_igv_junction)) +
  geom_line()


c2 <- big_data_dbscan_clean %>%
  dbscan::optics(minPts = 5)
c3 <- c2 %>% extractDBSCAN(eps_cl = 0.1)
c3
plot(big_data_dbscan_clean, col=c3$cluster+1, pch=20)
#plot(c1, big_data_dbscan_clean)


big_data_dbscan_long <- pivot_longer(big_data_dbscan_clean, cols = c("no_dox", "dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075"), names_to = "names", values_to = "values")
big_data_dbscan_long$names <- factor(big_data_dbscan_long$names, levels = c("no_dox", "dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075"))

###labels as kd level?



penguins <- big_data_dbscan %>% 
  drop_na() %>%
  #select(-year)%>%
  dplyr::mutate(ID=dplyr::row_number())

penguins_meta <- penguins %>%
  dplyr::select(ID, paste_into_igv_junction)

set.seed(142)
umap_fitt <- penguins %>%
  dplyr::select(where(is.numeric)) %>%
  tibble::column_to_rownames("ID") %>%
  scale()
#%>% 
umap_fit <- umap(umap_fitt, n_neighbors = 10, min_dist = 0.25)
umap_fit$config

umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  dplyr::rename(UMAP1="V1",
         UMAP2="V2") %>%
  dplyr::mutate(ID= dplyr::row_number())%>%
  dplyr::inner_join(penguins_meta, by="ID")

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2))+
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot")





#umap.defaults
big_data_umap_clean <- big_data_dbscan_clean %>% scale()
umap_c1 <- umap(big_data_umap_clean, n_neighbors = 15, min_dist = 0.1)
umap_c1$config

plot(umap_c1$layout)
dummy <- as.data.frame(umap_c1$layout) #%>%
dummy %>%
  ggplot(aes(x= V1, y = V2)) +
  geom_point() + theme_classic()
