##load targeted seq table and compute CE ratio / rename samples in agreement to UNC13A paper
maxc <- read.csv("~/Desktop/ce_bio/count_table_of_aligned_remove_200bp_plus_product_reads.csv") %>%
  mutate(fraction_cryptic = n_cryptic / n_counts_for_gene) %>%
  mutate(sample_name_ari = case_when(sample_name == "P64_11" ~ "HC1",
                                     sample_name == "P17_07" ~ "HC2",
                                     sample_name == "P47_11" ~ "HC3",
                                     sample_name == "P35_07" ~ "HC4",
                                     sample_name == "P56_13" ~ "FTD1",
                                     sample_name == "P45_15" ~ "FTD2",
                                     sample_name == "P63_05" ~ "FTD3",
                                     sample_name == "P28_07" ~ "FTD4",
                                     sample_name == "P40_04" ~ "FTD5",
                                     sample_name == "P07_15" ~ "FTD6",
                                     sample_name == "P16_09" ~ "FTD7",
                                     sample_name == "P86_08" ~ "FTD8",
                                     sample_name == "P13_13" ~ "FTD9",
                                     sample_name == "P11_07" ~ "FTE10",
                                     sample_name == "ddH2O" ~ "HO2O")) #%>%
  #pivot_wider(names_from = sample_name_ari, values_from = ratio_cryptic)
table(maxc$sample_name_ari)

maxc %>%
  #filter(gene_name == "UNC13A") %>% ###if to compare to unc13a paper
  pivot_longer(cols = c(n_counts_for_gene, n_cryptic)) %>%#+
  ggplot(aes(x = sample_name_ari, y = value, fill = name)) +  #, fill = gene_name)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_color_manual(values = "black") +
  facet_grid(facets = vars(name), scale = 'free_y') +
  theme_bw()

##filter targeted table (at least 5 counts) - DON'T THINK IS NEEDED
#maxc_filtered <- maxc %>%
#  filter(n_counts_for_gene > 5 & ratio_cryptic < 0.5) %>%
#  group_by(gene_name) %>% 
#  filter(n() > 10) #%>%
#resume2 <- as.data.frame(table(maxc_filtered$gene_name))
#table(maxc_filtered$sample_name_ari)

##use only TDP43path-specific genes (use list from AL) - this is exluding more than 2/3 of genes!!
al <- read.csv("~/Desktop/ce_bio/expression_by_pathology_updated.csv") %>%
  dplyr::filter(fraction_path > 0.01 & fraction_not_path < 0.005) %>% ####need to carefully check this with AL
  distinct(gene, .keep_all = T) ##check when binding with annotated validation table - should be ok since left_join
names(al)[4] <- "gene_name"
maxc_al <- maxc %>%
  left_join(al, by = "gene_name") %>%
  dplyr::filter(!is.na(paste_into_igv_junction))
table(maxc_al$sample_name_ari)
table(maxc_al$gene_name)

maxc_al %>%
  dplyr::filter(gene_name != "SETD5") %>%
  pivot_longer(cols = c(n_counts_for_gene, n_cryptic)) %>%#+
  #filter(gene_name != "SETD5") %>% #|
  #filter(gene_name == "UNC13A") %>% # |
  #         gene_name == "AARS") %>%
  ggplot(aes(x = sample_name_ari, y = value, fill = name)) +  #, fill = gene_name)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_color_manual(values = "black") +
  facet_grid(facets = vars(name), scale = 'free_y') +
  theme_bw()

##annotate manual_validation and bind to previous df
names(manual_validation_sy5y)[2] <- "gene_name"
namual <- manual_validation_sy5y %>%
  mutate(progression = ifelse(dox0125 > 0.2, "A-Early",
                              ifelse(dox0187 > 0.2, "Early",
                       ifelse(dox021 > 0.2 | dox025 > 0.2, "Intermediate",
                       "Late")))) %>%
  filter(is_cryptic == 1)
table(namual$progression)
#CELF5 - problem w paste_into_igv_junc

#use dz for these 2 genes - CEP290, SYNJ2
names(manual_validation_dz)[2] <- "gene_name"
namual_dz <- manual_validation_dz %>% #??????????? check wb levels
  mutate(progression_dz = ifelse(dox002 > 0.2 | dox003 > 0.2, "Early",
                       ifelse(dox004 > 0.2 | dox005 > 0.2, "Intermediate",
                              "Late"))) %>%
  filter(is_cryptic == 1)
table(namual_dz$progression_dz)
#-CEP290 - INTERMEDIATE??
#-SYNJ2 - INTERMEDIATE??

missing_genes <- c("ADGRL1", "ATP8A2", "ATXN1", "CAMK2B", "DLGAP1", "IGSF21", "KALRN", "LINGO1") #found KALRN and ADGRL1
al_missing <- al %>% filter(gene_name %in% missing_genes)
big_data_missing <- big_data %>% filter(paste_into_igv_junction %in% al_missing$paste_into_igv_junction)
big_delta_missing <- big_delta %>% filter(paste_into_igv_junction %in% al_missing$paste_into_igv_junction)


##check igv files for these genes
maxc_type <- maxc_al %>%
  select(-type) %>%
  mutate(gene_id = gene_name) %>%
  filter(gene_name != "SETD5") %>%
  left_join(umap_df2, by = "gene_id") %>%

  #left_join(namual) %>% #left_join(umap_df2) #%>%. #left_join(umap_df2_dz)
  mutate(progression = ifelse(gene_name == "CELF5" | gene_name == "KALRN", "Intermediate", #-CELF5 - problem w paste_into_igv_junc - intermediate?
                              ifelse(gene_name == "CEP290" | gene_name == "SYNJ2", "Early", #from dz curves - potentially early ?
                                     ifelse(gene_name == "ADGRL1", "Late", #from big data (as KALRN)
                                            ifelse(gene_name %in% c("ATP8A2", "ATXN1", "CAMK2B", "DLGAP1", "IGSF21", "LINGO1", "CEP290", "SYNJ2"), "NA",
                                     progression)))) %>%
  filter(!is.na(progression)) %>% 
  mutate(group = ifelse(grepl("FT", sample_name_ari), "FTD",
                        ifelse(grepl("HC", sample_name_ari), "HC",
                               "H2O"))) %>% 
  filter(cluster != "NA") %>%
  group_by(gene_name) %>%
  filter(sum(n_cryptic) > 0) %>%
  ungroup()
  #group_by(gene_name) %>%
  #filter(sum(n_cryptic) > 0) %>%
  #ungroup() %>%
  #group_by(gene_name, sample_name_ari) %>%
  #mutate(detected = ifelse(n_cryptic <= 4, "No", "Yes"))
table(maxc_type$gene_name)
table(maxc_type$sample_name_ari)
View(as.data.frame(table(maxc_type$gene_name, maxc_type$detected, maxc_type$group, maxc_type$progression)))

table(maxc_type$cluster, maxc_type$gene_name)
maxc_type %>%
  filter(gene_name != "PRELID3A" & gene_name != "IQCE") %>%
  mutate(gene_name = fct_reorder(gene_name, as.numeric(cluster))) %>%
  ggplot(aes(x = sample_name_ari, y = fraction_cryptic)) +
  geom_bar(aes(fill = gene_name), stat = "identity", position = "stack", color="black") +#position = position_dodge(width = 0.8)) +
  facet_grid(facets = vars(cluster)) +
  scale_fill_manual(values =colorspace::diverge_hcl(n = 15)) +
  #scale_fill_viridis_d() +
  #scale_fill_brewer(palette = "RdYlGn") +
 # scale_y_continuous(limits = c(0,0.3)) +
  #labs(x = "", y = "", fill = "Early") +
  #scale_fill_manual(values = randomcoloR::distinctColorPalette(k = 30)) +
  theme_bw()

write.csv(maxc_type, "~/Desktop/panel_annotated.csv")

maxc_type %>%
  filter(progression != "NA") %>%
  filter(gene_name != "SETD5") %>%
  filter(detected == "Yes") %>%
  filter(group != "H2O") %>%
  ggplot(aes(x = progression)) +
  geom_bar(aes(fill = progression)) +
  facet_wrap(facets = vars(group), ncol = 3) +
  scale_fill_brewer(palette = "Set1") +
  xlab("") +
  ylab("Count of cryptic exons detected by targeted panel") +
  labs(fill = "Category") +
  theme_bw()

maxc_early <- maxc_type %>%
  filter(progression == "Early") %>%
  filter(gene_name != "SETD5") %>% #|
  ggplot(aes(x = sample_name_ari, y = fraction_cryptic)) +
  geom_bar(aes(fill = gene_name), stat = "identity", position = "stack", color="black") +#position = position_dodge(width = 0.8)) +
  #facet_grid(facets = vars(gene_name)) +
  scale_fill_brewer(palette = "Reds") +
  scale_y_continuous(limits = c(0,0.3)) +
  labs(x = "", y = "", fill = "Early") +
  #scale_fill_manual(values = randomcoloR::distinctColorPalette(k = 30)) +
  theme_bw()
maxc_intermediate <- maxc_type %>%
  filter(progression == "Intermediate") %>%
  filter(gene_name != "SETD5") %>% #|
  ggplot(aes(x = sample_name_ari, y = fraction_cryptic)) +
  geom_bar(aes(fill = gene_name), stat = "identity", position = "stack", color="black") +#position = position_dodge(width = 0.8)) +
  #facet_grid(facets = vars(gene_name)) +
  scale_fill_brewer(palette = "Greens") +
  scale_y_continuous(limits = c(0,0.3)) +
  labs(x = "", y = "", fill = "Intermediate") +
  #scale_fill_manual(values = randomcoloR::distinctColorPalette(k = 30)) +
  theme_bw()
maxc_late <- maxc_type %>%
  filter(progression == "Late") %>%
  filter(gene_name != "SETD5") %>% #|
  ggplot(aes(x = sample_name_ari, y = fraction_cryptic)) +
  geom_bar(aes(fill = gene_name), stat = "identity", position = "stack", color="black") +#position = position_dodge(width = 0.8)) +
  #facet_grid(facets = vars(gene_name)) +
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(limits = c(0,0.3)) +
  labs(x = "", y = "", fill = "Late") +
  #scale_fill_manual(values = randomcoloR::distinctColorPalette(k = 30)) +
  theme_bw()
maxc_na <- maxc_type %>%
  filter(progression == "NA") %>%
  filter(gene_name != "SETD5") %>% #|
  ggplot(aes(x = sample_name_ari, y = fraction_cryptic)) +
  geom_bar(aes(fill = gene_name), stat = "identity", position = "stack", color="black") +#position = position_dodge(width = 0.8)) +
  #facet_grid(facets = vars(gene_name)) +
  scale_fill_grey() +
  scale_y_continuous(limits = c(0,0.3)) +
  labs(x = "", y = "", fill = "NA") +
  #scale_fill_manual(values = randomcoloR::distinctColorPalette(k = 30)) +
  theme_bw()
maxc_early / maxc_intermediate / maxc_late / maxc_na



maxc_type$progression <- as.factor(maxc_type$progression)
maxc_type <- arrange(maxc_type, progression)
maxc_type$gene_name <- factor(maxc_type$gene_name, levels = unique(maxc_type$gene_name))
maxc_type %>%
  filter(progression != "NA") %>%
  filter(gene_name != "SETD5") %>% #|
  ggplot(aes(x = reorder(gene_name, progression), y = fraction_cryptic)) +
  geom_bar(aes(fill = progression), stat = "identity", position = "stack", color="black") +#position = position_dodge(width = 0.8)) +
  facet_grid(facets = vars(sample_name_ari)) +
  #scale_fill_grey() +
  #scale_y_continuous(limits = c(0,0.3)) +
  #scale_fill_manual(values = randomcoloR::distinctColorPalette(k = 30)) +
  theme_bw()
reorder(maxc_type$gene_name, maxc_type$progression)
unique(maxc_type$gene_name)

maxc_type %>%
  filter(progression != "NA") %>%
  filter(gene_name != "SETD5") %>% #|
  ggplot(aes(x = sample_name_ari, y = fraction_cryptic)) +
  geom_col(aes(fill = gene_name), position = position_dodge(width = 0.8)) +
  facet_grid(facets = vars(progression)) +
  #scale_fill_brewer(palette = "Set1") +
  #scale_fill_manual(values = randomcoloR::distinctColorPalette(k = 30)) +
  theme_bw()

maxc_type %>%
  filter(progression != "NA") %>%
  filter(gene_name != "SETD5") %>%
  filter(group != "H2O") %>%
  ggplot(aes(x = group, y = fraction_cryptic)) +
  geom_boxplot(aes(fill = progression)) +
  facet_wrap(facets = vars(progression), ncol = 3) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw()

maxc_type %>%
  filter(progression != "NA") %>%
  filter(gene_name != "SETD5") %>%
  #filter(group != "H2O") %>%
  ggplot(aes(x = gene_name, y = fraction_cryptic)) +
  geom_col(aes(fill = progression)) +
  facet_wrap(facets = vars(sample_name_ari), ncol = 5) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

maxc_type %>%
  filter(progression != "NA") %>%
  filter(gene_name != "SETD5") %>%
  #filter(group != "H2O") %>%
  ggplot(aes(x = progression, y = fraction_cryptic)) +
  geom_boxplot(aes(fill = progression)) +
  facet_wrap(facets = vars(sample_name_ari), ncol = 5) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))


maxc_type %>%
  filter(group != "H2O") %>%
  ggplot(aes(x = sample_name_ari, y = gene_name, fill = fraction_cryptic)) +
  geom_tile() +
  #scale_fill_viridis_d() +
  scale_fill_gradient(low = "white", high = "black") +
  facet_grid(rows = vars(progression), cols = vars(group), 
             scales = "free", space = "free") +
  ylab("") + xlab("Sample name") + labs(fill = "CE fraction") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))



example_file <- "https://davetang.org/file/TagSeqExample.tab"
data <- read.delim(example_file, header = TRUE, row.names = "gene")
data_subset <- as.matrix(data[rowSums(data)>50000,])

maxc_matrix <- maxc_type %>%
  filter(progression != "NA") %>%
  filter(gene_name != "SETD5") %>%
  select(sample_name_ari, ratio_cryptic, gene_name) %>%
  arrange(sample_name_ari) %>%
  pivot_wider(names_from = sample_name_ari, values_from = ratio_cryptic) %>%
  mutate(across(where(anyNA), ~ replace_na(., 0))) %>%
  column_to_rownames("gene_name") %>%
  as.matrix() #
pheatmap::pheatmap(maxc_matrix)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(maxc_matrix, 1, cal_z_score))
pheatmap::pheatmap(data_subset_norm)

maxc_type %>%
  filter(grepl("FT", sample_name_ari)) %>%
  #filter(progression == "Early") %>%
  filter(gene_name != "SETD5") %>% #|
  ggplot(aes(x = sample_name_ari, y = ratio_cryptic)) +
  geom_bar(aes(fill = progression), stat = "identity", position = "stack", color="black") +#position = position_dodge(width = 0.8)) +
  facet_grid(facets = vars(progression)) +
  scale_fill_brewer(palette = "Set1") +
  #scale_fill_manual(values = randomcoloR::distinctColorPalette(k = 30)) +
  theme_bw()


maxc_type %>%
  #filter(progression != "NA") %>%
  filter(gene_name != "SETD5") %>% #|
  #filter(gene_name == "UNC13A") %>% # |
  #         gene_name == "AARS") %>%
  ggplot(aes(x = sample_name_ari, y = pc_cryptic, fill = gene_name, group = progression, color = "black")) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_color_manual(values = "black") +
  theme_bw()


panel_rpm_kd <- read.csv("~/Downloads/cryptic_rpm_per.csv")
panel_rpm_ctrl <- read.csv("~/Downloads/cryptic_rpm_ctrls.csv")
panel_list_broad <- cbind(panel_rpm_kd, panel_rpm_ctrl[,c(2:8)]) %>%
  mutate(total_kd = rowSums(panel_list_broad[,c(2:8)]),
         total_ctrl = rowSums(panel_list_broad[,c(9:15)])) %>%
  mutate(delta_tot = total_kd - total_ctrl) %>%
  arrange(desc(delta_tot)) %>%
  mutate(rank = 1:nrow(panel_list_broad))



panel_rpm_kd <- read.csv("~/Documents/phd/research_lines/ce_panel/cryptic_rpm_per.csv")
panel_rpm_ctrl <- read.csv("~/Documents/phd/research_lines/ce_panel/cryptic_rpm_ctrls.csv")
panel_list_broad <- cbind(panel_rpm_kd, panel_rpm_ctrl[,c(2:8)]) %>%
  mutate(total_kd = rowSums(panel_list_broad[,c(2:8)]),
         total_ctrl = rowSums(panel_list_broad[,c(9:15)])) %>%
  mutate(delta_tot = total_kd - total_ctrl) %>%
  arrange(desc(delta_tot)) %>%
  mutate(rank_i3n = 1:nrow(panel_list_broad)) %>%
  mutate(gene = gsub("_.*$", "", gene_name))

#panel_patient_specific <- read.csv("~/Documents/phd/research_lines/ce_panel/")
panel21 <- read.csv("~/Documents/phd/research_lines/ce_panel/CEs_label_Ver_2.1.csv")

table_cell_line <- read.csv("~/Desktop/cell_line_dataset.csv")
splicing_full_2 <- splicing_full %>%
  left_join(table_cell_line)

splicing_full_datasets_start <- splicing_full_2 %>% 
  arrange(desc(mean_dpsi_per_lsv_junction)) %>%
  group_by(paste_into_igv_junction, comparison) %>%
  distinct(.keep_all = T) %>% #16413 #17466 
  ungroup() %>% 
  group_by(paste_into_igv_junction) %>%
  mutate(datasets = n_distinct(comparison)) %>%
  # Add a column with the list of "comparison"
  mutate(list_comparison = list(unique(comparison))) %>%
  # Convert list column to character column with comma-separated values
  mutate(list_comparison = sapply(list_comparison, paste, collapse = ", ")) %>%
  mutate(cell_lines = n_distinct(cell_line)) %>%
  mutate(list_cell = list(unique(cell_line))) %>%
  mutate(list_cell = sapply(list_cell, paste, collapse = ", ")) %>%
  ungroup() %>%
  distinct(paste_into_igv_junction, .keep_all = T) %>%
  select(c(gene_name, paste_into_igv_junction, baseline_PSI, contrast_PSI, mean_dpsi_per_lsv_junction, probability_changing, junc_cat, cell_lines, list_cell, datasets, list_comparison)) %>%
  separate(paste_into_igv_junction, sep = ":", into = c("Seqname", "range")) %>% 
  separate(range, sep = "-", into = c("start", "end")) %>%
  mutate(#end = as.numeric(start),
    start = as.numeric(end)-1) 
splicing_full_datasets_start[35,3] <- 36573043

panel21_datasets_start <- panel21 %>% 
  left_join(splicing_full_datasets_start, by = c("Seqname", "start")) %>%
  distinct(gene, .keep_all = T) %>%
  filter(!is.na(gene_name)) %>% 
  rename(end = end.x) %>%
  select(-(end.y))

splicing_full_datasets_end <- splicing_full_2 %>% 
  arrange(desc(mean_dpsi_per_lsv_junction)) %>%
  group_by(paste_into_igv_junction, comparison) %>%
  distinct(.keep_all = T) %>% #16413 #17466 
  ungroup() %>% 
  group_by(paste_into_igv_junction) %>%
  mutate(datasets = n_distinct(comparison)) %>%
  # Add a column with the list of "comparison"
  mutate(list_comparison = list(unique(comparison))) %>%
  # Convert list column to character column with comma-separated values
  mutate(list_comparison = sapply(list_comparison, paste, collapse = ", ")) %>%
  mutate(cell_lines = n_distinct(cell_line)) %>%
  mutate(list_cell = list(unique(cell_line))) %>%
  mutate(list_cell = sapply(list_cell, paste, collapse = ", ")) %>%
  ungroup() %>%
  distinct(paste_into_igv_junction, .keep_all = T) %>%
  select(c(gene_name, paste_into_igv_junction, baseline_PSI, contrast_PSI, mean_dpsi_per_lsv_junction, probability_changing, junc_cat, cell_lines, list_cell, datasets, list_comparison)) %>%
  separate(paste_into_igv_junction, sep = ":", into = c("Seqname", "range")) %>% 
  separate(range, sep = "-", into = c("start", "end")) %>%
  mutate(end = as.numeric(start))
#start = as.numeric(end)-1)

panel21_datasets_end <- panel21 %>%
  left_join(splicing_full_datasets_end, by = c("Seqname", "end")) %>%
  distinct(gene, .keep_all = T) %>%
  filter(!is.na(gene_name)) %>%
  rename(start = start.x) %>%
  select(-(start.y))

panel21_datasets <- rbind(panel21_datasets_start, panel21_datasets_end) %>% 
  arrange(desc(datasets)) %>% 
  distinct(gene, .keep_all = T) %>% 
  left_join(panel_list_broad, by = "gene") %>%
  arrange(start) %>%
  arrange(Seqname)
  #mutate(chr = gsub("chr", "", Seqname)) %>%
  #arrange(chr) #%>%
  #mutate(Seqname = paste0("chr", chr))
panel21_datasets_check <- left_join(panel21_datasets, panel21) 

write.csv(panel21_datasets, "~/Documents/phd/research_lines/ce_panel/CE_panel_v2.1_annotated.csv", row.names = F)
write.table(panel21_datasets[,c(1:3)], "~/Desktop/CE_panel_v2_1.bed", sep = "\t", quote = F, row.names = F)

panel31 <- read.csv("~/Downloads/CEs_label_Ver_3.1.csv")

panel_change <- panel21_datasets %>% 
  left_join(panel31)

genes_nih <- readxl::read_xlsx("~/Desktop/Copy of CE for probe design_SS.xlsx")

ciaone <- genes_nih %>% left_join(panel21)
