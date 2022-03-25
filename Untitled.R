my_clean_reader = function(file){
  as.data.table(janitor::clean_names(fread(file)))
}
my_name_fixer = function(tbl){
  colnames(tbl) = gsub(colnames(tbl)[[9]],"psi",colnames(tbl))
  return(tbl)
}
file_path = "data/sj_tabs/normalized_annotated/"

suffix = "_normalized_annotated.csv"
psi_files = list.files(file_path,
                       pattern = suffix,
                       full.names = TRUE)
psi_list_full = purrr::map(psi_files,my_clean_reader)
samp_ids = base::tolower(purrr::simplify(purrr::map(psi_files, basename)))
samp_ids = gsub(suffix,"",samp_ids)
##add on the name as an additional column
psi_list_full = purrr::map2(psi_list_full, samp_ids, ~cbind(.x, SampleID = .y))
psi_list_full = purrr::map(psi_list_full,my_name_fixer)

psi_list_full = data.table::rbindlist(psi_list_full)




parent_cryptic <- psi_list_full %>%
  dplyr::filter(seqnames %in% exons_parent$seqnames & start > min(exons_parent$start) & end < max(exons_parent$end)) %>%
  #group_by(seqnames, start, end, strand_junction) %>%
  dplyr::mutate(sample_category = paste0(word(SampleID, 2, sep = "_"), "-", word(SampleID, 3, sep = "_"))) #%>%
  #dplyr::filter(sample_category == "dox-ctrl") %>% # | sample_category == "ctrl-ctrl") %>%
  #dplyr::filter(type != "annotated") %>%
  #dplyr::filter(psi > 0.1)


psi_list_unc <- psi_list_full %>%
  group_by(seqnames, start, end, strand_junction) %>%
  dplyr::filter(gene_id_junction == "ENSG00000130477.15") %>%
  ggplot(aes())






hedlund_axonal = fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/axon_seq_human_detected.csv")# doi: 10.1016/j.expneurol.2018.06.008.
maciel_2018 = fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/another_mn_transcriptome.csv")
maciel_2018 = maciel_2018 %>% 
  mutate(compart = case_when(padj<0.05 & log2FoldChange > 0 ~ "axonal",
                             padj<0.05 & log2FoldChange < 0 ~ "soma",
                             T ~ "ns"))
axonal_genes_either = unique(union(hedlund_axonal %>% pull(`Gene ID`),maciel_2018 %>% filter(compart == "axonal") %>% pull(`Gene Symbol`)))