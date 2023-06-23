new_model_b <- psi_list_full %>%
  dplyr::filter(seqnames == "chr19" & start > 17601328 & end < 17688365) %>%
  dplyr::mutate(sample_category = str_replace_all(paste0(word(SampleID, 2, sep = "_"), "-", word(SampleID, 3, sep = "_")),"[^[:alnum:]]", "")) %>%
  dplyr::filter(sample_category %in% contrast) %>% # | sample_category == "ctrl-ctrl") %>%
  dplyr::filter(type == "novel_donor" | type == "novel_acceptor") %>%
  dplyr::filter(psi > 0.1 & psi < 1)

new_model <- new_model_b %>%
  mutate(paste_into_igv_junction = paste0(seqnames, ":", start, "-", end)) %>%
  distinct(paste_into_igv_junction, .keep_all = T) %>%
  mutate(start = start - 1) %>%
  dplyr::select(c(1:4,7))
  
write.table(new_model, file="foo.bed", quote=F, sep="\t", row.names=F, col.names=F)
