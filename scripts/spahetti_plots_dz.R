####create a filtered table with deltas and normalized values
#write.csv(big_data[,c(1,2,5,9,11,17)], "~/Desktop/big_data.csv", row.names = F)

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
big_data_dz_filtereds <- big_data_dz[,c(1,2,5,9,11,17)] %>% #1901543/6 = 316923.8
  pivot_wider(names_from = .id, values_from = mean_psi_per_lsv_junction, values_fn = max) %>% ##195298
  filter(DZ_curves_0 < 0.06) %>%  #91339 ##remove this to include all splicing changes
  filter(DZ_curves_1 > 0.15) %>% #0.15 as in sy5y?
  mutate(delta101 = DZ_curves_1 - DZ_curves_01) %>%
  mutate(delta010075 = DZ_curves_01 - DZ_curves_0075) %>%
  mutate(delta0075005 = DZ_curves_0075 - DZ_curves_005) %>%
  mutate(delta005004 = DZ_curves_005 - DZ_curves_004) %>%
  mutate(delta004003 = DZ_curves_004 - DZ_curves_003) %>%
  mutate(delta003002 = DZ_curves_003 - DZ_curves_002) %>%
  mutate(delta0020 = DZ_curves_002 - DZ_curves_0) %>%
  filter((delta101 > - 0.1) &
           (delta010075 > - 0.1) &
           (delta0075005 > - 0.1) &
           (delta005004 > - 0.1) &
           (delta004003 > - 0.1) &
           (delta003002 > - 0.1) &
           (delta0020 > - 0.1)) #%>%
  #mutate(cryptic = ifelse(DZ_curves_0 < 0.06, "cryptic", "no_cryptic")) #915
big_data_dz_filtereds_dedup <- big_data_dz_filtereds[!duplicated(big_data_dz_filtereds$paste_into_igv_junction),] #901
table(big_data_dz_filtereds_dedup$cryptic)

big_data_dz_filtereds_early_late <- big_data_dz_filtereds_dedup %>%
  mutate(color_gene_name = as.factor(as.character(ifelse((DZ_curves_002-DZ_curves_0 > 0.1 | DZ_curves_003-DZ_curves_0 > 0.1) & 
                                                           DZ_curves_1-DZ_curves_0 > 0.2, "early", 
                                                         ifelse(DZ_curves_0075-DZ_curves_0 < 0.1 & 
                                                                  (DZ_curves_01-DZ_curves_0 > 0.1 | DZ_curves_1-DZ_curves_0 > 0.2), "late",
                                                                ifelse((DZ_curves_004-DZ_curves_0 > 0.1 | DZ_curves_005-DZ_curves_0 > 0.1 | DZ_curves_0075-DZ_curves_0 > 0.1) &
                                                                         DZ_curves_1-DZ_curves_0 > 0.2, "intermediate",
                                                                       "none")))))) %>%
  pivot_longer(cols = c("DZ_curves_0","DZ_curves_002", "DZ_curves_003", "DZ_curves_004", "DZ_curves_005", "DZ_curves_0075", "DZ_curves_01", "DZ_curves_1")) %>%
  mutate(name = factor(name, levels = c("DZ_curves_0","DZ_curves_002", "DZ_curves_003", "DZ_curves_004", "DZ_curves_005", "DZ_curves_0075", "DZ_curves_01", "DZ_curves_1")))
table(big_data_dz_filtereds_early_late$color_gene_name, big_data_dz_filtereds_early_late$cryptic)

plot_dz_early_late <- big_data_dz_filtereds_early_late %>%
  ggplot(mapping = aes(x = name, y = value, group = paste_into_igv_junction)) +
  geom_line(aes(color = color_gene_name, alpha = color_gene_name), show.legend = F) +
  #geom_text_repel(aes(label = label_junction), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  scale_color_manual(values = c("#FAAC5E","#D53E4F","#3288BD","gray")) +
  scale_alpha_manual(values = c(1,1,1,0.1)) +
  xlab("TDP-43 knockdown level") +
  ylab("Normalised PSI") +
  scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++", "++++++", "+++++++")) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  facet_wrap(cryptic ~ color_gene_name, nrow = 2) +
  #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_classic()
print(plot_dz_early_late)
#write.csv(big_data_filtereds_early_late, "~/Desktop/early_late_dz.csv", row.names = F)






big_data_filtereds_normal <- big_data_filtereds_dedup %>%
  mutate(DZ_curves_1normal = (DZ_curves_1 - DZ_curves_0)  / (DZ_curves_1 - DZ_curves_0)) %>%
  mutate(DZ_curves_01normal = (DZ_curves_01 - DZ_curves_0)  / (DZ_curves_1 - DZ_curves_0)) %>%
  mutate(DZ_curves_0075normal = (DZ_curves_0075 - DZ_curves_0)  / (DZ_curves_1 - DZ_curves_0)) %>%
  mutate(DZ_curves_005normal = (DZ_curves_005 - DZ_curves_0)  / (DZ_curves_1 - DZ_curves_0)) %>%
  mutate(DZ_curves_004normal = (DZ_curves_004 - DZ_curves_0)  / (DZ_curves_1 - DZ_curves_0)) %>%
  mutate(DZ_curves_003normal = (DZ_curves_003 - DZ_curves_0)  / (DZ_curves_1 - DZ_curves_0)) %>%
  mutate(DZ_curves_002normal = (DZ_curves_002 - DZ_curves_0)  / (DZ_curves_1 - DZ_curves_0)) %>%
  mutate(DZ_curves_0normal = (DZ_curves_0 - DZ_curves_0)  / (DZ_curves_1 - DZ_curves_0)) %>%
  mutate(delta101normal = DZ_curves_1normal - DZ_curves_01normal) %>%
  mutate(delta010075normal = DZ_curves_01normal - DZ_curves_0075normal) %>%
  mutate(delta0075005normal = DZ_curves_0075normal - DZ_curves_005normal) %>%
  mutate(delta005004normal = DZ_curves_005normal - DZ_curves_004normal) %>%
  mutate(delta004003normal = DZ_curves_004normal - DZ_curves_003normal) %>%
  mutate(delta003002normal = DZ_curves_003normal - DZ_curves_002normal) %>%
  mutate(delta0020normal = DZ_curves_002normal - DZ_curves_0normal) %>%
  filter((delta101normal   > - 0.1) &
           (delta010075normal > - 0.1) &
           (delta0075005normal > - 0.1) &
           (delta005004normal > - 0.1) &
           (delta004003normal > - 0.1) &
           (delta003002normal > - 0.1) &
           (delta0020normal > - 0.1)) #218     ###455
##play with different percentiles and report results for each
## 1st - 99th ## 5th - 95th ## 10th - 90th ## 20th - 80th ## 25th - 75th ## 33rd - 66th
percentile_top <- c(99, 95, 90, 80, 75, 67, 50)

###Method 0 - with names
list_gene <- c(#"CYFIP2", "SYNE1", 
  "AARS1")
big_data_filtereds_with_names <- big_data_filtereds_dedup %>%
  #mutate(alpha_gene_name = ifelse(gene_name %in% list_gene, 1,0.2)) %>%
  #dplyr::mutate(label_junction = case_when(.id == input_list[length(input_list) - 3] & (color_gene_name == "late" | color_gene_name == "early") ~ gene_name, T ~ ""))
  #as.factor(as.character(ifelse(dox_0187 > 0.3 & dox_075 > 0.3, "early", 
  #                                               ifelse(dox_025 < 0.1 & dox_075 > 0.3, "late", "none"))))) %>%
  pivot_longer(cols = c("DZ_curves_0","DZ_curves_002", "DZ_curves_003", "DZ_curves_004", "DZ_curves_005", "DZ_curves_0075", "DZ_curves_01", "DZ_curves_1")) %>%
  mutate(name = factor(name, levels = c("DZ_curves_0","DZ_curves_002", "DZ_curves_003", "DZ_curves_004", "DZ_curves_005", "DZ_curves_0075", "DZ_curves_01", "DZ_curves_1"))) %>%
  mutate(label_junction = case_when(name == "DZ_curves_1" &
                                       ((gene_name == "STMN2" & value > 0.5 ) | 
                                    #      (gene_name == "AARS1" & value > 0.12) |
                                         (gene_name == "UNC13A" & value > 0.51) | 
                                      (gene_name %in% list_gene & name == "DZ_curves_1")) ~ gene_name, T ~ "")) %>%
  mutate(color_gene_name = ifelse(gene_name %in% list_gene | 
                                    paste_into_igv_junction == "chr19:17641556-17642414" |
                                    paste_into_igv_junction == "chr19:17642541-17642845" |
                                    paste_into_igv_junction == "chr8:79613937-79616822", "3",
                                  #ifelse(paste_into_igv_junction == "chr8:79613937-79616822" |
                                  #         paste_into_igv_junction == "chr16:70271972-70272796" |
                                  #         paste_into_igv_junction == "chr19:17641556-17642414", "2",
                                  "1")) %>%
  mutate(alpha_gene_name = ifelse(color_gene_name == "3" | color_gene_name == "2", "2", "1"))         
table(big_data_filtereds_with_names$color_gene_name)


plot_early_late <- big_data_filtereds_with_names %>%
  ggplot(mapping = aes(x = name, y = value, group = paste_into_igv_junction)) +
  geom_line(aes(color = color_gene_name, alpha = alpha_gene_name), show.legend = F) +
  geom_text_repel(aes(label = label_junction), point.padding = 0.3,
                  nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  geom_hline(yintercept = 0.1, linetype = "dotted") +
  scale_color_manual(values = c("gray", "#3288BD", "#D53E4F")) +
  scale_alpha_manual(values = c(0.1,1)) +
  xlab("TDP-43 knockdown level") +
  ylab("PSI") +
  scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++", "++++++", "+++++++")) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_classic()
print(plot_early_late)


##### METHOD 1: early-late #####
#### raw PSI ####
big_data_filtereds_early_late <- big_data_filtereds_dedup %>%
  mutate(color_gene_name = as.factor(as.character(ifelse(DZ_curves_004 > 0.2 & DZ_curves_1 > 0.3, "early", 
                                                         ifelse(DZ_curves_01 < 0.1 & DZ_curves_1 > 0.3, "late", "none"))))) %>%
  pivot_longer(cols = c("DZ_curves_0","DZ_curves_002", "DZ_curves_003", "DZ_curves_004", "DZ_curves_005", "DZ_curves_0075", "DZ_curves_01", "DZ_curves_1")) %>%
  mutate(name = factor(name, levels = c("DZ_curves_0","DZ_curves_002", "DZ_curves_003", "DZ_curves_004", "DZ_curves_005", "DZ_curves_0075", "DZ_curves_01", "DZ_curves_1")))
table(big_data_filtereds_early_late$color_gene_name)

plot_early_late <- big_data_filtereds_early_late %>%
  ggplot(mapping = aes(x = name, y = value, group = paste_into_igv_junction)) +
  geom_line(aes(color = color_gene_name, alpha = color_gene_name), show.legend = F) +
  #geom_text_repel(aes(label = label_junction), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  scale_color_manual(values = c("#3288BD","#D53E4F","gray")) +
  scale_alpha_manual(values = c(1,1,0.2)) +
  xlab("TDP-43 knockdown level") +
  ylab("Normalised PSI") +
  scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++")) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_classic()
print(plot_early_late)
write.csv(big_data_filtereds_early_late, "~/Desktop/early_late_dz.csv", row.names = F)

#### using normalized PSI ####
big_data_filtereds_early_late_normal <- big_data_filtereds_normal %>%
  mutate(color_gene_name = as.factor(as.character(ifelse(DZ_curves_01normal < 0.1, "late", 
                                                         ifelse(DZ_curves_003normal > 0.3, "early", "none"))))) %>%
  pivot_longer(cols = c("DZ_curves_0normal","DZ_curves_002normal", "DZ_curves_003normal", "DZ_curves_004normal", 
                        "DZ_curves_005normal", "DZ_curves_0075normal", "DZ_curves_01normal", "DZ_curves_1normal")) %>%
  mutate(name = factor(name, levels = c("DZ_curves_0normal","DZ_curves_002normal", "DZ_curves_003normal", "DZ_curves_004normal", 
                                        "DZ_curves_005normal", "DZ_curves_0075normal", "DZ_curves_01normal", "DZ_curves_1normal")))
table(big_data_filtereds_early_late_normal$color_gene_name)

plot_early_late_normal <- big_data_filtereds_early_late_normal %>%
  ggplot(mapping = aes(x = name, y = value, group = paste_into_igv_junction)) +
  geom_line(aes(color = color_gene_name, alpha = color_gene_name), show.legend = F) +
  #geom_text_repel(aes(label = label_junction), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  scale_color_manual(values = c("#3288BD","#D53E4F","gray")) +
  scale_alpha_manual(values = c(1,1,0.2)) +
  xlab("TDP-43 knockdown level") +
  ylab("Normalised PSI") +
  scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++")) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_classic()
print(plot_early_late_normal)
write.csv(big_data_filtereds_early_late_normal, "~/Desktop/early_late_normal_dz.csv", row.names = F)



##### METHOD 2: using slope #####
#### raw PSI ####
for (i in percentile_top) {
  big_data_filtereds_slope <- big_data_filtereds_dedup %>%
    mutate(max_slope = pmax(delta101, delta010075, delta0075005, delta005004, delta004003, delta003002, delta0020)) %>%
    mutate(percentile = ntile(max_slope, 100)) %>%
    mutate(color_gene_name = as.factor(as.character(ifelse(percentile > i, "high_slope",
                                                           ifelse(percentile <= (100-i), "low_slope", "mid_slope"))))) %>%
    pivot_longer(cols = c("DZ_curves_0","DZ_curves_002", "DZ_curves_003", "DZ_curves_004", "DZ_curves_005", "DZ_curves_0075", "DZ_curves_01", "DZ_curves_1")) %>%
    mutate(name = factor(name, levels = c("DZ_curves_0","DZ_curves_002", "DZ_curves_003", "DZ_curves_004", "DZ_curves_005", "DZ_curves_0075", "DZ_curves_01", "DZ_curves_1"))) %>%
    dplyr::mutate(alpha_gene_name = as.factor(as.character(ifelse(color_gene_name == "high_slope" | color_gene_name == "low_slope", 1, 2)))) #%>%
  #dplyr::mutate(label_junction = case_when(name =="dox_0187" & (color_gene_name == "late" | color_gene_name == "early") ~ gene_name, T ~ ""))
  print(table(big_data_filtereds_slope$color_gene_name))
  
  plot <- big_data_filtereds_slope %>%
    ggplot(mapping = aes(x = name, y = value, group = paste_into_igv_junction)) +
    geom_line(aes(color = color_gene_name, alpha = alpha_gene_name), show.legend = F) +
    #geom_text_repel(aes(label = label_junction), point.padding = 0.3,
    #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
    scale_color_manual(values = c("#046C9A", 
                                  #"#FAAE61", 
                                  "#D53E4F", 
                                  "gray")) +
    scale_alpha_manual(values = c(1, 0.2)) +
    #geom_vline(linetype = "dotted", xintercept = 2) +
    #geom_vline(linetype = "dotted", xintercept = 3) +
    #geom_vline(linetype = "dotted", xintercept = 4) +
    #geom_vline(linetype = "dotted", xintercept = 5) +
    xlab("TDP-43 knockdown level") +
    ylab("PSI") +
    #scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++")) +
    scale_x_discrete(labels = c("100%","77%", "25%", "8%", "4%", "0%")) +
    scale_y_continuous(breaks=seq(0,1,0.2)) +
    #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
    theme_classic()
  print(plot)
  write.csv(big_data_filtereds_slope, file = paste0("~/Desktop/slope_", i,"_dz.csv"), row.names = F)
}

#### normalized PSI ####
for (i in percentile_top) {
  big_data_filtereds_slope_normal <- big_data_filtereds_normal %>%
    mutate(max_slope = pmax(delta101normal, delta010075normal, delta0075005normal, delta005004normal, delta004003normal, delta003002normal, delta0020normal)) %>%
    mutate(percentile = ntile(max_slope, 100)) %>%
    mutate(color_gene_name = as.factor(as.character(ifelse(percentile > i, "high_slope",
                                                           ifelse(percentile <= (100-i), "low_slope", "mid_slope"))))) %>%
    pivot_longer(cols = c("DZ_curves_0normal","DZ_curves_002normal", "DZ_curves_003normal", "DZ_curves_004normal", 
                          "DZ_curves_005normal", "DZ_curves_0075normal", "DZ_curves_01normal", "DZ_curves_1normal")) %>%
    mutate(name = factor(name, levels = c("DZ_curves_0normal","DZ_curves_002normal", "DZ_curves_003normal", "DZ_curves_004normal", 
                                          "DZ_curves_005normal", "DZ_curves_0075normal", "DZ_curves_01normal", "DZ_curves_1normal"))) %>%
  dplyr::mutate(alpha_gene_name = as.factor(as.character(ifelse(color_gene_name == "high_slope" | color_gene_name == "low_slope", 1, 2)))) 
  print(table(big_data_filtereds_slope_normal$color_gene_name))
  
  plot <- big_data_filtereds_slope_normal %>%
    ggplot(mapping = aes(x = name, y = value, group = paste_into_igv_junction)) +
    geom_line(aes(color = color_gene_name, alpha = color_gene_name), show.legend = F) +
    #geom_text_repel(aes(label = label_junction), point.padding = 0.3,
    #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
    scale_color_manual(values = c("#3288BD","#D53E4F","gray")) +
    scale_alpha_manual(values = c(1,1,0.2)) +
    xlab("TDP-43 knockdown level") +
    ylab("Normalised PSI") +
    scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++")) +
    scale_y_continuous(breaks=seq(0,1,0.2)) +
    #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
    theme_classic()
  print(plot)
  write.csv(big_data_filtereds_slope_normal, file = paste0("~/Desktop/slope_normal_", i,"_dz.csv"), row.names = F)
}



#### METHOD 3: raw PSI when TDP43 is 0% ####
for (i in percentile_top) {
  big_data_filtereds_075 <- big_data_filtereds_dedup %>%
    mutate(percentile = ntile(DZ_curves_1, 100)) %>%
    mutate(color_gene_name = as.factor(as.character(ifelse(percentile > i, "high_psi",
                                                           ifelse(percentile <= (100-i), "low_psi", "mid_psi"))))) %>%
    pivot_longer(cols = c("DZ_curves_0","DZ_curves_002", "DZ_curves_003", "DZ_curves_004", "DZ_curves_005", "DZ_curves_0075", "DZ_curves_01", "DZ_curves_1")) %>%
    mutate(name = factor(name, levels = c("DZ_curves_0","DZ_curves_002", "DZ_curves_003", "DZ_curves_004", "DZ_curves_005", "DZ_curves_0075", "DZ_curves_01", "DZ_curves_1"))) %>%
    dplyr::mutate(alpha_gene_name = as.factor(as.character(ifelse(color_gene_name == "high_psi" | color_gene_name == "low_psi", 1, 2))))
  print(table(big_data_filtereds_075$color_gene_name))
  
  plot <- big_data_filtereds_075 %>%
    ggplot(mapping = aes(x = name, y = value, group = paste_into_igv_junction)) +
    geom_line(aes(color = color_gene_name, alpha = alpha_gene_name), show.legend = F) +
    #geom_text_repel(aes(label = label_junction), point.padding = 0.3,
    #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
    scale_color_manual(values = c("#C93312", "#046C9A","#899DA4")) +
    scale_alpha_manual(values = c(1,0.2)) +
    xlab("TDP-43 knockdown level") +
    ylab("PSI") +
    scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++")) +
    scale_y_continuous(breaks=seq(0,1,0.2)) +
    #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
    theme_classic()
  print(plot)
  write.csv(big_data_filtereds_075, file = paste0("~/Desktop/max_psi_", i,"_dz.csv"), row.names = F)
}


##### METHOD: extra - IGNORE FOR NOW #####
big_data_filtereds_combo <- big_data_filtereds_dedup %>%
  mutate(max_slope = pmax(delta075025, delta025021, delta021018, delta018012, delta012000)) %>%
  mutate(percentile = ntile(max_slope, 100)) %>%
  mutate(color_gene_name = as.factor(as.character(#ifelse(gene_name == "STMN2" | gene_name == "UNC13A", 0,
    ifelse(percentile > 80 & (dox_0125 > 0.2 | dox_0187 > 0.2), 0, 
           ifelse(percentile > 80 & dox_021 > 0.2 | dox_025 > 0.2, 1, 
                  ifelse(percentile > 80, 2, 3)))))) %>%
  pivot_longer(cols = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075")) %>%
  mutate(name = factor(name, levels = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075"))) %>%
  dplyr::mutate(alpha_gene_name = as.factor(as.character(ifelse(color_gene_name == 0 | color_gene_name == 1 | color_gene_name == 2, 1, 2)))) %>%
  dplyr::mutate(label_junction = case_when(name =="dox_0187" & (color_gene_name == "late" | color_gene_name == "early") ~ gene_name, T ~ ""))




#### then re-include all junctions (based on p value) - IGNORE FOR NOW ####

#names(big_delta)[c(10,11,13)] <- c("dpsi", "probability", "psi_nodox")

#big_delta_filtereds <- big_delta[,c(1,7,10,11,13,18,22,36)] %>% #74067/6 = 12344.5
#  pivot_wider(names_from = .id, values_from = c(dpsi, probability), values_fn = max) %>% #41250
#  filter(probability_dox0075 > 0.8) %>% #2737
#  #filter(psi_nodox < 0.05) #244
#  filter(!is.na(probability_dox0025) & !is.na(probability_dox0021) & !is.na(probability_dox00187) & !is.na(probability_dox00125)) %>%#249
#  mutate(psi_dox00125 = dpsi_dox00125 + psi_nodox) %>%
#  mutate(psi_dox00187 = dpsi_dox00187 + psi_nodox) %>%
#  mutate(psi_dox0021 = dpsi_dox0021 +  psi_nodox) %>%
#  mutate(psi_dox0025 = dpsi_dox0025 +  psi_nodox) %>%
#  mutate(psi_dox0075 = dpsi_dox0075 +  psi_nodox)
#big_delta_filtereds_dedup <- big_delta_filtereds[!duplicated(big_delta_filtereds$paste_into_igv_junction),] #2375


#big_delta_filtereds_early_late <- big_delta_filtereds_dedup %>%
#  #dplyr::filter(psi_dox0075 > psi_nodox) %>%
#  mutate(color_gene_name = as.character(ifelse(psi_dox0075 < psi_nodox, "late", "early"))) %>%
#  pivot_longer(cols = c("psi_nodox","psi_dox00125", "psi_dox00187", "psi_dox0021", "psi_dox0025", "psi_dox0075")) %>%
#  mutate(name = factor(name, levels = c("psi_nodox","psi_dox00125", "psi_dox00187", "psi_dox0021", "psi_dox0025", "psi_dox0075")))

#plot_early_late <- big_delta_filtereds_early_late %>%
#  ggplot(mapping = aes(x = name, y = value, group = paste_into_igv_junction)) +
#  geom_line(aes(color = color_gene_name), show.legend = F) +
#  #geom_text_repel(aes(label = label_junction), point.padding = 0.3,
#  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
#  scale_color_manual(values = c("#3288BD","#D53E4F","gray")) +
#  scale_alpha_manual(values = c(1,1,0.2)) +
#  xlab("TDP-43 knockdown level") +
#  ylab("PSI") +
#  scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++")) +
#  #scale_y_continuous(breaks=seq(0,1,0.2)) +
#  #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
#  theme_classic()
#print(plot_early_late)

big_data_filtereds2 <- big_data[,c(1,2,5,9,11,17)] %>% #1701012 /6 = 283502
  pivot_wider(names_from = .id, values_from = mean_psi_per_lsv_junction, values_fn = max) %>% ##232472
  filter(dox_075 - no_dox > 0.1) %>% #6383
  mutate(delta075025 = dox_075 - dox_025) %>%
  mutate(delta025021 = dox_025 - dox_021) %>%
  mutate(delta021018 = dox_021 - dox_0187) %>%
  mutate(delta018012 = dox_0187 - dox_0125) %>%
  mutate(delta012000 = dox_0125 - no_dox) %>%
  filter((delta075025   > - 0.1) &
           (delta025021 > - 0.1) &
           (delta021018 > - 0.1) &
           (delta018012 > - 0.1) &
           (delta012000 > - 0.1)) %>% #4417
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
           (delta012000normal > - 0.1)) ###418   #1027
big_data_filtereds2_dedup <- big_data_filtereds2[!duplicated(big_data_filtereds2$paste_into_igv_junction),]

big_data_filtereds2_early_late <- big_data_filtereds2 %>%
  mutate(color_gene_name = as.factor(as.character(ifelse(dox_0187 > 0.5, "early", "none")))) %>%
  pivot_longer(cols = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075")) %>%
  mutate(name = factor(name, levels = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075")))

plot2_early_late <- big_data_filtereds2_early_late %>%
  ggplot(mapping = aes(x = name, y = value, group = paste_into_igv_junction)) +
  geom_line(aes(color = color_gene_name), show.legend = F) +
  #geom_text_repel(aes(label = label_junction), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  scale_color_manual(values = c("#3288BD","#D53E4F","gray")) +
  scale_alpha_manual(values = c(1,1,0.2)) +
  xlab("TDP-43 knockdown level") +
  ylab("PSI") +
  scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++")) +
  #scale_y_continuous(breaks=seq(0,1,0.2)) +
  #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_classic()
print(plot2_early_late)



#### ####