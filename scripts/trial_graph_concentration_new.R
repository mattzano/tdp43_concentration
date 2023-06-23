####create a filtered table with deltas and normalized values
#write.csv(big_data[,c(1,2,5,9,11,17)], "~/Desktop/big_data.csv", row.names = F)

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
big_data_filtereds <- big_data[,c(1,2,5,9,11,17)] %>% #1701012/6 = 283502
  pivot_wider(names_from = .id, values_from = mean_psi_per_lsv_junction, values_fn = max) %>% ##232472
  #filter(no_dox < 0.05) %>%  #114616  ####to include all splicing events comment this line
  filter(dox_075-no_dox > 0.1) %>% #1669
  mutate(delta075025 = dox_075 - dox_025) %>%
  mutate(delta025021 = dox_025 - dox_021) %>%
  mutate(delta021018 = dox_021 - dox_0187) %>%
  mutate(delta018012 = dox_0187 - dox_0125) %>%
  mutate(delta012000 = dox_0125 - no_dox) %>%
  filter((delta075025   > - 0.1) &
           (delta025021 > - 0.1) &
           (delta021018 > - 0.1) &
           (delta018012 > - 0.1) &
           (delta012000 > - 0.1)) %>%
  mutate(cryptic = ifelse(no_dox < 0.05 & de_novo_junctions == 0, "cryptic", "no_cryptic"))#1346
big_data_filtereds_dedup <- big_data_filtereds[!duplicated(big_data_filtereds$paste_into_igv_junction),] #1325

table(big_data_filtereds$cryptic)
big_data_filtereds_dedup %>% filter(no_dox < 0.05)
873/4350


big_data_filtereds_normal <- big_data_filtereds_dedup %>%
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
           (delta012000normal > - 0.1)) ###410
#big_data_filtereds_normal_dedup <- big_data_filtereds_normal[!duplicated(big_data_filtereds_normal$paste_into_igv_junction),] #410

##play with different percentiles and report results for each
## 1st - 99th ## 5th - 95th ## 10th - 90th ## 20th - 80th ## 25th - 75th ## 33rd - 66th
percentile_top <- c(#99, 95, 90, 80, 75, 
  67) #, 50) ###I think 67 or 75 is a good place to start for motif etc


#####method?
big_data_filtereds_early_late_new <- big_data_filtereds_dedup %>%
  mutate(color_gene_name = as.factor(as.character(ifelse((dox_0125-no_dox > 0.1 | dox_0187-no_dox > 0.1) & 
                                                           dox_075-no_dox > 0.2, "early", 
                                                                ifelse(dox_025-no_dox < 0.1 & 
                                                                         dox_075-no_dox > 0.2, "late",
                                                                       ifelse((dox_021-no_dox > 0.1 | dox_025-no_dox > 0.1) & 
                                                                                dox_075-no_dox > 0.2, "intermediate",
                                                                              "none")))))) %>%
  pivot_longer(cols = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075")) %>%
  mutate(name = factor(name, levels = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075")))
#data_plot2 <- as.data.frame(table(big_data_filtereds_early_late_new$color_gene_name, big_data_filtereds_early_late_new$cryptic))

plot_early_late <- big_data_filtereds_early_late_new %>%
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


###Method 0 - with names
list_gene <- c("ATG4B")
big_data_filtereds_with_names <- big_data_filtereds_dedup %>%
  filter(no_dox < 0.05) %>%
  #mutate(alpha_gene_name = ifelse(gene_name %in% list_gene, 1,0.2)) %>%
  #dplyr::mutate(label_junction = case_when(.id == input_list[length(input_list) - 3] & (color_gene_name == "late" | color_gene_name == "early") ~ gene_name, T ~ ""))
           #as.factor(as.character(ifelse(dox_0187 > 0.3 & dox_075 > 0.3, "early", 
          #                                               ifelse(dox_025 < 0.1 & dox_075 > 0.3, "late", "none"))))) %>%
  pivot_longer(cols = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075")) %>%
  mutate(name = factor(name, levels = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075"))) %>%
  mutate(label_junction = case_when(name == "dox_025" &
                                     #(gene_name == "STMN2" | 
                                        gene_name %in% list_gene #|
                                        #(gene_name == "HDGFL2" & value > 0.40) |
                                        #(gene_name == "UNC13A" & value > 0.4))) #| (gene_name %in% list_gene & name == "dox_075") 
                                    ~ gene_name, T ~ "")) %>%
  mutate(color_gene_name = ifelse(gene_name %in% list_gene, "2", "1")) %>%
                                  #ifelse(paste_into_igv_junction == "chr8:79613937-79616822" |
                                  #       paste_into_igv_junction == "chr19:4492152-4493703" |
                                  #       paste_into_igv_junction == "chr19:17641556-17642414", "2","1")) %>%
  mutate(alpha_gene_name = ifelse(color_gene_name == "2", "2", "1"))         
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
  scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++")) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_classic()
print(plot_early_late)
#ggsave(filename = "~/Desktop/spaghetti_plot_with_name_KCNQ2.png")


##### METHOD 1: early-late #####
#### raw PSI ####
big_data_filtereds_early_late <- big_data_filtereds_dedup %>%
  mutate(color_gene_name = as.factor(as.character(ifelse(paste_into_igv_junction == "chr8:79613937-79616822" |
                                                           paste_into_igv_junction == "chr19:17641556-17642414", "special",
                                                         ifelse(dox_0187 > 0.4 & dox_075 > 0.4, "early", 
                                                         ifelse(dox_025 < 0.1 & dox_075 > 0.4, "late", "none")))))) %>%
  pivot_longer(cols = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075")) %>%
  mutate(name = factor(name, levels = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075")))
table(big_data_filtereds_early_late$color_gene_name)

plot_early_late <- big_data_filtereds_early_late %>%
  ggplot(mapping = aes(x = name, y = value, group = paste_into_igv_junction)) +
  geom_line(aes(color = color_gene_name, alpha = color_gene_name), show.legend = F) +
  #geom_text_repel(aes(label = label_junction), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  scale_color_manual(values = c("#FAAC5E","#D53E4F","gray","#3288BD")) +
  scale_alpha_manual(values = c(1,1,0.2,1)) +
  geom_hline(yintercept = 0.1, linetype = "dotted") +
  xlab("TDP-43 knockdown level") +
  ylab(expression(paste(Psi))) +
  scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++")) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_classic() +
  theme(text = element_text(size = 22))
print(plot_early_late)
#write.csv(big_data_filtereds_early_late, "~/Desktop/early_late.csv", row.names = F)

#### using normalized PSI ####
big_data_filtereds_early_late_normal <- big_data_filtereds_normal %>%
  mutate(color_gene_name = as.factor(as.character(ifelse(dox_025normal < 0.1, "late", 
                                                         ifelse(dox_0125normal > 0.3, "early", "none"))))) %>%
  pivot_longer(cols = c("no_doxnormal","dox_0125normal", "dox_0187normal", "dox_021normal", "dox_025normal", "dox_075normal")) %>%
  mutate(name = factor(name, levels = c("no_doxnormal","dox_0125normal", "dox_0187normal", "dox_021normal", "dox_025normal", "dox_075normal")))
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
#write.csv(big_data_filtereds_early_late_normal, "~/Desktop/early_late_normal.csv", row.names = F)



##### METHOD 2: using slope #####
#### raw PSI ####
for (i in percentile_top) {
big_data_filtereds_slope <- big_data_filtereds_dedup %>%
  mutate(max_slope = pmax(delta075025, delta025021, delta021018, delta018012, delta012000)) %>%
  mutate(percentile = ntile(max_slope, 100)) %>%
  mutate(color_gene_name = as.factor(as.character(ifelse(percentile > i, "high_slope",
                                                         ifelse(percentile <= (100-i), "low_slope", "mid_slope"))))) %>%
  pivot_longer(cols = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075")) %>%
  mutate(name = factor(name, levels = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075"))) %>%
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
#write.csv(big_data_filtereds_slope, file = paste0("~/Desktop/slope_", i,".csv"), row.names = F)
}

#### normalized PSI ####
for (i in percentile_top) {
big_data_filtereds_slope_normal <- big_data_filtereds_normal %>%
  mutate(max_slope = pmax(delta075025normal, delta025021normal, delta021018normal, delta018012normal, delta012000normal)) %>%
  mutate(percentile = ntile(max_slope, 100)) %>%
  mutate(color_gene_name = as.factor(as.character(ifelse(percentile > i, "high_slope",
                                                         ifelse(percentile <= (100-i), "low_slope", "mid_slope"))))) %>%
  pivot_longer(cols = c("no_doxnormal","dox_0125normal", "dox_0187normal", "dox_021normal", "dox_025normal", "dox_075normal")) %>%
  mutate(name = factor(name, levels = c("no_doxnormal","dox_0125normal", "dox_0187normal", "dox_021normal", "dox_025normal", "dox_075normal"))) %>%
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
#write.csv(big_data_filtereds_slope_normal, file = paste0("~/Desktop/slope_normal_", i,".csv"), row.names = F)
}



#### METHOD 3: raw PSI when TDP43 is 0% ####
for (i in percentile_top) {
big_data_filtereds_075 <- big_data_filtereds_dedup %>%
  mutate(color_gene_name = as.factor(as.character(ifelse(paste_into_igv_junction == "chr8:79613937-79616822" |
                                                           paste_into_igv_junction == "chr19:17641556-17642414", "special",
                                                         ifelse(dox_075 > 0.5, "early", 
                                                                ifelse(dox_075 < 0.2, "late", "none")))))) %>%
  
  #mutate(percentile = ntile(dox_075, 100)) %>%
  #mutate(color_gene_name = as.factor(as.character(ifelse(percentile > i, "high_psi",
  #                                                       ifelse(percentile <= (100-i), "low_psi", "mid_psi"))))) %>%
  pivot_longer(cols = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075")) %>%
  mutate(name = factor(name, levels = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075"))) %>%
  dplyr::mutate(alpha_gene_name = as.factor(as.character(ifelse(color_gene_name == "high_psi" | color_gene_name == "low_psi", 1, 2))))
print(table(big_data_filtereds_075$color_gene_name))

plot <- big_data_filtereds_075 %>%
  ggplot(mapping = aes(x = name, y = value, group = paste_into_igv_junction)) +
  geom_line(aes(color = color_gene_name, alpha = alpha_gene_name), show.legend = F) +
  #geom_text_repel(aes(label = label_junction), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  scale_color_manual(values = c("#FAAC5E","#D53E4F","gray","#3288BD")) +
  scale_alpha_manual(values = c(1,1,0.2,1)) +
  geom_hline(yintercept = 0.1, linetype = "dotted") +
  xlab("TDP-43 knockdown level") +
  ylab(expression(paste(Psi))) +
  scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++")) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_classic() +
  theme(text = element_text(size = 22))
print(plot)
#write.csv(big_data_filtereds_075, file = paste0("~/Desktop/max_psi_", i,".csv"), row.names = F)
}


##### METHOD: extra - IGNORE FOR NOW #####
#big_data_filtereds_combo <- big_data_filtereds_dedup %>%
#  mutate(max_slope = pmax(delta075025, delta025021, delta021018, delta018012, delta012000)) %>%
#  mutate(percentile = ntile(max_slope, 100)) %>%
#  mutate(color_gene_name = as.factor(as.character(#ifelse(gene_name == "STMN2" | gene_name == "UNC13A", 0,
#    ifelse(percentile > 80 & (dox_0125 > 0.2 | dox_0187 > 0.2), 0, 
#           ifelse(percentile > 80 & dox_021 > 0.2 | dox_025 > 0.2, 1, 
#                  ifelse(percentile > 80, 2, 3)))))) %>%
#  pivot_longer(cols = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075")) %>%
#  mutate(name = factor(name, levels = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075"))) %>%
#  dplyr::mutate(alpha_gene_name = as.factor(as.character(ifelse(color_gene_name == 0 | color_gene_name == 1 | color_gene_name == 2, 1, 2)))) %>%
#  dplyr::mutate(label_junction = case_when(name =="dox_0187" & (color_gene_name == "late" | color_gene_name == "early") ~ gene_name, T ~ ""))




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
  
#big_data_filtereds2 <- big_data[,c(1,2,5,9,11,17)] %>% #1701012 /6 = 283502
#  pivot_wider(names_from = .id, values_from = mean_psi_per_lsv_junction, values_fn = max) %>% ##232472
#  filter(dox_075 - no_dox > 0.1) %>% #6383
#  mutate(delta075025 = dox_075 - dox_025) %>%
#  mutate(delta025021 = dox_025 - dox_021) %>%
#  mutate(delta021018 = dox_021 - dox_0187) %>%
#  mutate(delta018012 = dox_0187 - dox_0125) %>%
#  mutate(delta012000 = dox_0125 - no_dox) %>%
#  filter((delta075025   > - 0.1) &
#           (delta025021 > - 0.1) &
#           (delta021018 > - 0.1) &
#           (delta018012 > - 0.1) &
#           (delta012000 > - 0.1)) %>% #4417
#  mutate(dox_075normal = (dox_075 - no_dox)  / (dox_075 - no_dox)) %>%
#  mutate(dox_025normal = (dox_025 - no_dox)  / (dox_075 - no_dox)) %>%
#  mutate(dox_021normal = (dox_021 - no_dox)  / (dox_075 - no_dox)) %>%
#  mutate(dox_0187normal = (dox_0187 - no_dox)  / (dox_075 - no_dox)) %>%
#  mutate(dox_0125normal = (dox_0125 - no_dox)  / (dox_075 - no_dox)) %>%
#  mutate(no_doxnormal = (no_dox - no_dox)  / (dox_075 - no_dox)) %>%
#  mutate(delta075025normal = dox_075normal - dox_025normal) %>%
#  mutate(delta025021normal = dox_025normal - dox_021normal) %>%
#  mutate(delta021018normal = dox_021normal - dox_0187normal) %>%
#  mutate(delta018012normal = dox_0187normal - dox_0125normal) %>%
#  mutate(delta012000normal = dox_0125normal - no_doxnormal) %>%
#  filter((delta075025normal   > - 0.1) &
#           (delta025021normal > - 0.1) &
#           (delta021018normal > - 0.1) &
#           (delta018012normal > - 0.1) &
#           (delta012000normal > - 0.1)) ###418   #1027
#big_data_filtereds2_dedup <- big_data_filtereds2[!duplicated(big_data_filtereds2$paste_into_igv_junction),]

#big_data_filtereds2_early_late <- big_data_filtereds2 %>%
#  mutate(color_gene_name = as.factor(as.character(ifelse(dox_0187 > 0.5, "early", "none")))) %>%
#  pivot_longer(cols = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075")) %>%
#  mutate(name = factor(name, levels = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075")))

#plot2_early_late <- big_data_filtereds2_early_late %>%
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
#print(plot2_early_late)



#### ####