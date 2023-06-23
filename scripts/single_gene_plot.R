library(rstatix)
library(tidyverse)

###try and use log2 fc
results_table_sy_merged <- results_table_sy5y_curves_0.075[,c(9,3,7)] %>%
  filter(symbol == "TARDBP"  | symbol == "UNC13A" | 
           symbol == "STMN2" | symbol == "RPTOR" |
           symbol == "CYFIP2" | symbol == "SYNE1" |
           symbol == "AARS" | symbol == "HDGFL2" |
           symbol == "PTPRN2") %>%
  left_join(results_table_sy5y_curves_0.025[,c(3,9,7)], by = "symbol",  suffix = c("_0.075", "_0.025")) %>%
  left_join(results_table_sy5y_curves_0.021[,c(3,9,7)], by = "symbol",  suffix = c("", "_0.021")) %>%
  left_join(results_table_sy5y_curves_0.0187[,c(3,9,7)], by = "symbol", suffix = c("", "_0.0187")) %>%
  left_join(results_table_sy5y_curves_0.0125[,c(3,9,7)], by = "symbol", suffix = c("", "_0.0125")) %>%
  mutate(log2FoldChange_0 = 0,
         padj_0 = 1)
names(results_table_sy_merged)[c(6,7)] <- c("log2FoldChange_0.021", "padj_0.021")

results_table_sy_merged_longer <- results_table_sy_merged %>%
  pivot_longer(cols = c(2,4,6,8,10,12)) %>%
  #mutate(symbol = factor(symbol, levels = c("TARDBP", "STMN2", "UNC13A", "CYFIP2", "RPTOR", "SYNE1", "HDGFL2", "AARS"))) %>%
  mutate(name = gsub(".*_", "", name))

results_table_sy_merged_longer_tdp <- results_table_sy_merged_longer %>%
  dplyr::filter(symbol == "TARDBP")

results_table_sy_merged_longer_target <- results_table_sy_merged_longer %>%
  dplyr::filter(symbol == "PTPRN2")
  #dplyr::filter(symbol != "TARDBP" & symbol != "STMN2")

results_table_sy_merged_longer_target %>%
  #dplyr::filter(symbol == "PTPRN2") %>%
  ggplot() +
  geom_col(data = results_table_sy_merged_longer_tdp, aes(x = name, y = value, fill = symbol), alpha = 0.2, show.legend = T) +
  geom_line(aes(x = name, y = value, color = symbol, group = symbol), size = 1.2) +
  geom_point(aes(x = name, y = value, color = symbol, group = symbol), size = 2) + #
  scale_color_manual(values = c("#FEE08B")) + #, "#66BD63", "#1A9850", "#006837", "#4575B4", "#313695")) +
  scale_fill_manual(values = "#AA0026") + 
  theme_classic() +
  labs(x ="Doxycycline concentration",
       y = "Log2 Fold Change",
       fill = "Gene name",
       color = "")
#ggsave("~/Desktop/downreg_dark_ce.png")


###log2fc in dzs
results_table_dz_merged <- results_table_dz_curves_1[,c(9,3,7)] %>%
  filter(symbol == "TARDBP"  | symbol == "UNC13A" | 
           symbol == "STMN2" | symbol == "RPTOR" |
           symbol == "CYFIP2" | symbol == "SYNE1" |
           symbol == "AARS" | symbol == "HDGFL2" |
           symbol == "PTPRN2") %>%
  left_join(results_table_dz_curves_0.1[,c(3,9,7)], by = "symbol",  suffix = c("_1", "_0.1")) %>%
  left_join(results_table_dz_curves_0.075[,c(3,9,7)], by = "symbol",  suffix = c("", "_0.075")) %>%
  left_join(results_table_dz_curves_0.05[,c(3,9,7)], by = "symbol", suffix = c("", "_0.05")) %>%
  left_join(results_table_dz_curves_0.04[,c(3,9,7)], by = "symbol", suffix = c("", "_0.04")) %>%
  left_join(results_table_dz_curves_0.03[,c(3,9,7)], by = "symbol", suffix = c("", "_0.03")) %>%
  left_join(results_table_dz_curves_0.02[,c(3,9,7)], by = "symbol", suffix = c("", "_0.02")) %>%
  mutate(log2FoldChange_0 = 0,
         padj_0 = 1)
names(results_table_dz_merged)[c(6,7)] <- c("log2FoldChange_0.075", "padj_0.075")

results_table_dz_merged_longer <- results_table_dz_merged %>%
  pivot_longer(cols = c(2,4,6,8,10,12,14,16)) %>%
  mutate(symbol = factor(symbol, levels = c("TARDBP", "STMN2", "UNC13A", "CYFIP2", "RPTOR", "SYNE1", "HDGFL2", "AARS", "PTPRN2"))) %>%
  mutate(name = gsub(".*_", "", name))

results_table_dz_merged_longer_tdp <- results_table_dz_merged_longer %>%
  dplyr::filter(symbol == "TARDBP")

results_table_dz_merged_longer_target <- results_table_dz_merged_longer %>%
  dplyr::filter(symbol != "TARDBP") #& symbol != "STMN2")


results_table_dz_merged_stats <- results_table_dz_merged %>%
  dplyr::filter(symbol != "TARDBP") %>%
  select(c(1,2,3)) %>%
  mutate(symbol = factor(symbol, levels = c("TARDBP", "STMN2", "UNC13A", "CYFIP2", "RPTOR", "SYNE1", "HDGFL2", "AARS"))) %>%
  #mutate(padj_1 = gsub(".*_", "", padj_1)) %>%
  mutate(label_plot = ifelse(log2FoldChange_1 > 0, "",
                             ifelse(padj_1 < 0.0001, "****",
                                    ifelse(padj_1 < 0.001, "***",
                                           ifelse(padj_1 < 0.01, "**",
                                                  ifelse(padj_1 < 0.05, "*", ""))))))
         

results_table_dz_merged_longer_target %>%
  dplyr::filter(symbol == "PTPRN2") %>%
  ggplot() +
  geom_col(data = results_table_dz_merged_longer_tdp, aes(x = name, y = value, fill = symbol), alpha = 0.2, show.legend = T) +
  geom_line(aes(x = name, y = value, color = symbol, group = symbol), size = 1.2) +
  geom_point(aes(x = name, y = value, color = symbol, group = symbol), size = 2) + #
  #geom_text(data = results_table_dz_merged_stats, aes(x = 8.3, y = log2FoldChange_1,label = label_plot)) +
  scale_color_manual(values = c("#FDAE61")) + #, "#FEE08B", "#66BD63", "#1A9850", "#006837", "#4575B4", "#313695")) +
  scale_fill_manual(values = "#AA0026") + 
  theme_classic() +
  labs(x ="Doxycycline concentration",
       y = "Log2 Fold Change",
       fill = "Gene name",
       color = "")
#ggsave("~/Desktop/dark_ce_curves.png")










for_single_gene_plot <- normed_counts_long %>%
  filter(symbol == "UNC13A" | 
           symbol == "STMN2" | symbol == "RPTOR" |
           symbol == "CYFIP2" | symbol == "SYNE1" |
           symbol == "AARS" | symbol == "HDGFL2"
  ) %>%
  group_by(symbol, group) %>%
  summarise(mean(value)) %>%
  mutate(symbol = factor(symbol, levels = c("STMN2", "UNC13A", "CYFIP2", "RPTOR", "SYNE1", "HDGFL2", "AARS")))
#stattest <- for_single_gene_plot %>%
#  group_by(symbol) %>%
#  t_test(value ~ group, ref.group = "0") %>%
#  add_significance() %>% 
#  add_xy_position(x="group")
#stattest

tardbp_levels <- normed_counts_long %>%
  dplyr::filter(symbol == "TARDBP") %>%
  group_by(group, symbol) %>%
  summarise(mean(value))

for_single_gene_plot %>%
  ggplot() + #aes(x = group, y = `mean(value)`, color = symbol, group = symbol)) +
  geom_col(data = tardbp_levels, aes(x = group, y = `mean(value)`, fill = "#AA0026"), alpha = 0.2, show.legend = F) +
  geom_line(aes(x = group, y = `mean(value)`, color = symbol, group = symbol), size = 1.2) +
  geom_point(aes(x = group, y = `mean(value)`, color = symbol, group = symbol), size = 2) + #position = position_dodge2(width = 1)) +
  #stat_pvalue_manual(stattest, hide.ns = T) + #[c(1,2),]) +
  scale_y_log10(limits = c(300,30000), oob = rescale_none) +
  scale_color_manual(values = c("#FDAE61", "#FEE08B", "#66BD63", "#1A9850", "#006837", "#4575B4", "#313695")) +
  #scale_color_brewer(palette = "RdYlGn") +
  #facet_wrap(facets = vars(symbol), scales = "free_y") +
  theme_classic() +
  labs(x ="Doxycycline concentration",
       y = "Normalised gene counts",
       color = "Gene name")
  #theme(axis.text.x = element_text(angle = 90))


#compare chx and upf1
AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), source = x)
  }))
}

normed_long_both <- AppendMe(c("normed_counts_chx_long", "normed_counts_UPF1_long")) 

normed_long_both %>%
  dplyr::filter(symbol == "TARDBP" | symbol == "UPF1" | symbol == "CYFIP2" | symbol == "INSR") %>%
  ggplot(aes(x = group, y = value, color = source)) +
  geom_boxplot(alpha = 0.1, show.legend = F) +
  geom_point(position = position_dodge2(width = 1), show.legend = F) +
  #stat_pvalue_manual(stattest, hide.ns = T) + #[c(1,2),]) +
  #scale_y_log10() +
  facet_wrap(facets = vars(symbol), scales = "free_y") +
  scale_color_brewer(palette = "Set1") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

###see if more of a general trend
results_table_ctrl_both <- left_join(results_table_cont, results_table_UPF1_ctrl, by = "symbol")

results_table_ctrl_both %>%
  ggplot(aes(x = `log2FoldChange.x`, y = `log2FoldChange.y`)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  scale_color_brewer(palette = "Set1") +
  theme_classic()



results_table_sy_merged_last <- results_table_sy5y_curves_0.075[,c(9,3,7)] %>%
 left_join(results_table_sy5y_curves_0.025[,c(3,9,7)], by = "symbol",  suffix = c("_0.075", "_0.025")) %>%
  #mutate(log2FoldChange_0 = 0,
  #       padj_0 = 1) %>%
  dplyr::filter(!is.na(`log2FoldChange_0.025`) &
                  !is.na(`log2FoldChange_0.075`)) %>%
  dplyr::filter(`log2FoldChange_0.075` > `log2FoldChange_0.025` &
                  `padj_0.025` < 0.05 &
                  `log2FoldChange_0.025` < -0.5 &
                  !is.na(symbol)) #%>%
  dplyr::mutate(big_change = ifelse(`log2FoldChange_0.075` - `log2FoldChange_0.025` > 0.5, "rescue > 0.5", 
                                    ifelse(`log2FoldChange_0.075` - `log2FoldChange_0.025` < 0, "no_rescue", "intermediate")))
#names(results_table_sy_merged)[c(6,7)] <- c("log2FoldChange_0.021", "padj_0.021")

results_table_dz_merged_last <- results_table_dz_curves_1[,c(9,3,7)] %>%
  left_join(results_table_dz_curves_0.1[,c(3,9,7)], by = "symbol",  suffix = c("_1", "_0.1")) %>%
  #left_join(results_table_dz_curves_0.075[,c(3,9,7)], by = "symbol",  suffix = c("", "_0.075")) %>%
  #mutate(log2FoldChange_0 = 0,
  #       padj_0 = 1) %>%
  dplyr::filter(!is.na(`log2FoldChange_1`) &
                                    !is.na(`log2FoldChange_0.1`)) %>%
  dplyr::filter(#`log2FoldChange_1` > `log2FoldChange_0.1` &
                  #`padj_0.1` < 0.05 &
                  #`log2FoldChange_0.1` < -0.5 &
                  !is.na(symbol))

#names(results_table_dz_merged)[c(6,7)] <- c("log2FoldChange_0.075", "padj_0.075")

results_merged_sy_dz <- left_join(results_table_sy_merged_last, results_table_dz_merged_last, by = "symbol") %>%
  dplyr::filter(!is.na(`log2FoldChange_0.025`) &
                  !is.na(`log2FoldChange_0.075`) &
                  !is.na(`log2FoldChange_1`) &
                  !is.na(`log2FoldChange_0.1`))
write.table(results_merged_sy_dz, "~/Desktop/resuced_sy_no_dz.csv", sep = ",", row.names = F)


results_table_sy_merged_last %>%
  pivot_longer(cols = c(2,4,6)) %>%
  #mutate(symbol = factor(symbol, levels = c("TARDBP", "STMN2", "UNC13A", "CYFIP2", "RPTOR", "SYNE1", "HDGFL2", "AARS"))) %>%
  mutate(name = gsub(".*_", "", name)) %>%
  ggplot() +
  geom_line(aes(x = name, y = value, group = symbol)) +
  geom_point(aes(x = name, y = value), size = 2) + #
  #scale_color_manual(values = c("#FEE08B")) + #, "#66BD63", "#1A9850", "#006837", "#4575B4", "#313695")) +
  #scale_fill_manual(values = "#AA0026") + 
  facet_wrap(vars(big_change)) +
  theme_classic() +
  labs(x ="Doxycycline concentration",
       y = "Log2 Fold Change",
       fill = "Gene name",
       color = "")

table_hs_rbp <- read.table("~/Desktop/rbp_bed/Table_HS_RBP.txt", header = T) %>%
  mutate(symbol = gsub("*_HUMAN", "", Entry_Name))

results_table_sy_merged_last_annot <- results_table_sy_merged_last %>%
  left_join(table_hs_rbp, by = "symbol")


