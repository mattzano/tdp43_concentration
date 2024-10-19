#write.csv(big_delta_de, "~/Desktop/big_delta_de.csv")
#write.csv(big_delta_de_dz, "~/Desktop/big_delta_de_dz.csv")

for_plot <- big_delta_de[,c(3,7,9,10)] %>%
  dplyr::filter(symbol == "STMN2" | symbol == "AARS1" | symbol == "HDGFL2" | symbol == "UNC13A" |
                  symbol == "SYNE1" | symbol == "CYFIP2")
for_plot_0 <- for_plot %>% 
  distinct(symbol, .keep_all = T) %>% 
  mutate(log2FoldChange = 0, padj = 0, source = 0)
for_plot <- rbind(for_plot_0,for_plot) %>% 
  mutate(symbol = factor(symbol, levels = c("AARS1", "HDGFL2", "SYNE1", "CYFIP2", "UNC13A", "STMN2")))

a <- for_plot %>%   
  ggplot() +
  geom_col(aes(x = source, y = log2FoldChange, fill = symbol), data = big_delta_de_tdp, alpha = 0.3, show.legend = F) +
  geom_point(aes(x = source, y = log2FoldChange, color = symbol, group = symbol), show.legend = F) +
  geom_line(aes(x = source, y = log2FoldChange, color = symbol, group = symbol), show.legend = F) +
  scale_fill_manual(values = "black") +
  scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#FF7F00", "#FDBF6F", "#FB9A99", "#E31A1C"))+ # "#B2DF8A" "#33A02C" "#CAB2D6" "#6A3D9A")) +
  scale_y_continuous(limits = c(-6.2,0.6)) +
  labs(x = "Doxycycline concentration", y = "Log2 Fold Change compared to control", color = "", fill = "Gene") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
ggsave("~/Desktop/curves_nmd.png", width = 6, height = 4)


big_delta_de_dz_tdp <- big_delta_de_dz[,c(3,7,9,10)] %>% 
  dplyr::filter(symbol == "TARDBP")
for_plot_dz_0_tdp <- big_delta_de_dz_tdp[1,] %>% 
  mutate(log2FoldChange = 0, padj = 0, source = 0)
big_delta_de_dz_tdp <- rbind(for_plot_dz_0_tdp,big_delta_de_dz_tdp)

for_plot_dz <- big_delta_de_dz[,c(3,7,9,10)] %>%
  dplyr::filter(symbol == "STMN2" | symbol == "AARS1" | symbol == "HDGFL2" | symbol == "UNC13A" |
                  symbol == "SYNE1" | symbol == "CYFIP2")
for_plot_dz_0 <- for_plot_dz %>% 
  distinct(symbol, .keep_all = T) %>% 
  mutate(log2FoldChange = 0, padj = 0, source = 0)
for_plot_dz <- rbind(for_plot_dz_0,for_plot_dz) %>% 
  mutate(symbol = factor(symbol, levels = c("AARS1", "HDGFL2", "SYNE1", "CYFIP2", "UNC13A", "STMN2")))

b <- for_plot_dz %>%   
  ggplot() +
  geom_col(aes(x = source, y = log2FoldChange, fill = symbol), data = big_delta_de_dz_tdp, alpha = 0.3) +
  geom_point(aes(x = source, y = log2FoldChange, color = symbol, group = symbol)) +
  geom_line(aes(x = source, y = log2FoldChange, color = symbol, group = symbol)) +
  scale_fill_manual(values = "black") +
  scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#FF7F00", "#FDBF6F", "#FB9A99", "#E31A1C"))+ # "#B2DF8A" "#33A02C" "#CAB2D6" "#6A3D9A")) +
  scale_y_continuous(limits = c(-5.2,0.5)) +
  labs(x = "Doxycycline concentration", y = "Log2 Fold Change compared to control", color = "", fill = "Gene") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
ggsave("~/Desktop/curves_nmd_dz.png", width = 6, height = 4)

library(patchwork)
a+b -> caio
ggsave("~/Desktop/caio.png", caio, width = 8, height = 4)
