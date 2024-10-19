library(tidyverse)


##load sy5y annot

data_puja <- data.frame(c(rep("STMN2",6), rep("INSR",6), rep("UNC13A",6)),
                     rep(c("-", "+", "++", "+++", "++++", "+++++"),3),
                     c("no",rep("yes",5),"no","no",rep("yes",4),"no","no",rep("yes",4)))
names(data_puja) <- c("Gene", "TDP-43 knockdown levels", "presence")
data_puja$Gene <- factor(data_puja$Gene, levels = c("STMN2", "INSR", "UNC13A"))

hello <- data_puja %>% 
  ggplot(aes(x = `TDP-43 knockdown levels`, y = Gene, color = presence, group = Gene)) +
  geom_line(#linewidth = 1.5, 
    show.legend = F) +
  geom_point(#size = 6, 
    show.legend = F) +
  scale_color_manual(values = c("white", "black")) +
  scale_y_discrete(limits=rev) +
  theme_classic()
print(hello)


list_gene <- c("STMN2", "INSR", "UNC13A")
plot_early_late <- sy5y_full_annot %>%
  #group_by(gene_id) %>% 
  arrange(desc(dox075)) %>% 
  distinct(gene_id, .keep_all = T) %>% 
  filter(cryptic == "yes") %>% 
  pivot_longer(cols = c(8:13)) %>% 
  mutate(label_junction = case_when(name == "dox025" &
                                      gene_id %in% list_gene
                                    ~ gene_id, T ~ "")) %>%
  mutate(color_gene_name = ifelse(gene_id %in% list_gene, "2", "1")) %>%
  mutate(alpha_gene_name = ifelse(color_gene_name == "2", "2", "1")) %>% 
  ggplot(mapping = aes(x = name, y = value, group = paste_into_igv_junction)) +
  geom_line(aes(color = color_gene_name, alpha = alpha_gene_name), show.legend = F) +
  geom_text_repel(aes(label = label_junction), point.padding = 0.2,
                  nudge_y = 0.2, min.segment.length = 0.2, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  geom_hline(yintercept = 0.1, linetype = "dotted") +
  scale_color_manual(values = c("gray", "#3288BD", "#D53E4F")) +
  scale_alpha_manual(values = c(0.1,1)) +
  xlab("") +
  ylab("PSI") +
  scale_x_discrete(labels = c("-", "+", "++", "+++", "++++", "+++++")) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  #  labels = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_classic()
print(plot_early_late)

heya <- plot_early_late / hello + patchwork::plot_layout(heights = c(3, 1))
ggsave("/Users/matteozanovello/Desktop/curves_puja.png", heya, width = 12)
          
           
