

genedata <- read.csv("data/genefinder_ELAVL3.csv")
#genedata$plotby <- c(2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,1,1,1)
genedata$level <- as.factor(genedata$level)

plot_elavl3 <- ggplot(genedata,aes_string(x="level", y="count", fill="level")) +
  geom_boxplot(aes_string(group = "level"), outlier.shape = NA,alpha = 0.7) + theme_bw() +
  scale_y_log10(name="Normalized counts - log10 scale") +
  labs(title=paste0("Normalized counts for ELAVL3")) +
  scale_x_discrete(name = "Doxycycline concentration") +
  geom_jitter(aes_string(group = "level"),position = position_jitter(width = 0.0001)) +
  scale_fill_discrete(name = "Experimental\nconditions")
plot(plot_elavl3)
ggsave("elavl3_counts.pdf", plot_celf5)

