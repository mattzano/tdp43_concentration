#roundel plot
manual_validation_dz <- read.table("~/Downloads/validation_curves_dz.csv", sep = ",", header = T)
table(manual_validation_dz$cryptic)

manual_validation_sy <- read_xlsx("~/Documents/phd/research_lines/tdp-43 curves/manual_validation_curves/validation_curves.xlsx")
  #mutate(paste_into_igv_junction = paste0(chr, ":", start, "-", end)) %>%
df <- manual_validation_sy %>% 
  mutate(appearance = Type) %>%
  filter(validation == "yes" & appearance != "0") %>%# & cryptic == "yes")
  mutate(cryptic_clean = ifelse(is.na(cryptic), "no", 
                                ifelse(cryptic == "no", "no", "yes"))) %>%
  select(c("appearance", "cryptic_clean", "junc_cat")) %>%
  arrange(appearance)


ciao1 <- as.data.frame(table(df$appearance, df$cryptic_clean, df$junc_cat)) 
ciao2 <- as.data.frame(table(df$appearance, df$cryptic_clean)) %>%
  mutate(Var3 = NA)
ciao3 <- as.data.frame(table(df$appearance)) %>%
  mutate(Var3 = NA,
         Var2 = NA)

ciao <- rbind(ciao1, ciao2, ciao3) %>%
  mutate(value = ifelse(is.na(Var3) & is.na(Var2), as.character(Var1),
                        ifelse(is.na(Var3), as.character(Var2), as.character(Var3)))) %>%
  mutate(ymid = as.numeric(ifelse(is.na(Var3) & is.na(Var2), 1,
                                  ifelse(is.na(Var3), 2, 3))),
         ymax = ymid + 0.5, 
         ymin = ymid - 0.5) %>%
  #arrange(Var1) %>%
  group_by(ymid) %>%
  arrange(Var1, Var2) %>%
  mutate(xmin = c(0, head(cumsum(Freq), -1)),
         xmax = cumsum(Freq),
         xmid = (xmax + xmin) / 2)
ciao %>%
  mutate(value_label = paste0(Var1,Var2,Var3)) %>%
    #group_by(ymid) %>%
  ggplot(aes(xmid, ymid, fill = Var1)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  alpha = as.character(ymid), 
                  color = Var1)) +
  #geom_text(aes(y = ymid + 0.25, label = value_label), label = ciao$value) +
  #geomtextpath::geom_textpath(aes(y = ymid + 0.25, label = c(value), group = value)) +
  geomtextpath::geom_textpath(aes(y = ymid + 0.25, label = value, group = value_label)) +
  scale_alpha_manual(values = c(1, 0.5, 0.2)) +
  scale_fill_manual(values = c("#cd9900", "#00817e", "darkred")) +
  scale_colour_manual(values = c("#cd9900", "#00817e", "darkred")) +
  scale_y_continuous(limits = c(-0.5, 4)) +
    #geomtextpath::coord_curvedpolar(theta = "x") +
  coord_polar() +
  theme_void() +
  theme(legend.position = "none")





#####


ciao <- df %>%
  #pivot_longer(1:3) %>%
  #group_by(name, value) %>%
  mutate(width = n()) %>%
  unique() %>%
  arrange(name) %>%
  group_by(name) %>%
  mutate(ymid = as.numeric(case_when(name == "appearance" ~ 1,
                                    name == "cryptic_clean" ~ 2,
                                    name == "junc_cat" ~ 3)),
         ymax = ymid + 0.5, 
         ymin = ymid - 0.5,
         xmin = c(0, head(cumsum(width), -1)),
         xmax = cumsum(width),
         xmid = (xmax + xmin) / 2) %>%
  #mutate(width2 = ifelse(ymid == 3, width/sum(width),
  #                       ifelse( ymid == 2, width/sum(width), width)))
  ggplot(aes(xmid, ymid, fill = top_level)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                alpha = name, color = top_level)) +
  geom_text(aes(y = ymid + 0.25, label = value, group = value)) +
  #geomtextpath::geom_textpath(aes(y = ymid + 0.25, label = value, 
  #                                group = c(name))) +
  scale_alpha_manual(values = c(1, 0.5, 0.2)) +
  scale_fill_manual(values = c("#cd9900", "#00817e", "darkred")) +
  scale_colour_manual(values = c("#cd9900", "#00817e", "darkred")) +
  scale_y_continuous(limits = c(-0.5, 4)) +
  coord_polar() +
  theme_void() #+
  theme(legend.position = "none")

