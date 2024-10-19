### CE compendium
library(tidyverse)


## start wide
splicing_full <- read.csv("~/Desktop/ce_bio/splicing_full_delta_psi_tables.csv") %>% ##458890
  #filter(comparison != "tdp43ctrl-tdp43upf1" & comparison != "controltdp43kd-cycloheximidetdp43kd") %>%
  #comparison != "controlliufacsneurons-tdp43kdliufacsneurons" &
  #comparison != "controlfergusonhela-tdp43kdfergusonhela") %>%
  filter(contrast_PSI > 0.1 & baseline_PSI < 0.05) #%>% #21967 #17058 #18111 
#table(splicing_full$comparison)
#table(splicing_full$junc_cat)
#reorder so that you keep only highest psi
splicing_full <- splicing_full[order(splicing_full$mean_dpsi_per_lsv_junction, decreasing = T),]


events_per_dataset <- table(splicing_full$paste_into_igv_junction) %>%
  as.data.frame() %>%
  filter(Freq > 1) #1) #to see more than one
sum(table(events_per_dataset$Freq)) #13120 -3382   #8905 - 3040   #9688 - 3078
#filter
splicing_full_unique <- splicing_full %>%
  group_by(paste_into_igv_junction, comparison) %>%
  distinct(.keep_all = T) %>% #16413 #17466 
  ungroup() %>% 
  group_by(paste_into_igv_junction) %>%
  mutate(datasets = n_distinct(comparison)) %>%
  ungroup() %>%
  #distinct(paste_into_igv_junction, .keep_all = T) %>%
  select(c(gene_name, paste_into_igv_junction, baseline_PSI, contrast_PSI, mean_dpsi_per_lsv_junction, probability_changing, datasets, junc_cat))
table(splicing_full_unique$datasets)
  
events_per_dataset_filtered <- table(splicing_full_unique$paste_into_igv_junction) %>%
  as.data.frame() %>%
  filter(Freq > 2)
sum(table(events_per_dataset_filtered$Freq)) #8905 2896   #9688 of which 2934 in more than 1 dataset

#strict conditions
splicing_strict <- splicing_full_unique %>%
  filter(mean_dpsi_per_lsv_junction > 0.2 & 
           probability_changing > 0.9)
events_per_dataset_strict <- table(splicing_strict$paste_into_igv_junction) %>%
  as.data.frame() %>%
  filter(Freq > 1)
sum(table(events_per_dataset_strict$Freq))  #1307 669 (649)  #1440 of which 657 in more than 1 dataset
#create a simplified one for merging
splicing_strict_unique <- splicing_strict %>% 
  distinct(paste_into_igv_junction, .keep_all = T)

splicing_strict_unique %>% distinct(gene_name, .keep_all = T) -> gene_table


###cross with patient specific
nygc_specific <- read.csv("~/Desktop/ce_bio/expression_by_pathology_updated.csv") %>%
  mutate(gene_name = gene,
         junc_cat = type) %>%
  select(-gene, -type) %>%
  filter(fraction_path > 0.01 & fraction_not_path < 0.005)
nygc_specific %>% distinct(gene_name) %>% pull() %>% length() #287
##using strict CE
splicing_strict_selective <- splicing_strict_unique %>% 
  left_join(nygc_specific) %>%
  filter(!is.na(strand))
splicing_strict_selective %>% distinct(gene_name) %>% pull() #71 #74 genes
#using broad list
splicing_full_selective <- splicing_full_unique %>% 
  left_join(nygc_specific) %>%
  filter(!is.na(strand))
splicing_full_selective %>% distinct(gene_name) %>% pull() #71 #74 genes
events_per_dataset_strict <- table(splicing_full_selective$paste_into_igv_junction) %>%
  as.data.frame() %>%
  filter(Freq > 2)
sum(table(events_per_dataset_strict$Freq))

View(splicing_strict_selective)
#write.csv(splicing_strict_selective, "~/Desktop/panel_patient_specific_broad.csv")

al_bed <- read.table("~/Downloads/ce_multiple_datasets.bed")
al_bed_gene <- al_bed %>%
  separate(V4, into = c("gene_name", "cat"), sep = "_") %>%
  mutate(length = V3-V2) %>%
  left_join(splicing_strict_selective) %>%
  separate(paste_into_igv_junction, into = c("chr", "range"), sep = ":") %>%
  separate(range, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end)) %>%
  mutate(end_minus_one = end-1) %>%
  mutate(condition_start = ifelse(V2 == end_minus_one, "yes", "no"),
         condition_end = ifelse(V3 == start, "yes", "no")) %>% 
  mutate(paste_coord_exon = paste0(V1, ":", V2, "-", V3)) %>%
  distinct(paste_coord_exon, .keep_all = T)
  
  group_by(gene_name) %>%
  #slice(which.min(length)) %>%
  ungroup()

check <- splicing_strict_selective %>%
  right_join(al_bed_gene)
check %>% distinct(gene_name) %>% pull()




pulling <- splicing_strict_selective %>%
  separate(paste_into_igv_junction, into = c("chr", "range"), sep = ":") %>%
  separate(range, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end)) %>%
  mutate(end_minus_one = end-1)

ciaone <- al_bed %>%
  separate(V4, into = c("gene_name", "cat"), sep = "_") %>% 
  group_by(gene_name) %>% 
  filter(n() > 1) %>%
  #mutate(ifelse(V3 %in% pulling$end, "yes", "no")) %>%
  mutate(condition_start = ifelse(V2 %in% pulling$end_minus_one, "yes", "no")) %>%
  mutate(condition_end = ifelse(V3 %in% pulling$start, "yes", "no")) #>%
  #ungroup() %>%
  #filter(condition_start == "yes" |
  #         condition_end == "yes")


for_checking <- splicing_full %>%
  separate(paste_into_igv_junction, into = c("chr", "range"), sep = ":") %>%
  separate(range, into = c("start", "end"), sep = "-") %>% 
  mutate(start = as.numeric(start),
         end = as.numeric(end)) %>%
  mutate(end_minus_one = end-1)
  
splicing_full_al <- for_checking %>%
  mutate(ifelse(start %in% al_bed$V3, "yes", "no")) %>%
  mutate(ifelse(end_minus_one %in% al_bed$V2, "yes", "no"))



splicing_full_selective %>%
  select(2) %>%
  separate(paste_into_igv_junction, into = c("chr", "range"), sep = ":") %>%
  separate(range, into = c("start", "end"), sep = "-") -> bed_selective

write.table(bed_selective, "~/Desktop/selective2.bed", quote = F, sep = "\t", row.names = F)



###cross also with curves
#curves_sy <- read.csv("~/Documents/phd/research_lines/tdp-43 curves/manual_validation_curves/shsy_validated_events.ann.csv") %>% 
#  mutate(gene_name = gene_id,
#         paste_into_igv_junction = paste0(chr, ":", start, "-", end)) %>%
#  select(-gene_id, -junc_cat) %>%
#  mutate(progression = ifelse())
#where is namual?
namuel <- namual %>%
  select(-junc_cat)

splicing_strict_curves <- splicing_strict_unique %>% 
  full_join(namuel) %>%
  filter(!is.na(Type))
splicing_strict_curves %>% distinct(gene_name) %>% pull()

splicing_full_curves <- splicing_full_unique %>% 
  full_join(namuel) %>%
  filter(!is.na(Type))
splicing_full_curves %>% distinct(gene_name) %>% pull()


splicing_full_curves_selective <- splicing_full_selective %>% 
  full_join(splicing_full_curves)

###cross with nmd data
nmd_data <- read.csv("~/Desktop/ce_bio/nmd_chx.csv")



