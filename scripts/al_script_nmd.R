affected_by_chx = merged_nmd_ce |> 
  select(gene_name.x,paste_into_igv_junction,Control_Control:Cycloheximide_TDP43KD,color_gene_name.x) |> 
  filter(Control_Control < 0.05 & 
           color_gene_name.x == "Delta PSI > 0.05") |> 
  pull(paste_into_igv_junction)

unaffected_by_chx = merged_nmd_ce |> 
  select(gene_name.x,paste_into_igv_junction,Control_Control:Cycloheximide_TDP43KD,color_gene_name.x) |> 
  filter(Control_Control < 0.05 & 
           color_gene_name.x == "Delta PSI < 0.05") |> 
  pull(paste_into_igv_junction)

affected_by_upf1 = merged_nmd_ce |> 
  select(gene_name.y,paste_into_igv_junction,ctrl_ctrl:TDP43_UPF1,color_gene_name.y) |> 
  filter(ctrl_ctrl < 0.05 & 
           color_gene_name.y == "Delta PSI > 0.05") |> 
  pull(paste_into_igv_junction)

unaffected_by_upf1 = merged_nmd_ce |> 
  select(gene_name.y,paste_into_igv_junction,ctrl_ctrl:TDP43_UPF1,color_gene_name.y) |> 
  filter(ctrl_ctrl < 0.05 & 
           color_gene_name.y == "Delta PSI < 0.05") |> 
  pull(paste_into_igv_junction)




late_sushi = fread("data/shsy5y_curve_cate.csv") |> mutate(paste_into_igv_junction = glue::glue("{chr}:{start}-{end}")) |> janitor::clean_names()





# Create contingency table
upf_tbl = late_sushi[is_cryptic == 'Cryptic'] |> 
  filter(paste_into_igv_junction %in% c(affected_by_upf1,unaffected_by_upf1)) |> 
  mutate(nmd_effect = case_when(paste_into_igv_junction %in% affected_by_upf1 ~ "increased - UPF1 KD", 
                                paste_into_igv_junction %in% unaffected_by_upf1 ~ "unaffected - UPF1 KD",
                                T ~'oops')) |> 
  select(nmd_effect,type) |> table()

chx_tbl = late_sushi[is_cryptic == 'Cryptic'] |> 
  filter(paste_into_igv_junction %in% c(affected_by_chx,unaffected_by_chx)) |> 
  mutate(nmd_effect = case_when(paste_into_igv_junction %in% affected_by_chx ~ "increased - CHX", 
                                paste_into_igv_junction %in% unaffected_by_chx ~ "unaffected - CHX",
                                T ~'oops')) |> 
  select(type,nmd_effect) |> table()








