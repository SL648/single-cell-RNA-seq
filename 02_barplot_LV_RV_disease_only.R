# Barplot focused on disease contrasts only.
# Keeps LV DCM vs donor and RV DCM vs donor.
# Region contrasts (donor LV vs RV, DCM LV vs RV) are intentionally excluded from the main barplot.

read_pb_by_celltype <- function(celltype) {
  pb <- read_pseudobulk(celltype)
  validate_pb(pb, celltype)
  pb$meta <- add_disease2(pb$meta)
  pb
}

collect_region_disease_deg <- function(celltype, region) {
  pb <- read_pb_by_celltype(celltype)
  keep <- pb$meta$Region_x == region & pb$meta$disease2 %in% c("donor", "DCM")
  pb2 <- subset_pb(pb, keep)
  res <- run_two_group_contrast(pb2, rep(TRUE, nrow(pb2$meta)), "disease2", c("donor", "DCM"), "disease2DCM")
  if (is.null(res)) return(NULL)

  df <- prep_gene_df(res$deg, PDE_GENES,
                     cell_type = ifelse(celltype == "CM", "Cardiomyocytes", "Endothelial cells"),
                     contrast = paste0(region, "_DCM_vs_donor"))
  df$region <- region
  df
}

df_list <- list(
  collect_region_disease_deg("CM", "LV"),
  collect_region_disease_deg("CM", "RV"),
  collect_region_disease_deg("EC", "LV"),
  collect_region_disease_deg("EC", "RV")
)

bar_df <- do.call(rbind, df_list)
write.csv(bar_df, "Fig_barplot_PDE_LV_RV_DCM_vs_donor.csv", row.names = FALSE)

p_bar <- plot_bar_disease_by_region(
  bar_df,
  title = "PDE differential expression by region (DCM vs donor)"
)

save_plot_png("Fig_barplot_PDE_LV_RV_DCM_vs_donor.png", p_bar, 8, 6)

message("Barplot colors: LV = #4DBBD5, RV = #E64B35")
message("02_barplot_LV_RV_disease_only.R done")
