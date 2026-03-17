# Main volcano script
# Outputs PNG only.
# Notes:
# - Region/disease volcanoes are generated for donor LV vs RV, DCM LV vs RV,
#   LV DCM vs donor, and RV DCM vs donor.
# - Genotype volcanoes are generated vs donor (LV only by default), including PVneg vs donor.
# - Optional supplementary genotype vs PVneg volcanoes can also be turned on.

GENOTYPE_VOLCANO_REGION <- "LV"
MAKE_GENOTYPE_VS_PVNEG_SUPP <- FALSE

read_pb_by_celltype <- function(celltype) {
  pb <- read_pseudobulk(celltype)
  validate_pb(pb, celltype)
  pb$meta <- add_disease2(pb$meta)
  pb$meta$Primary.Genetic.Diagnosis <- clean_genotype(pb$meta$Primary.Genetic.Diagnosis)
  pb
}

run_region_and_disease_volcanoes <- function(celltype) {
  pb <- read_pb_by_celltype(celltype)

  # donor: LV vs RV
  keep <- pb$meta$disease2 == "donor" & pb$meta$Region_x %in% c("LV", "RV")
  pb1 <- subset_pb(pb, keep)
  pb1$meta$region2 <- factor(pb1$meta$Region_x, levels = c("RV", "LV"))
  res <- run_two_group_contrast(pb1, rep(TRUE, nrow(pb1$meta)), "region2", c("RV", "LV"), "region2LV")
  if (!is.null(res)) {
    write.csv(res$deg, paste0(celltype, "_donor_LV_vs_RV_edgeR_all.csv"), row.names = FALSE)
    p <- make_volcano(res$deg,
      title = paste0(celltype, ": donor LV vs RV"),
      pde_genes = PDE_GENES, mech_genes = MECH_GENES)
    save_plot_png(paste0(celltype, "_donor_LV_vs_RV_volcano.png"), p, 7, 5)
  }

  # DCM: LV vs RV
  keep <- pb$meta$disease2 == "DCM" & pb$meta$Region_x %in% c("LV", "RV")
  pb2 <- subset_pb(pb, keep)
  pb2$meta$region2 <- factor(pb2$meta$Region_x, levels = c("RV", "LV"))
  res <- run_two_group_contrast(pb2, rep(TRUE, nrow(pb2$meta)), "region2", c("RV", "LV"), "region2LV")
  if (!is.null(res)) {
    write.csv(res$deg, paste0(celltype, "_DCM_LV_vs_RV_edgeR_all.csv"), row.names = FALSE)
    p <- make_volcano(res$deg,
      title = paste0(celltype, ": DCM LV vs RV"),
      pde_genes = PDE_GENES, mech_genes = MECH_GENES)
    save_plot_png(paste0(celltype, "_DCM_LV_vs_RV_volcano.png"), p, 7, 5)
  }

  # LV: DCM vs donor
  keep <- pb$meta$Region_x == "LV" & pb$meta$disease2 %in% c("donor", "DCM")
  pb3 <- subset_pb(pb, keep)
  res <- run_two_group_contrast(pb3, rep(TRUE, nrow(pb3$meta)), "disease2", c("donor", "DCM"), "disease2DCM")
  if (!is.null(res)) {
    write.csv(res$deg, paste0(celltype, "_LV_DCM_vs_donor_edgeR_all.csv"), row.names = FALSE)
    p <- make_volcano(res$deg,
      title = paste0(celltype, ": LV DCM vs donor"),
      pde_genes = PDE_GENES, mech_genes = MECH_GENES)
    save_plot_png(paste0(celltype, "_LV_DCM_vs_donor_volcano.png"), p, 7, 5)
  }

  # RV: DCM vs donor
  keep <- pb$meta$Region_x == "RV" & pb$meta$disease2 %in% c("donor", "DCM")
  pb4 <- subset_pb(pb, keep)
  res <- run_two_group_contrast(pb4, rep(TRUE, nrow(pb4$meta)), "disease2", c("donor", "DCM"), "disease2DCM")
  if (!is.null(res)) {
    write.csv(res$deg, paste0(celltype, "_RV_DCM_vs_donor_edgeR_all.csv"), row.names = FALSE)
    p <- make_volcano(res$deg,
      title = paste0(celltype, ": RV DCM vs donor"),
      pde_genes = PDE_GENES, mech_genes = MECH_GENES)
    save_plot_png(paste0(celltype, "_RV_DCM_vs_donor_volcano.png"), p, 7, 5)
  }
}

run_genotype_vs_donor_volcanoes <- function(celltype, region = "LV") {
  pb <- read_pb_by_celltype(celltype)
  for (g in GENOTYPE_GROUPS) {
    keep <- pb$meta$Region_x == region & (
      pb$meta$disease == "normal" |
      (pb$meta$disease == "dilated cardiomyopathy" & pb$meta$Primary.Genetic.Diagnosis == g)
    )
    pb2 <- subset_pb(pb, keep)
    pb2$meta$geno2 <- ifelse(pb2$meta$disease == "normal", "donor", g)
    res <- run_two_group_contrast(pb2, rep(TRUE, nrow(pb2$meta)), "geno2", c("donor", g), paste0("geno2", g))
    if (is.null(res)) next

    write.csv(res$deg, paste0(celltype, "_", region, "_", g, "_vs_donor_edgeR_all.csv"), row.names = FALSE)
    p <- make_volcano(res$deg,
      title = paste0(celltype, ": ", region, " ", g, " vs donor"),
      pde_genes = PDE_GENES, mech_genes = MECH_GENES)
    save_plot_png(paste0(celltype, "_", region, "_", g, "_vs_donor_volcano.png"), p, 7, 5)
  }
}

run_genotype_vs_pvneg_supp <- function(celltype, region = "LV") {
  pb <- read_pb_by_celltype(celltype)
  for (g in setdiff(GENOTYPE_GROUPS, "PVneg")) {
    keep <- pb$meta$Region_x == region & pb$meta$disease == "dilated cardiomyopathy" &
      pb$meta$Primary.Genetic.Diagnosis %in% c("PVneg", g)
    pb2 <- subset_pb(pb, keep)
    pb2$meta$geno2 <- factor(pb2$meta$Primary.Genetic.Diagnosis, levels = c("PVneg", g))
    res <- run_two_group_contrast(pb2, rep(TRUE, nrow(pb2$meta)), "geno2", c("PVneg", g), paste0("geno2", g))
    if (is.null(res)) next

    write.csv(res$deg, paste0(celltype, "_", region, "_", g, "_vs_PVneg_edgeR_all.csv"), row.names = FALSE)
    p <- make_volcano(res$deg,
      title = paste0(celltype, ": ", region, " ", g, " vs PVneg (supp)"),
      pde_genes = PDE_GENES, mech_genes = MECH_GENES)
    save_plot_png(paste0(celltype, "_", region, "_", g, "_vs_PVneg_volcano.png"), p, 7, 5)
  }
}

for (ct in CELLTYPES) {
  message("Running volcano analyses for ", ct)
  run_region_and_disease_volcanoes(ct)
  run_genotype_vs_donor_volcanoes(ct, region = GENOTYPE_VOLCANO_REGION)
  if (MAKE_GENOTYPE_VS_PVNEG_SUPP) run_genotype_vs_pvneg_supp(ct, region = GENOTYPE_VOLCANO_REGION)
}

message("01_volcano_LV_RV_plus_genotype.R done")
