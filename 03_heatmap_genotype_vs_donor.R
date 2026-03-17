# Genotype heatmaps: main analysis compares each genotype to donor, separately for LV and RV.
# Heatmaps are split into PDE-only and mechanosensitive-only panels.
# Mechanosensitive genes are manually ordered by Type from the Excel sheet (no clustering).
# Optional supplementary heatmaps compare genotype vs PVneg.

MAKE_PVNEG_SUPP_HEATMAP <- TRUE

read_pb_by_celltype <- function(celltype) {
  pb <- read_pseudobulk(celltype)
  validate_pb(pb, celltype)
  pb$meta <- add_disease2(pb$meta)
  pb$meta$Primary.Genetic.Diagnosis <- clean_genotype(pb$meta$Primary.Genetic.Diagnosis)
  pb
}

run_genotype_vs_donor <- function(pb, genotype, region) {
  keep <- pb$meta$Region_x == region & (
    pb$meta$disease == "normal" |
      (pb$meta$disease == "dilated cardiomyopathy" & pb$meta$Primary.Genetic.Diagnosis == genotype)
  )
  pb2 <- subset_pb(pb, keep)
  pb2$meta$geno2 <- ifelse(pb2$meta$disease == "normal", "donor", genotype)
  run_two_group_contrast(pb2, rep(TRUE, nrow(pb2$meta)), "geno2", c("donor", genotype), paste0("geno2", genotype))
}

run_genotype_vs_pvneg <- function(pb, genotype, region) {
  keep <- pb$meta$Region_x == region & pb$meta$disease == "dilated cardiomyopathy" &
    pb$meta$Primary.Genetic.Diagnosis %in% c("PVneg", genotype)
  pb2 <- subset_pb(pb, keep)
  pb2$meta$geno2 <- factor(pb2$meta$Primary.Genetic.Diagnosis, levels = c("PVneg", genotype))
  run_two_group_contrast(pb2, rep(TRUE, nrow(pb2$meta)), "geno2", c("PVneg", genotype), paste0("geno2", genotype))
}

add_mech_type <- function(df) {
  df$type <- unname(MECH_TYPE_MAP[as.character(df$symbol)])
  df
}

plot_heatmap_tile_mech <- function(df, title) {
  df$symbol <- factor(df$symbol, levels = unique(df$symbol))
  contrast_levels <- unique(as.character(df$contrast))
  df$contrast_num <- match(as.character(df$contrast), contrast_levels)

  type_breaks <- aggregate(as.numeric(df$symbol), by = list(type = df$type), FUN = min)
  type_breaks <- type_breaks[order(type_breaks$x), , drop = FALSE]

  annot_df <- unique(df[, c("symbol", "type")])
  annot_df$annot_x <- 0

  type_levels <- unique(as.character(annot_df$type))
  type_levels <- type_levels[!is.na(type_levels) & type_levels != ""]
  type_cols <- setNames(grDevices::hcl.colors(length(type_levels), palette = "Set 3"), type_levels)

  ggplot() +
    geom_tile(data = df, aes(x = contrast_num, y = symbol, fill = logFC), color = "white") +
    geom_text(data = df, aes(x = contrast_num, y = symbol, label = star), size = 3) +
    geom_point(data = annot_df,
               aes(x = annot_x, y = symbol, color = type),
               shape = 15, size = 4.2, show.legend = TRUE) +
    scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35", midpoint = 0,
                         na.value = "grey90") +
    scale_color_manual(values = type_cols, na.translate = FALSE) +
    scale_x_continuous(breaks = c(0, seq_along(contrast_levels)),
                       labels = c("Type", contrast_levels),
                       expand = expansion(mult = c(0.02, 0.02))) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          legend.box = "vertical") +
    labs(title = title, x = NULL, y = NULL, fill = "logFC", color = "Type") +
    geom_hline(yintercept = type_breaks$x - 0.5, color = "grey70", linewidth = 0.5)
}

write_split_heatmaps <- function(results_named_list, contrast_order, ct, region, cell_label, suffix) {
  # PDE panel
  pde_df <- build_heatmap_df(results_named_list, PDE_GENES, contrast_order)
  if (!is.null(pde_df) && nrow(pde_df) > 0) {
    write.csv(pde_df, paste0(ct, "_", region, "_heatmap_PDE_", suffix, ".csv"), row.names = FALSE)
    p_pde <- plot_heatmap_tile(pde_df, title = paste0(cell_label, " ", region, ": PDE (", suffix, ")"))
    save_plot_png(paste0(ct, "_", region, "_heatmap_PDE_", suffix, ".png"), p_pde, 6.5, 6.5)
  }

  # Mechanosensitive panel
  mech_df <- build_heatmap_df(results_named_list, MECH_GENES_ORDERED, contrast_order)
  if (!is.null(mech_df) && nrow(mech_df) > 0) {
    mech_df <- add_mech_type(mech_df)
    write.csv(mech_df, paste0(ct, "_", region, "_heatmap_mechanosensitive_", suffix, ".csv"), row.names = FALSE)
    p_mech <- plot_heatmap_tile_mech(mech_df, title = paste0(cell_label, " ", region, ": mechanosensitive proteins (", suffix, ")"))
    save_plot_png(paste0(ct, "_", region, "_heatmap_mechanosensitive_", suffix, ".png"), p_mech, 8.5, 12)
  }
}

for (ct in CELLTYPES) {
  pb <- read_pb_by_celltype(ct)
  cell_label <- ifelse(ct == "CM", "Cardiomyocytes", "Endothelial cells")

  for (region in c("LV", "RV")) {
    donor_results <- list()
    for (g in GENOTYPE_GROUPS) {
      donor_results[[g]] <- run_genotype_vs_donor(pb, g, region)
      if (!is.null(donor_results[[g]])) {
        write.csv(donor_results[[g]]$deg,
                  paste0(ct, "_", region, "_", g, "_vs_donor_edgeR_all.csv"),
                  row.names = FALSE)
      }
    }
    write_split_heatmaps(donor_results, GENOTYPE_GROUPS, ct, region, cell_label, "genotype_vs_donor")

    if (MAKE_PVNEG_SUPP_HEATMAP) {
      pv_results <- list()
      pv_order <- setdiff(GENOTYPE_GROUPS, "PVneg")
      for (g in pv_order) {
        pv_results[[g]] <- run_genotype_vs_pvneg(pb, g, region)
        if (!is.null(pv_results[[g]])) {
          write.csv(pv_results[[g]]$deg,
                    paste0(ct, "_", region, "_", g, "_vs_PVneg_edgeR_all.csv"),
                    row.names = FALSE)
        }
      }
      write_split_heatmaps(pv_results, pv_order, ct, region, cell_label, "genotype_vs_PVneg")
    }
  }
}

message("03_heatmap_genotype_vs_donor.R done")
