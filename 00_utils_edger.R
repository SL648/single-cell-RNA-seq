suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(ggrepel)
  library(readxl)
})

read_pseudobulk <- function(prefix) {
  mtx <- readMM(paste0(prefix, ".counts.mtx"))
  genes <- read.delim(paste0(prefix, ".genes.tsv"), header = FALSE, stringsAsFactors = FALSE)
  samples <- read.delim(paste0(prefix, ".samples.tsv"), header = FALSE, stringsAsFactors = FALSE)
  meta <- read.csv(paste0(prefix, ".meta.csv"), row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

  rownames(mtx) <- genes$V1
  colnames(mtx) <- samples$V1
  meta <- meta[colnames(mtx), , drop = FALSE]
  stopifnot(all(colnames(mtx) == rownames(meta)))

  list(counts = mtx, meta = meta)
}

validate_pb <- function(pb, label = deparse(substitute(pb))) {
  message("Checking ", label, " ...")
  stopifnot(all(colnames(pb$counts) == rownames(pb$meta)))
  invisible(TRUE)
}

subset_pb <- function(pb, keep) {
  keep <- as.logical(keep)
  keep[is.na(keep)] <- FALSE

  meta2 <- pb$meta[keep, , drop = FALSE]
  idx <- match(rownames(meta2), colnames(pb$counts))

  if (any(is.na(idx))) {
    bad <- head(rownames(meta2)[is.na(idx)], 10)
    stop("Some meta rownames are not found in counts colnames: ", paste(bad, collapse = ", "))
  }

  counts2 <- pb$counts[, idx, drop = FALSE]
  colnames(counts2) <- colnames(pb$counts)[idx]
  stopifnot(all(colnames(counts2) == rownames(meta2)))

  list(counts = counts2, meta = meta2)
}

add_symbol_from_ensembl <- function(df, ensembl_col = "ensembl") {
  df$ensembl_clean <- sub("\\..*$", "", df[[ensembl_col]])
  mp <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(df$ensembl_clean),
    keytype = "ENSEMBL",
    columns = c("SYMBOL")
  )
  # Keep first mapping per ENSEMBL to avoid accidental length inflation / ambiguity.
  mp <- mp[!duplicated(mp$ENSEMBL), , drop = FALSE]
  df$symbol <- mp$SYMBOL[match(df$ensembl_clean, mp$ENSEMBL)]
  df
}

add_disease2 <- function(meta) {
  meta$disease2 <- NA_character_
  meta$disease2[meta$disease == "normal"] <- "donor"
  meta$disease2[meta$disease == "dilated cardiomyopathy"] <- "DCM"
  meta
}

clean_genotype <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "<NA>")] <- NA
  x
}

p_to_star <- function(fdr) {
  out <- rep("", length(fdr))
  out[!is.na(fdr) & fdr < 0.05] <- "*"
  out[!is.na(fdr) & fdr < 0.01] <- "**"
  out[!is.na(fdr) & fdr < 0.001] <- "***"
  out
}

save_plot_png <- function(filename, plot, width, height, dpi = 300) {
  ggsave(filename, plot, width = width, height = height, dpi = dpi)
}

run_edger_generic <- function(counts, meta, group_col, group_levels,
                              covariates = c("sex"), coef_name = NULL) {
  stopifnot(all(colnames(counts) == rownames(meta)))

  meta <- as.data.frame(meta, stringsAsFactors = FALSE)

  keep_samples <- !is.na(meta[[group_col]]) & meta[[group_col]] %in% group_levels
  present_covariates <- covariates[covariates %in% colnames(meta)]
  for (cv in present_covariates) {
    keep_samples <- keep_samples & !is.na(meta[[cv]])
  }

  meta <- meta[keep_samples, , drop = FALSE]
  counts <- counts[, rownames(meta), drop = FALSE]
  stopifnot(all(colnames(counts) == rownames(meta)))

  meta[[group_col]] <- factor(meta[[group_col]], levels = group_levels)
  present_groups <- levels(droplevels(meta[[group_col]]))
  if (length(present_groups) < 2) {
    stop("Need at least two groups in ", group_col)
  }

  usable_covariates <- character(0)
  dropped_covariates <- setdiff(covariates, present_covariates)

  for (cv in present_covariates) {
    vals <- meta[[cv]]
    vals2 <- vals[!is.na(vals)]
    if (length(unique(vals2)) < 2) {
      dropped_covariates <- c(dropped_covariates, cv)
      next
    }
    meta[[cv]] <- factor(meta[[cv]])
    usable_covariates <- c(usable_covariates, cv)
  }

  build_design <- function(use_covariates) {
    rhs <- c(group_col, use_covariates)
    fml <- as.formula(paste("~", paste(rhs, collapse = " + ")))
    model.matrix(fml, data = meta)
  }

  design <- build_design(usable_covariates)

  # If any covariate causes rank deficiency (common in tiny genotype subsets),
  # fall back to the main group-only model instead of failing the whole script.
  if (qr(design)$rank < ncol(design) && length(usable_covariates) > 0) {
    message("Design not full rank for ", group_col,
            "; retrying without covariates: ",
            paste(usable_covariates, collapse = ", "))
    dropped_covariates <- unique(c(dropped_covariates, usable_covariates))
    usable_covariates <- character(0)
    design <- build_design(usable_covariates)
  }

  if (qr(design)$rank < ncol(design)) {
    stop("Design matrix still not full rank even after dropping covariates.")
  }

  if (length(dropped_covariates) > 0) {
    message("Dropped covariates: ", paste(unique(dropped_covariates), collapse = ", "))
  }
  if (length(usable_covariates) > 0) {
    message("Using covariates: ", paste(usable_covariates, collapse = ", "))
  } else {
    message("Using no covariates; model is ~ ", group_col)
  }

  y <- DGEList(counts = as.matrix(counts))
  keep_genes <- filterByExpr(y, group = meta[[group_col]])
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  design_cols <- colnames(design)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)

  if (is.null(coef_name)) {
    if (ncol(design) < 2) stop("Design has fewer than 2 columns; cannot use default coef = 2")
    res <- glmQLFTest(fit, coef = 2)
  } else {
    idx <- which(design_cols == coef_name)
    if (length(idx) != 1) {
      stop("coef_name not found uniquely in design: ", coef_name,
           " | available: ", paste(design_cols, collapse = ", "))
    }
    res <- glmQLFTest(fit, coef = idx)
  }

  deg <- topTags(res, n = Inf)$table
  deg$ensembl <- rownames(deg)
  deg <- add_symbol_from_ensembl(deg, "ensembl")

  list(meta = meta, design = design, deg = deg, res = res)
}

run_two_group_contrast <- function(pb, keep, group_col, group_levels, coef_name, covariates = c("sex")) {
  pb2 <- subset_pb(pb, keep)
  n_tab <- table(pb2$meta[[group_col]], useNA = "ifany")
  needed <- group_levels
  if (!all(needed %in% names(n_tab)) || any(n_tab[needed] < 2)) {
    message("Skip contrast ", coef_name, ": insufficient samples. Counts = ",
            paste(names(n_tab), n_tab, sep = ":", collapse = ", "))
    return(NULL)
  }
  run_edger_generic(pb2$counts, pb2$meta, group_col = group_col,
                    group_levels = group_levels, covariates = covariates,
                    coef_name = coef_name)
}

make_volcano <- function(deg,
                         title,
                         fdr_cut = 0.05,
                         lfc_cut = 1,
                         pde_genes = character(),
                         mech_genes = character()) {
  deg$neglog10FDR <- -log10(pmax(deg$FDR, .Machine$double.xmin))
  deg$status <- "NS"
  deg$status[deg$FDR < fdr_cut & deg$logFC >  lfc_cut] <- "Up"
  deg$status[deg$FDR < fdr_cut & deg$logFC < -lfc_cut] <- "Down"

  pde_set  <- unique(na.omit(as.character(pde_genes)))
  mech_set <- unique(na.omit(as.character(mech_genes)))

  deg$to_label <- !is.na(deg$symbol) & (
    (deg$symbol %in% pde_set) |
    ((deg$symbol %in% mech_set) & deg$FDR < fdr_cut)
  )

  ggplot(deg, aes(x = logFC, y = neglog10FDR)) +
    geom_point(aes(color = status), size = 1, alpha = 0.75) +
    scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "grey80")) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed") +
    ggrepel::geom_text_repel(
      data = deg[deg$to_label, , drop = FALSE],
      aes(label = symbol), size = 3,
      max.overlaps = Inf, box.padding = 0.3, point.padding = 0.2
    ) +
    theme_classic() +
    labs(title = title, x = "log2 Fold Change", y = "-log10(FDR)")
}

prep_gene_df <- function(deg, genes, cell_type, contrast) {
  out <- deg[deg$symbol %in% genes, c("symbol", "logFC", "FDR"), drop = FALSE]
  out$cell_type <- cell_type
  out$contrast <- contrast
  out
}

plot_bar_disease_by_region <- function(df, title) {
  df$region <- factor(df$region, levels = c("LV", "RV"))
  df$star <- p_to_star(df$FDR)
  ggplot(df, aes(x = symbol, y = logFC, fill = region)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    geom_text(aes(label = star, group = region),
              position = position_dodge(width = 0.8),
              vjust = ifelse(df$logFC >= 0, -0.2, 1.1), size = 3) +
    facet_wrap(~ cell_type, ncol = 1) +
    scale_fill_manual(values = c("LV" = "#4DBBD5", "RV" = "#E64B35")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, x = "PDE gene", y = "log2 Fold Change (DCM vs donor)", fill = "Region")
}

build_heatmap_df <- function(results_named_list, genes, contrast_order) {
  out <- list()
  for (nm in names(results_named_list)) {
    res <- results_named_list[[nm]]
    if (is.null(res)) next
    tmp <- res$deg[match(genes, res$deg$symbol), c("symbol", "logFC", "FDR"), drop = FALSE]
    tmp$symbol <- genes
    tmp$contrast <- nm
    out[[nm]] <- tmp
  }
  if (length(out) == 0) {
    return(data.frame(
      symbol = factor(character(), levels = rev(genes)),
      logFC = numeric(),
      FDR = numeric(),
      contrast = factor(character(), levels = contrast_order),
      star = character(),
      stringsAsFactors = FALSE
    ))
  }
  df <- do.call(rbind, out)
  rownames(df) <- NULL
  df$symbol <- factor(df$symbol, levels = rev(genes))
  df$contrast <- factor(df$contrast, levels = contrast_order)
  df$star <- p_to_star(df$FDR)
  df
}

plot_heatmap_tile <- function(df, title) {
  ggplot(df, aes(x = contrast, y = symbol, fill = logFC)) +
    geom_tile(color = "white") +
    geom_text(aes(label = star), size = 3) +
    scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35", midpoint = 0,
                         na.value = "grey90") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()) +
    labs(title = title, x = NULL, y = NULL, fill = "logFC")
}

message("Loaded 00_utils_edger.R")
