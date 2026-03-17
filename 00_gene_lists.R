PDE_GENES <- c(
  "PDE1A","PDE1B","PDE1C","PDE2A","PDE3A","PDE3B",
  "PDE4A","PDE4B","PDE4C","PDE4D","PDE5A","PDE7A","PDE8A","PDE9A"
)

mech_tbl <- readxl::read_excel("Mechanosensitive_Protein_list_updated.xlsx")
mech_tbl$gene <- toupper(trimws(mech_tbl[["Gene Name"]]))
mech_tbl$type <- trimws(as.character(mech_tbl[["Type"]]))
mech_tbl <- mech_tbl[!is.na(mech_tbl$gene) & mech_tbl$gene != "", c("gene", "type")]
mech_tbl <- mech_tbl[!duplicated(mech_tbl$gene), , drop = FALSE]

MECH_GENES <- mech_tbl$gene
MECH_TYPE_MAP <- setNames(mech_tbl$type, mech_tbl$gene)
MECH_TYPE_LEVELS <- unique(na.omit(mech_tbl$type))
MECH_GENES_ORDERED <- unlist(lapply(MECH_TYPE_LEVELS, function(tp) {
  sort(mech_tbl$gene[mech_tbl$type == tp])
}), use.names = FALSE)

CORE_HEATMAP_GENES <- unique(c(PDE_GENES, MECH_GENES))
GENOTYPE_GROUPS <- c("TTN", "RBM20", "LMNA", "PVneg")
CELLTYPES <- c("CM", "EC")

message("Loaded 00_gene_lists.R | PDE genes: ", length(PDE_GENES),
        " | Mechanosensitive genes: ", length(MECH_GENES),
        " | Types: ", paste(MECH_TYPE_LEVELS, collapse = ", "))
