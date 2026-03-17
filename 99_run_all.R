# Run this from the pseudobulk working directory containing:
#   CM.counts.mtx / CM.genes.tsv / CM.samples.tsv / CM.meta.csv
#   EC.counts.mtx / EC.genes.tsv / EC.samples.tsv / EC.meta.csv
#   Mechanosensitive_Protein_list_updated.xlsx
# and this script folder.
c("#7FC97F","#BEAED4","#FDC086","#386CB0")
#

message("Working directory: ", getwd())
# setwd("E:/PhD/Science_single_cell/CM_EC")
# Preflight file existence checks
stopifnot(
  file.exists("CM.counts.mtx"), file.exists("CM.genes.tsv"),
  file.exists("CM.samples.tsv"), file.exists("CM.meta.csv"),
  file.exists("EC.counts.mtx"), file.exists("EC.genes.tsv"),
  file.exists("EC.samples.tsv"), file.exists("EC.meta.csv"),
  file.exists("Mechanosensitive_Protein_list_updated.xlsx")
)

source("00_utils_edger.R")
source("00_gene_lists.R")
source("01_volcano_LV_RV_plus_genotype.R")
source("02_barplot_LV_RV_disease_only.R")
source("03_heatmap_genotype_vs_donor.R")
message("All scripts finished.")