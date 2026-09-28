# Purpose:  Partition expression variation across the 12 matched pseudobulk
#           profiles (4 animals x VisiumHD / MERSCOPE / Xenium) into Platform,
#           Group (WT vs KO) and Residual with PERMANOVA (vegan::adonis2), for
#           the Figure 2 variance-explained panel.
#
#           Distance: Euclidean on TMM log-CPM (prior.count = 2, all common
#           genes), the same matrix and distance as the Fig 2 pseudobulk MDS
#           (fig2_gene_comparison.R).
#           Model: dist ~ Platform + Group, sequential (type I) SS. The design
#           is balanced, so Platform and Group are orthogonal and sequential R2
#           equals marginal R2. Residual = between-animal (within-group)
#           variation plus the Animal x Platform interaction.
#           Permutations:
#             - Platform: permuted within Animal (blocks); (3!)^4 = 1296
#               permutations, complete enumeration, so P is exact.
#             - Group: constant within Animal, so not testable by within-animal
#               permutation; permuting whole animals gives only 3 distinct
#               splits (min P = 1/3) and the wrong error stratum. Group P is
#               set to NA; Group R2 is descriptive.
# Inputs:   results/03_benchmarking/gene_comparison/dge.rds
# Outputs:  results/03_benchmarking/variance_partition/permanova_results.rds
#               list: table (term, Df, SumOfSqs, R2, pct_var, F, p_value),
#                     adonis2, dist, design, distance, ss_type, term_order,
#                     permutation, n_perm
#           results/03_benchmarking/variance_partition/variance_explained.rds
#               tibble: Factor (Platform / Group / Residual), pct_var

suppressPackageStartupMessages({
  library(edgeR)
  library(vegan)
  library(permute)
  library(dplyr)
  library(tibble)
  library(optparse)
})

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--dge"),     type = "character",
              default = "results/03_benchmarking/gene_comparison/dge.rds",
              help    = "DGEList from gene_comparison.R [default: %default]"),
  make_option(c("--out_dir"), type = "character",
              default = "results/03_benchmarking/variance_partition",
              help    = "Output directory [default: %default]"),
  make_option(c("--seed"),    type = "integer", default = 42L,
              help    = "RNG seed (only used if permutations are sampled) [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!file.exists(opt$dge)) stop("Not found: ", opt$dge)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(opt$seed)

# ---------------------------------------------------------------------------
# Expression matrix and distance
# ---------------------------------------------------------------------------
message("Loading DGEList...")
dge      <- readRDS(opt$dge) |> calcNormFactors()
logcpm   <- cpm(dge, log = TRUE, prior.count = 2)
dist_mat <- dist(t(logcpm), method = "euclidean")

# Animal identifies the matched profiles across platforms
design_df <- as.data.frame(dge$samples) |>
  rownames_to_column("ColName") |>
  transmute(
    ColName,
    Platform = factor(Platform),
    Group    = droplevels(factor(Type)),
    Animal   = factor(IDnum)
  )
stopifnot(identical(design_df$ColName, labels(dist_mat)))

# Check the matched design: each animal on every platform, one Group per animal
stopifnot(all(table(design_df$Animal, design_df$Platform) == 1))
stopifnot(all(rowSums(table(design_df$Animal, design_df$Group) > 0) == 1))
message("  ", nrow(logcpm), " genes x ", ncol(logcpm), " samples; ",
        nlevels(design_df$Platform), " platforms, ",
        nlevels(design_df$Animal), " animals, groups: ",
        paste(levels(design_df$Group), collapse = "/"))

# ---------------------------------------------------------------------------
# PERMANOVA
# ---------------------------------------------------------------------------
# Permute platform labels within each animal (valid for Platform only)
perm_ctrl <- how(blocks = design_df$Animal, nperm = 9999)
n_perm    <- nrow(shuffleSet(nrow(design_df), control = perm_ctrl))
message("Running adonis2 (", n_perm, " within-animal permutations)...")

permanova <- adonis2(
  dist_mat ~ Platform + Group,
  data         = design_df,
  by           = "terms",
  permutations = perm_ctrl
)
# Withhold the Group P value (also in the stored adonis2 object): Group SS is
# invariant under within-animal permutations, so its "P" only reflects
# changes in the residual and is not a valid test of WT vs KO
permanova[["Pr(>F)"]][rownames(permanova) == "Group"] <- NA_real_

permanova_tbl <- as.data.frame(permanova) |>
  rownames_to_column("term") |>
  as_tibble() |>
  rename(F = `F`, p_value = `Pr(>F)`) |>
  mutate(
    pct_var = R2 * 100) |>
  select(term, Df, SumOfSqs, R2, pct_var, F, p_value)
print(as.data.frame(permanova_tbl))

# ---------------------------------------------------------------------------
# Variance explained
# ---------------------------------------------------------------------------
factor_levels <- c("Platform", "Group", "Residual")
variance_explained <- permanova_tbl |>
  filter(term %in% factor_levels) |>
  transmute(Factor = factor(term, levels = factor_levels), pct_var) |>
  arrange(Factor)
stopifnot(nrow(variance_explained) == 3,
          abs(sum(variance_explained$pct_var) - 100) < 1e-6)

message("  PERMANOVA R2 (%): ",
        paste(variance_explained$Factor, sprintf("%.1f", variance_explained$pct_var),
              sep = " = ", collapse = ", "))
message("  Platform permutation P (within-animal, exact): ",
        signif(permanova_tbl$p_value[permanova_tbl$term == "Platform"], 3),
        "; Group P not reported (not testable under this design)")

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
permanova_results <- list(
  table       = permanova_tbl,
  adonis2     = permanova,
  dist        = dist_mat,
  design      = design_df,
  distance    = "euclidean on TMM log-CPM (prior.count = 2), all common genes",
  ss_type     = "sequential (type I; adonis2 by = 'terms')",
  term_order  = c("Platform", "Group"),
  permutation = paste("Platform: labels permuted within Animal (blocks),",
                      "complete enumeration; Group: not testable, P = NA"),
  n_perm      = n_perm
)

saveRDS(permanova_results,  file.path(opt$out_dir, "permanova_results.rds"))
saveRDS(variance_explained, file.path(opt$out_dir, "variance_explained.rds"))
message("Saved: permanova_results.rds, variance_explained.rds to ", opt$out_dir)
message("Done.")
