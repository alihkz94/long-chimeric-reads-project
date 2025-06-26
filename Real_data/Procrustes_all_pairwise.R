## Compare All Chimera‐Removal Methods Pairwise with Procrustes Analysis

setwd("~/xml_denovo/procrusters")
## 1. Load libraries
library(data.table)
library(vegan)
library(plyr)
library(ggplot2)
library(openxlsx)
library(janitor)
library(RColorBrewer)
library(grid)        # for unit()

## 2. Load and clean sample metadata
META <- fread("Jamy_metadata.txt")
setnames(META,
         old = colnames(META),
         new = make_clean_names(colnames(META)))

## 2a. Only keep samples you trust
valid_samps <- META$run_accession

## 3. Discover method directories containing OTU_table.txt
all_dirs       <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
is_method_dir  <- sapply(all_dirs, function(d) file.exists(file.path(d, "OTU_table.txt")))
method_dirs    <- all_dirs[is_method_dir]
method_names   <- make_clean_names(basename(method_dirs))
names(method_dirs) <- method_names

message("Detected methods:")
print(method_names)

## 4. Build all unique pairwise comparisons of methods
method_pairs <- t(combn(method_names, 2))
pair_dt      <- data.table(method1 = method_pairs[,1],
                           method2 = method_pairs[,2])
message("Will compare ", nrow(pair_dt), " pairs of methods.")

## 5. Load each OTU table once and filter to trusted samples
load_otu <- function(meth) {
  dt <- fread(file.path(method_dirs[[meth]], "OTU_table.txt"))
  setnames(dt, "OTU", "otu")
  
  # keep only trusted samples
  sample_cols <- setdiff(colnames(dt), "otu")
  keep_cols   <- intersect(sample_cols, valid_samps)
  if (length(keep_cols) < length(sample_cols)) {
    dropped <- setdiff(sample_cols, keep_cols)
    message("  → dropping ", length(dropped),
            " samples from ", meth, ": ", paste(dropped, collapse = ", "))
  }
  dt[, c("otu", keep_cols), with = FALSE]
}

OTU_tables <- llply(unique(c(pair_dt$method1, pair_dt$method2)),
                    load_otu, .progress = "none")
names(OTU_tables) <- unique(c(pair_dt$method1, pair_dt$method2))

## 6. Skip any subsampling—use raw count tables directly
TABS <- OTU_tables

## 7. Ordination helper functions (PA + Bray–Curtis)
get_matrix <- function(tab) {
  m <- as.matrix(tab[, -1, with = FALSE])
  rownames(m) <- tab$otu
  t(m)
}
get_ord <- function(mat) {
  mat_pa <- decostand(mat, method = "pa")     # presence/absence
  d      <- vegdist(mat_pa, method = "bray")  # Bray–Curtis on PA
  monoMDS(d, k = 2)
}

## 8. Run Procrustes + permutation test for every pair
RES_list <- alply(seq_len(nrow(pair_dt)), .margins = 1, function(i) {
  m1 <- pair_dt[i, method1]
  m2 <- pair_dt[i, method2]
  dt1 <- TABS[[m1]]; dt2 <- TABS[[m2]]
  shared_samps <- sort(intersect(colnames(dt1)[-1], colnames(dt2)[-1]))
  dt1 <- dt1[, c("otu", shared_samps), with = FALSE]
  dt2 <- dt2[, c("otu", shared_samps), with = FALSE]
  
  ord1 <- get_ord(get_matrix(dt1))
  ord2 <- get_ord(get_matrix(dt2))
  
  prc <- procrustes(ord1, ord2)
  pt  <- protest(ord1, ord2, permutations = 999)
  
  res_vec <- residuals(prc)
  data.table(
    method1   = m1,
    method2   = m2,
    sample    = names(res_vec),
    residual  = as.numeric(res_vec),
    R_value   = pt$t0,
    p_value   = pt$signif
  )
}, .progress = "text")

RES <- rbindlist(RES_list)
RES[, Comparison := paste(method1, "vs", method2)]

## 9. Merge with sample metadata
keep_meta <- c("run_accession", "isolation_source")
RES <- merge(RES,
             META[, ..keep_meta],
             by.x = "sample", by.y = "run_accession",
             all.x = TRUE)

## 10. DROP all rows where isolation_source is NA
RES <- RES[!is.na(isolation_source)]

## 11. Drop isolation sources with only one sample
src_counts <- RES[, .N, by = isolation_source]
valid_srcs <- src_counts[N > 1, isolation_source]
RES <- RES[isolation_source %in% valid_srcs]

## 12. Export full residual table AND summary stats
# 12a. Summary table: one row per method pair
summary_stats <- RES[, .(
  R_value = unique(R_value),
  p_value = unique(p_value)
), by = .(method1, method2, Comparison)]

# 12b. Write both sheets into one workbook
write.xlsx(
  x = list(
    Residuals = RES,
    Summary   = summary_stats
  ),
  file      = "Procrustes_residuals_all_pairs_with_stats.xlsx",
  colNames  = TRUE
)

## 13. Build display names for plotting
get_display_name <- function(method_name) {
  if (method_name == "chimeras_denovo_custom") {
    "chimeras_denovo adjusted+FP-FN"
  } else if (method_name == "uchime_denovo_custom") {
    "uchime_denovo adjusted+FP-FN"
  } else if (method_name == "remove_bimera_denovo_custom") {
    "removeBimeraDenovo+FP-FN"
  } else if (method_name == "chimeras_denovo_default") {
    "chimeras_denovo default"
  } else if (method_name == "uchime_denovo_default") {
    "uchime_denovo default"
  } else if (method_name == "remove_bimera_denovo_default") {
    "removeBimeraDenovo default"
  } else {
    gsub("_", " ", method_name)
  }
}

RES[, method1_display := sapply(method1, get_display_name)]
RES[, method2_display := sapply(method2, get_display_name)]
RES[, ComparisonLabel := sprintf("%s vs %s", method1_display, method2_display)]

## 14. Palette for facets
cmp_levels <- unique(RES$ComparisonLabel)
n_cmp      <- length(cmp_levels)
base_cols  <- brewer.pal(min(n_cmp, 12), "Paired")

cmp_cols <- if (n_cmp > 12) {
  setNames(colorRampPalette(base_cols)(n_cmp), cmp_levels)
} else {
  setNames(base_cols[1:n_cmp], cmp_levels)
}

## 15. Plot: Boxplots by isolation source, faceted by ComparisonLabel
plt_box <- ggplot(RES, aes(
  x    = isolation_source,
  y    = residual,
  fill = ComparisonLabel
)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(
    name   = "Comparison",
    values = cmp_cols
  ) +
  facet_wrap(~ ComparisonLabel, ncol = 3) +
  labs(x = "Isolation source",
       y = "residuals",
       title = "Residual distributions by isolation source\nfor each chimera-removal comparison") +
  theme_light() +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, face = "bold"),
    strip.text        = element_text(size = 8),
    legend.position   = "none"
  )

print(plt_box)

ggsave("boxplots_all_pairs_by_source.svg", plt_box,
       width  = 16, height = 10, device = "svg")
ggsave("boxplots_all_pairs_by_source.pdf", plt_box,
       width  = 16, height = 10, device = "pdf")





# extract R per ComparisonLabel from RES

# assume RES has a column ComparisonLabel and R_value
summary_R <- RES[, .(R_value = unique(R_value)), by = ComparisonLabel]

# 2. build the label vector
label_r <- setNames(
  sprintf("%s (R=%.2f)",
          summary_R$ComparisonLabel,
          summary_R$R_value),
  summary_R$ComparisonLabel
)

# 3. redraw the boxplot with R in facet strips
library(ggplot2)

plt_box_with_R <- ggplot(RES, aes(
  x    = isolation_source,
  y    = residual,
  fill = ComparisonLabel
)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(name   = "Comparison",
                    values = cmp_cols) +
  facet_wrap(
    ~ ComparisonLabel,
    ncol    = 3,
    labeller = labeller(ComparisonLabel = label_r)
  ) +
  labs(
    x     = "isolation source",
    y     = "residual",
    title = "Residual distributions by isolation source for each chimera removal comaprison"
  ) +
  theme_light() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, face = "bold"),
    strip.text      = element_text(size = 8),
    legend.position = "none"
  )

print(plt_box_with_R)

# 4. save if desired
ggsave("boxplots_all_pairs_with_R.svg", plt_box_with_R,
       width = 16, height = 10, device = "svg")
ggsave("boxplots_all_pairs_with_R.pdf", plt_box_with_R,
       width = 16, height = 10, device = "pdf")
