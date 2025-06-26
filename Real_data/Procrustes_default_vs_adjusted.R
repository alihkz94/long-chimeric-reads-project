## Compare Default vs Custom Chimera‐Removal with Procrustes Analysis 

library(data.table)
library(vegan)
library(plyr)
library(ggplot2)
library(openxlsx)
library(janitor)

#--- read and clean metadata
META <- fread("Jamy_metadata.txt")
setnames(META,
         old = colnames(META),
         new = make_clean_names(colnames(META)))

#--- which samples to keep
valid_samps <- META$run_accession

#--- find all method folders with an OTU table
all_dirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
is_method_dir <- sapply(all_dirs, function(d) {
  file.exists(file.path(d, "OTU_table.txt"))
})
method_dirs_full <- all_dirs[is_method_dir]
method_names <- trimws(basename(method_dirs_full))
names(method_dirs_full) <- method_names

#--- pair up “default” vs “custom” methods
defaults <- grep(" default$", method_names, value = TRUE)
pair_list <- lapply(defaults, function(dn) {
  cn <- sub(" default$", "_custom", dn)
  if (cn %in% method_names) {
    return(c(default = dn, custom = cn))
  } else {
    warning("No matching custom folder for default: ", dn)
    return(NULL)
  }
})
pair_list <- Filter(Negate(is.null), pair_list)
if (length(pair_list) == 0) {
  stop("No valid default–custom pairs found. Check folder names.")
}
pair_matrix <- do.call(rbind, pair_list)
message("## Detected pairs for analysis:")
print(pair_matrix)

#--- function to load each OTU table
load_otu <- function(name) {
  dt <- fread(file.path(method_dirs_full[[name]], "OTU_table.txt"))
  sample_cols <- setdiff(colnames(dt), "OTU")
  keep_samps  <- intersect(sample_cols, valid_samps)
  if (length(keep_samps) < length(sample_cols)) {
    dropped <- setdiff(sample_cols, keep_samps)
    message("Dropping ", length(dropped), 
            " samples from ", name, ": ", paste(dropped, collapse =", "))
  }
  dt <- dt[, c("OTU", keep_samps), with = FALSE]
  return(dt)
}

#--- load all the OTU tables, raw counts
unique_methods <- unique(as.vector(pair_matrix))
TABS_unstd <- llply(unique_methods, load_otu)
names(TABS_unstd) <- unique_methods

#--- no subsampling: just use raw tables
TABS <- TABS_unstd

#--- helper to get sample-by-OTU matrix
get_matrix <- function(tab) {
  m <- as.matrix(tab[, -1, with = FALSE])
  rownames(m) <- tab$OTU
  t(m)
}

#--- ordination now uses presence/absence + Bray–Curtis
get_ordination <- function(mat) {
  mat_pa <- decostand(mat, method = "pa")    # presence / absence
  d <- vegdist(mat_pa, method = "bray")     # Bray–Curtis on PA
  monoMDS(d, k = 2)
}

#--- Procrustes + permutation test for each default/custom pair
RES_list <- alply(seq_len(nrow(pair_matrix)), .margins = 1, function(i) {
  dname <- pair_matrix[i, "default"]
  cname <- pair_matrix[i, "custom"]
  df1 <- TABS[[dname]]; df2 <- TABS[[cname]]
  shared <- sort(intersect(colnames(df1)[-1], colnames(df2)[-1]))
  df1 <- df1[, c("OTU", shared), with = FALSE]
  df2 <- df2[, c("OTU", shared), with = FALSE]
  
  ord1 <- get_ordination(get_matrix(df1))
  ord2 <- get_ordination(get_matrix(df2))
  
  prc <- procrustes(ord1, ord2)
  pt  <- protest(ord1, ord2, permutations = 999)
  
  res <- residuals(prc)
  data.table(
    default   = dname,
    custom    = cname,
    sample    = names(res),
    residual  = as.numeric(res),
    R_value   = pt$t0,
    p_value   = pt$signif
  )
}, .progress = "text")

RES <- rbindlist(RES_list)
RES[, Comparison := paste(default, "vs", custom)]

#--- add metadata
meta_cols <- c("run_accession", "sample_accession", "isolation_source",
               "environment_biome", "environment_feature",
               "environment_material", "geographic_location_country_and_or_sea",
               "geographic_location_depth")
RES <- merge(RES, META[, ..meta_cols],
             by.x = "sample", by.y = "run_accession", all.x = TRUE)

#--- write results to Excel
write.xlsx(list("Residuals" = RES),
           file = "Procrustes_residuals_default_vs_custom.xlsx",
           colNames = TRUE)

#--- set up labels for plotting
static_labels <- c(
  "chimeras_denovo_custom"    = "chimeras_denovo adjusted+FP-FN",
  "uchime_denovo_custom"      = "uchime_denovo adjusted+FP-FN",
  "removeBimeraDenovo_custom" = "removeBimeraDenovo+FP-FN"
)
RES[, method_label := static_labels[custom]]
RES[, method_label := factor(
  method_label,
  levels = c(
    "removeBimeraDenovo+FP-FN",
    "chimeras_denovo adjusted+FP-FN",
    "uchime_denovo adjusted+FP-FN"
  )
)]

#--- final boxplot of residuals by isolation source
plt <- ggplot(RES, aes(x = isolation_source,
                       y = residual,
                       fill = method_label)) +
  geom_boxplot() +
  scale_fill_manual(
    values = c(
      "uchime_denovo adjusted+FP-FN"   = "#5177B8",
      "chimeras_denovo adjusted+FP-FN" = "#E97169",
      "removeBimeraDenovo+FP-FN"       = "#0A9F37"
    ),
    breaks = c(
      "removeBimeraDenovo+FP-FN",
      "chimeras_denovo adjusted+FP-FN",
      "uchime_denovo adjusted+FP-FN"
    )
  ) +
  labs(
    x    = "isolation source",
    y    = "residual",
    fill = "Chimera removal method"
  ) +
  ggtitle("Procrustes residuals (Default vs adjusted chimera filtering)") +
  theme_light() +
  theme(axis.text.x = element_text(hjust = 1))

ggsave(
  filename = "Procrustes_residuals_default_vs_custom.svg",
  plot     = plt,
  width    = 12,
  height   = 7,
  device   = "svg"
)
