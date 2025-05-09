## Compare Default vs Custom Chimera‐Removal with Procrustes Analysis 

library(data.table)
library(vegan)
library(SRS)
library(plyr)
library(ggplot2)
library(openxlsx)
library(janitor)

## 1. Load and clean sample metadata
META <- fread("Jamy_metadata.txt")
setnames(META,
         old = colnames(META),
         new = make_clean_names(colnames(META)))

## 2. Discover method directories containing OTU_table.txt
all_dirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
is_method_dir <- sapply(all_dirs, function(d) {
  file.exists(file.path(d, "OTU_table.txt"))
})
method_dirs_full <- all_dirs[is_method_dir]
method_names <- trimws(basename(method_dirs_full))
names(method_dirs_full) <- method_names

## 3. Match each <prefix>_default to <prefix>_custom
defaults <- grep("_default$", method_names, value = TRUE)
pair_list <- lapply(defaults, function(dn) {
  cn <- sub("_default$", "_custom", dn)
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

## 4. Load OTU tables (uncompressed)
load_otu <- function(name) {
  fread(file.path(method_dirs_full[[name]], "OTU_table.txt"))
}
unique_methods <- unique(as.vector(pair_matrix))
TABS_unstd <- llply(unique_methods, load_otu)
names(TABS_unstd) <- unique_methods

## 5. Compute minimum sequencing depth across all methods
min_depth <- min(sapply(TABS_unstd, function(dt) {
  min(colSums(dt[ , -1, with = FALSE]))
}))
message("Minimum library size for SRS: ", min_depth)

## 6. Standardize depths via SRS
TABS <- llply(TABS_unstd, function(dt) {
  df <- as.data.frame(dt)
  rownames(df) <- df$OTU; df$OTU <- NULL
  srs_mat <- SRS(df, Cmin = min_depth, seed = 42)
  res <- as.data.table(srs_mat, keep.rownames = "OTU")
  setcolorder(res, c("OTU", setdiff(colnames(res), "OTU")))
  return(res)
}, .progress = "none")

## 7. Ordination helper functions
get_matrix <- function(tab) {
  m <- as.matrix(tab[, -1, with = FALSE])
  rownames(m) <- tab$OTU
  t(m)
}
get_ordination <- function(mat) {
  h <- decostand(mat, method = "hellinger")
  d <- vegdist(h, method = "bray")
  monoMDS(d, k = 2)
}

## 8. Compute Procrustes residuals for each pair
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
  res <- residuals(prc)
  data.table(default = dname,
             custom  = cname,
             sample  = names(res),
             residual = as.numeric(res))
}, .progress = "text")
RES <- rbindlist(RES_list)
RES[, Comparison := paste(default, "vs", custom)]

## 9. Merge with sample metadata
meta_cols <- c("run_accession", "sample_accession", "isolation_source",
               "environment_biome", "environment_feature",
               "environment_material", "geographic_location_country_and_or_sea",
               "geographic_location_depth")
RES <- merge(RES, META[, ..meta_cols],
             by.x = "sample", by.y = "run_accession", all.x = TRUE)

## 10. Export to Excel
write.xlsx(list("Residuals" = RES),
           file = "Procrustes_residuals_Default_vs_Custom_debugged.xlsx",
           colNames = TRUE)

## 11. Relabel custom methods for plotting
label_map <- c(
  "chimeras_denovo_custom"    = "chimeras_denovo adjusted+FP-FN",
  "uchime_denovo_custom"      = "uchime_denovo adjusted+FP-FN",
  "removeBimeraDenovo_custom" = "removeBimeraDenovo adjusted+FP-FN"
)
RES[, method_label := label_map[custom]]

## 12. Order methods by mean residual and set factor levels
order_dt <- RES[, .(MeanResid = mean(residual)), by = method_label][order(-MeanResid)]
RES[, method_label := factor(method_label, levels = order_dt$method_label)]

## 13. Plot boxplots of residuals by isolation source with new legend labels
p <- ggplot(RES, aes(x = isolation_source, y = residual, fill = method_label)) +
  geom_boxplot() +
  labs(x    = "isolation source",
       y    = "residual",
       fill = "Chimera removal method") +
  ggtitle("Procrustes residuals (Default vs adjusted chimera filtering)") +
  theme_light() +
  theme(axis.text.x = element_text(hjust = 1))

## 14. Save as SVG
ggsave("Procrustes_residuals_default_vs_custom.svg",
       plot   = p,
       width  = 10,
       height = 6,
       device = "svg")
