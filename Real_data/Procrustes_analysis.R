## Estimate Procrustes residuals between various chimera removal methods and raw data
## The higher the residual, the more different the two datasets are (in terms of community composition)

## Expected input data structure:
# ├── working_directory
# │   ├── Jamy_metadata.txt
# │   └── OTUs
# │       ├── chimeras_denovo_custom
# │       │   ├── OTUs.fasta.gz
# │       │   └── OTU_table.txt.gz
# │       ├── chimeras_denovo_default
# │       │   ├── OTUs.fasta.gz
# │       │   └── OTU_table.txt.gz
# │       ├── Raw
# │       │   ├── OTUs.fasta.gz
# │       │   └── OTU_table.txt.gz
# │       ├── removeBiemraDenovo_default
# │       │   ├── OTUs.fasta.gz
# │       │   └── OTU_table.txt.gz
# │       ├── removeBimeraDenovo_custom
# │       │   ├── OTUs.fasta.gz
# │       │   └── OTU_table.txt.gz
# │       ├── uchime_denovo_custom
# │       │   ├── OTUs.fasta.gz
# │       │   └── OTU_table.txt.gz
# │       └── uchime_denovo_default
# │           ├── OTUs.fasta.gz
# │           └── OTU_table.txt.gz


library(data.table)
library(vegan)
library(SRS)
library(plyr)
library(ggplot2)
library(openxlsx)

## Load sample metadata
META <- fread("Jamy_metadata.txt")

## Clean column names
setnames(META,
  old = colnames(META),
  new = janitor::make_clean_names(colnames(META)))

## Load OTU tables
dirs <- list.dirs("OTUs", recursive = FALSE)
TABS_unstandardized <- alply(
  .data = dirs,
  .margins = 1,
  .fun = function(dir) { fread(file.path(dir, "OTU_table.txt.gz")) })
names(TABS_unstandardized) <- basename(dirs)

## Get minimum sequencing depth
min_depth <- ldply(.data = TABS_unstandardized, .fun = function(x){
  data.table(MinDepth = min(colSums(x[ , -1])))
})
min_depth$Method <- names(TABS_unstandardized)

min(min_depth$MinDepth)   # 41094


## Standardize sequencing depth
TABS <- llply(.data = TABS_unstandardized, .fun = function(x){
  setDF(x)
  rownames(x) <- x$OTU
  x$OTU <- NULL
  res <- SRS(data = x, Cmin = 41000, seed = 42)
  res$OTU <- rownames(res)
  setDT(res)
  setcolorder(res, "OTU")
  return(res)
}, .progress = "text")


## Estimate sepcies richness and diversity (in raw data), add to metadata
DIV <- data.table(
  SampleID = colnames(TABS$Raw)[-1],
  Richness = specnumber(x = t(TABS$Raw[,-1])),
  Shannon  = diversity(x = t(TABS$Raw[,-1]))
)

## Merge with metadata
META <- merge(x = META, y = DIV,
  by.x = "run_accession", by.y = "SampleID",
  all.x = TRUE)


## Function to prepare matrix for ordination
get_matrix <- function(tab) {
    tb <- as.matrix(tab[, -1])
    rownames(tb) <- tab$OTU
    tb <- t(tb)
    return(tb)
}

## Function to get ordination scores
get_ordination <- function(tb) {
    tb <- decostand(tb, method = "hellinger")
    dd <- vegdist(tb, method = "bray")
    mds <- monoMDS(dd, k = 2)
    return(mds)
}

## Run Procrsutes test and estimate residuals between samples

## All pairwise comparisons
# cmb <- combn(names(TABS), 2)

## Comparisons with raw data (dataset without chimera removal)
c2 <- names(TABS)[ ! names(TABS) %in% "Raw"]
c1 <- rep("Raw", times = length(c2))
cmb <- matrix(c(c1, c2), nrow = 2, byrow = TRUE)

RES <- alply(
    .data = 1:ncol(cmb),
    .margins = 1,
    .fun = function(i) {

        ## Get data
        df1 <- TABS[[ cmb[1, i] ]]
        df2 <- TABS[[ cmb[2, i] ]]

        ## Subset data to the same set of samples
        smp <- sort(intersect(colnames(df1)[-1], colnames(df2)[-1]))
        unq <- setdiff(colnames(df1)[-1], colnames(df2)[-1])
        if(length(unq) > 0) {
            cat("Initial sample sets are not identical. Missing samples: ", unq, 
                "; cmb: ", cmb[1, i], " vs ", cmb[2, i], "\n")
        }

        ## Subset data to the same set of samples
        clz <- c("OTU", smp)
        df1 <- df1[ , ..clz ]
        df2 <- df2[ , ..clz ]

        ## Prepare matrices for ordination
        tb1 <- get_matrix(df1)
        tb2 <- get_matrix(df2)

        ## Get ordination scores
        ord1 <- get_ordination(tb1)
        ord2 <- get_ordination(tb2)

        ## Run procrsutes test
        prc <- procrustes(ord1, ord2)

        ## Extract residuals
        res <- residuals(prc)
        res <- data.table(
            table1   = cmb[1, i],
            table2   = cmb[2, i],
            sample   = names(res),
            residual = res)
        return(res)
}, .progress = "text")

RES <- rbindlist(RES)
RES[ , Comparsion := paste(table1, table2, sep = " vs ") ]

## Add sample metadata
meta_clz <- c("run_accession", "sample_accession", "isolation_source",
  "environment_biome", "environment_feature", "environment_material", 
  "geographic_location_country_and_or_sea", "geographic_location_depth",
  "Richness", "Shannon")

if(any(!unique(RES$sample) %in% META$run_accession)) {
  cat("Some samples are not present in the metadata. Missing samples: ", 
  setdiff(unique(RES$sample), META$run_accession), "\n")
}

RES <- merge(
    x = RES, 
    y = META[, ..meta_clz],
    by.x = "sample", by.y = "run_accession", all.x = TRUE)


## Export results
write.xlsx(list(
  "Residuals" = RES
  ),
  file = "Jamy__Procrustes_residuals_VsRaw.xlsx",
  colNames = TRUE)


####### Compare residuals between various isolation sources

## Reorder by residual
smr <- RES[ , .(MeanResid = mean(residual)), by = "table2" ]
setorder(smr, -MeanResid)

RES$table2 <- factor(RES$table2, levels = smr$table2)

## Plot
pp <- ggplot(RES, aes(x = isolation_source, y = residual)) +
    geom_boxplot(aes(fill = table2)) + 
    labs(fill = "Chimera removal method") +
    ggtitle("Procrustes residuals (Raw data vs Various chimera removal methods)") +
    theme_light()

## Save plot
ggsave(pp,
  filename = "Jamy__Procrustes_residuals_VsRaw.pdf",
  width = 10, height = 6)




###### OTHER METHODOLOGY (PRESENCE/ABSENCE) #########
#––– libraries
library(data.table)
library(vegan)
library(plyr)
library(ggplot2)
library(openxlsx)
library(janitor)

#––– 1. Load & clean metadata
META <- fread("Jamy_metadata.txt")
setnames(META,
         old = colnames(META),
         new = make_clean_names(colnames(META))
)

#––– 1a. Only keep samples you trust
valid_samps <- META$run_accession

#––– 2. Discover method directories containing OTU_table.txt
all_dirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
is_method_dir <- sapply(all_dirs, function(d) {
  file.exists(file.path(d, "OTU_table.txt"))
})
method_dirs_full <- all_dirs[is_method_dir]
method_names <- trimws(basename(method_dirs_full))
names(method_dirs_full) <- method_names

if (length(method_dirs_full) == 0) {
  stop("No directories with OTU_table.txt found. Check your directory structure.")
}

message("Found method directories:")
print(method_names)

#––– 3. Load OTU tables and drop unwanted samples
load_otu <- function(name) {
  dt <- fread(file.path(method_dirs_full[[name]], "OTU_table.txt"))
  sample_cols <- setdiff(colnames(dt), "OTU")
  keep_cols   <- intersect(sample_cols, valid_samps)
  if (length(keep_cols) < length(sample_cols)) {
    dropped <- setdiff(sample_cols, keep_cols)
    message("  → dropping ", length(dropped),
            " samples from ", name, ": ", paste(dropped, collapse = ", "))
  }
  dt[, c("OTU", keep_cols), with = FALSE]
}

TABS_unstd <- llply(method_names, load_otu)
names(TABS_unstd) <- method_names

#––– 3a. Canonicalize method names: replace spaces with underscores, rename "reference" → "Raw"
canonical_names <- gsub(" ", "_", method_names)
canonical_names[canonical_names == "reference"] <- "Raw"
names(TABS_unstd) <- canonical_names

message("Loaded methods with canonical names:")
for (i in seq_along(method_names)) {
  message("  '", method_names[i], "' -> '", canonical_names[i], "'")
}

#––– 3b. Check that Raw (reference) is present
if (!"Raw" %in% names(TABS_unstd)) {
  stop("Raw data not found. Available methods: ", paste(names(TABS_unstd), collapse = ", "))
}

#––– 4. Skip SRS subsampling; use raw counts directly
TABS <- TABS_unstd

#––– 5. Compute alpha diversity on raw counts
raw_mat <- as.matrix(TABS$Raw[, -1, with = FALSE])
rownames(raw_mat) <- TABS$Raw$OTU

DIV <- data.table(
  SampleID = colnames(raw_mat),
  Richness = specnumber(t(raw_mat)),
  Shannon  = diversity(t(raw_mat))
)
META <- merge(
  META, DIV,
  by.x = "run_accession", by.y = "SampleID",
  all.x = TRUE
)

#––– 6. Ordination helpers (PA + Bray–Curtis)
get_matrix <- function(tab) {
  m <- as.matrix(tab[, -1, with = FALSE])
  rownames(m) <- tab$OTU
  t(m)
}
get_ordination <- function(mat) {
  mat_pa <- decostand(mat, "pa")     # presence / absence
  d      <- vegdist(mat_pa, "bray")  # Bray–Curtis on PA
  monoMDS(d, k = 2)
}

#––– 7. Procrustes + protest (R & p) vs Raw
methods <- setdiff(names(TABS), "Raw")
cmb <- rbind(rep("Raw", length(methods)), methods)

RES_list <- alply(seq_len(ncol(cmb)), 1, function(i) {
  m1 <- cmb[1, i]; m2 <- cmb[2, i]
  df1 <- TABS[[m1]]; df2 <- TABS[[m2]]
  smp <- intersect(colnames(df1)[-1], colnames(df2)[-1])
  df1 <- df1[, c("OTU", smp), with = FALSE]
  df2 <- df2[, c("OTU", smp), with = FALSE]
  
  ord1 <- get_ordination(get_matrix(df1))
  ord2 <- get_ordination(get_matrix(df2))
  
  pr <- procrustes(ord1, ord2)
  pt <- protest(ord1, ord2, permutations = 999)
  
  r_vec <- residuals(pr)
  data.table(
    table1   = m1,
    table2   = m2,
    sample   = names(r_vec),
    residual = as.numeric(r_vec),
    R_value  = pt$t0,
    p_value  = pt$signif
  )
}, .progress = "text")

RES <- rbindlist(RES_list)
RES[, Comparison := paste(table1, table2, sep = " vs ")]

#––– 8. Merge sample metadata
keep_cols <- c("run_accession", "sample_accession", "isolation_source",
               "environment_biome", "environment_feature",
               "environment_material", "geographic_location_country_and_or_sea",
               "geographic_location_depth", "Richness", "Shannon")
RES <- merge(RES, META[, ..keep_cols],
             by.x = "sample", by.y = "run_accession", all.x = TRUE)

#––– 9. Order factor levels of table2 by mean residual
ord2 <- RES[, .(MeanResid = mean(residual)), by = "table2"]
setorder(ord2, -MeanResid)
RES$table2 <- factor(RES$table2, levels = ord2$table2)

#––– 10. Define colors & labels
method_colors <- c(
  chimeras_denovo_custom     = "#E97169",
  removeBimeraDenovo_custom  = "#0A9F37",
  uchime_denovo_custom       = "#5177B8",
  chimeras_denovo_default    = "#FF8C00",
  removeBimeraDenovo_default = "#9370DB",
  uchime_denovo_default      = "#FFD700"
)
method_labels <- c(
  chimeras_denovo_custom     = "chimeras_denovo adjusted+FP-FN",
  uchime_denovo_custom       = "uchime_denovo adjusted+FP-FN",
  removeBimeraDenovo_custom  = "removeBimeraDenovo+FP-FN",
  chimeras_denovo_default    = "chimeras_denovo default",
  uchime_denovo_default      = "uchime_denovo default",
  removeBimeraDenovo_default = "removeBimeraDenovo default"
)

#––– 11. Plot (print + save PDF & SVG)
pp <- ggplot(RES, aes(isolation_source, residual, fill = table2)) +
  geom_boxplot() +
  scale_fill_manual(
    name   = "Chimera removal method",
    values = method_colors,
    labels = method_labels
  ) +
  labs(
    x     = "isolation source",
    y     = "residual",
    title = "Procrustes residuals: Raw vs chimera removal methods"
  ) +
  theme_light() +
  theme(axis.text.x = element_text(hjust = 1))

print(pp)

ggsave(pp,
       filename = "Jamy__Procrustes_residuals_VsRaw_chat.pdf",
       width = 12, height = 7)
ggsave(pp,
       filename = "Jamy__Procrustes_residuals_VsRaw_chat.svg",
       width = 12, height = 7, device = "svg")

#––– 12. Export results (including R & p) to Excel
write.xlsx(
  list(Residuals = RES),
  file     = "Jamy__Procrustes_residuals_VsRaw_with_stats_chat.xlsx",
  colNames = TRUE
)

message("Analysis complete! Files saved:")
message("  - Jamy__Procrustes_residuals_VsRaw_chat.pdf")
message("  - Jamy__Procrustes_residuals_VsRaw_chat.svg")
message("  - Jamy__Procrustes_residuals_VsRaw_with_stats.xlsx")
