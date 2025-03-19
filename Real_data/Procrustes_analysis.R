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
