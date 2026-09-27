# ============================================================
# Sepsis scRNA-seq Analysis — Full Pipeline Script
# GSE163668
# ============================================================

library(GEOquery)
library(Seurat)
library(dplyr)

# ============================================================
# 1. DOWNLOAD DATA
# ============================================================

base_path <- "E:/Sepsis Project/Sepsis.Data/GSE163668/extracted"

getGEOSuppFiles(
  "GSE163668",
  makeDirectory = TRUE,
  baseDir = "E:/Sepsis Project/Sepsis.Data/"
)

untar(
  "E:/Sepsis Project/Sepsis.Data/GSE163668/GSE163668_RAW.tar",
  exdir = base_path
)

gse_info <- getGEO("GSE163668", GSEMatrix = TRUE)
meta <- pData(gse_info[[1]])

# ============================================================
# 2. PROCESS POOLED SAMPLES (DEMULTIPLEXING)
# ============================================================

pooled_samples <- c(
  "GSM4995425_MVIR1-POOL-SCG1",
  "GSM4995426_MVIR1-POOL-SCG2",
  "GSM4995427_MVIR1-POOL-SCG3",
  "GSM4995428_MVIR1-POOL-SCG4",
  "GSM4995429_MVIR1-POOL-SCG8",
  "GSM4995430_MVIR1-POOL-SCG9"
)

process_pool <- function(sample_id, base_path) {
  message("Processing: ", sample_id)

  cluster_assign <- read.table(
    file.path(base_path, paste0(sample_id, "_cluster_assignments.tsv.gz")),
    header = TRUE, stringsAsFactors = FALSE
  )
  colnames(cluster_assign) <- c("cluster_id", "patient")

  samples_file <- read.table(
    file.path(base_path, paste0(sample_id, ".clust1.samples.txt.gz")),
    header = TRUE, sep = "\t", stringsAsFactors = FALSE
  )

  samples_clean <- samples_file %>%
    filter(DROPLET.TYPE == "SNG")

  samples_clean$cluster_id <- as.integer(
    sapply(strsplit(samples_clean$BEST.GUESS, ","), `[`, 1)
  )

  samples_clean <- samples_clean %>%
    left_join(cluster_assign, by = "cluster_id")

  sample_dir <- file.path(base_path, sample_id)
  dir.create(sample_dir, showWarnings = FALSE)

  file.copy(file.path(base_path, paste0(sample_id, "_barcodes.tsv.gz")),
            file.path(sample_dir, "barcodes.tsv.gz"), overwrite = TRUE)
  file.copy(file.path(base_path, paste0(sample_id, "_features.tsv.gz")),
            file.path(sample_dir, "features.tsv.gz"), overwrite = TRUE)
  file.copy(file.path(base_path, paste0(sample_id, "_matrix.mtx.gz")),
            file.path(sample_dir, "matrix.mtx.gz"), overwrite = TRUE)

  counts <- Read10X(sample_dir)

  valid_barcodes <- intersect(colnames(counts), samples_clean$BARCODE)
  counts_filtered <- counts[, valid_barcodes]

  pool_obj <- CreateSeuratObject(counts = counts_filtered, project = sample_id)

  patient_map <- setNames(samples_clean$patient, samples_clean$BARCODE)
  patient_vec <- patient_map[colnames(pool_obj)]
  names(patient_vec) <- colnames(pool_obj)
  pool_obj <- AddMetaData(pool_obj, metadata = patient_vec, col.name = "patient")

  pool_obj <- subset(pool_obj, subset = !is.na(patient))

  return(pool_obj)
}

pool_objects <- lapply(pooled_samples, process_pool, base_path = base_path)
names(pool_objects) <- pooled_samples

combined_pooled <- merge(
  pool_objects[[1]],
  y = pool_objects[-1],
  add.cell.ids = pooled_samples
)

saveRDS(combined_pooled, "E:/Sepsis Project/Sepsis.Data/combined_pooled_14patients.rds")

# ============================================================
# 3. PROCESS INDIVIDUAL SAMPLES (32 samples)
# ============================================================

individual_rds_dir <- "E:/Sepsis Project/Sepsis.Data/individual_rds"
dir.create(individual_rds_dir, showWarnings = FALSE)

individual_samples <- c(
  "GSM4995431_MVIR1-HS12-D0BLD1-SCG1", "GSM4995432_MVIR1-HS14-D0BLD1-SCG1",
  "GSM4995433_MVIR1-HS16-D0BLD1-SCG1", "GSM4995434_MVIR1-HS25-D0BLD1-SCG1",
  "GSM4995435_MVIR1-HS26-D0BLD1-SCG1", "GSM4995436_MVIR1-HS29-D0BLD1-SCG1",
  "GSM4995437_MVIR1-HS31-D0BLD1-SCG1", "GSM4995438_MVIR1-HS38-D0BLD1-SCG1",
  "GSM4995439_MVIR1-HS44-D0BLD1-SCG1", "GSM4995440_MVIR1-HS45-D0BLD1-SCG1",
  "GSM4995441_MVIR1-HS46-D0BLD1-SCG1", "GSM4995442_MVIR1-HS47-D0BLD1-SCG1",
  "GSM4995443_MVIR1-HS49-D0BLD1-SCG1", "GSM4995444_MVIR1-HS52-D0BLD1-SCG1",
  "GSM4995445_MVIR1-HS54-D0BLD1-SCG1", "GSM4995446_MVIR1-HS55-D0BLD1-SCG1",
  "GSM4995447_MVIR1-HS60-D0BLD1-SCG1", "GSM4995448_MVIR1-HS62-D0BLD1-SCG1",
  "GSM4995449_XHLT1-HS3-BLD1-SCG1",  "GSM4995450_XHLT1-HS4-BLD1-SCG1",
  "GSM4995451_XHLT1-HS5-BLD1-SCG1",  "GSM4995452_XHLT1-HS6-BLD2-SCG1",
  "GSM4995453_XHLT1-HS7-BLD2-SCG1",  "GSM4995454_XHLT1-HS14-BLD1-SCG1",
  "GSM4995455_XHLT1-HS19-BLD1-SCG1", "GSM4995456_XHLT1-HS20-BLD1-SCG1",
  "GSM4995457_XHLT1-HS23-BLD1-SCG1", "GSM4995458_XHLT1-HS29-BLD1-SCG1",
  "GSM4995459_XHLT1-HS30-BLD1-SCG1", "GSM4995460_XHLT1-HS31-BLD1-SCG1",
  "GSM4995461_XHLT1-HS32-BLD1-SCG1", "GSM4995462_XHLT1-HS33-BLD1-SCG1"
)

process_and_save_individual_v2 <- function(sample_id, base_path, meta, out_dir) {
  message("Processing: ", sample_id)
  gsm_id <- sub("_.*", "", sample_id)
  sample_dir <- file.path(base_path, sample_id)

  counts <- Read10X(sample_dir)

  obj <- CreateSeuratObject(counts = counts, project = sample_id,
                             min.cells = 3, min.features = 200)

  obj$patient <- sub("^[^_]*_", "", sample_id)
  obj$group <- ifelse(grepl("MVIR1", sample_id), "Patient", "Control")
  obj$covid_status <- meta[gsm_id, "covid_status:ch1"]
  obj$severity <- meta[gsm_id, "phenotype:ch1"]

  out_path <- file.path(out_dir, paste0(sample_id, ".rds"))
  saveRDS(obj, out_path)

  rm(counts, obj)
  gc()
  return(out_path)
}

for (s in individual_samples) {
  process_and_save_individual_v2(s, base_path, meta, individual_rds_dir)
}

# merge individual samples (memory-efficient cbind approach)
rds_files <- list.files(individual_rds_dir, pattern = "\\.rds$", full.names = TRUE)

counts_list <- list()
meta_list <- list()
for (i in seq_along(rds_files)) {
  obj <- readRDS(rds_files[i])
  counts_list[[i]] <- GetAssayData(obj, layer = "counts")
  meta_list[[i]] <- obj@meta.data
  rm(obj); gc()
}
names(counts_list) <- individual_samples
names(meta_list) <- individual_samples

common_genes <- Reduce(intersect, lapply(counts_list, rownames))
counts_list <- lapply(counts_list, function(m) m[common_genes, ])

for (i in seq_along(counts_list)) {
  colnames(counts_list[[i]]) <- paste0(names(counts_list)[i], "_", colnames(counts_list[[i]]))
}
for (i in seq_along(meta_list)) {
  rownames(meta_list[[i]]) <- paste0(names(meta_list)[i], "_", rownames(meta_list[[i]]))
}

combined_counts <- do.call(cbind, counts_list)
combined_meta <- do.call(rbind, meta_list)
rownames(combined_meta) <- colnames(combined_counts)

combined_individual <- CreateSeuratObject(counts = combined_counts, meta.data = combined_meta)
rm(combined_counts, counts_list); gc()

saveRDS(combined_individual, "E:/Sepsis Project/Sepsis.Data/combined_individual_32samples.rds")

# ============================================================
# 4. MERGE POOLED + INDIVIDUAL INTO FINAL OBJECT
# ============================================================

combined_pooled <- JoinLayers(combined_pooled)

pooled_counts <- GetAssayData(combined_pooled, layer = "counts")
pooled_meta   <- combined_pooled@meta.data
individual_counts <- GetAssayData(combined_individual, layer = "counts")
individual_meta   <- combined_individual@meta.data

pooled_meta$group <- "Patient"
pooled_meta$severity <- NA
pooled_meta$covid_status <- NA

common_cols <- c("patient", "group", "severity", "covid_status")
pooled_meta <- pooled_meta[, common_cols]
individual_meta <- individual_meta[, common_cols]

common_genes_final <- intersect(rownames(pooled_counts), rownames(individual_counts))
pooled_counts <- pooled_counts[common_genes_final, ]
individual_counts <- individual_counts[common_genes_final, ]

final_counts <- cbind(pooled_counts, individual_counts)
final_meta <- rbind(pooled_meta, individual_meta)
rownames(final_meta) <- colnames(final_counts)

sepsis_full <- CreateSeuratObject(counts = final_counts, meta.data = final_meta)
rm(pooled_counts, individual_counts, final_counts); gc()

saveRDS(sepsis_full, "E:/Sepsis Project/Sepsis.Data/sepsis_full_46people.rds")

# ============================================================
# 5. QUALITY CONTROL
# ============================================================

sepsis_full[["percent.mt"]] <- PercentageFeatureSet(sepsis_full, pattern = "^MT-")

sepsis_full <- subset(
  sepsis_full,
  subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10
)

library(SingleCellExperiment)
library(scDblFinder)

sce_full <- as.SingleCellExperiment(sepsis_full)
set.seed(100)
sce_full <- scDblFinder(sce_full, samples = "patient")

sepsis_full$doublet_score <- colData(sce_full)$scDblFinder.score
sepsis_full$doublet_class <- colData(sce_full)$scDblFinder.class
rm(sce_full); gc()

sepsis_full <- subset(sepsis_full, subset = doublet_class == "singlet")

saveRDS(sepsis_full, "E:/Sepsis Project/Sepsis.Data/sepsis_full_QC_complete.rds")

# ============================================================
# 6. NORMALIZATION, FEATURE SELECTION, PCA
# ============================================================

sepsis_full <- NormalizeData(sepsis_full, normalization.method = "LogNormalize", scale.factor = 10000)
sepsis_full <- FindVariableFeatures(sepsis_full, selection.method = "vst", nfeatures = 2000)
sepsis_full <- ScaleData(sepsis_full)
sepsis_full <- RunPCA(sepsis_full, features = VariableFeatures(sepsis_full))

ElbowPlot(sepsis_full, ndims = 50)

# ============================================================
# 7. INTEGRATION (HARMONY)
# ============================================================

library(harmony)
sepsis_full <- RunHarmony(sepsis_full, group.by.vars = "patient")

saveRDS(sepsis_full, "E:/Sepsis Project/Sepsis.Data/sepsis_full_integrated.rds")

# ============================================================
# 8. UMAP + CLUSTERING
# ============================================================

sepsis_full <- RunUMAP(sepsis_full, reduction = "harmony", dims = 1:15)
sepsis_full <- FindNeighbors(sepsis_full, reduction = "harmony", dims = 1:15)
sepsis_full <- FindClusters(sepsis_full, resolution = 0.5)

saveRDS(sepsis_full, "E:/Sepsis Project/Sepsis.Data/sepsis_full_clustered.rds")

DimPlot(sepsis_full, reduction = "umap", label = TRUE) + NoLegend()
DimPlot(sepsis_full, reduction = "umap", group.by = "patient") + NoLegend()

# ============================================================
# 9. CLUSTER MARKERS + ANNOTATION
# ============================================================

cluster_markers <- FindAllMarkers(
  sepsis_full, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25
)

top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

library(writexl)
write_xlsx(top_markers, "E:/Sepsis Project/Sepsis.Data/top_markers_per_cluster_integrated.xlsx")

# remove low-quality / ambiguous clusters
sepsis_full <- subset(
  sepsis_full,
  subset = seurat_clusters %in% c("8", "12", "14", "17"),
  invert = TRUE
)

cluster_names <- c(
  "0"  = "Neutrophils",
  "1"  = "Naive T cells",
  "2"  = "Classical Monocytes",
  "3"  = "Memory CD4+ T cells",
  "4"  = "CD8+ Cytotoxic T cells",
  "5"  = "Monocytes",
  "6"  = "Platelets",
  "7"  = "NK cells",
  "9"  = "B cells",
  "10" = "Erythrocytes",
  "11" = "Non-classical Monocytes",
  "13" = "Plasma cells",
  "15" = "Proliferating cells",
  "16" = "Dendritic cells"
)

sepsis_full$cell_type <- unname(cluster_names[as.character(sepsis_full$seurat_clusters)])

saveRDS(sepsis_full, "E:/Sepsis Project/Sepsis.Data/sepsis_full_annotated_clean.rds")

DimPlot(sepsis_full, reduction = "umap", group.by = "cell_type", label = TRUE) + NoLegend()

# ============================================================
# 10. CELL COMPOSITION: PATIENT vs CONTROL
# ============================================================

cell_composition <- sepsis_full@meta.data %>%
  group_by(group, cell_type) %>%
  summarise(n = n()) %>%
  mutate(percentage = n / sum(n) * 100)

library(ggplot2)
ggplot(cell_composition, aes(x = cell_type, y = percentage, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Immune Cell Composition: Patient vs Control",
       y = "Percentage of Cells (%)", x = "Cell Type")

composition_table <- table(sepsis_full$cell_type, sepsis_full$group)
chisq_result <- chisq.test(composition_table)
chisq_result

# ============================================================
# 11. DIFFERENTIAL EXPRESSION PER CELL TYPE (Patient vs Control)
# ============================================================

Idents(sepsis_full) <- "cell_type"

cell_types_to_test <- c(
  "Neutrophils", "Classical Monocytes", "Monocytes", "Naive T cells",
  "Memory CD4+ T cells", "CD8+ Cytotoxic T cells", "NK cells", "B cells"
)

run_de <- function(ct) {
  de <- tryCatch({
    FindMarkers(
      sepsis_full, ident.1 = "Patient", ident.2 = "Control",
      group.by = "group", subset.ident = ct,
      min.pct = 0.25, logfc.threshold = 0.25
    )
  }, error = function(e) NULL)
  if (is.null(de)) return(NULL)
  de <- cbind(gene = rownames(de), de)
  rownames(de) <- NULL
  return(de)
}

de_results <- lapply(cell_types_to_test, run_de)
names(de_results) <- cell_types_to_test
de_results <- de_results[!sapply(de_results, is.null)]

write_xlsx(de_results, "E:/Sepsis Project/Sepsis.Data/DE_all_celltypes_Patient_vs_Control.xlsx")

# ============================================================
# 12. PATHWAY (GO) ENRICHMENT
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Classical Monocytes
monocyte_de <- de_results[["Classical Monocytes"]]
monocyte_genes_up <- monocyte_de %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>% pull(gene)

go_results <- enrichGO(
  gene = monocyte_genes_up, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05
)

dotplot(go_results, showCategory = 15) +
  ggtitle("Pathway Enrichment: Classical Monocytes (Patient vs Control)")
ggsave("E:/Sepsis Project/Sepsis.Data/pathway_monocytes_dotplot.png", width = 8, height = 6, dpi = 300)
write_xlsx(as.data.frame(go_results), "E:/Sepsis Project/Sepsis.Data/GO_pathways_Monocytes.xlsx")

# Neutrophils
neutrophil_de <- de_results[["Neutrophils"]]
neutrophil_genes_up <- neutrophil_de %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>% pull(gene)

go_results_neutro <- enrichGO(
  gene = neutrophil_genes_up, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05
)

dotplot(go_results_neutro, showCategory = 15) +
  ggtitle("Pathway Enrichment: Neutrophils (Patient vs Control)")
ggsave("E:/Sepsis Project/Sepsis.Data/pathway_neutrophils_dotplot.png", width = 8, height = 6, dpi = 300)
write_xlsx(as.data.frame(go_results_neutro), "E:/Sepsis Project/Sepsis.Data/GO_pathways_Neutrophils.xlsx")

# ============================================================
# 13. SEVERITY ANALYSIS (Mild/Moderate vs Severe)
# ============================================================

severity_composition <- sepsis_full@meta.data %>%
  filter(!is.na(severity), severity != "CONTROL") %>%
  group_by(severity, cell_type) %>%
  summarise(n = n()) %>%
  mutate(percentage = n / sum(n) * 100)

ggplot(severity_composition, aes(x = cell_type, y = percentage, fill = severity)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Immune Cell Composition: Mild/Moderate vs Severe",
       y = "Percentage of Cells (%)", x = "Cell Type")

severity_table <- table(
  sepsis_full$cell_type[!is.na(sepsis_full$severity) & sepsis_full$severity != "CONTROL"],
  sepsis_full$severity[!is.na(sepsis_full$severity) & sepsis_full$severity != "CONTROL"]
)
chisq_severity <- chisq.test(severity_table)
chisq_severity

# ============================================================
# 14. FINAL SAVE
# ============================================================

saveRDS(sepsis_full, "E:/Sepsis Project/Sepsis.Data/sepsis_full_FINAL.rds")

write_xlsx(
  list(
    "Composition_Patient_vs_Control" = cell_composition,
    "Composition_Severity" = severity_composition,
    "Chisq_Patient_Control" = data.frame(statistic = chisq_result$statistic, p_value = chisq_result$p.value),
    "Chisq_Severity" = data.frame(statistic = chisq_severity$statistic, p_value = chisq_severity$p.value)
  ),
  "E:/Sepsis Project/Sepsis.Data/Final_Statistical_Results.xlsx"
)
