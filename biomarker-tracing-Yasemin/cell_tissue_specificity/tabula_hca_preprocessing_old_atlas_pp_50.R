library("seismicGWAS")
library("plyr")
library("Seurat")
library("qs")
library("scRNAseq")
library("dplyr")
library("data.table")
options(future.globals.maxSize = 10000 * 1024^2)

obj = qread("/hpc/users/sucuy01/Yasemin_projects/gnpc_biomarker/gnpc_run_2/gnpc_updated_pipeline/tabula_sapiens_cortical_cells_atlas_all-genes_downsamp.qs")

######################### SINGLE-CELL OBJECT PREPROCESSING #########################
obj[[]]$tissue = gsub("[_ ]", ".", obj[[]]$tissue)
obj[[]]$cell_type = gsub("[_, /():-]", ".", obj[[]]$cell_type)
obj[[]] <- obj[[]] %>%
  mutate(cell_tissue = paste(tissue, cell_type, sep = "_"))
Idents(obj) = "cell_tissue"

obj = subset(x = obj, downsample = 50)

DefaultAssay(obj) = "RNA"
obj = NormalizeData(obj)

plasma_proteins = read.csv("/hpc/users/sucuy01/Yasemin_projects/gnpc_biomarker/updated_biomarker/biomarker-tracing/gnpc_sumstats/AD.csv")$gene
var_genes = intersect(plasma_proteins, rownames(obj))
cat("Plasma proteins found in atlas:", length(var_genes), "\n")
arr = LayerData(obj, assay = "RNA", layer = "counts")[var_genes,]
obj_smal = CreateSeuratObject(arr, assay = "RNA", meta.data = obj[[]])
obj_smal = NormalizeData(obj_smal)

obj_smal[[]] <- obj_smal[[]] %>% dplyr::mutate(cell_tissue = gsub("[- ]", ".", cell_tissue))

save_path <- "/hpc/users/sucuy01/Yasemin_projects/gnpc_biomarker/updated_biomarker/biomarker-tracing/tutorial/atlas_outputs"
dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
saveRDS(obj_smal, file.path(save_path, "tabula_sapiens_old_atlas_pp_50_small.rds"))
