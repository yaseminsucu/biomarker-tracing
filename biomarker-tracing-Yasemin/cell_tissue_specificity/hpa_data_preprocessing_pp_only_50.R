# Install these packages first before running the code
library("seismicGWAS")
library("plyr")
library("Seurat")
library("qs")
library("scRNAseq")
library("dplyr")
library("data.table")
options(future.globals.maxSize = 10000 * 1024^2)


# Load the data (replace the path to the path of the human protein atlas data)
path = "/hpc/users/sucuy01/Yasemin_projects/gnpc_biomarker/updated_biomarker/biomarker-tracing/tutorial/atlas_inputs/merged_human_protein_atlas_single_cell.qs"
obj = qread(path)

######################### SINGLE-CELL OBJECT PREPROCESSING #########################
# Create cell_tissue ident
obj[[]]$tissue = gsub("_", ".", obj[[]]$tissue)
obj[[]]$cell_type = gsub("_", ".", obj[[]]$cell_type)
obj[[]] <- obj[[]] %>%
  mutate(
    cell_tissue = paste(tissue, cell_type, sep = "_")
  )
ident = "cell_tissue"
Idents(obj) = ident

# Subset the dataset to reduce data size
obj = subset(x = obj, downsample = 50)

# Find variable genes
DefaultAssay(obj) = "RNA"
obj = NormalizeData(obj)
obj = FindVariableFeatures(obj, nfeatures = 10000)

# Only retain variable genes and subset the dataset to reduce data size
plasma_proteins = read.csv("/hpc/users/sucuy01/Yasemin_projects/gnpc_biomarker/updated_biomarker/biomarker-tracing/gnpc_sumstats/AD.csv")$gene
var_genes = union(VariableFeatures(obj), plasma_proteins)
var_genes = intersect(var_genes, rownames(obj))
arr = LayerData(obj, assay = "RNA", layer = "counts")[var_genes,]
obj_smal = CreateSeuratObject(
    arr, assay = "RNA", meta.data = obj[[]]
)
obj_smal = NormalizeData(obj_smal)

# Reformat the string
obj_smal[[]] <- obj_smal[[]] %>% dplyr::mutate(cell_tissue = gsub("[- ]", ".", cell_tissue))

# Save the object
save_path <- "/hpc/users/sucuy01/Yasemin_projects/gnpc_biomarker/updated_biomarker/biomarker-tracing/tutorial/atlas_outputs"
dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
saveRDS(obj_smal, file.path(save_path, "human_protein_atlas_single_cell_pp_50_small.rds"))
