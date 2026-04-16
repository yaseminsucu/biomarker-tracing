#Scripts and Results from the March 12th run


This branch contains the minimally modified scripts of the biomarker-tracing pipeline with the run results included. 

Data used: Pheneome Proteome Atlas and GNPC data for Alzheimer's Disease 
Atlasses used: Human Protein Cell Atlas and Tabula Sapiens merged with Human Cortical Cell Atlas 

The focus is on the **pp_50 configuration**, where plasma protein genes are incorporated into atlas preprocessing.


---

#Folders 

### `atlas_files/`

Processed atlas objects and corresponding cell-type specificity matrices.

Includes two atlas sources:

**Human Protein Atlas (HPA)**
- `human_protein_atlas_single_cell_pp_50_small.rds` (**pp refers to plasma proteins**, 50 refers to downsampling size)
- `human_protein_atlas_seismic_cell_spec_matrix_pp_50.tsv`

**Tabula Sapiens + Human Cortical Atlas (HCA)**
- `tabula_sapiens_old_atlas_pp_50_small.rds`
- `tabula_sapiens_old_atlas_pp_50_seismic_cell_spec_matrix.tsv`

---

### `cell_tissue_specificity/`

These scripts implement the **pp_50 modification**.

---

### `pipeline_yml/`

pipeline runs:

- `univar_multivar_hpa_pp_50.yml` → HPA-based analysis
- `univar_multivar_ts_hca_pp_50.yml` → Tabula + HCA analysis


---

### `results/`

Contains final outputs from the AD analysis.

Includes:

- `pipeline_univar_multivar_hpa_pp_50/`
- `pipeline_univar_multivar_ts_hca_pp_50/`

Each pipeline directory contains:
- univariate association results
- elastic net models
- random forest outputs
- stability analyses (skip this one since it doesn't generated the results due to packing issue)

Also includes summary plots:
- `hpa_default_vs_pp_comparison.png`
- `hpa_50_vs_200_comparison.png` (comparison of the downsampling size effect with 50 vs 200, and went with 50 downsampling like the default version)
- `ts_hca_default_vs_pp_comparison.png`

---

### `notebooks/`

Run result exploration and comparison.

- `results_gnpc_phenome_p_HPA_HPC.ipynb`  
  → analysis comparing default pipeline vs plasma protein–enhanced run

---

## Notes

- This package corresponds specifically to the **March 12 AD-only run**
- Core pipeline logic remains unchanged. Minor adjustments were made for environment compatibility (e.g., file paths and data loading). The primary modification is in atlas preprocessing, where plasma protein genes were incorporated into variable feature selection.

---
