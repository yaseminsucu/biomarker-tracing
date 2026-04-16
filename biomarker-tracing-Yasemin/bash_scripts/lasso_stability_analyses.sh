#!/usr/bin/env bash
set -euo pipefail

# Hard-pin python from your conda env (works in LSF non-interactive shells)
PY=/hpc/users/sucuy01/.conda/envs/scanpy-env/bin/python

# (Optional but highly recommended) log what we're using
echo "Python being used:"
$PY -c "import sys, sklearn; print(sys.executable); print('sklearn', sklearn.__version__)"

# Get the names of params
script_dir="python_main_cell_type_spec_method/lasso_stability_analyses.py"
placeholder=$1
atlas_smal_path=$2
prot_data_path=$3
save_path=$4
disease=$5
output_label=$6

num_alphas=$7
subset_size=$8
gene_weight=$9
num_iter=${10}
samp_method=${11}
optim_alpha_mode=${12}
ztransform_type=${13}

# Version of the code when splitting the data by cell-tissue pair
$PY $script_dir \
--atlas_smal_path $atlas_smal_path \
--prot_data_path $prot_data_path \
--save_path $save_path \
--disease $disease \
--num_alphas $num_alphas \
--subset_size $subset_size \
--num_iter $num_iter \
--samp_method $samp_method \
--optim_alpha_mode $optim_alpha_mode \
--gene_weight $gene_weight \
--output_label $output_label \
--ztransform_type $ztransform_type
