#!/usr/bin/env bash
set -euo pipefail

# Hard-pin python from your conda env (works in LSF non-interactive shells)
PY=/hpc/users/sucuy01/.conda/envs/scanpy-env/bin/python

# (Optional but highly recommended) log what we're using
echo "Python being used:"
$PY -c "import sys, sklearn; print(sys.executable); print('sklearn', sklearn.__version__)"

# Get the names of params
script_dir="python_main_cell_type_spec_method/univar_association_testing.py"
placeholder=$1
atlas_smal_path=$2
prot_data_path=$3
save_path=$4
disease=$5
abs_hr=$6
output_label=$7
covar_df=$8
covar_gini=$9
ztransform_type=${10}

# Version of the code when splitting the data by cell-tissue pair
$PY $script_dir \
--atlas_smal_path $atlas_smal_path \
--prot_data_path $prot_data_path \
--save_path $save_path \
--disease $disease \
--abs_hr $abs_hr \
--output_label $output_label \
--covar_df $covar_df \
--covar_gini $covar_gini \
--ztransform_type $ztransform_type
