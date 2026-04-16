#!/usr/bin/env bash
set -euo pipefail

# Hard-pin python from your conda env (works in LSF non-interactive shells)
PY=/hpc/users/sucuy01/.conda/envs/scanpy-env/bin/python

# (Optional but highly recommended) log what we're using
echo "Python being used:"
$PY -c "import sys, sklearn; print(sys.executable); print('sklearn', sklearn.__version__)"

# Get the names of params
script_dir="python_main_cell_type_spec_method/tree_based_methods.py"
placeholder=$1
atlas_smal_path=$2
prot_data_path=$3
save_path=$4
disease=$5
output_label=$6

num_trees=$7
min_samples_split=$8
min_samples_leaf=$9
max_samples=${10}

kfold_n=${11}
n_permute_repeat=${12}
param_search=${13}
abs_hr=${14}
ztransform_type=${15}

# Version of the code when splitting the data by cell-tissue pair
$PY $script_dir \
--atlas_smal_path $atlas_smal_path \
--prot_data_path $prot_data_path \
--save_path $save_path \
--disease $disease \
--param_search $param_search \
--num_trees $num_trees \
--min_samples_split $min_samples_split \
--min_samples_leaf $min_samples_leaf \
--max_samples $max_samples \
--kfold_n $kfold_n \
--n_permute_repeat $n_permute_repeat \
--output_label $output_label \
--abs_hr $abs_hr \
--ztransform_type $ztransform_type
