#!/bin/bash
# this is intended for running DPS jobs; the input directory is where four files have been pulled because download=TRUE in the algorithm_config.yaml file
# a tar file of biomass models, a data table csv, and two raster stack geotiff files

#source activate icesat2_boreal
basedir=$( cd "$(dirname "$0")" ; pwd -P )

#unset PROJ_LIB

#pip install --user -r ${basedir}/requirements.txt

mkdir output

# Note: the numbered args are fed in with the in_param_dict in the Run DPS chunk of 3.4_dps.ipynb
ATL08_tindex_master_fn=${1}
TOPO_TIF=${2}
HLS_TIF=${3}
LC_TIF=${4}
DO_SLOPE_VALID_MASK=${5}
ATL08_SAMPLE_CSV=${6}
in_tile_num=${7}
in_tile_fn=${8}
in_tile_field=${9}
iters=${10}
calculate_uncertainty=${11}
minDOY=${12}
maxDOY=${13}
max_sol_el=${14}
expand_training=${15}
local_train_perc=${16}
min_n=${17}
boreal_vect_fn=${18}
predict_var=${19}
max_n=${20}
pred_vars=${21}
bio_models_tar_fn=${22}
SAR_TIF=${23}
year=${24}

source activate icesat2_boreal

tar -xf ${bio_models_tar_fn}

# This PWD is wherever the job is run (where the .sh is called from) 
OUTPUTDIR="${PWD}/output"

# Get the output merged CSV of filtered ATL08 for the input tile and its neighbors
cmd="python ${basedir}/../../lib/merge_neighbors_atl08.py -in_tile_num ${in_tile_num} -in_tile_fn ${in_tile_fn} -in_tile_field ${in_tile_field} -csv_list_fn ${ATL08_tindex_master_fn} -out_dir ${OUTPUTDIR}"

echo $cmd
eval $cmd

# Set the output merged CSV name to a var
MERGED_ATL08_CSV=$(ls ${OUTPUTDIR}/*_merge_neighbors_*.csv | head -1)

echo $MERGED_ATL08_CSV
echo $ATL08_SAMPLE_CSV

source activate r

# required arguments
args=(--atl08_path "${MERGED_ATL08_CSV}")
args+=(--broad_path "${6}")
args+=(--topo_path "${2}")
args+=(--hls_path "${3}")
args+=(--lc_path "${4}")
args+=(--boreal_vector_path "${18}")
args+=(--year "${24}")

# optional arguments
[[ -n "${23}" ]] && args+=(--sar_path "${23}")
[[ -n "${5}" ]] && args+=(--mask "${5}")
[[ -n "${10}" ]] && args+=(--uncertainty_iterations "${10}")
[[ -n "${11}" ]] && args+=(--calculate_uncertainty "${11}")
[[ -n "${12}" ]] && args+=(--minDOY "${12}")
[[ -n "${13}" ]] && args+=(--maxDOY "${13}")
[[ -n "${14}" ]] && args+=(--max_sol_el "${14}")
[[ -n "${15}" ]] && args+=(--expand_training "${15}")
[[ -n "${16}" ]] && args+=(--local_train_perc "${16}")
[[ -n "${17}" ]] && args+=(--min_samples "${17}")
[[ -n "${20}" ]] && args+=(--max_samples "${20}")
[[ -n "${19}" ]] && args+=(--predict_var "${19}")
[[ -n "${21}" ]] && args+=(--pred_vars "${21}")

mapBoreal_R="Rscript ${basedir}/../../lib/mapBoreal_simple.R"
command=("${mapBoreal_R}" "${args[@]}")
echo "${command}"
"${command[@]}"
