#!/bin/bash
# this is intended for running DPS jobs; the input directory is where four files have been pulled because download=TRUE in the algorithm_config.yaml file
# a tar file of biomass models, a data table csv, and two raster stack geotiff files

#source activate icesat2_boreal
basedir=$( cd "$(dirname "$0")" ; pwd -P )
libdir=$(dirname "$(dirname "${basedir}")")/lib

mkdir output
source activate icesat2_boreal
echo "basedir=${basedir}, bio_models=${22}, outdir=${OUTPUTDIR}, pwd=${PWD}, home=${HOME}"

# This PWD is wherever the job is run (where the .sh is called from) 
OUTPUTDIR="${PWD}/output"

# required arguments
args=(--in_tile_num "${7}")
args+=(--in_tile_fn "${8}")
args+=(--out_dir "${OUTPUTDIR}")

# optional arguments
[[ -n "${1}" ]] && args+=(--csv_list_fn "${1}")
[[ -n "${9}" ]] && args+=(--in_tile_field "${9}")
[[ -n "${29}" ]] && args+=(--atl08_year_list "${29}")

command=(python "${libdir}/merge_neighbors_atl08.py" "${args[@]}")
echo "${command[@]}"
"${command[@]}"

echo $command
eval $command

# Set the output merged CSV name to a var
MERGED_ATL08_CSV=$(ls ${OUTPUTDIR}/*_merge_neighbors_*.csv | head -1)

echo $MERGED_ATL08_CSV
echo $ATL08_SAMPLE_CSV

source activate r

# required arguments
args=(--atl08_path "${MERGED_ATL08_CSV}")
args+=(--broad_path "${6}") # should be optional
args+=(--topo_path "${2}")
args+=(--hls_path "${3}")
args+=(--lc_path "${4}")
args+=(--boreal_vector_path "${17}") # should be optional
args+=(--ecoregions_path "${18}") # should be optional
args+=(--countries_path "${30}") # should be optional
args+=(--biomass_models_path "${22}")
args+=(--year "${24}")

# optional arguments
[[ -n "${5}" ]] && args+=(--mask "${5}")
[[ -n "${10}" ]] && args+=(--n_iters "${10}")
[[ -n "${11}" ]] && args+=(--minDOY "${11}")
[[ -n "${12}" ]] && args+=(--maxDOY "${12}")
[[ -n "${13}" ]] && args+=(--max_sol_el "${13}")
[[ -n "${14}" ]] && args+=(--expand_training "${14}")
[[ -n "${15}" ]] && args+=(--local_train_perc "${15}")
[[ -n "${16}" ]] && args+=(--min_samples "${16}")
[[ -n "${19}" ]] && args+=(--predict_var "${19}")
[[ -n "${20}" ]] && args+=(--max_samples "${20}")
[[ -n "${21}" ]] && args+=(--pred_vars "${21}")
[[ -n "${23}" ]] && args+=(--sar_path "${23}")
[[ -n "${25}" ]] && args+=(--cores "${25}")
[[ -n "${26}" ]] && args+=(--ntree "${26}")
[[ -n "${27}" ]] && args+=(--zero_short_veg_height "${27}")
[[ -n "${28}" ]] && args+=(--slope_thresh "${28}")

command=(Rscript "${libdir}/mapBoreal_simple.R" "${args[@]}")
echo "${command[@]}"
"${command[@]}"
