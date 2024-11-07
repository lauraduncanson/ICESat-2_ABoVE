# This script was developed by Laura Duncanson to produce tiled boreal biomass 30-m maps with inputs from Carlos A Silva, Alex Mandel, and Paul Montesano.
# Inputs are stacks of Landsat and Copernicus DEM data, and tables of linked 30-m ATL08 data to coincident stack attrbitutes

#----------------------------------------------#

# 3.4 ICESat-2 biomass 30m ATL08

#### i) Description
#### Algorithm for tile-based AGBD mapping with ICESat-2, Copernicus DEM, and Landsat/HLS composites. This work was funded by NASA's ABoVE program, PI Duncanson, Co-Is Paul Montesano, Amy Neuenschwander, and NASA's ICESat-2 Science Team (PI Neuenschwader, Co-I Duncanson). Inputs and code from Nathan Thomas, Carlos Silva, Eric Guenther, Alex Mandel, and implementation assistance from NASA MAAP team George Change, Sujen Shah, Brian Satorius

#### ii) How the algorithm works?
#### Data tables of linked 30-m ATL08 RH metrics and covariate stack metrics are imported (outputs from tile_atl08.py)
#### AGBD models (externally developed) are loaded in R and applied over the ICESat-2 30m ATL08 data to data tables
#### A set of random forest models are fit to predict 30m ICESat-2 biomass as a function of covariates per tile
#### Raster stacks of covariates are sub-tiled to reduce memory usage
#### The suite of rf models are applied to each tile, mean and SD are output
#### Subtiles are recombined and written to disc as cogs

#### iii) Inputs
####  - rds_models: list of ICESat-2 simulation-derived AGB model paths
####  - stack: a combined raster stack of Landsat and Copernicus DEM data
####  - ice2_30_atl08: list containing the path to the data tables
####  - offset: offset applied in the model

#### iii) Outputs
####  COGs of predicted 30-m AGBD, SD AGBD

get_height_column_names <- function(in_data){
  return(
    names(in_data)[grep('^RH_[0-9]{2}$', names(in_data))]
  )
}

rename_height_columns_to_match_pretrained_models <- function(in_data){
  return(
    in_data |>
      rename_with(~gsub('\\brh([0-9]{2})\\b', 'RH_\\1', .x), matches='^rh[0-9]{2}$') |>
      rename(RH_98=h_canopy)
  )
}

offset_RH_columns <- function(all_train_data, offset){
  RH_columns <- get_height_column_names(all_train_data)

  return(
    all_train_data |>
      mutate(across(RH_columns, ~. + offset))
  )
}

set_model_id_for_AGB_prediction <- function(in_data, offset){
  # TODO: uncomment correct model ids once tested against old results
  return(
    in_data |>
      mutate(model_id = case_when(
        segment_landcover %in% c(111, 113, 121, 123) ~ "m8",# "m3",
        segment_landcover %in% c(112, 114, 122, 124) ~ "m8", # "m1",
        TRUE ~ "m8"
      )
      )
  )
}

randomize_AGB_model <- function(model){
  # modify coeffients through sampling variance covariance matrix
  model_coeffs <- mvrnorm(n=50, mu=model$coefficients, Sigma=vcov(model))
  model$coefficients <- model_coeffs[1,]

  return(model)
}

GEDI2AT08AGB<-function(rds_models, df, randomize=FALSE, max_n=10000, sample=TRUE){
  if (sample && nrow(df) > max_n)
    df <- reduce_sample_size(df, max_n)

  df$AGB <- NA
  df$SE <- NA

  ids<-unique(df$model_id)
  n_models <- length(ids)

  for (i in ids){
    model_i <- rds_models[[i]]

    # Modify coeffients through sampling variance covariance matrix
    if(randomize)
      model_i <- randomize_AGB_model(model_i)

    # Predict AGB and SE
    df$AGB[df$model_id==i] <- predict(model_i, newdata=df[df$model_id==i,])
    df$SE[df$model_id==i] <- summary(model_i)$sigma^2

    df$AGB[df$AGB < 0] <- 0.0

    # Calculate Correction Factor C
    C <- mean(model_i$fitted.values^2)/mean(model_i$model$`sqrt(AGBD)`^2)

    # Bias correction in case there is a systematic over or under estimation in the model
    df$AGB[df$model_id==i] <- C*(df$AGB[df$model_id==i]^2)
  }

  # Apply slopemask, validmask and landcover masks
  bad_lc <- c(0, 60, 80, 200, 50, 70)
  df$AGB[df$slopemask == 0 |
                df$ValidMask == 0 |
                df$segment_landcover %in% bad_lc] <- 0.0
  return(df)
}

partial_sd <- function(arr){
  partial_sd_arr <- rep(0, length(arr) - 1)
  for(i in 2:length(arr)){
    partial_sd_arr[i-1] <- sd(arr[1:i], na.rm = T)
  }
  return(partial_sd_arr)
}

sd_change_relative_to_baseline <- function(arr, last_n){
  partial_sd_arr <- partial_sd(arr)
  paritial_std_arr_last_n_out <- head(partial_sd_arr, max(1, length(arr) - last_n))

  baseline_sd <- mean(paritial_std_arr_last_n_out, na.rm=T)
  full_sd <-  mean(partial_sd_arr, na.rm=T)

  if (baseline_sd)
    relative_sd_change <-  abs(full_sd - baseline_sd) / baseline_sd
  else
    relative_sd_change <-  Inf

  return(relative_sd_change)
}

DOY_and_solar_filter <- function(tile_data, start_DOY, end_DOY, solar_elevation){
  filter <- which(tile_data$doy >= start_DOY & tile_data$doy <= end_DOY & tile_data$solar_elevation <= solar_elevation)
  return(filter)
}

late_season_filter <- function(tile_data, minDOY, maxDOY, min_icesat2_samples, max_sol_el){

  maxDOY_in_data = max(tile_data$doy)

  for(late_months in 0:3){

    late_season_DOY = maxDOY + 30 * late_months
    filter <- DOY_and_solar_filter(tile_data, minDOY, late_season_DOY, max_sol_el)

    if(late_season_DOY >= maxDOY_in_data | length(filter) >= min_icesat2_samples)
      break
  }

  return(list(filter=filter, late_season_DOY=late_season_DOY))
}

early_and_late_season_filter <- function(tile_data, minDOY, late_season_DOY, min_icesat2_samples, max_sol_el){
  minDOY_in_data = min(tile_data$doy)

  for(early_months in 0:3){

    early_season_DOY = minDOY - 30 * early_months
    filter <- DOY_and_solar_filter(tile_data, early_season_DOY, late_season_DOY, max_sol_el)

    if(early_season_DOY <= minDOY_in_data | length(filter) >= min_icesat2_samples)
      break
  }

  return(list(filter=filter, early_season_DOY=early_season_DOY))
}

expand_training_around_growing_season <- function(tile_data, minDOY, maxDOY, max_sol_el, min_icesat2_samples){

  # first try with no solar elevation
  filter <- DOY_and_solar_filter(tile_data, minDOY, maxDOY, 0)
  if(length(filter) >= min_icesat2_samples){
    print('returning enough data with min and max DOY and 0 elevation...')
    return(tile_data[filter,])
  }

  # next try with max solar elevation
  filter <- DOY_and_solar_filter(tile_data, minDOY, maxDOY, max_sol_el)
  if(length(filter) >= min_icesat2_samples)
    return(tile_data[filter,])

  # next try expanding 1 month later in growing season, iteratively, up to 3 months
  late_season_filter_and_doy <- late_season_filter(tile_data, minDOY, maxDOY, min_icesat2_samples, max_sol_el)
  if(length(late_season_filter_and_doy$filter) >= min_icesat2_samples)
    return(tile_data[late_season_filter_and_doy$filter,])

  # next try expanding 1 month earlier in growing season, iteratively, up to 3 months
  # Note that the upper window might be later in the growing season from the previous call

  early_season_filter_and_doy <- early_and_late_season_filter(tile_data, minDOY, late_season_filter_and_doy$late_season_DOY, min_icesat2_samples, max_sol_el)
  if(length(early_season_filter_and_doy$filter) >= min_icesat2_samples)
    return(tile_data[early_season_filter_and_doy$filter,])

  print("WARNING: min_icesat2_samples condition was not met, applying the extended filter to tile_data anyways")
  return(tile_data[early_season_filter_and_doy$filter,])
}

reduce_sample_size <- function(df, sample_size){

  sample_ids <- seq(1, nrow(df))
  df_sample_ids <- sample(sample_ids, sample_size, replace=FALSE)

  return(df[df_sample_ids,])
}

remove_stale_columns <- function(df, column_names) {

  columns_to_remove <- intersect(names(df), column_names)
  df <- df[, !names(df) %in% columns_to_remove, drop=FALSE]

  return(df)
}

sample_broad_data_within_latitude <- function(broad_data, lat, threshold){

  broad_within_lat <- which(broad_data$lat > (lat-threshold) & broad_data$lat < (lat+threshold))
  broad_data <- broad_data[broad_within_lat,]
  broad_samp_ids <- seq(1, nrow(broad_data))
  broad_sample_ids <- sample(broad_samp_ids, nrow(broad_data), replace=FALSE)
  return(broad_data[broad_sample_ids,])

}

expand_training_with_broad_data <- function(broad_data, tile_data){

  broad_data <- sample_broad_data_within_latitude(broad_data, min(tile_data$lat), 5)
  return(rbind(tile_data, broad_data))
}

remove_height_outliers <- function(all_train_data){
  # remove height outliers based on more than 3SD from the landcover mean

  return(
  all_train_data |>
    group_by(segment_landcover) |>
    summarise(thresh=mean(h_canopy, na.rm=T) + 3 * sd(h_canopy, na.rm=T)) |>
    right_join(all_train_data, by='segment_landcover') |>
    filter(h_canopy <= thresh)
  )
}

set_output_file_names <- function(predict_var, tile_num){
  key <- if (predict_var == 'AGB') 'agb' else 'ht'
  out_fn_stem = paste(
    paste0('output/boreal_', key), format(Sys.time(),"%Y%m%d%s"), str_pad(tile_num, 7, pad = "0"),
    sep="_"
  )

  fn_suffix <- c('.tif', '_summary.csv', '_train_data.csv', '_stats.Rds', '_model.Rds')
  names <- c('map', 'summary', 'train', 'stats', 'model')

  output_file_names <- paste0(out_fn_stem, fn_suffix)
  names(output_file_names) <- names

  return(output_file_names)
}

write_ATL08_table <- function(target, df, out_file_path){
    out_columns <- if(target=='AGB') c('lon', 'lat', 'AGB', 'SE') else c('lon', 'lat', 'h_canopy')
    write.csv(df[, out_columns], file=out_file_path, row.names=FALSE)
}

write_single_model_summary <- function(model, df, target, out_fns){
  target <- if(target == 'AGB') df$AGB else df$RH_98
  local_model <- lm(model$predicted ~ target, na.rm=TRUE)
  saveRDS(model, file=out_fns['model'])

  rsq <- max(model$rsq, na.rm=T)
  cat('rsq_model: ', rsq, '\n')

  rsq_local <- summary(local_model)$r.squared
  cat('rsq_local: ', rsq_local, '\n')

  na_data <- which(is.na(local_model$predicted==TRUE))

  if(length(na_data) == 0)
    rmse_local <- sqrt(mean(local_model$residuals^2))

  cat('rmse_local: ', rmse_local, '\n')

  imp_vars <- model$importance
  out_accuracy <- list(rsq_local, rmse_local, imp_vars)
  saveRDS(out_accuracy, file=out_fns['stats'])
}

write_output_raster_map <- function(maps, output_fn){
  maps <- subst(maps, -9999, NA)
  if (nlyr(maps) > 1){
    out_sd <- app(maps, sd)
    out_sd <- subst(out_sd, -9999, NA)
    maps <- c(subset(maps, 1), out_sd)
  }
  NAflag(maps)
  options <- c("COMPRESS=LZW", overwrite=TRUE, gdal=c("COMPRESS=LZW", "OVERVIEW_RESAMPLING=AVERAGE"))
  writeRaster(maps, filename=output_fn, filetype="COG", gdal=options)
}

write_output_summaries <-function(tile_summaries, boreal_summaries, target, output_fn){
  if (target == 'AGB')
    column_names <- c('tile_total', 'boreal_total')
  else
    column_names <- c('tile_mean', 'boreal_mean')

  df <- data.frame(tile_summaries, boreal_summaries)
  names(df) <- column_names

  write.csv(df, output_fn, row.names=FALSE)
}


prepare_training_data <- function(ice2_30_atl08_path, ice2_30_sample_path, expand_training, minDOY, maxDOY, max_sol_el, min_icesat2_samples, local_train_perc, offset){
  tile_data <- read.csv(ice2_30_atl08_path)

  if(expand_training)
    tile_data <- expand_training_around_growing_season(tile_data, minDOY, maxDOY, max_sol_el, min_icesat2_samples)

  cat('n_avail training:', nrow(tile_data), '\n')
  tile_data <- remove_stale_columns(tile_data, c("binsize", "num_bins"))
  broad_data <- read.csv(ice2_30_sample_path)
  broad_data <- remove_stale_columns(broad_data, c("X__index_level_0__", "geometry"))

  # take proportion of broad data we want based on local_train_perc
  sample_local <- nrow(tile_data) * (local_train_perc/100)
  cat('sample_local:', sample_local, '\n')

  if (sample_local < min_icesat2_samples)
    tile_data <- reduce_sample_size(tile_data, sample_local)

  # sample from broad data to complete sample size
  # this will work if either there aren't enough local samples for n_min OR if there is forced broad sampling
  n_broad <- min_icesat2_samples - nrow(tile_data)
  if(n_broad > 1)
    all_train_data <- expand_training_with_broad_data(broad_data, tile_data)
  else
    all_train_data <- tile_data

  all_train_data <- remove_height_outliers(all_train_data)
  all_train_data <- na.omit(as.data.frame(all_train_data))
  all_train_data <- rename_height_columns_to_match_pretrained_models(all_train_data)
  all_train_data$h_canopy <- all_train_data$RH_98
  all_train_data <- offset_RH_columns(all_train_data, offset)
  all_train_data <- set_model_id_for_AGB_prediction(all_train_data)
  str(all_train_data)
  cat('table for model training generated with ', nrow(all_train_data), ' observations\n')
  return(all_train_data)
}

get_rds_models <- function(){
  rds_model_fns <- list.files(path='~/dps_output/', pattern='*.rds', full.names = TRUE)
  rds_models <- lapply(rds_model_fns, readRDS)
  names(rds_models) <- paste0("m",1:length(rds_models))
  return(rds_models)
}

predict_stack <- function(model, stack){
  stack <- na.omit(stack)
  map <- predict(stack, model, na.rm=TRUE)
  # set slope and valid mask to zero
  # TODO maybe mask can skip over these pixels by default?
  map <- mask(map, stack$slopemask, maskvalues=0, updatevalue=0)
  map <- mask(map, stack$ValidMask, maskvalues=0, updatevalue=0)

  return(map)
}

tile_and_boreal_summary <- function(map, predict_var, boreal_poly, summary_and_convert_functions){
  convert_fun <- summary_and_convert_functions[['convert_fun']]
  summary_fun <- summary_and_convert_functions[['summary_fun']]

  map_conv <- if (is.null(convert_fun)) map else app(map, convert_fun)

  tile_summary <- global(map_conv, summary_fun, na.rm=TRUE)[[summary_fun]]

  boreal <- extract(map_conv, boreal_poly, fun=summary_fun, na.rm=TRUE)
  boreal_summary <- if(summary_fun=='sum') sum(boreal$lyr.1, na.rm=TRUE) else boreal$lyr1[1]

  return(list(tile_summary=tile_summary, boreal_summary=boreal_summary))
}

fit_model <- function(model, model_config, train_df, pred_vars, predict_var){
  y_fit <- if (predict_var == 'Ht') train_df$h_canopy else train_df$AGB
  x_fit <- train_df[pred_vars]

  model_fit <- do.call(model, modifyList(model_config, list(y=y_fit, x=x_fit)))
  return(model_fit)
}

run_modeling_pipeline <-function(rds_models, all_train_data, boreal_poly,
                                 model, model_config, randomize,
                                 summary_and_convert_functions,
                                 max_n, sample, pred_vars, predict_var, stack){
  t1 <- Sys.time()
  print('creating AGB traing data frame.')
  train_df <- GEDI2AT08AGB(rds_models, all_train_data, randomize, max_n, sample)
  print('fitting model')
  model <- fit_model(model, model_config, train_df, pred_vars, predict_var)
  print('predicting biomass map')
  map <- predict_stack(model, stack)
  print('calculating tile and boreal summaries')
  summary <- tile_and_boreal_summary(map, predict_var, boreal_poly, summary_and_convert_functions)
  t2 <- Sys.time()
  cat('pipeline runtime:', t2 - t1, ' (s)\n')
  return(list(
    train_df=train_df, model=model, map=map,
    tile_summary=summary[['tile_summary']], boreal_summary=summary[['boreal_summary']]
  ))
}

get_summary_and_convert_functions <- function(predict_var){
  if (predict_var == 'AGB'){
    summary_fun <- 'sum'
    convert_fun <- function(x){(x*0.09)/1000000000}
  }
  else{
    summary_fun <- 'mean'
    convert_fun <- NULL
  }
  return(list(summary_fun=summary_fun, convert_fun=convert_fun))
}

adjust_sd_thresh <- function(n_models, default_sd_thresh=0.05){
  if(n_models > 75)
    return(0.06)
  else if (n_models > 100)
    return(0.08)
  else if (n_models > 200)
    return(0.1)
  return(default_sd_thresh)
}

mapBoreal<-function(ice2_30_atl08_path,
                    ice2_30_sample_path,
                    offset=100,
                    rep=10,
                    stack=stack,
                    minDOY=121,
                    maxDOY=273,
                    max_sol_el=0,
                    expand_training=TRUE,
                    calculate_uncertainty=TRUE,
                    local_train_perc=100,
                    min_icesat2_samples=5000,
                    boreal_poly=boreal_poly,
                    predict_var,
                    max_n=10000,
                    pred_vars=c('elev', 'slope')){

  tile_num = tail(unlist(strsplit(path_ext_remove(ice2_30_atl08_path), "_")), n=1)
  cat("Modelling and mapping boreal AGB tile: ", tile_num, "\n")

  all_train_data <- prepare_training_data(
    ice2_30_atl08_path, ice2_30_sample_path, minDOY, maxDOY, max_sol_el,
    expand_training, min_icesat2_samples, local_train_perc, offset
  )

  fixed_modeling_pipeline_params <- list(
    rds_models=get_rds_models(), all_train_data=all_train_data, boreal_poly=boreal_poly,
    pred_vars=pred_vars, predict_var=predict_var, stack=stack,
    summary_and_convert_functions=get_summary_and_convert_functions(predict_var),
    model=randomForest, sample=TRUE
  )

  results <- do.call(run_modeling_pipeline, modifyList(
    fixed_modeling_pipeline_params,
    list(max_n=max_n, randomize=FALSE, model_config=list(ntree=NTREE, mtry=6))
  ))

  output_fns <- set_output_file_names(predict_var, tile_num)

  write_ATL08_table(predict_var, results[['train_df']], output_fns[['train']])
  write_single_model_summary(results[['model']], results[['train_df']],  predict_var, output_fns)

  maps <- c(results[['map']])
  tile_summaries <- c(results[['tile_summary']])
  boreal_summaries <- c(results[['boreal_summary']])

  print('First Prediction Results:')
  cat('tile_summary:', results[['tile_summary']])
  cat(' boreal_summary:', results[['tile_summary']], '\n')

  if (calculate_uncertainty) {
    # loop variables
    sd_thresh <- 0.05
    last_n <- 9
    # can be initialized to anything bigger than sd_thresh
    sd_diff <- sd_thresh + 1
    this_rep <- 1
    # Here sample size is reduced from 10K to 1K, and by setting randmize=TRUE
    # the linear AGB models are randomized, also, the default random forest mtry is used
    params <- modifyList(
      fixed_modeling_pipeline_params,
      list(max_n=1000, randomize=TRUE, model_config=list(ntree=NTREE))
    )

    while(sd_diff > sd_thresh && this_rep < rep){
      cat('Uncertainty loop, iteration:', this_rep, '\n')
      results <- do.call(run_modeling_pipeline, params)
      cat('tile_summary:', results[['tile_summary']])
      cat(' boreal_summary:', results[['tile_summary']], '\n')
      maps <- c(maps, results[['map']])
      tile_summaries <- c(tile_summaries, results[['tile_summary']])
      boreal_summaries <- c(boreal_summaries, results[['boreal_summary']])
      if (this_rep > last_n){
        sd_diff <- sd_change_relative_to_baseline(tile_summaries, last_n=last_n)
        cat('sd_diff:', sd_diff, '\n')
      }
      sd_thresh <- adjust_sd_thresh(this_rep)
      this_rep <- this_rep + 1
    }

    # report_convergence(sd_diff, sd_thresh, this_rep, max_iters)
  }
  print('AGB successfully predicted!')
  write_output_summaries(tile_summaries, boreal_summaries, predict_var,  output_fns[['summary']])
  write_output_raster_map(maps, output_fns[['map']])

}

# ####################### Run code ##############################

# Get command line args
# args = commandArgs(trailingOnly=TRUE)
# #rds_filelist <- args[1]
# data_table_file <- args[1]
# topo_stack_file <- args[2]
# l8_stack_file <- args[3]
# LC_mask_file <- args[4]
# DO_MASK_WITH_STACK_VARS <- args[5]
# data_sample_file <- args[6]
# iters <- args[7]
# ppside <- args[8]
# minDOY <- args[9]
# maxDOY <- args[10]
# max_sol_el <- args[11]
# expand_training <- args[12]
# local_train_perc <- args[13]
# min_n <- args[14]
# boreal_vect <- args[15]
# predict_var <- args[16]
# max_n <- args[17]
# pred_vars <- args[18]
# print(pred_vars)
# print('max_n:')
# print(max_n)
pred_vars = '~/Downloads/dps_output/pred_vars.txt'
pred_vars <- as.character(read.table(pred_vars, header=FALSE, sep=' ')[1,])
setwd('~/')
# print('pred_vars:')
# print(pred_vars)

#for debugging replace args with hard paths
#data_table_file <- '/projects/my-private-bucket/dps_output/run_tile_atl08_ubuntu/tile_atl08/2022/11/30/19/22/04/120959/atl08_005_30m_filt_topo_landsat_20221130_1216.csv'
#topo_stack_file <- '/projects/shared-buckets/nathanmthomas/alg_34_testing/Copernicus_1216_covars_cog_topo_stack.tif'
#l8_stack_file <- '/projects/shared-buckets/nathanmthomas/alg_34_testing/HLS_1216_06-15_09-01_2019_2021.tif'
#LC_mask_file <- '/projects/shared-buckets/nathanmthomas/alg_34_testing/esa_worldcover_v100_2020_1216_cog.tif'

## data_table_file <- '~/Downloads/dps_output/atl08_006_030m_2020_2020_06_09_filt_covars_merge_neighbors_034673.csv'
## topo_stack_file <- '~/Downloads/dps_output/CopernicusGLO30_34673_cog_topo_stack.tif'
## l8_stack_file <- '~/Downloads/dps_output/HLS_34673_07-01_08-31_2023_2023.tif'
## LC_mask_file <- '~/Downloads/dps_output/esa_worldcover_v100_2020_34673_cog.tif'
## data_sample_file <- '~/Downloads/dps_output/boreal_train_data_2020_n3.csv'
## boreal_vect <- '~/Downloads/dps_output/wwf_circumboreal_Dissolve.geojson'

## data_table_file <- '~/Downloads/inputs/atl08_006_030m_2020_2020_06_09_filt_covars_merge_neighbors_363400.csv'
## topo_stack_file <- '~/Downloads/inputs/CopernicusGLO30_363400_cog_topo_stack.tif'
## LC_mask_file <- '~/Downloads/inputs/esa_worldcover_v100_2020_363400_cog.tif'
## l8_stack_file <- '~/Downloads/inputs/HLS_363400_07-01_08-31_2020_2020.tif'
## SAR_stack_file <- '~/Downloads/inputs/SAR_S1_2020_363400_cog.tif'

data_sample_file <- '~/Downloads/dps_outputs/boreal_train_data_2020_n3.csv'
boreal_vect <- '~/Downloads/dps_output/wwf_circumboreal_Dissolve.geojson'

data_table_file <- '~/Downloads/inputs/atl08_006_030m_2020_2020_06_09_filt_covars_merge_neighbors_001613.csv'
topo_stack_file <- '~/Downloads/inputs/CopernicusGLO30_1613_cog_topo_stack.tif'
LC_mask_file <- '~/Downloads/inputs/esa_worldcover_v100_2020_1613_cog.tif'
l8_stack_file <- '~/Downloads/inputs/HLS_1613_07-01_08-31_2020_2020.tif'
SAR_stack_file <- '~/Downloads/inputs/SAR_S1_2020_1613_cog.tif'

## data_table_file <- '~/Downloads/dps_outputs/atl08_006_030m_2020_2020_06_09_filt_covars_merge_neighbors_004104.csv'
## topo_stack_file <- '~/Downloads/dps_outputs/CopernicusGLO30_4104_cog_topo_stack.tif'
## LC_mask_file <- '~/Downloads/dps_outputs/esa_worldcover_v100_2020_4104_cog.tif'
## l8_stack_file <- '~/Downloads/dps_outputs/HLS_4104_07-01_08-31_2020_2020.tif'
## data_sample_file <- '~/Downloads/dps_outputs/boreal_train_data_2020_n3.csv'
## boreal_vect <- '~/Downloads/dps_outputs/wwf_circumboreal_Dissolve.geojson'

mask_stack <- 'TRUE'
iters <- 30
ppside <- 1
minDOY <- 130
maxDOY <- 250
max_sol_el <- 5
expand_training <- 'TRUE'
local_train_perc <- 100
min_icesat2_samples <- 5000
max_n <- 10000
calculate_uncertainty <- TRUE
predict_var <- 'AGB'
#predict_var <- 'Ht'

ppside <- as.double(ppside)
minDOY <- as.double(minDOY)
maxDOY <- as.double(maxDOY)
max_sol_el <- as.double(max_sol_el)
local_train_perc <- as.double(local_train_perc)
print(paste0("Do mask? ", mask_stack))

library(randomForest)
library(dplyr)
# library(rockchalk)
library(terra)


resample_or_reproject_inputs <- function(){
  topo <- rast(topo_stack_file)
  l8 <- rast(l8_stack_file)
  lc <- rast(LC_mask_file)
  #sar <- rast(SAR_stack_file)

  if (nrow(topo) != nrow(l8) || ncol(topo) != ncol(l8)){
    topo <- resample(topo, l8, method = 'near')
    ext(topo) <- ext(l8)
  }
  if (nrow(lc) != nrow(l8) || ncol(lc) != ncol(l8)){
    lc <- resample(lc, l8, method = 'near')
    ext(lc) <- ext(l8)
  }
  ## if (nrow(sar) != nrow(l8) || ncol(sar) != ncol(l8)){
  ##   sar <- resample(sar, l8, method = 'near')
  ##   ext(sar) <- ext(l8)
  ## }

  stack <- c(l8, topo, lc)
  boreal_poly <- project(vect(boreal_vect), crs(l8))

  return(list("stack" = stack, "boreal_poly" = boreal_poly))
}

mask_input_stack <- function(stack){
  MASK_LYR_NAMES = c('slopemask', 'ValidMask')
  MASK_LANDCOVER_NAMES = c(50, 70, 80, 100)

  print("Masking stack...")
  # Bricking the stack will make the masking faster (i think)
  # brick = rast(stack)

  for(LYR_NAME in MASK_LYR_NAMES){
    m <- terra::subset(stack, grep(LYR_NAME, names(stack), value = T))
    stack <- mask(stack, m == 0, maskvalue=TRUE)
  }

  for(LC_NAME in MASK_LANDCOVER_NAMES){
    n <- terra::subset(stack, grep('esa_worldcover_v100_2020', names(stack), value=LC_NAME))
    stack <- mask(stack, n == LC_NAME, maskvalue=TRUE)
  }

  return(stack)
}

stack_poly <- resample_or_reproject_inputs()
stack <- if (mask_stack) mask_input_stack(stack_poly[['stack']]) else stack_poly[['stack']]
boreal_poly <- stack_poly[['boreal_poly']]


print("modelling begins")

set.seed(123)
NTREE = 30
maps<-mapBoreal(ice2_30_atl08_path=data_table_file,
                ice2_30_sample=data_sample_file,
                offset=100.0,
                rep=iters,
                stack=stack,
                minDOY=minDOY,
                maxDOY=maxDOY,
                max_sol_el=max_sol_el,
                expand_training=expand_training,
                calculate_uncertainty=calculate_uncertainty,
                local_train_perc=local_train_perc,
                min_icesat2_samples=min_icesat2_samples,
                boreal_poly=boreal_poly,
                predict_var=predict_var,
                max_n=max_n,
                pred_vars=pred_vars)
