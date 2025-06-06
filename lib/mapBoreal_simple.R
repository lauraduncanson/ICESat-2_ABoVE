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

library(optparse)
library(randomForest)
library(dplyr)
library(fs)
library(stringr)
library(rockchalk)
library(terra)
library(parallel)

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
      mutate(across(all_of(RH_columns), ~. + offset))
  )
}

set_model_id_for_AGB_prediction <- function(in_data, offset){
  # TODO: uncomment correct model ids once tested against old results
  return(
    in_data |>
      mutate(model_id = case_when(
        segment_landcover %in% c(111, 113, 121, 123) ~ "m3", # needle leaf
        segment_landcover %in% c(112, 114, 122, 124) ~ "m1", # broad leaf
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
  filter <- which((tile_data$doy >= start_DOY) &
                    (tile_data$doy <= end_DOY) &
                    (tile_data$solar_elevation < solar_elevation)
                  )
  return(filter)
}

late_season_filter <- function(tile_data, minDOY, maxDOY,
                               default_maxDOY, min_icesat2_samples, max_sol_el){
  n_late <- 0
  for(late_months in 0:3) {
    if(n_late < min_icesat2_samples) {

      default_maxDOY <- default_maxDOY + 30 * late_months

      if(default_maxDOY < maxDOY){
        filter <- DOY_and_solar_filter(minDOY, default_maxDOY, max_sol_el)
        n_late <- length(filter)
      }
    }
  }
  return(list(filter=filter, default_maxDOY=default_maxDOY))
}

early_and_late_season_filter <- function(tile_data, minDOY,
                                         default_minDOY, default_maxDOY,
                                         min_icesat2_samples, max_sol_el){
  n_early <- 0
  for(early_months in 0:3){
    if(n_early < min_icesat2_samples){
      default_minDOY <- default_minDOY - 30 * early_months

      if(default_minDOY > minDOY){
        filter <- DOY_and_solar_filter(default_minDOY, default_maxDOY, max_sol_el)
        n_early <- length(filter)
      }
    }

  }
  return(list(filter=filter, default_minDOY=default_minDOY))
}

expand_training_around_season <- function(tile_data, minDOY, maxDOY,
                                          default_minDOY, default_maxDOY,
                                          max_sol_el, min_icesat2_samples){
  filter <- DOY_and_solar_filter(tile_data, minDOY, maxDOY, max_sol_el)
  if(length(filter) >= min_icesat2_samples){
    cat('Found nough data with max_solar_elevation:', max_sol_el, '\n')
    return(tile_data[filter,])
  }
  # next try expanding 1 month later in growing season, iteratively, up to 3 months
  filter <- late_season_filter(
    tile_data, minDOY, maxDOY, default_maxDOY, min_icesat2_samples, max_sol_el
  )
  if(length(filter$filter) >= min_icesat2_samples){
    cat('Found enough data when expanding into late season DOY:', filter$default_maxDOY, '\n')
    return(tile_data[filter$filter,])
  }

  # next try expanding 1 month earlier in growing season, iteratively, up to 3 months
  # Note that the upper window might be later in the growing season from the previous call
  current_default_maxDOY <- filter$default_maxDOY
  filter <- early_and_late_season_filter(
    tile_data, minDOY, default_minDOY, current_default_maxDOY, min_icesat2_samples, max_sol_el
  )
  if(length(filter$filter) >= min_icesat2_samples){
    cat(
      'Found enough data when expanding into early and late season DOY:[',
      filter$default_minDOY, ' ', default_minDOY,  ']\n'
    )
    return(tile_data[filter$filter,])
  }

  print('Search into late and early season did not return enough data')
  print('applying basic filter')
  tile_data <- tile_data[DOY_and_solar_filter(tile_data, default_minDOY, default_maxDOY, 0),]
  return(tile_data)
}

reduce_sample_size <- function(df, sample_size){
  return(df[sample(row.names(df), sample_size, replace=FALSE),])
}

remove_stale_columns <- function(df, column_names) {
  columns_to_remove <- intersect(names(df), column_names)
  df <- df[, !names(df) %in% columns_to_remove, drop=FALSE]

  return(df)
}

sample_broad_data_within_latitude <- function(broad_data, lat, threshold, samples_needed){
  broad_within_lat <- which(broad_data$lat > (lat-threshold) & broad_data$lat < (lat+threshold))
  broad_data <- broad_data[broad_within_lat,]
  return(broad_data[sample(row.names(broad_data), samples_needed, replace=TRUE), ])
}

expand_training_with_broad_data <- function(broad_data, tile_data, samples_needed){
  broad_data <- sample_broad_data_within_latitude(broad_data, min(tile_data$lat), 5, samples_needed)
  # TODO this check may no longer be needed, ask Paul
  if (!setequal(names(broad_data), names(tile_data))) {
    only_in_broad <- setdiff(names(broad_data), names(tile_data))
    only_in_local <- setdiff(names(tile_data), names(broad_data))
    diff <- union(only_in_broad, only_in_local)
    print('Warning!')
    print('Boreal wide training data and local data have non matching columns!')
    print('Will drop the following non matching columns and continue:')
    print(diff)
    common <- intersect(names(broad_data), names(tile_data))
    tile_data <- tile_data[common]
    broad_data <- broad_data[common]
    print(common)
  }

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

set_short_veg_height_to_zero <- function(df, slope_thresh){
  height_columns <- c(names(df)[grep('^rh[0-9]{2}$', names(df))],
                      "h_canopy","h_min_canopy", "h_max_canopy", "h_mean_canopy")
  cond <- df$segment_landcover == 100
  cond <- cond | ((df$segment_landcover %in% c(20, 30, 60, 100)) & (df$slope > slope_thresh))
  cond[is.na(cond)] <- FALSE
  df[cond, height_columns] <- 0.0
  return(df)
}

short_veg_filter <- function(df, slope_thresh){
  cond <- (df$segment_landcover %in% c(20, 30, 60, 100)) & (df$slope > slope_thresh)
  cond[is.na(cond)] <- FALSE
  df <- df[!cond, ,drop=FALSE]
  return(df)
}

set_output_file_names <- function(predict_var, tile_num, year){
  key <- if (predict_var == 'AGB') 'agb' else 'ht'
  out_fn_stem = paste(
    paste0('output/boreal_', key, '_', year), format(Sys.time(),"%Y%m%d%s"), str_pad(tile_num, 7, pad = "0"),
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
  local_model <- lm(model$predicted ~ target)
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
  out_accuracy <- list(RSQ=rsq_local, RMSE=rmse_local, importance=imp_vars)
  saveRDS(out_accuracy, file=out_fns['stats'])
}

write_output_raster_map <- function(maps, std = NULL, output_fn) {
  # Set NA flag for primary maps
  if (is.list(maps) || nlyr(maps) > 0) {
    for (i in 1:nlyr(maps)) {
      NAflag(maps[[i]]) <- -9999
    }
  } else {
    NAflag(maps) <- -9999
  }

  if (!is.null(std)) {
    # Uncertainty case: combine mean and std
    NAflag(std) <- -9999
    output_maps <- c(maps, std)
  } else if (nlyr(maps) > 1) {
    # Multiple layers case: calculate mean and sd
    output_maps <- c(app(maps, mean), app(maps, sd))
  } else {
    # Single layer case
    output_maps <- maps
  }

  raster_options <- c("COMPRESS=LZW", overwrite = TRUE,
                      gdal = c("COMPRESS=LZW", "OVERVIEW_RESAMPLING=AVERAGE"))
  writeRaster(output_maps, filename = output_fn, filetype = "COG",
              gdal = raster_options, NAflag = -9999)
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

read_and_filter_training_data <- function(ice2_30_atl08_path, expand_training, min_icesat2_samples, minDOY, maxDOY, max_sol_el){
  default_maxDOY <- 273
  default_minDOY <- 121
  tile_data <- read.csv(ice2_30_atl08_path)

  night_time_in_season <- DOY_and_solar_filter(tile_data, default_minDOY, default_maxDOY, 0)
  cat('train data size before any filtering:', nrow(tile_data), '\n')
  cat('length(night_time_in_season):', length(night_time_in_season), 'expand_training:', expand_training, ' min_n:', min_icesat2_samples, '\n')

  if (length(night_time_in_season) < min_icesat2_samples && expand_training){
    cat('running expansion:', length(night_time_in_season), '<', min_icesat2_samples, '\n')
    tile_data <- expand_training_around_season(tile_data, minDOY, maxDOY, default_minDOY, default_maxDOY, max_sol_el, min_icesat2_samples)
  }
  else{
    print('night time filter only')
    tile_data <- tile_data[night_time_in_season, ]
  }
  cat('training data size after filtering:', nrow(tile_data), '\n')
  tile_data <- remove_stale_columns(tile_data, c("binsize", "num_bins"))
  return(tile_data)
}

augment_training_data_with_broad_data <- function(tile_data, ice2_30_sample_path, local_train_perc, min_icesat2_samples){
  broad_data <- read.csv(ice2_30_sample_path)
  broad_data <- remove_stale_columns(broad_data, c("X__index_level_0__", "geometry"))

  # take proportion of broad data we want based on local_train_perc
  sample_local <- ceiling(nrow(tile_data) * local_train_perc / 100)
  cat('sample_local:', sample_local, '\n')

  if (sample_local < min_icesat2_samples){
    cat('reducing sample size to', sample_local, ' from ', nrow(tile_data), 'to complete with broad data \n')
    tile_data <- reduce_sample_size(tile_data, sample_local)
  }

  # sample from broad data to complete sample size
  # this will work if either there aren't enough local samples for n_min OR if there is forced broad sampling
  n_broad <- min_icesat2_samples - nrow(tile_data)
  if(n_broad > 1){
    tile_data <- expand_training_with_broad_data(broad_data, tile_data, n_broad)
    cat('training data size after augmenting with broad data:', nrow(tile_data), '\n')
  }
  return(tile_data)
}

reformat_training_data_for_AGB_modeling <- function(tile_data, offset){
  tile_data <- rename_height_columns_to_match_pretrained_models(tile_data)
  tile_data$h_canopy <- tile_data$RH_98
  tile_data <- offset_RH_columns(tile_data, offset)
  tile_data <- set_model_id_for_AGB_prediction(tile_data)
  return(tile_data)
}

prepare_training_data <- function(ice2_30_atl08_path, ice2_30_sample_path,
                                  expand_training, minDOY, maxDOY, max_sol_el,
                                  min_icesat2_samples, local_train_perc, offset, stack_vars,
                                  remove_short_veg, zero_short_veg_height, slope_thresh){
  print('preparing training data ...')

  tile_data <- read_and_filter_training_data(
    ice2_30_atl08_path, expand_training,
    min_icesat2_samples, minDOY, maxDOY, max_sol_el
  )

  tile_data <- augment_training_data_with_broad_data(
    tile_data, ice2_30_sample_path, local_train_perc, min_icesat2_samples
  )

  if (remove_short_veg)
    tile_data <- short_veg_filter(tile_data, slope_thresh)

  if (zero_short_veg_height)
    tile_data <- set_short_veg_height_to_zero(tile_data, slope_thresh)

  needed_cols <- union(
    c('lat', 'lon', 'segment_landcover', 'h_canopy', 'rh25', 'rh50', 'rh60',
      'rh70', 'rh75', 'rh80', 'rh85', 'rh90', 'rh95'),
    stack_vars
  )

  tile_data <- tile_data |> select(all_of(needed_cols))

  tile_data <- reformat_training_data_for_AGB_modeling(tile_data, offset)

  tile_data <- remove_height_outliers(tile_data)
  cat('training data size after removing height outliers:', nrow(tile_data), '\n')

  tile_data <- tile_data |> filter(if_all(everything(), ~ !is.na(.x) & .x != -9999))
  cat('training data size after removing NAs:', nrow(tile_data), '\n')

  str(tile_data)
  cat('table for model training generated with ', nrow(tile_data), ' observations\n')

  if (nrow(tile_data) <= 1) {
    # TODO another option could be to drop SAR columns and continue
    stop('No traing data available, likley due to SAR being all -9999')
  }
  return(tile_data)
}

get_rds_models <- function(){
  rds_model_fns <- list.files(pattern='*.rds')
  rds_models <- lapply(rds_model_fns, readRDS)
  names(rds_models) <- paste0("m",1:length(rds_models))
  print(rds_models)
  return(rds_models)
}

create_predict_function <- function(cores){
  if (cores == 1) {
    predict_stack <- function(model, stack) {
      stack <- na.omit(stack)
      map <- predict(stack, model, na.rm = TRUE)
      map <- mask(map, stack$slopemask, maskvalues = 0, updatevalue = 0)
      map <- mask(map, stack$ValidMask, maskvalues = 0, updatevalue = 0)
      return(map)
    }
  }
  else {
    predict_stack_manual_chunks <- function(model, stack_path) {
      n_chunks <- cores
      stack <- rast(stack_path)
      chunk_size <- ceiling(nrow(stack) / n_chunks)
      chunks <- vector("list", n_chunks)

      for (i in seq_len(n_chunks)) {
        row_start <- (i - 1) * chunk_size + 1
        row_end   <- min(i * chunk_size, nrow(stack))
        chunks[[i]] <- list(id = i, row_start = row_start, row_end = row_end)
      }

      process_chunk <- function(ch) {
        stack <- rast(stack_path)
        n_rows <- ch$row_end - ch$row_start + 1
        vals <- terra::values(stack, row = ch$row_start, nrows = n_rows)
        valid <- complete.cases(vals)
        preds <- rep(NA_real_, nrow(vals))
        if (any(valid)){
          preds[valid] <- predict(model, vals[valid, , drop = FALSE])
        }
        list(preds = preds, row_start = ch$row_start, row_end = ch$row_end)
      }

      message("Running on cluster with ", n_chunks, " workers")
      results <- parallel::mclapply(chunks, process_chunk,
                                    mc.cores = n_chunks, mc.preschedule = FALSE)
      # Reassemble

      all_preds <- numeric(nrow(stack) * ncol(stack))
      for (res in results) {
        offset <- (res$row_start - 1) * ncol(stack) + 1
        all_preds[offset:(offset + length(res$preds) - 1)] <- res$preds
      }

      rast_pred <- rast(matrix(all_preds, nrow = nrow(stack), ncol = ncol(stack), byrow = TRUE))
      ext(rast_pred) <- ext(stack)
      crs(rast_pred) <- crs(stack)
      rast_pred
    }
  }
}

tile_and_boreal_summary <- function(map, predict_var, boreal_poly, summary_and_convert_functions){
  convert_fun <- summary_and_convert_functions[['convert_fun']]
  summary_fun <- summary_and_convert_functions[['summary_fun']]

  tile_summary <- convert_fun(global(map, summary_fun, na.rm=TRUE)[[summary_fun]])
  boreal_extract <- extract(map, boreal_poly, fun=summary_fun, na.rm=TRUE, touches=TRUE)
  boreal_summary <- convert_fun(sum(boreal_extract$lyr.1, na.rm=TRUE))

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
                                 summary_and_convert_functions, predict_function,
                                 max_samples, sample, pred_vars, predict_var, stack){
  t1 <- Sys.time()
  print('creating AGB traing data frame.')
  train_df <- GEDI2AT08AGB(rds_models, all_train_data, randomize, max_samples, sample)

  print('fitting model')
  model <- fit_model(model, model_config, train_df, pred_vars, predict_var)

  print('predicting biomass map')
  map <- predict_function(model, stack)

  print('calculating tile and boreal summaries')
  summary <- tile_and_boreal_summary(map, predict_var, boreal_poly, summary_and_convert_functions)
  cat('tile_summary:', summary$tile_summary, ' boreal summary:', summary$boreal_summary, '\n')

  t2 <- Sys.time()
  cat('pipeline runtime:', difftime(t2, t1, units="mins"), ' (m)\n')

  return(list(
    train_df=train_df, model=model, map=map,
    tile_summary=summary[['tile_summary']], boreal_summary=summary[['boreal_summary']]
  ))
}

get_summary_and_convert_functions <- function(predict_var){
  if (predict_var == 'AGB'){
    summary_fun <- 'sum'
    convert_fun <- function(x){x * 0.09 * 1e-9}
  }
  else{
    summary_fun <- 'mean'
    # no conversion needed for Ht
    convert_fun <- function(x){x}
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

welford_update <- function(count, mu, M2, new_value){
    count <- count + 1
    delta = new_value - mu
    mu <- mu + delta / count
    delta2 = new_value - mu
    M2 <- M2 + delta * delta2
    return (list(count=count, mu=mu, M2=M2))
}

run_uncertainty_calculation <- function(fixed_modeling_pipeline_params, max_iters, min_iters, results){
  sd_thresh <- 0.05
  last_n <- 9 # kind of arbitrary
  # sd_diff can be initialized to anything bigger than sd_thresh
  sd_diff <- sd_thresh + 1
  this_iter <- 1

  mu <- c(results[['map']])
  tile_summary <- c(results[['tile_summary']])
  boreal_summary <- c(results[['boreal_summary']])

  params <- modifyList(
    fixed_modeling_pipeline_params,
    list(max_samples=1000, randomize=TRUE)
  )

  # initializing to 0, with crs of mu
  M2 <- mu
  values(M2) <- 0.0

  while(((sd_diff > sd_thresh) && (this_iter < max_iters)) || (this_iter < min_iters)){
    cat('Uncertainty loop, iteration:', this_iter, '\n')
    results <- do.call(run_modeling_pipeline, params)
    res <- welford_update(this_iter, mu, M2, results[['map']])
    tile_summary <- c(tile_summary, results[['tile_summary']])
    boreal_summary <- c(boreal_summary, results[['boreal_summary']])

    if (this_iter > last_n){
      sd_diff <- sd_change_relative_to_baseline(tile_summary, last_n=last_n)
      cat('sd_diff:', sd_diff, '\n')
    }

    sd_thresh <- adjust_sd_thresh(this_iter)
    this_iter <- this_iter + 1
    mu <- res[['mu']]
    M2 <- res[['M2']]
  }

  return(list(map=mu, std=sqrt(M2/(this_iter - 1)), n_iters=this_iter,
              tile_summary=tile_summary, boreal_summary=boreal_summary))
}

resample_if_needed <- function(src, des){
  if (nrow(src) != nrow(des) || ncol(src) != ncol(des)){
    src <- resample(src, des, method = 'near')
    ext(src) <- ext(des)
  }
  return(src)
}

prepare_raster <- function(path, subset_bands=NULL, extra_bands=NULL, dest_raster=NULL){
  raster <- rast(path)
  raster_bands <- names(raster)

  if (!is.null(subset_bands))
    raster_bands <- intersect(raster_bands, subset_bands)

  if (!is.null(extra_bands))
    raster_bands <- c(raster_bands, extra_bands)

  raster <- subset(raster, raster_bands)

  if (!is.null(dest_raster))
    raster <- resample_if_needed(raster, dest_raster)

  return(raster)
}

resample_reproject_and_mask <- function(topo_path, hls_path, lc_path, pred_vars, mask, sar_path=NULL){
  hls <- prepare_raster(hls_path, subset_bands=pred_vars, extra_bands='ValidMask')
  topo <- prepare_raster(topo_path, subset_bands=pred_vars, extra_bands='slopemask', dest_raster=hls)
  lc <- prepare_raster(lc_path, dest_raster=hls)

  if (!is.null(sar_path)) {
    sar_path <- sub("^s3://", "/vsis3/", sar_path)
    sar <- prepare_raster(sar_path, subset_bands=pred_vars, dest_raster=hls)
    stack <- c(hls, sar, topo, lc)
  }
  else {
    stack <- c(hls, topo, lc)
  }

  if(mask)
    stack <- mask_input_stack(stack)

  # landcover mask is not needed anymore
  # TODO I think masks are not probably needed once we leave this function
  stack <- subset(stack, names(stack) != 'esa_worldcover_v100_2020')
  return(stack)
}

mask_input_stack <- function(stack){
  MASK_LYR_NAMES = c('slopemask', 'ValidMask')
  MASK_LANDCOVER_NAMES = c(50, 60, 70, 80)

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

parse_pred_vars <- function(pred_vars, remove_sar){
  pred_vars <- unlist(strsplit(pred_vars, split = " "))

  if(remove_sar){
    print('Removing default SAR variables from pred_vars')
    sar_vars <- list(
      "vv_median_frozen", "vh_median_frozen",
      "vv_median_summer", "vh_median_summer",
      "vv_median_shoulder", "vh_median_shoulder"
    )
    pred_vars <- pred_vars[!pred_vars %in% sar_vars]
  }

  return(pred_vars)
}

mapBoreal<-function(atl08_path, broad_path, hls_path, topo_path, lc_path, boreal_vector_path, year,
                    sar_path=NULL, mask=TRUE, max_sol_el=5, offset=100, minDOY=130, maxDOY=250,
                    expand_training=TRUE, calculate_uncertainty=TRUE, max_iters=30, min_iters=0,
                    local_train_perc=100, min_samples=5000, max_samples=10000, cores=1, ntree=100,
                    predict_var='AGB', pred_vars=c('elevation', 'slope', 'NDVI'),
                    remove_short_veg=FALSE, zero_short_veg_height=FALSE, slope_thresh=15){

  tile_num = tail(unlist(strsplit(path_ext_remove(atl08_path), "_")), n=1)
  cat("Modelling and mapping boreal AGB tile: ", tile_num, "\n")

  pred_vars <- parse_pred_vars(pred_vars, remove_sar=is.null(sar_path))
  print(pred_vars)

  stack <- resample_reproject_and_mask(topo_path, hls_path, lc_path, pred_vars, mask, sar_path)
  boreal_poly <- project(vect(boreal_vector_path), crs(stack))

  all_train_data <- prepare_training_data(
    atl08_path, broad_path, expand_training, minDOY,
    maxDOY, max_sol_el, min_samples, local_train_perc, offset, names(stack),
    remove_short_veg, zero_short_veg_height, slope_thresh
  )

  if (cores > 1) {
    # saving stack and passing path instead because otherwise when parallel section runs
    # the memory overhead multiplies by number of nodes
    stack_path <- './stack.tif'
    writeRaster(stack, stack_path, overwrite = TRUE)
    rm(stack); gc()
    stack <- stack_path
  }

  fixed_modeling_pipeline_params <- list(
    rds_models=get_rds_models(), all_train_data=all_train_data, boreal_poly=boreal_poly,
    pred_vars=pred_vars, predict_var=predict_var, stack=stack,
    summary_and_convert_functions=get_summary_and_convert_functions(predict_var),
    model=randomForest, model_config=list(ntree=ntree), sample=TRUE,
    predict_function=create_predict_function(cores=cores)
  )

  results <- do.call(run_modeling_pipeline, modifyList(
    fixed_modeling_pipeline_params,
    list(max_samples=max_samples, randomize=FALSE)
  ))

  output_fns <- set_output_file_names(predict_var, tile_num, year)

  write_ATL08_table(predict_var, results[['train_df']], output_fns[['train']])
  write_single_model_summary(results[['model']], results[['train_df']],  predict_var, output_fns)

  if (calculate_uncertainty) {
    results <- run_uncertainty_calculation(fixed_modeling_pipeline_params, max_iters, min_iters, results)
  }
  print('AGB successfully predicted!')
  write_output_summaries(results[['tile_summary']], results[['boreal_summary']], predict_var,  output_fns[['summary']])

  if (calculate_uncertainty) {
    write_output_raster_map(results[['map']], results[['std']], output_fns[['map']])
  } else {
    write_output_raster_map(results[['map']], output_fn = output_fns[['map']])
  }

}

option_list <- list(
  make_option(
    c("-a", "--atl08_path"), type = "character",
    help = "Path to the atl08 training data"
  ),
  make_option(
    c("-b", "--broad_path"), type = "character",
    help = "Path to the boreal wide training data"
  ),
  make_option(
    c("-t", "--topo_path"), type = "character",
    help = "Path to the topo stack file"
  ),
  make_option(
    c("-h", "--hls_path"), type = "character",
    help = "Path to the HLS stack file"
  ),
  make_option(
    c("-l", "--lc_path"), type = "character",
    help = "Path to the land cover mask file"
  ),
  make_option(
    c("-s", "--sar_path"), type = "character", default = NULL,
    help = "Path to the land cover mask file"
  ),
  make_option(
    c("-v", "--boreal_vector_path"), type = "character",
    help = "Path to the boreal vector file",
    ),
  make_option(
    c("-y", "--year"), type = "character",
    help = "Year of the input HLS imagery"
  ),
  make_option(
    c("-m", "--mask"), type = "logical", default = TRUE,
    help = "Whether to mask imagery [default: %default]"
  ),
  make_option(
    c("--max_sol_el"), type = "numeric", default = 5,
    help = "Maximum solar elevation degree allowed in training data [default: %default]"
  ),
  make_option(
    c("--minDOY"), type = "integer", default = 130,
    help = "Minimum day of year allowed in training data [default: %default]"
  ),
  make_option(
    c("--maxDOY"), type = "integer", default = 250,
    help = "Maximum day of year allowed in training data [default: %default]"
  ),
  make_option(
    c("--min_samples"), type = "integer", default = 5000,
    help = "Minimum number of samples to avoid augmenting the training data with broad data [default: %default]"
  ),
  make_option(
    c("--max_samples"), type = "integer", default = 10000,
    help = "Maximum number of samples used for training [default: %default]"
  ),
  make_option(
    c("-e", "--expand_training"), type = "logical", default = TRUE,
    help = "Whether to expand training around the season [default: %default]"
  ),
  make_option(
    c("-u", "--calculate_uncertainty"), type = "logical", default = TRUE,
    help = "Whether to calculate uncertainty [default: %default]"
  ),
  make_option(
    c("--max_iters"), type = "integer", default = 30,
    help = "Max number of uncertainty iterations [default: %default]"
  ),
  make_option(
    c("--min_iters"), type = "integer", default = 0,
    help = "Min number of uncertainty iterations [default: %default]"
  ),
  make_option(
    c("-c", "--cores"), type = "integer", default = 1,
    help = "Number of cores used for parallel prediction steps [default: %default]"
  ),
 make_option(
    c ("--ntree"), type = "integer", default = 100,
    help = "Number of random forest trees [default: %default]"
  ),
 make_option(
    c ("--remove_short_veg"), type = "logical", default = FALSE,
    help = "removes shrubs, herbaceous, moss/lichen and bare/spare veg classes from training dataset [default: %default]"
  ),
 make_option(
    c ("--zero_short_veg_height"), type = "logical", default = FALSE,
    help = "Sets the RH metrics of shrubs, herbaceous, moss/lichen and bare/spare veg classes from training dataset to zero [default: %default]"
  ),
 make_option(
    c ("--slope_thresh"), type = "numeric", default = 15,
    help = "slope threshold beyond which short veg height is set to zero [default: %default]"
  ),
 make_option(
    c("-p", "--local_train_perc"), type = "integer", default = 100,
    help = "Percent of atl08 data to be used in case it is augmented with broad data [default: %default]"
  ),
  make_option(
    c("--predict_var"), type = "character", default = "AGB",
    help = "Variable to predict, it can be either AGB or Ht [default: %default]"
  ),
  make_option(
    c("--pred_vars"), type = "character",
    default = paste(
      "Red Green elevation slope tsri tpi NIR SWIR SWIR2 NDVI",
      "SAVI MSAVI NDMI EVI NBR NBR2 TCB TCG TCW",
      "vv_median_frozen vh_median_frozen vv_median_summer",
      "vh_median_summer vv_median_shoulder vh_median_shoulder"
    ),
    help = paste(
      "List of predictor variables, must be a subset from the default options",
      "seperated by space, e.g, NDVI slope\n [default: %default]"
    )
  ),
  make_option(
    c("--help"), action = "store_true",
    help = "Show help message"
  )
)

opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
opt <- parse_args(opt_parser)

cat("Parsed arguments:\n")
print(opt)

if (!is.null(opt$help)) {
  print_help(opt_parser)
}
do.call(mapBoreal, opt)
