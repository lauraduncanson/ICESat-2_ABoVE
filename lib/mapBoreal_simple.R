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
library(arrow)

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

GEDI2AT08AGB<-function(biomass_models, df, iter=1, max_n=10000, sample=TRUE){
  if (sample && nrow(df) > max_n)
    df <- reduce_sample_size(df, max_n)

  df$AGB <- NA
  df$SE <- NA

  ids<-unique(df$model_id)
  n_models <- length(ids)

  for (i in ids){
    model_id_iter <- paste0(i, '_', iter)
    model_i <- biomass_models[[model_id_iter]]
    cat('\niter:', iter, 'model_id:', model_id_iter, '\n')
    # Predict AGB and SE
    df$AGB[df$model_id==i] <- predict(model_i, newdata=df[df$model_id==i,])
    df$SE[df$model_id==i] <- summary(model_i)$sigma^2

    df$AGB[df$AGB < 0] <- 0.0

    # Calculate Correction Factor C
    C <- mean(model_i$model$`sqrt(AGBD)`^2) / mean(model_i$fitted.values^2)

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
    paste0('output/boreal_', key, '_', year),
    format(Sys.time(),"%Y%m%d%s"),
    str_pad(tile_num, 7, pad = "0"),
    sep="_"
  )

  fn_suffix <- c('.tif', '_overall.csv', '_north.csv', '_boreal.csv', '_boreal_eco.csv',
                 '_by_lc.csv', '_by_slope.csv','_by_ecoregion.csv', '_by_country.csv',
                 '_model_stats.csv', '_ensemble_stats.csv', '_val.csv', '_train.parquet')
  names <- c('map', 'overall', 'north', 'boreal', 'boreal_eco',
             'by_lc', 'by_slope', 'by_ecoregion', 'by_country',
             'model_stats', 'ensemble_stats', 'validation', 'training')

  output_file_names <- paste0(out_fn_stem, fn_suffix)
  names(output_file_names) <- names

  return(output_file_names)
}

get_model_stats <- function(model){
  rsq <- tail(model$rsq, 1)

  mse <- tail(model$mse, 1)
  rmse <- if(is.null(mse) || is.na(mse)) NA else sqrt(mse)

  summary <- as.data.frame(t(model$importance))
  summary$OOB_MSE <- mse
  summary$OOB_RMSE <- rmse
  summary$OOB_R2 <- rsq

  return(summary)
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

  tmp_output_fn <- str_replace(output_fn, '.tif', '_temp.tif')
  writeRaster(
    output_maps,
    filename = tmp_output_fn,
    filetype = "GTiff",
    gdal = c("COMPRESS=DEFLATE", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
    NAflag = -9999
  )
  # adding custom overviews for faster visualization
  system(paste("gdaladdo",
               "-r average",
               "--config GDAL_TIFF_OVR_BLOCKSIZE 256",
               tmp_output_fn,
               "2 4 8 16",
               collapse = " "))
  # wrapping the tmp tiff in COG and forcing the same overviews and blocksize
  system(paste("gdal_translate",
               tmp_output_fn,
               output_fn,
               "-of COG",
               "-co BLOCKSIZE=256",
               "-co OVERVIEWS=FORCE_USE_EXISTING",
               "-co COMPRESS=DEFLATE",
               collapse = " "))
  file.remove(tmp_output_fn)
}

write_output_summaries_and_stats <-function(summaries, model_stats, output_fns){

  for(k in names(summaries)){
    write.csv(summaries[[k]], output_fns[[k]], row.names=FALSE)
  }

  all_model_stats <- bind_rows(model_stats)
  row.names(all_model_stats) <- 1:nrow(all_model_stats)
  write.csv(all_model_stats, output_fns[['model_stats']])
}

read_and_filter_training_data <- function(atl08_path, expand_training, min_samples, minDOY, maxDOY, max_sol_el){
  default_maxDOY <- 273
  default_minDOY <- 121
  tile_data <- read.csv(atl08_path)

  night_time_in_season <- DOY_and_solar_filter(tile_data, default_minDOY, default_maxDOY, 0)
  cat('train data size before any filtering:', nrow(tile_data), '\n')
  cat('length(night_time_in_season):', length(night_time_in_season), 'expand_training:', expand_training, ' min_n:', min_samples, '\n')

  if (length(night_time_in_season) < min_samples && expand_training){
    cat('running expansion:', length(night_time_in_season), '<', min_samples, '\n')
    tile_data <- expand_training_around_season(tile_data, minDOY, maxDOY, default_minDOY, default_maxDOY, max_sol_el, min_samples)
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
                                  zero_short_veg_height, slope_thresh,
                                  year, val_thresh, val_frac, biomass_models, out_train_data_fn){
  print('preparing training data ...')

  tile_data <- read_and_filter_training_data(
    ice2_30_atl08_path, expand_training,
    min_icesat2_samples, minDOY, maxDOY, max_sol_el
  )

  tile_data <- augment_training_data_with_broad_data(
    tile_data, ice2_30_sample_path, local_train_perc, min_icesat2_samples
  )

  if (zero_short_veg_height)
    tile_data <- set_short_veg_height_to_zero(tile_data, slope_thresh)

  needed_cols <- setdiff(union(
    c('y', 'lat', 'lon', 'segment_landcover', 'h_canopy', 'rh25', 'rh50', 'rh60',
      'rh70', 'rh75', 'rh80', 'rh85', 'rh90', 'rh95'),
    stack_vars
  ), c('esa_worldcover_v100_2020'))

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

  # to get AGB and SE and publish the tile_data
  tile_data <- GEDI2AT08AGB(biomass_models, tile_data, iter=1, sample=FALSE)
  # reset the offset before saving
  tile_data <- offset_RH_columns(tile_data, -1 * offset)
  write_parquet(tile_data[c('lon', 'lat', 'segment_landcover',
                            get_height_column_names(tile_data),
                            'AGB', 'SE', 'model_id')], out_train_data_fn)
  # add the offset back and drop AGB and SE
  tile_data <- offset_RH_columns(tile_data, offset)
  tile_data <- tile_data[,!(names(tile_data) %in% c('AGB', 'SE'))]

  tile_data <- train_test_split_if_enough_data(tile_data, year, val_thresh, val_frac)

  return(tile_data)
}

train_test_split_if_enough_data <- function(df, year, val_thresh, val_frac){
  df_y <- df[df$y == year,]
  if (nrow(df_y) > val_thresh && val_frac > 0) {
    val_sz <- floor(val_frac * nrow(df_y))
    val_idx <- sample(row.names(df_y), val_sz, replace=FALSE)
    val_data <- df_y[val_idx,]
    train_data <- df[!row.names(df) %in% val_idx,]
  }
  else {
    val_data <- NULL
    train_data <- df
  }
  return(list(val_data=val_data, train_data=train_data))
}

read_randomized_biomass_models <- function(biomass_models_path, n){
  base_dir <- dirname(biomass_models_path)
  untar(biomass_models_path, exdir=base_dir)

  models <- vector(mode='list', length = n*3)
  model_ids <- c('m1', 'm3', 'm8')
  names <- c()
  for (model_id in model_ids){
    names <- c(names, paste0(model_id, '_', 1:n))
  }
  names(models) <- names

  for (model_name in names){
    models[[model_name]] <- readRDS(file.path(base_dir, paste0(model_name, '.rds')))
  }
  return(models)
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

classify_slope <- function(slope_raster, layer_name){
  m <- c(-Inf, 0, 0,
         0, 5, 1,
         5, 10, 2,
         10, 15, 3,
         15, 20, 4,
         20, 30, 5,
         30, 40, 6,
         40, 50, 7,
         50, 90, 8,
         90, Inf, 9)

  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  slope_class <- classify(slope_raster, rclmat, include.lowest=TRUE, right=FALSE)
  names(slope_class) <- layer_name
  return(slope_class)
}

rasterize_boundaries <- function(template, poly, layer_name, field=NULL){
  # rasterize poly using template raster to set dim, crs, res, ...
  # assumes poly is in EPSG 4326
  bbox_4326 <- terra::buffer(project(as.polygons(ext(template), crs(template)), 'EPSG:4326'), 1000)
  poly_cropped <- project(crop(poly, bbox_4326), crs(template))

  if (nrow(poly_cropped) == 0) {
    # No overlap case return all NA raster
    r <- rast(template)
    values(r) <- NA
  }
  else if (!is.null(field)) {
    r <- rasterize(poly_cropped, template, field=field, touches=FALSE)
  }
  else {
    r <- rasterize(poly_cropped, template, touches=TRUE)
  }
  names(r) <- layer_name
  return(r)
}

clip_to_north_lat <- function(template, north_lat=51.6){
  poly <- as.polygons(ext(template), crs(template))
  poly_4326 <- project(poly, 'EPSG:4326')
  poly_4326_ext <- ext(poly_4326)

  if (ymin(poly_4326_ext) <= north_lat && ymax(poly_4326_ext) >= north_lat){
    relative_to_north_lat <- 'intersects'
    ymin(poly_4326_ext) <- north_lat
    poly_4326_north <- as.polygons(poly_4326_ext, 'EPSG:4326')
  }
  else if (ymax(poly_4326_ext) < north_lat){
    relative_to_north_lat <- 'south'
    poly_4326_north <- NULL
  }
  else {
    relative_to_north_lat <- 'north'
    poly_4326_north <- NULL
  }

  return(list(poly_4326_north=poly_4326_north,
              relative_to_north_lat=relative_to_north_lat))
}

prep_summary_layers <- function(slope_raster, lc_raster, ecoregions, boreal_poly, countries){
  slope_lyr <- classify_slope(slope_raster, 'slope')
  names(lc_raster) <- 'lc'
  zones <- c(slope_lyr, lc_raster)
  zones_info <- list(has_eco=FALSE, has_boreal=FALSE, has_country=FALSE,
                     has_boreal_eco=FALSE, relative_to_north_lat=NULL)
  # slope_raster is only used as a template raster for the rasterizer
  # e.g, dims, res, crs, ...
  ecoregions_lyr <- rasterize_boundaries(slope_raster, ecoregions, 'eco', field='ECO_ID')
  if (any(!is.na(values(ecoregions_lyr)))){
    zones <- c(zones, ecoregions_lyr)
    zones_info$has_eco <- TRUE
  }

  countries_lyr <- rasterize_boundaries(slope_raster, countries, 'country', field='iso3_code')
  if (any(!is.na(values(countries_lyr)))){
    zones <- c(zones, countries_lyr)
    zones_info$has_country <- TRUE
  }

  boreal_lyr <- rasterize_boundaries(slope_raster, boreal_poly, 'boreal')
  if (any(!is.na(values(boreal_lyr)))){
    zones <- c(zones, boreal_lyr)
    zones_info$has_boreal <- TRUE
  }

  if (zones_info$has_eco && zones_info$has_boreal){
    boreal_eco_lyr <- rasterize_boundaries(
      slope_raster,
      terra::intersect(ecoregions, boreal_poly),
      'boreal_eco',
      field='ECO_ID'
    )
    if (any(!is.na(values(boreal_eco_lyr)))){
      zones <- c(zones, boreal_eco_lyr)
      zones_info$has_boreal_eco <- TRUE
    }
  }

  north_result <- clip_to_north_lat(slope_raster)
  zones_info$relative_to_north_lat <- north_result[['relative_to_north_lat']]
  if (north_result[['relative_to_north_lat']] == 'intersects'){
    north_lyr <- rasterize_boundaries(slope_raster, north_result[['poly_4326_north']], 'north')
    zones <- c(zones, north_lyr)
  }

  return(list(zones=zones, zones_info=zones_info))
}

calculate_zonal_summary <- function(map, zones, zones_info, agg_fun, cores, map_name, result_name) {

  summary_fun <-function(lyr){
    aggregates <- zonal(map, zones[[lyr]], agg_fun, na.rm = TRUE)
    names(aggregates)[names(aggregates) == map_name] <- result_name

    counts <- zonal(map, zones[[lyr]], "notNA", na.rm = TRUE)
    names(counts)[names(counts) == map_name] <- 'count'

    full_join(aggregates, counts, by=lyr)
  }

  tasks <- list(
    by_lc = function(){
      summary_fun('lc')
    },
    by_slope = function(){
      summary_fun('slope')
    },
    by_ecoregion = function(){
      if(zones_info$has_eco)
        summary_fun('eco')
      else NULL
    },
    by_country = function(){
      if(zones_info$has_country)
        summary_fun('country')
      else NULL
    },
    boreal = function(){
      if(zones_info$has_boreal)
        summary_fun('boreal')
      else NULL
    },
    boreal_eco = function(){
      if(zones_info$has_boreal_eco)
        summary_fun('boreal_eco')
      else NULL
    },
    north = function(){
      if(zones_info$relative_to_north_lat == 'intersects')
        summary_fun('north')
      else if(zones_info$relative_to_north_lat == 'south')
        NULL
    }
  )

  results_list <- mclapply(tasks, function(f) f(), mc.cores = cores)
  agg_fun <- match.fun(agg_fun)
  results_list$overall <- results_list$by_slope |>
    summarise(!!result_name := agg_fun(.data[[result_name]]), count=sum(count))

  if (zones_info$relative_to_north_lat == 'north') {
    results_list$north <- results_list$overall
  }
  return(results_list)
}


fit_model <- function(model, model_config, train_df, pred_vars, predict_var){
  y_fit <- if (predict_var == 'Ht') train_df$h_canopy else train_df$AGB
  x_fit <- train_df[pred_vars]

  model_fit <- do.call(model, modifyList(model_config, list(y=y_fit, x=x_fit)))
  return(model_fit)
}

run_modeling_pipeline <-function(biomass_models, all_train_data, zones, zones_info,
                                 model, model_config, iter,
                                 predict_function, cores, agg_fun, map_name, summary_column_name,
                                 max_samples, sample, pred_vars, predict_var, stack){
  t1 <- Sys.time()
  print('creating AGB traing data frame.')
  train_df <- GEDI2AT08AGB(biomass_models, all_train_data, iter=iter, max_samples, sample)

  print('fitting model')
  model <- fit_model(model, model_config, train_df, pred_vars, predict_var)

  gc(verbose=TRUE)

  print('predicting biomass map')
  preds <- predict_function(model, stack)
  names(preds) <- map_name

  print('calculating zonal summaries')
  zonal_summary <- calculate_zonal_summary(preds, zones, zones_info, agg_fun,
                                           cores, map_name, summary_column_name)

  model_stats <- get_model_stats(model)
  print('model stats:')
  print(model_stats[c('OOB_MSE', 'OOB_RMSE', 'OOB_R2')])

  t2 <- Sys.time()
  cat('pipeline runtime:', difftime(t2, t1, units="mins"), ' (m)\n')

  return(list(model_stats=model_stats, map=preds, zonal_summary=zonal_summary))
}

welford_update <- function(count, mu, M2, new_value){
    count <- count + 1
    delta = new_value - mu
    mu <- mu + delta / count
    delta2 = new_value - mu
    M2 <- M2 + delta * delta2
    return (list(count=count, mu=mu, M2=M2))
}

run_uncertainty_calculation <- function(fixed_modeling_pipeline_params, n_iters){
  results <- do.call(run_modeling_pipeline, modifyList(
    fixed_modeling_pipeline_params,
    list(iter=1)
  ))

  summary_keys <- names(results$zonal_summary)
  valid_keys <- summary_keys[sapply(results$zonal_summary[summary_keys], Negate(is.null))]
  this_iter <- 1
  for (k in valid_keys){
    results$zonal_summary[[k]]$iter <- this_iter
  }

  model_stats <- list(results[['model_stats']])

  # initializing to 0, with crs of mu
  mu <- c(results[['map']])
  M2 <- mu
  values(M2) <- 0.0

  while(this_iter < n_iters){
    cat('Uncertainty loop, iteration:', this_iter, '\n')
    params <- modifyList(
      fixed_modeling_pipeline_params,
      list(iter=this_iter+1)
    )
    new_results <- do.call(run_modeling_pipeline, params)

    updated <- welford_update(this_iter, mu, M2, new_results[['map']])
    mu <- updated[['mu']]
    M2 <- updated[['M2']]

    this_iter <- this_iter + 1

    for (k in valid_keys){
      new_results$zonal_summary[[k]]$iter <- this_iter
      results$zonal_summary[[k]] <- rbind(
        results$zonal_summary[[k]],
        new_results$zonal_summary[[k]]
      )
    }

    model_stats[[this_iter]] <- new_results[['model_stats']]
  }

  s <- calculate_zonal_summary(mu, params$zones, params$zones_info,
                               params$agg_fun, params$cores, params$map_name,
                               params$summary_column_name)
  for (k in valid_keys){
    s[[k]]$iter <- -1
    results$zonal_summary[[k]] <- rbind(results$zonal_summary[[k]], s[[k]])
  }

  std <- sqrt(M2 / (this_iter - 1))
  names(std) <- paste0('std_', if (params$predict_var == 'AGB') 'agbd' else 'ht')

  return(list(
    map=mu, std=std,
    zonal_summary=results$zonal_summary[valid_keys],
    model_stats=model_stats
  ))
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

write_ensemble_stats <- function(val_df, out_fn=NULL){
  val_df$res <- val_df$y_true - val_df$y_pred
  lm_ <- lm(y_true ~ y_pred, data=val_df)

  stats <- data.frame(
    MAE=mean(abs(val_df$res)),
    MSE=mean(val_df$res^2),
    RMSE=sqrt(mean(val_df$res^2)),
    R2=1-sum(val_df$res^2)/sum((val_df$y_true-mean(val_df$y_true))^2),
    lm_R2=summary(lm_)$r.squared,
    lm_RMSE=sqrt(mean((val_df$y_true-predict(lm_))^2))
  )

  if (!is.null(out_fn)){
    write.csv(stats, out_fn)
  }
  return(stats)
}

sample_map_at_lidar_points <-function(df, map, biomass_models, year, predict_var, out_fn=NULL){
  df <- df[df$y == as.integer(year), ]
  df <- GEDI2AT08AGB(biomass_models, df, iter=1, sample=FALSE)

  points_vect <- vect(df, geom = c('lon', 'lat'), crs = 'EPSG:4326')

  if (crs(points_vect) != crs(map)) {
    points_vect <- project(points_vect, crs(map))
  }

  # Extract mean Ht or AGB (pred_var) values from first layer
  extracted_values <- extract(map[[1]], points_vect)

  target <- if (predict_var == 'AGB') 'AGB' else 'h_canopy'

  val_df <- data.frame(
    lat = df[['lat']],
    lon = df[['lon']],
    y_pred = extracted_values[, 2],
    y_true = df[[target]]
  )
  val_df <- na.omit(val_df)

  if (!is.null(out_fn) && nrow(val_df) > 0){
    write.csv(val_df, out_fn)
  }

  return(val_df)
}

convert_AGBD_Mg_to_AGB_Pg <- function(results, AGBD_column_name){
  # assuming a pixel is 900 m^2
  Mg_per_ha_to_Pg <- 0.09 * 1e-9

  for (k in seq_along(results)){
    results[[k]][[AGBD_column_name]] <- results[[k]][[AGBD_column_name]] * Mg_per_ha_to_Pg
    names(results[[k]])[names(results[[k]]) == AGBD_column_name] <- 'AGB_Pg'
  }

  return(results)
}

mapBoreal <- function(atl08_path,
                      broad_path,
                      hls_path,
                      topo_path,
                      lc_path,
                      boreal_vector_path,
                      ecoregions_path,
                      countries_path,
                      biomass_models_path,
                      year,
                      sar_path=NULL,
                      mask=TRUE,
                      max_sol_el=5,
                      offset=100,
                      minDOY=130,
                      maxDOY=250,
                      expand_training=TRUE,
                      n_iters=30,
                      local_train_perc=100,
                      min_samples=5000,
                      max_samples=10000,
                      cores=2,
                      ntree=50,
                      predict_var='AGB',
                      pred_vars=c('elevation', 'slope', 'NDVI'),
                      zero_short_veg_height=FALSE,
                      slope_thresh=15,
                      val_thresh=11000,
                      val_frac=0.10
                      )
{

  tile_num = tail(unlist(strsplit(path_ext_remove(atl08_path), "_")), n=1)
  cat("Modelling and mapping boreal AGB tile: ", tile_num, "\n")

  pred_vars <- parse_pred_vars(pred_vars, remove_sar=is.null(sar_path))
  print(pred_vars)

  stack <- resample_reproject_and_mask(topo_path, hls_path, lc_path, pred_vars, mask, sar_path)
  boreal_poly <- vect(boreal_vector_path)
  ecoregions <- vect(ecoregions_path)
  countries <- vect(countries_path)
  zones <- prep_summary_layers(
    stack[['slope']],
    stack[['esa_worldcover_v100_2020']],
    ecoregions,
    boreal_poly,
    countries
  )
  # there are rge objects and now rasterized in zones layers above
  rm(boreal_poly)
  rm(ecoregions)
  rm(countries)
  biomass_models <- read_randomized_biomass_models(biomass_models_path, n_iters)
  output_fns <- set_output_file_names(predict_var, tile_num, year)

  all_data <- prepare_training_data(
    atl08_path, broad_path, expand_training, minDOY,
    maxDOY, max_sol_el, min_samples, local_train_perc, offset, names(stack),
    zero_short_veg_height, slope_thresh,
    as.integer(year), val_thresh, val_frac, biomass_models, output_fns[['training']]
  )

  if (cores > 1) {
    # saving stack and passing path instead because otherwise when parallel section runs
    # the memory overhead multiplies by number of nodes
    stack_path <- './stack.tif'
    writeRaster(stack, stack_path, overwrite = TRUE)
    rm(stack); gc()
    stack <- stack_path
  }

  map_name <- if(predict_var=='AGB') 'mean_agbd' else 'mean_ht'
  agg_fun <- if(predict_var=='AGB') 'sum' else 'mean'

  fixed_modeling_pipeline_params <- list(
    biomass_models=biomass_models, all_train_data=all_data[['train_data']],
    pred_vars=pred_vars, predict_var=predict_var, stack=stack, zones=zones[['zones']],
    zones_info=zones[['zones_info']], cores=cores, max_samples=max_samples,
    map_name=map_name,
    agg_fun=agg_fun,
    summary_column_name=paste0(map_name, '_', agg_fun),
    model=randomForest, model_config=list(ntree=ntree), sample=TRUE,
    predict_function=create_predict_function(cores=cores)
  )

  results <- run_uncertainty_calculation(fixed_modeling_pipeline_params, n_iters)
  cat(predict_var,  'successfully predicted!\n')

  if (predict_var == 'AGB'){
    converted_summary <- convert_AGBD_Mg_to_AGB_Pg(
      results[['zonal_summary']],
      fixed_modeling_pipeline_params$summary_column_name
    )
  }
  else {
    converted_summary <- results[['zonal_summary']]
  }
  write_output_summaries_and_stats(
    converted_summary,
    results[['model_stats']],
    output_fns
  )

  write_output_raster_map(results[['map']], results[['std']], output_fns[['map']])

  if (!is.null(all_data[['val_data']])){
    val_df <- sample_map_at_lidar_points(
      all_data[['val_data']],
      results[['map']],
      fixed_modeling_pipeline_params[['biomass_models']],
      year,
      predict_var,
      output_fns[['validation']]
    )
    if (nrow(val_df) > 0) {
      stats <- write_ensemble_stats(val_df, output_fns[['ensemble_stats']])
      print(stats)
    }
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
    c("--ecoregions_path"), type = "character",
    help = "Path to the ecoregions vector file",
    ),
  make_option(
    c("--countries_path"), type = "character",
    help = "Path to the world admin boundary (at country level) vector file with a iso3_code attrbitute",
    ),
  make_option(
    c("--biomass_models_path"), type = "character",
    help = "Path to the tarbarr of biomass rds models",
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
    c("--n_iters"), type = "integer", default = 30,
    help = "number of bootstrap iterations, must be >= 2 [default: %default]"
  ),
  make_option(
    c("-c", "--cores"), type = "integer", default = 2,
    help = "Number of cores used for parallel prediction steps [default: %default]"
  ),
  make_option(
    c ("--ntree"), type = "integer", default = 50,
    help = "Number of random forest trees [default: %default]"
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
    c ("--val_thresh"), type = "numeric", default = 11000,
    help = "min number of atl08 samples (after all filtering) needed to validate the model [default: %default]"
  ),
  make_option(
    c ("--val_frac"), type = "numeric", default = 0.10,
    help = "Fraction of the current year's atl08 samples to use for validation if val_thresh has reached [default: %default]"
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
for(arg in c('atl08_path', 'broad_path', 'topo_path', 'hls_path', 'lc_path',
             'boreal_vector_path', 'ecoregions_path', 'biomass_models_path',
             'countries_path')){
  if (is.null(opt[[arg]])){
    # TODO: some of these args should actually be optional.
    stop(paste0("ERROR: --",arg, " is required."))
  }
}

set.seed(27182)
do.call(mapBoreal, opt)
