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
####  - models_id: models id
####  - stack: a combined raster stack of Landsat and Copernicus DEM data
####  - ice2_30_atl08: list containing the path to the data tables
####  - offset: offset applied in the model

#### iii) Outputs
####  COGs of predicted 30-m AGBD, SD AGBD

#----------------------------------------------#
############# functions ########################
#----------------------------------------------#

#applyModels takes a list of models and applies them to a raster stack, calculates per tile total (mean for height). 
applyModels <- function(models=models,
                           stack=stack,
                           pred_vars=pred_vars,
                           predict_var=predict_var,
                           tile_num){
        rem <- length(models)
        if(rem>1){
            models <- models[-rem]
        }

        #split stack into list of files
        if(ppside > 1){
            tile_list <- SplitRas(raster=stack,ppside=ppside,save=FALSE)
    
        print(paste0('tiles successfully split into ', length(tile_list), ' tiles'))

        #run mapping over each tile in a loop, create a list of tiled rasters for each layer

        n_subtiles <- length(tile_list)
        print(paste0('nsubtiles:', n_subtiles))
            if(exists('out_map')==TRUE){rm(out_map)}
            
         for (tile in 1:n_subtiles){
            tile_stack <- tile_list[[tile]]
            #for a subtile that is all NA, combine it with the next subtile
             print('tile number:')
             print(tile)
              if(predict_var=='AGB'){
                maps<-agb_mapping(
                         model_list=models,
                         tile_num=tile_num,
                         stack=tile_stack,
                         boreal_poly=boreal_poly)
                  
                if((exists('out_map')==FALSE) | tile==1){
                     if(length(maps)>1){
                         out_map <- maps[[1]]
                         tile_total <- maps[[2]]
                     }
                 } 

                 if((exists('out_map')==TRUE) & (length(maps)>1) & (tile>1)){
                     out_map <- mosaic(maps[[1]], out_map, fun="max")
                     if(exists('tile_total')==FALSE){tile_total <- 0.0}
                     print('tile total:')
                     tile_total <- tile_total + maps[[2]]
                     rm(maps)
                 }
             }
             
             if(predict_var=='Ht'){
                 maps<-ht_mapping(
                         model_list=models,
                         tile_num=tile_num,
                         stack=tile_stack,
                         boreal_poly=boreal_poly)
             
             if((exists('out_map')==FALSE) | tile==1){
                     if(length(maps)>1){
                         out_map <- maps[[1]]
                         tile_mean <- maps[[2]]$Tile_Mean
                         print(tile_mean)
                     }
                 } 

                 if((exists('out_map')==TRUE) & (length(maps)>1) & (tile>1)){
                     out_map <- mosaic(maps[[1]], out_map, fun="max")
                     if(exists('tile_mean')==FALSE){tile_mean <- maps[[2]]$Tile_Mean}
                     print('tile mean:')
                     tile_mean <- mean(c(tile_mean, maps[[2]]$Tile_Mean))
                     print(tile_mean)
                     rm(maps)
                 }
             }   
            }
           }
        if (ppside == 1){
          mapping_fun <- if (predict_var == 'AGB') agb_mapping else ht_mapping
          out_map <- mapping_fun(model_list=models, tile_num=tile_num, stack=stack, boreal_poly=boreal_poly)
        }
          return(list(out_map[[1]], NULL))
    }

combine_temp_files <- function(target, tile_num){
  if (target == 'AGB'){
    pattern <- '_total.csv'
    summary_fun <- sum
    names <- c('tile_total', 'tile_boreal_total')
    output_prefix <- 'output/boreal_agb'
    output_suffix <- '_total_all.csv'
  }
  else if (target == 'Ht'){
    pattern <- '_mean.csv'
    summary_fun <- mean
    names <- c('tile_mean', 'tile_boreal_mean')
    output_prefix <- 'output/boreal_ht'
    output_suffix <- '_mean_all.csv'
  }
  else {stop('Target should be one of AGB or Ht')}

  csv_files <- list.files(path='output', pattern=pattern, full.names=TRUE)
  all_data <- bind_cols(lapply(csv_files, read.csv))
  file.remove(csv_files)

  combined_boreal <- all_data |> select(matches('Boreal')) |> apply(1, summary_fun, na.rm = TRUE)
  combined_tile <- all_data |> select(matches('Tile')) |> apply(1, summary_fun, na.rm = TRUE)
  combined <- data.frame(cbind(combined_tile, combined_boreal))
  names(combined) <- names

  out_fn_stem = paste(output_prefix, format(Sys.time(),"%Y%m%d%s"), str_pad(tile_num, 4, pad = "0"), sep="_")
  out_fn <- paste0(out_fn_stem, output_suffix)
  write.csv(file=out_fn, combined, row.names=FALSE)
  return(combined_tile)
}

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

reformat_for_AGB_prediction <- function(in_data, offset){
  RH_columns <- get_height_column_names(in_data)

  in_data <- in_data |>
    select(c(RH_columns, segment_landcover)) |>
    mutate(across(RH_columns, ~. + offset)) |>
    mutate(model_id = case_when(
      segment_landcover %in% c(111, 113, 121, 123) ~ "m8",# "m3",
      segment_landcover %in% c(112, 114, 122, 124) ~ "m8", # "m1",
      TRUE ~ "m8"
    )
    ) |>
    select(-segment_landcover)

  return(in_data)
}

randomize <- function(model_i){
  # modify coeffients through sampling variance covariance matrix
  model_varcov <- vcov(model_i)
  coeffs <- model_i$coefficients

  mod.coeffs <- mvrnorm(n = 50, mu=coeffs, Sigma = model_varcov)
  model_i$coefficients <- mod.coeffs[1,]

  return(model_i)
}

GEDI2AT08AGB<-function(rds_models, models_id, in_data, offset=100, DO_MASK=FALSE, one_model=TRUE, max_n=5000.0, sample=TRUE){
  # TODO the three steps below should actually be done only once outside in the caller
  in_data <- na.omit(as.data.frame(in_data))
  if (sample && nrow(in_data) > max_n)
    in_data <- reduce_sample_size(in_data, max_n)
  in_data <- rename_height_columns_to_match_pretrained_models(in_data)

  df <- reformat_for_AGB_prediction(in_data, offset)

  df$AGB<-NA
  df$SE<-NA

  ids<-unique(df$model_id)
  n_models <- length(ids)

  for (i in ids){
    model_i<-readRDS(rds_models[names(rds_models)==i])

    # Modify coeffients through sampling variance covariance matrix
    if(!one_model)
      model_i <- randomize(model_i)

    # Predict AGB and SE
    df$AGB[df$model_id==i] <- predict(model_i, newdata=df[df$model_id==i,])
    df$SE[df$model_id==i] <- summary(model_i)$sigma^2

    df$AGB[df$AGB < 0] <- 0.0

    # Calculate Correction Factor C
    C <- mean(model_i$fitted.values^2)/mean(model_i$model$`sqrt(AGBD)`^2)

    # Bias correction in case there is a systematic over or under estimation in the model
    df$AGB[df$model_id==i] <- C*(df$AGB[df$model_id==i]^2)
  }
  in_data <- in_data |> bind_cols(AGB=df$AGB, SE=df$SE)

  # Apply slopemask, validmask and landcover masks
  bad_lc <- c(0, 60, 80, 200, 50, 70)
  in_data$AGB[in_data$slopemask == 0 |
                in_data$ValidMask == 0 |
                in_data$segment_landcover %in% bad_lc] <- 0.0

  return(in_data)
}

# stats
StatModel <- function( obs, est){
  xy<-na.omit(cbind(obs,est))
  obs<-xy[,1]
  est<-xy[,2]
  rmse <- sqrt( sum(( est - obs )^2)/length(obs) ) # Root mean square error
  bias <- mean( est - obs ) # bias
  rmseR <- 100 * sqrt( sum(( est - obs )^2)/length(obs) ) / mean( obs )
  biasR <- 100 * mean( est - obs ) / mean( obs )
  r <- cor(est,obs)
  r2<-summary(lm(obs~est))$r.squared
  Stats<-data.frame( Stat=c("rmse","rmseR","bias","biasR","r","r2"),
                     Values=round(c(rmse,rmseR,bias,biasR,r,r2),2)) 
  return(Stats)
}

stratRandomSample<-function(agb=y,breaks, p){
  #require(data.table)
  n<-length(agb)
  ids<-1:n
  s<-round(n*p)
  agb[agb==0]<-0.0000000001
  ids_cut<-cut(agb,breaks=breaks, labels=F)
  df<-cbind(agb,ids,ids_cut)
  df<-data.table(df[!is.na(df[,1]),])
  number_sample<-ceiling(s/(length(breaks)-1))
  sel_all<-df[,.SD[sample(.N, min(number_sample,.N), replace = T)],by=ids_cut]
  return(ids_selected=sel_all$ids)
}

agbModeling<-function(rds_models, models_id, in_data, pred_vars, offset=100, DO_MASK, rep=100, predict_var){
  model_list <- list()
  rep <- rep + 1
  for (j in 1:rep){
    # reduce max_n for faster modeling
    current_max_n <- if (j == 1) max_n else 1000
    one_model <- j == 1

    xtable <- GEDI2AT08AGB(rds_models=rds_models,
                       models_id=models_id,
                       in_data=in_data,
                       offset=offset,
                       DO_MASK=DO_MASK, # not needed
                       one_model=one_model,
                       max_n=current_max_n,
                       sample=TRUE)
    if (j == 1)
      AGB_training_table <- xtable

    y_fit <- if (predict_var == 'Ht') xtable$RH_98 else xtable$AGB
    x_fit <- xtable[pred_vars]

    # TODO need to pass these models as a param not hard coded RF
    if (j == 1)
      rf_model <- randomForest(y=y_fit, x=x_fit, ntree=NTREE, mtry=6)
    else
      rf_model <- randomForest(y=y_fit, x=x_fit, ntree=NTREE)

    model_list <- list.append(model_list, rf_model)
  }

  return(list(AGB_training_table=AGB_training_table, model_list=model_list))
}

#split raster into subtiles, run mapping, recombine

# The function spatially aggregates the original raster
# it turns each aggregated cell into a polygon
# then the extent of each polygon is used to crop
# the original raster.
# The function returns a list with all the pieces
# in case you want to keep them in the memory. 
# it saves and plots each piece
# The arguments are:
# raster = raster to be chopped            (raster object)
# ppside = pieces per side                 (integer)

SplitRas <- function(raster,ppside,save){
  h        <- ceiling(ncol(raster)/ppside)
  v        <- ceiling(nrow(raster)/ppside)
  agg      <- aggregate(raster,fact=c(h,v))
  agg[]    <- 1:ncell(agg)
  agg_poly <- as.polygons(agg)
  names(agg_poly) <- "polis"
  r_list <- list()
  for(i in 1:ncell(agg)){
    e1          <- ext(agg_poly[agg_poly$polis==i,])
    r_list[[i]] <- crop(raster,e1)
  }
  if(save==T){
    for(i in 1:length(r_list)){
      writeRaster(r_list[[i]],filename=paste("SplitRas",i,sep=""),
                  format="GTiff",datatype="FLT4S",overwrite=TRUE)  
    }
  }
  return(r_list)
}

agb_mapping <- function(model_list=model_list, tile_num=tile_num, stack=stack, boreal_poly=boreal_poly) {
  return(generic_mapping(
    model_list=model_list, tile_num=tile_num, stack=stack, boreal_poly=boreal_poly,
    summary_fun='sum', convert_fun=function(x){(x*0.09)/1000000000},
    out_df_names=c('Tile_Total', 'Boreal_Total'), output_type='abg', output_csv_suffix='_total'
    )
  )
}

ht_mapping <- function(model_list=model_list, tile_num=tile_num, stack=stack, boreal_poly=boreal_poly) {
  return(generic_mapping(
    model_list=model_list, tile_num=tile_num, stack=stack, boreal_poly=boreal_poly,
    summary_fun='mean', convert_fun=NULL, out_df_names=c('Tile_mean', 'Boreal_mean'),
    output_type='ht', output_csv_suffix='_mean'
    )
  )
}

generic_mapping <-function(model_list, tile_num, stack, boreal_poly, summary_fun, convert_fun,
                           out_df_names, output_type, output_csv_suffix) {
  pred_stack <- na.omit(stack)
  pred_map = c()
  tile_summary = c()
  boreal_summary = c()
  n_models <- length(model_list)

  for (model_i in model_list){
    print('generic mapping iter')
    pred_map_i <- predict(pred_stack, model_i, na.rm=TRUE)
    # set slope and valid mask to zero
    pred_map_i <- mask(pred_map_i, pred_stack$slopemask, maskvalues=0, updatevalue=0)
    pred_map_i <- mask(pred_map_i, pred_stack$ValidMask, maskvalues=0, updatevalue=0)
    pred_map <- c(pred_map, pred_map_i)
    print('map pred done')
    pred_map_conv_i <- if (is.null(convert_fun)) pred_map_i else app(pred_map_i, convert_fun)
    tile_summary_i <- global(pred_map_conv_i, summary_fun, na.rm=TRUE)[[summary_fun]]
    tile_summary <- c(tile_summary, tile_summary_i)
    print('tile summary done')
    # repeat for just boreal
    boreal_i <- extract(pred_map_conv_i, boreal_poly, fun=summary_fun, na.rm=TRUE)
    boreal_summary_i <- if(summary_fun=='sum') sum(boreal_i$lyr.1, na.rm=TRUE) else boreal_i$lyr1[1]
    boreal_summary <- c(boreal_summary, boreal_summary_i)
    print('boreal summary done')
  }
  pred_map <- rast(pred_map)
  if (n_models > 1)
    sd_map <- app(pred_map, sd)

  tile_and_boreal_summary_df <- as.data.frame(cbind(tile_summary, boreal_summary))
  names(tile_and_boreal_summary_df) <- out_df_names
  saveDataFrame(tile_and_boreal_summary_df, output_type, tile_num, output_csv_suffix)

  if(n_models>1)
    maps <- list(c(pred_map[[1]], sd_map, pred_map), tile_and_boreal_summary_df)
  else
    maps <- list(c(pred_map[[1]]), tile_and_boreal_summary_df)

  return(maps)
}

saveDataFrame <- function(df, product, tile_num, suffix){
  out_fn_stem = paste("output/boreal", product, format(Sys.time(),"%Y%m%d%s"), str_pad(tile_num, 4, pad = "0"), sep="_")
  out_fn <- paste0(out_fn_stem, suffix, '.csv')
  write.csv(file=out_fn, df, row.names=FALSE)
}

partial_sd <- function(arr){
  partial_sd_arr <- rep(0, length(arr) - 1)
  for(i in 2:length(arr)){
    partial_sd_arr[i-1] <- sd(arr[1:i], na.rm = T)
  }
  return(partial_sd_arr)
}

sd_change_relative_to_baseline <- function(arr, n){
  partial_sd_arr <- partial_sd(arr)
  paritial_std_arr_last_n_out <- head(partial_sd_arr, max(1, length(arr) - n))

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

combine_csv_outpus <- function(target, tile_num){
  if (target == 'AGB'){
    pattern <- '_total_all.csv'
    names <- c('tile_total', 'tile_boreal_total')
    output_prefix <- 'output/boreal_agb'
    output_suffix <- '_total_iters.csv'
  }
  else if (target == 'Ht'){
    pattern <- '_mean_all.csv'
    names <- c('tile_mean', 'tile_boreal_mean')
    output_prefix <- 'output/boreal_ht'
    output_suffix <- '_mean_iters.csv'
  }
  else {stop('Target should be one of AGB or Ht')}

  csv_files <- list.files(path='output', pattern=pattern, full.names=TRUE)

  out_df <- bind_rows(lapply(csv_files, read.csv))
  names(out_df) <- names
  file.remove(csv_files)

  out_fn_stem = paste(output_prefix, format(Sys.time(),"%Y%m%d%s"), str_pad(tile_num, 4, pad = "0"), sep="_")
  output_fn = paste0(out_fn_stem, output_suffix)
  write.csv(file=output_fn, out_df, row.names=FALSE)

  return(out_fn_stem)
}

set_output_file_names <- function(out_fn_stem){
  fn_types <- c('tmp.tif', '.tif', '.csv', '_train_data.csv', '_stats.Rds', '_model.Rds')
  output_file_names <- paste0(out_fn_stem, fn_types)

  names <- c('tif', 'cog', 'csv', 'train', 'stats', 'model')
  names(output_file_names) <- names

  return(output_file_names)
}


write_ATL08_table <- function(target, df, out_file_path){
    out_columns <- if(target=='AGB') c('lon', 'lat', 'AGB', 'SE') else c('lon', 'lat', 'RH_98')
    write.csv(df[, out_columns], file=out_file_path, row.names=FALSE)
}

write_single_model_summary <- function(df, target, pred_vars, out_fns, nrow_tile){
    target <- if(target == 'AGB') df$AGB else df$RH_98

    rf_single <- randomForest(y=target, x=df[pred_vars], ntree=NTREE, importance=TRUE, mtry=6)
    local_model <- lm(rf_single$predicted[1:nrow_tile] ~ target[1:nrow_tile], na.rm=TRUE)
    saveRDS(rf_single, file=out_fns['model'])

    rsq <- max(rf_single$rsq, na.rm=T)
    cat('rsq: ', rsq, '\n')

    rsq_local <- summary(local_model)$r.squared
    cat('rmax_iters <- sq_local: ', rsq_local, '\n')

    na_data <- which(is.na(local_model$predicted==TRUE))

    if(length(na_data) == 0)
      rmse_local <- sqrt(mean(local_model$residuals^2))

    cat('rmse_local: ', rmse_local, '\n')

    imp_vars <- rf_single$importance
    out_accuracy <- list(rsq_local, rmse_local, imp_vars)
    saveRDS(out_accuracy, file=out_fns['stats'])
}

write_output_raster_map <- function(out_map, out_map_all, output_fn){
  # change -9999 to NA
  out_map <- subst(out_map, -9999, NA)
  out_sd <- app(out_map_all, sd)
  out_sd <- subst(out_sd, -9999, NA)
  out_map <- c(out_map, out_sd)
  NAflag(out_map)
  options <- c("COMPRESS=LZW", overwrite=TRUE, gdal=c("COMPRESS=LZW", "OVERVIEW_RESAMPLING=AVERAGE"))
  writeRaster(out_map, filename=output_fn, filetype="COG", gdal=options)
  cat("Write COG tif: ", output_fn, '\n')
}


mapBoreal<-function(rds_models,
                    models_id,
                    ice2_30_atl08_path, 
                    ice2_30_sample_path,
                    offset=100,
                    s_train=70, 
                    rep=10,
                    ppside=2,
                    stack=stack,
                    strat_random=TRUE,
                    output,
                    minDOY=121,
                    maxDOY=273,
                    max_sol_el=0,
                    expand_training=TRUE,
                    local_train_perc=100,
                    min_icesat2_samples=3000,
                    DO_MASK=FALSE,
                    boreal_poly=boreal_poly,
                    predict_var,
                    max_n=3000,
                    pred_vars=c('elev', 'slope')){

    # Get tile num
    tile_num = tail(unlist(strsplit(path_ext_remove(ice2_30_atl08_path), "_")), n=1)
    print("Modelling and mapping boreal AGB")
    cat("tile: ", tile_num, '\n')
    cat("ATL08 input: ", ice2_30_atl08_path, '\n')

    tile_data <- read.csv(ice2_30_atl08_path)

    if(expand_training)
      tile_data <- expand_training_around_growing_season(tile_data, minDOY, maxDOY, max_sol_el, min_icesat2_samples)

    cat('n_avail training:', nrow(tile_data), '\n')
    tile_data <- remove_stale_columns(tile_data, c("binsize", "num_bins"))
    broad_data <- read.csv(ice2_30_sample_path)
    broad_data <- remove_stale_columns(broad_data, c("X__index_level_0__", "geometry"))
    print(ice2_30_sample_path)

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

    str(all_train_data)
    all_train_data <- remove_height_outliers(all_train_data)

    tile_data_output <- tile_data # probably not needed
    cat('table for model training generated with ', nrow(all_train_data), ' observations\n')

    models<-agbModeling(rds_models=rds_models,
                            models_id=models_id,
                            in_data=all_train_data,
                            pred_vars=pred_vars,
                            offset=offset,
                            DO_MASK=DO_MASK,
                            rep=rep,
                            predict_var=predict_var)
    
    print('model fitting complete!')

    final_map <- applyModels(models[['model_list']], stack, pred_vars, predict_var, tile_num)

    xtable <- models[['AGB_training_table']]
    combined_totals <- combine_temp_files(predict_var, tile_num)

    #subset out the iteration bands
    out_map_all <- subset(final_map[[1]], 3:nlyr(final_map[[1]]))
    
    #just pull the mean for out_map, sd will be added later
    out_map <- subset(final_map[[1]], 1)
    
    rm(final_map)

    #set the variance threshold - 0.05 = 5%
    var_thresh <- 0.05
    
    if(rep>1){
        var_diff <- sd_change_relative_to_baseline(combined_totals, 9)
        cat('var_diff:', var_diff, '\n')

        #if larger difference, need more models and more iterations
        #save(combined_totals, file='/projects/lduncanson/testing/test_totals.Rdata')
        #set some maximum number of iterations
        max_iters <- 100
        # this if statement seems to have no effect. Initially length(combined_totals) = 2
        # which is < 100 (= max_iters), was it meant to be part of the while loop condition?
        if(length(combined_totals)<max_iters){
            while(var_diff > var_thresh){
            print('Adding more interations...')
            new_models <- agbModeling(rds_models=rds_models,
                            models_id=models_id,
                            in_data=all_train_data,
                            pred_vars=pred_vars,
                            offset=offset,
                            DO_MASK=DO_MASK,
                            rep=10,
                            predict_var=predict_var)
                
            new_final_map <- applyModels(new_models[['model_list']], stack, pred_vars, predict_var, tile_num)
            combined_totals_new <- combine_temp_files(predict_var, tile_num)

            #combine original map with new iterations map
            out_map_all <- c(out_map_all, subset(new_final_map[[1]], 3:nlyr(new_final_map[[1]])))
            rm(new_final_map)
            combined_totals <- c(combined_totals, combined_totals_new)
            var_diff <- sd_change_relative_to_baseline(combined_totals, 9)

            if(length(combined_totals)>75){
                var_thresh <- 0.06
                }
            if(length(combined_totals)>100){
                var_thresh <- 0.08
                }
            if(length(combined_totals)>200){
                var_thresh <- 0.1
                }
            }
        }
    }

    out_fn_stem <- combine_csv_outpus(predict_var, tile_num)
    print('AGB successfully predicted!')
    print('mosaics completed!')

    out_fns <- set_output_file_names(out_fn_stem)

    write_output_raster_map(out_map, out_map_all, out_fns['cog'])
    write_ATL08_table(predict_var, xtable, out_fns['train'])
    write_single_model_summary(xtable, predict_var, pred_vars, out_fns, nrow(tile_data_output))

    print("Returning names of COG and CSV...")
    return(list(out_fns['cog'], out_fns['csv']))
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
# print('pred_vars:')
# print(pred_vars)

#for debugging replace args with hard paths
#data_table_file <- '/projects/my-private-bucket/dps_output/run_tile_atl08_ubuntu/tile_atl08/2022/11/30/19/22/04/120959/atl08_005_30m_filt_topo_landsat_20221130_1216.csv'
#topo_stack_file <- '/projects/shared-buckets/nathanmthomas/alg_34_testing/Copernicus_1216_covars_cog_topo_stack.tif'
#l8_stack_file <- '/projects/shared-buckets/nathanmthomas/alg_34_testing/HLS_1216_06-15_09-01_2019_2021.tif'
#LC_mask_file <- '/projects/shared-buckets/nathanmthomas/alg_34_testing/esa_worldcover_v100_2020_1216_cog.tif'

data_table_file <- '~/Downloads/dps_output/atl08_006_030m_2020_2020_06_09_filt_covars_merge_neighbors_034673.csv'
topo_stack_file <- '~/Downloads/dps_output/CopernicusGLO30_34673_cog_topo_stack.tif'
l8_stack_file <- '~/Downloads/dps_output/HLS_34673_07-01_08-31_2023_2023.tif'
LC_mask_file <- '~/Downloads/dps_output/esa_worldcover_v100_2020_34673_cog.tif'
data_sample_file <- '~/Downloads/dps_output/boreal_train_data_2020_n3.csv'
boreal_vect <- '~/Downloads/dps_output/wwf_circumboreal_Dissolve.geojson'

# data_table_file <- '~/Downloads/dps_outputs/atl08_006_030m_2020_2020_06_09_filt_covars_merge_neighbors_004104.csv'
# topo_stack_file <- '~/Downloads/dps_outputs/CopernicusGLO30_4104_cog_topo_stack.tif'
# LC_mask_file <- '~/Downloads/dps_outputs/esa_worldcover_v100_2020_4104_cog.tif'
# l8_stack_file <- '~/Downloads/dps_outputs/HLS_4104_07-01_08-31_2020_2020.tif'
# data_sample_file <- '~/Downloads/dps_outputs/boreal_train_data_2020_n3.csv'
# boreal_vect <- '~/Downloads/dps_outputs/wwf_circumboreal_Dissolve.geojson'

DO_MASK_WITH_STACK_VARS <- 'TRUE'
#data_sample_file <- '/projects/my-private-bucket/boreal_train_data_v11.csv'
iters <- 30
ppside <- 1
minDOY <- 130
maxDOY <- 250
max_sol_el <- 5
expand_training <- 'TRUE'
local_train_perc <- 100
min_icesat2_samples <- 5000
max_n <- 10000

#boreal_vect <- '/projects/shared-buckets/nathanmthomas/boreal_tiles_v003.gpkg'
predict_var <- 'AGB'
#predict_var <- 'Ht'

ppside <- as.double(ppside)
minDOY <- as.double(minDOY)
maxDOY <- as.double(maxDOY)
max_sol_el <- as.double(max_sol_el)
local_train_perc <- as.double(local_train_perc)

MASK_LYR_NAMES = c('slopemask', 'ValidMask')

#MASK_LANDCOVER_NAMES = c(0,13,15,16)
MASK_LANDCOVER_NAMES = c(50,70,80,100)

print(paste0("Do mask? ", DO_MASK_WITH_STACK_VARS))

# loading packages and functions
#----------------------------------------------#
library(randomForest)
#library(rgdal)
library(data.table)
library(ggplot2)
library(dplyr)
library(rlist)
library(fs)
library(stringr)
#library(gdalUtils)
library(rockchalk)
library(terra)
# run code
# adding model ids
rds_models <- list.files(path='~/dps_output/', pattern='*.rds', full.names = TRUE)
models_id<-names(rds_models)<-paste0("m",1:length(rds_models))

# make sure data are linked properly
#check extents
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

stack_poly <- resample_or_reproject_inputs()
stack <- stack_poly[['stack']]
boreal_poly <- stack_poly[['boreal_poly']]


if(DO_MASK_WITH_STACK_VARS){
    print("Masking stack...")
    # Bricking the stack will make the masking faster (i think)
    #brick = rast(stack)
    for(LYR_NAME in MASK_LYR_NAMES){
        m <- terra::subset(stack, grep(LYR_NAME, names(stack), value = T))

        stack <- mask(stack, m == 0, maskvalue=TRUE)

    }
    for(LC_NAME in MASK_LANDCOVER_NAMES){
        n <- terra::subset(stack, grep('esa_worldcover_v100_2020', names(stack), value=LC_NAME))
        stack <- mask(stack, n == LC_NAME, maskvalue=TRUE)

    }
    rm(m)
}


print("modelling begins")

print('file name:')
print(data_sample_file)
set.seed(123)
NTREE = 30
maps<-mapBoreal(rds_models=rds_models,
                models_id=models_id,
                ice2_30_atl08_path=data_table_file,
                ice2_30_sample=data_sample_file,
                offset=100.0,
                s_train=70,
                rep=iters,
                ppside=ppside,
                stack=stack,
                strat_random=FALSE,
                output=out_fn,
                minDOY=minDOY,
                maxDOY=maxDOY,
                max_sol_el=max_sol_el,
                expand_training=expand_training,
                local_train_perc=local_train_perc,
                min_icesat2_samples=min_icesat2_samples,
                DO_MASK=DO_MASK_WITH_STACK_VARS,
                boreal_poly=boreal_poly,
                predict_var=predict_var,
                max_n=max_n,
                pred_vars=pred_vars)
