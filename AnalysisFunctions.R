# required packages
necessary.packages<-c("tidyverse", "here", "raster", "rasterVis", "rgdal", "rgeos", "maptools", "dplyr", "ggplot2", "sf", "foreach", "spatialEco", "Metrics")
already.installed<-necessary.packages%in%installed.packages()[, 'Package']    #asks if the necessary packages are already installed in the library?
if (length(necessary.packages[!already.installed])>=1) {    #if any are NOT installed, download them now.
  install.packages(necessary.packages[!already.installed],dep=T)   #are the dependencies really necessary (there are lots!)?
}
sapply(necessary.packages,function(p) {require(p,quietly=T,character.only=T)})

#set working file
setwd("E:/Gemma/DataProcessing")

#allow identification of sandbanks by corresponding satellite tile
cefas_mapping <- list()
cefas_mapping[["31UCV"]] <- list(
  list("NR_SA", "../CefasData/Cefas_acoustic_2016/20160608_NR_SA_1d0.tiff"),
  list("ID_CSA", "../CefasData/Cefas_acoustic_2016/20160606_ID_CSA_1d0.tiff")
)
cefas_mapping[["31UDV"]] <- list(
  list("IDFB", "../CefasData/Cefas_acoustic_2016/NNSB_20160601_MBFP_1d0_UTM31N_IDFBcubeB.tiff"),
  list("WCT001B", "../CefasData/Cefas_acoustic_2016/NNSB_20160603_MBFP_1d0_UTM31N_WCT001B.tiff"),
  list("WCT002B", "../CefasData/Cefas_acoustic_2016/NNSB_20160603_MBFP_1d0_UTM31N_WCT002B.tiff"),
  list("WCT003B", "../CefasData/Cefas_acoustic_2016/NNSB_20160603_MBFP_1d0_UTM31N_WCT003B.tiff"),
  list("LMBKBcross", "../CefasData/Cefas_acoustic_2016/NNSB_20160613_MBFP_1d0_UTM31N__LMBKBcross.tiff"),
  list("LMBKB", "../CefasData/Cefas_acoustic_2016/NNSB_20160613_MBFP_1d0_UTM31N_LMBKB.tiff"),
  list("WCT004B", "../CefasData/Cefas_acoustic_2016/NNSB_20160615_MBFP_1d0_UTM31N_WCT004B.tiff"),
  list("WCT005B", "../CefasData/Cefas_acoustic_2016/NNSB_20160615_MBFP_1d0_UTM31N_WCT005B.tiff"),
  list("WCT009B", "../CefasData/Cefas_acoustic_2016/NNSB_20160615_MBFP_1d0_UTM31N_WCT009B.tiff"),
  list("WCT0010B", "../CefasData/Cefas_acoustic_2016/NNSB_20160615_MBFP_1d0_UTM31N_WCT0010B.tiff"),
  list("WCT0011B", "../CefasData/Cefas_acoustic_2016/NNSB_20160615_MBFP_1d0_UTM31N_WCT011B_NEW.tiff"),
  list("WMCSB", "../CefasData/Cefas_acoustic_2016/NNSB_20160617_MBFP_1d0_UTM31N__WMCSB.tiff"),
  list("WCT006B", "../CefasData/Cefas_acoustic_2016/NNSB_20160615_MBFP_1d0_UTM31N_WCT006B.tiff")
)
cefas_mapping[["31UDU"]] <- list(
  list("HHW_SA", "../CefasData/Cefas_acoustic_2016/20160613_HHW_SA_1d0.tiff"),
  list("SWH_SA", "../CefasData/Cefas_acoustic_2016/20160614_SWH_SA_1d0.tiff"),
  list("HASA_SA", "../CefasData/Cefas_acoustic_2016/20160615_HHW_HASA_SA_1d0.tiff"),
  list("WCT006", "../CefasData/Cefas_acoustic_2016/NNSB_20160615_MBFP_1d0_UTM31N_WCT006B.tiff"),
  list("WCT007B", "../CefasData/Cefas_acoustic_2016/NNSB_20160615_MBFP_1d0_UTM31N_WCT007B.tiff"),
  list("WCT008B", "../CefasData/Cefas_acoustic_2016/NNSB_20160615_MBFP_1d0_UTM31N_WCT008B.tiff")
)

#use as a constant reference for remapping ARCSI data.
crs_cefas <- crs(raster(cefas_mapping[['31UCV']][[1]][[2]]))

#unique element of the processing for ACOLITE
do_acolite <- function(date, tilename) {
  #allows for auotmatic file retrieval and confirms files presence
  year <- substr(date, nchar(date)-3, nchar(date))
  input_dirname <- paste0("../SentinelData/Acolite/", year, "/", date, "/", tilename, "/",  collapse = "\\")
  tiff_files_in_input_dir <- dir(input_dirname, recursive=TRUE, full.names = TRUE, pattern="\\.tif$")
  if(length(tiff_files_in_input_dir) == 0)
    return()
  stopifnot(length(tiff_files_in_input_dir) == 1)
  input_filename <- tiff_files_in_input_dir[1]
  print(paste0('acolite ', tilename, ' ', date))
  #load acolite data and extract required bands
  raw_acolite <- brick(input_filename)

  coastal_layer <- raster(raw_acolite, layer = 'rhow0442')
  green_layer <- raster(raw_acolite, layer = 'rhow0559')
  
  #creates log values and generates SDB for whole tile.
  lncoast <- log(coastal_layer)
  lngreen <- log(green_layer)
  sat_depth <- lncoast/lngreen

  print(paste0(tilename, ' stage 1 done'))
  produce_results_csv(date = date, analysis_name = 'ACOLITE', tilename = tilename, sat_depth = sat_depth)
}
#unique element of processing for ARCSI
do_arcsi <- function(date, tilename) {
  #allows for auotmatic file retrieval and confirms files presence
  year <- substr(date, nchar(date)-3, nchar(date))
  input_dirname <- paste0("../SentinelData/Arcsi/", year, "/", date, "/", tilename, "/",  collapse = "\\")
  tiff_files_in_input_dir <- dir(input_dirname, recursive=TRUE, full.names = TRUE, pattern="\\.tif$")
  if(length(tiff_files_in_input_dir) == 0)
    return()
  stopifnot(length(tiff_files_in_input_dir) == 2)
  input_file <- tiff_files_in_input_dir[!endsWith(tiff_files_in_input_dir, "clouds.tif")]
  cloud_file <- tiff_files_in_input_dir[endsWith(tiff_files_in_input_dir, "clouds.tif")]

  stopifnot(length(input_file) == 1)
  stopifnot(length(cloud_file) == 1)
  
  print(paste0('arcsi ', tilename, date))
  #reads in ARCSI data
  raw_arcsi <- brick(input_file)
  cloud_layer <- raster(cloud_file)
  
  print('reraster starting')
  #edits co-ordinate reference system as ARCSI is tmerc, not utm.
  raw_arcsi <- projectRaster(raw_arcsi, crs = crs_cefas)
  cloud_layer <- projectRaster(cloud_layer, crs = crs_cefas)

  print('reraster finished')
  #generates log values and sdb
  green_layer <- raster(raw_arci, layer = 'Green')
  green_layer <- mask(green_layer, cloud_layer, inverse=TRUE, maskvalue = NA)
  blue_layer <- raster(raw_arcsi, layer = 'Blue')
  blue_layer <- mask(blue_layer, cloud_layer, inverse=TRUE, maskvalue = NA)

  lngreen <- log(green_layer)
  lnblue <- log(blue_layer)
  sat_depth <- lngreen/lnblue

  rm(raw_arcsi, cloud_later, green_layer, blue_layer, lngreen, lnblue)
  
  produce_results_csv(date = date, analysis_name = 'ARCSI', tilename = tilename, sat_depth = sat_depth)
}

produce_results_csv <- function(date = date, analysis_name = analysis_name,  tilename = tilename, sat_depth = sat_depth) {
#reads in Cefas MBES data
  for (cefas_info in cefas_mapping[[tilename]]) {
    sand_ID <- cefas_info[[1]]
    print(sand_ID)
    file_cefas <- cefas_info[[2]]
    cefas_filename <- basename(file_cefas)
    cefas_2016 <- raster(file_cefas)

    #crop - use an initially larger extent, and then crop smaller to ensure extents are statistically similar
    #extracts data from sdb which matches Cefas MBES   
    print('cropping')
    shape <- trim(cefas_2016, padding = 20)
    print('trimmed cefas')
    shape <- rasterToPolygons(shape)
    print('converted to polygons')
    filepath <- paste0('../CefasData/ShapeFile/multistep_', sand_ID, collapse="\\")
    raster::shapefile(shape,filepath, overwrite=TRUE) # be aware of size of file - can cause failures if to large
    rm(shape)
    filepath <- paste0('../CefasData/ShapeFile/multistep_', sand_ID, '.shp', collapse="\\")
    multi_step <- readOGR(filepath)
    crop1 <- crop(x = sat_depth, y = multi_step)
    shape2 <- rasterToPolygons(crop1)
    filepath <- paste0('../CefasData/ShapeFile/multistep_', sand_ID, '_2', collapse="\\")
    raster::shapefile(shape2, filepath, overwrite=TRUE) # be aware of file size
    rm(shape2)
    filepath <- paste0(filepath, '.shp', collapse = "\\")
    step_2 <- readOGR(filepath)
    crop2 <- crop(x = cefas_2016, y = step_2)
    crop2.1 <- crop(x = sat_depth, y = step_2)

    #create zonal statistics for the cefas data
    filepath <- paste0('../CefasData/Statistics/', sand_ID,'/', cefas_filename, '_', date,  '_', analysis_name, '_mean', collapse="\\")
    cefas_mean <- aggregate(crop2, fact=10, fun=mean)
    writeRaster(cefas_mean, filepath, overwrite=TRUE)
    filepath <- paste0('../CefasData/Statistics/', sand_ID,'/',cefas_filename, '_', date,  '_', analysis_name, '_min', collapse="\\")
    cefas_min <- aggregate(crop2, fact=10, fun=min)
    writeRaster(cefas_min, filepath, overwrite=TRUE)
    filepath <- paste0('../CefasData/Statistics/', sand_ID,'/',cefas_filename, '_', date,  '_', analysis_name, '_max', collapse="\\")
    cefas_max <- aggregate(crop2, fact=10, fun=max)
    writeRaster(cefas_max, filepath, overwrite=TRUE)

    satellite_data <- crop2.1

    #create data frame from stats
    cefas_points <- rasterToPoints(cefas_mean, spatial=TRUE)
    file <- raster::extract(satellite_data, cefas_points)
    file <- cbind(cefas_points,file)
    file <- as.data.frame(file)
    file <- file[names(file)!="geometry"]
    file <- file[complete.cases(file),]

    #ensures naming is accurate for processing
    col_names <- colnames(file)
    names(file)[names(file) == col_names[1]] <- "Depth"
    names(file)[names(file) == col_names[2]] <- "Lg_Ratio"

    #save data sets for additional analysis
    filepath <- paste0('DataSets/', analysis_name, '_', date,'_', sand_ID, '_dataset.csv', collapse = "\\")
    write.csv(file, file=filepath)

    #aggregate dataset, using Cefas MBES to 1DP
    file[,'Depth']=round(file[,'Depth'],1)
    file <- aggregate(file$Lg_Ratio, by=list(Depth=file$Depth), FUN=mean)

    #save aggregated data sets for additional analysis 
    filepath_agg <- paste0('DataSets/', analysis_name, '_', date, '_', sand_ID, '_dataset_aggregated.csv', collapse = "\\")
    write.csv(file, file=filepath_agg)
    
    do_complete_graphs(date = date, tilename = tilename, analysis = analysis_name, sand_ID = sand_ID)
    do_aggregated_graphs(date = date, tilename = tilename, analysis = analysis_name, sand_ID = sand_ID)
    
  }

}

do_complete_graphs <- function(date, tilename, analysis_name, sand_ID) {
  #creates graphic analysis
  analysis <- analysis_name
  filename <- paste0('./DataSets/', analysis, '_', date, '_', sand_ID, '_dataset.csv', collapse="\\")
  print('Complete Graphs')
  file <- read.csv(filename)
  #confirms names
  col_names <- colnames(file)
  names(file)[names(file) == col_names[2]] <- "Depth"
  names(file)[names(file) == col_names[3]] <- "Lg_Ratio"
  #autogenerates key required data
  title <- paste0(sand_ID, ' - ', date, " Sandbank - ", analysis, " Data")
  xmin <- (min(file$Lg_Ratio))
  xmax <- (max(file$Lg_Ratio))
  if (analysis == "ACOLITE") {
    xlabel <- "Lg(Coastal Aerosol) / Lg(Green)" }
    else { (analysis == "ARCSI")
    xlabel <- "Lg(Green) / Lg(Blue)"
  }

  ymin <- (min(file$Depth))
  ymax <- (max(file$Depth))
  fit <- lm(file$Depth~file$Lg_Ratio)
  ggplot(file, aes(x=Lg_Ratio, y=Depth)) + geom_point() + xlim(c(xmin, xmax)) +
  scale_y_continuous(c(ymin, ymax)) +  geom_smooth(method=lm, col="red") +
  ggtitle(title) + theme(plot.title = element_text(hjust=.5)) + labs(x=xlabel, y = "Depth (m)") +
  theme(panel.border=element_rect(colour="black", fill=NA), panel.background = element_blank(),
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12)) +
  geom_label(aes(x=xmin+0.05, y=ymin+3), hjust = 0, label = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared,4),
                                                         "\nIntercept = ",signif(fit$coef[[1]],4),
                                                         "\nSlope =",signif(fit$coef[[2]],4)))

  filepath <- paste0('./Graphs/', analysis, '_', sand_ID, '_', date, '_full.png')
  ggsave(filepath)

  do_statistics(date = date, tilename = tilename, analysis = analysis_name,  sand_ID = sand_ID, file = file)
}


do_aggregated_graphs <- function(date, tilename, analysis_name, sand_ID) {
  #creates graphical analysis
  analysis <- analysis_name
  filename <- paste0('./DataSets/', analysis, '_', date, '_', sand_ID, '_dataset_aggregated.csv', collapse="\\")
  print('Aggregated Graphs')
  file <- read.csv(filename)
  #confirms column names
  col_names <- colnames(file)
  names(file)[names(file) == col_names[2]] <- "Depth"
  names(file)[names(file) == col_names[3]] <- "Lg_Ratio"
  #autogenerates all key information
  title <- paste0(sand_ID, ' - ', date, " Sandbank - Aggregated", analysis, " Data")
  xmin <- (min(file$Lg_Ratio))
  xmax <- (max(file$Lg_Ratio))
  ymin <- (min(file$Depth))
  ymax <- (max(file$Depth))
   if (analysis == "ACOLITE") {
    xlabel <- "Lg(Coastal Aerosol) / Lg(Green)" }
    else { (analysis == "ARCSI")
    xlabel <- "Lg(Green) / Lg(Blue)"
  }
  fit <- lm(file$Depth~file$Lg_Ratio)
  ggplot(file, aes(x=Lg_Ratio, y=Depth)) + geom_point() + xlim(c(xmin, xmax)) +
  scale_y_continuous(c(ymin, ymax)) +  geom_smooth(method=lm, col="red") +
  ggtitle(title) + theme(plot.title = element_text(hjust=.5)) + labs(x=xlabel, y = "Depth (m)") +
  theme(panel.border=element_rect(colour="black", fill=NA), panel.background = element_blank(),
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12)) +
  geom_label(aes(x=xmin+0.05, y=ymin+3), hjust = 0, label = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared,4),
                                                         "\nIntercept = ",signif(fit$coef[[1]],4),
                                                         "\nSlope =",signif(fit$coef[[2]],4)))

  filepath <- paste0('./Graphs/', analysis, '_', sand_ID, '_', date, '_aggregated.png')
  ggsave(filepath)

  do_statistics(date = date, tilename = tilename, analysis = analysis, sand_ID = sand_ID, file = file)
}

do_statistics <- function(date, tilename, analysis, sand_ID, file) {
  results <- read.csv('./DataSets/results.csv')
  #generates linear regression data
  fit <- lm(file$Depth~file$Lg_Ratio)
  #extracts R^2 value
  R <- summary(fit)$adj.r.squared
  #allows for data to be extracted and written externally for future uses.
  fit <- fit$coef
  equation <- (fit[2]*(file$Lg_Ratio)) + fit[1]
  print('equation')
  RMSE <- rmse(file$Depth, equation)
  MAE <- mae(file$Depth, equation)
  fit[1]=round(fit[1],2)
  fit[2]=round(fit[2],2)
  fit_eq <- paste0('y = ', fit[2], 'x', fit[1] )
  x <- paste0('analysis', ',', date, ',', tilename, ',', sand_ID,',', fit_eq,',', RMSE,',', MAE)
  x <- write.table(x, 'results.csv', append=TRUE)
  write.csv2('./DataSets/results.csv')
}

process_everything_for_date_and_tile <- function(date, tilename) {
  do_acolite(date, tilename)
  #do_arcsi(date, tilename)
}

process_everything_for_date <- function(date) {
  process_everything_for_date_and_tile(date, "31UCV")
  #process_everything_for_date_and_tile(date, "31UDU")
  #process_everything_for_date_and_tile(date, "31UDV")
}

main <- function() {
  process_everything_for_date("28-08-2015")
  process_everything_for_date("03-07-2016")
  process_everything_for_date("12-08-2016")
  process_everything_for_date("18-06-2017")
  process_everything_for_date("28-06-2018")
  process_everything_for_date("23-07-2019")
  process_everything_for_date("27-08-2019")
  process_everything_for_date("21-09-2019")
}

main()
