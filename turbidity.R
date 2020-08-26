# required packages
necessary.packages<-c("tidyverse", "here", "raster", "rasterVis", "rgdal", "rgeos", "rgdal", "rgeos", "maptools", "dplyr", "ggplot2", "sf", "bnspatial", "Metrics", "waterquality")
already.installed<-necessary.packages%in%installed.packages()[, 'Package']    #asks if the necessary packages are already installed in the library?
if (length(necessary.packages[!already.installed])>=1) {    #if any are NOT installed, download them now.
  install.packages(necessary.packages[!already.installed],dep=T)   #are the dependencies really necessary (there are lots!)?
}
sapply(necessary.packages,function(p) {require(p,quietly=T,character.only=T)})

#set working file
setwd("E:/Gemma/DataProcessing")


cefas_data <- dir("../CefasData/Cefas_acoustic_2016", recursive=TRUE, full.names=TRUE, pattern="\\.tiff$")
cefas_filename <- str_match(cefas_data,".*/((NNSB_)?(\\d{4})(\\d{2})(\\d{2}).*)\\.tiff$")
foreach(x=1:20) %do% {
  cefas_filename[x] <- paste0('cefas.',cefas_filename[x,2], collapse="\\")
  print(cefas_filename[x])
}
cefas_filename <-cefas_filename[,1]

#sentinel-2 data that matches the NR_SA selection

file_2015 <- stack("../SentinelData/Acolite/2015/28-08-2015/31UCV/S2A_rhow_acolite_L3_T31UCV_0d_20150828_110046.tif")
file_J2016 <- stack("../SentinelData/Acolite/2016/03-07-2016/31UCV/S2A_rhow_acolite_L3_T31UCV_0d_20160703_105622.tif")
file_A2016 <- stack("../SentinelData/Acolite/2016/12-08-2016/31UCV/S2A_rhow_acolite_L3_T31UCV_0d_20160812_105622.tif")
file_2017 <- stack("../SentinelData/Acolite/2017/18-06-2017/31UCV/S2A_rhow_acolite_L3_T31UCV_0d_20170618_105621.tif")
file_J2019 <- stack("../SentinelData/Acolite/2019/23-07-2019/31UCV/S2B_rhow_acolite_L3_T31UCV_0d_20190723_105629.tif")
file_A2019 <- stack("../SentinelData/Acolite/2019/25-08-2019/31UCV/S2B_rhow_acolite_L3_T31UCV_0d_20190825_110629.tif")
file_S2019 <- stack("../SentinelData/Acolite/2019/21-09-2019/31UCV/S2B_rhow_acolite_L3_T31UCV_0d_20190921_105739.tif")
file_c <- raster(cefas_data[[2]])
sand_ID <- 'NR_SA'

#turbidity NIR over green
coastal_green2015 <- wq_calc(raster_stack = file_2015, alg = 'TurbChip09NIROverGreen', sat = 'sentinel2')
coastal_greenJ2016 <- wq_calc(raster_stack = file_J2016, alg = 'TurbChip09NIROverGreen', sat = 'sentinel2')
coastal_greenA2016 <- wq_calc(raster_stack = file_A2016, alg = 'TurbChip09NIROverGreen', sat = 'sentinel2')
coastal_green2017 <- wq_calc(raster_stack = file_2017, alg = 'TurbChip09NIROverGreen', sat = 'sentinel2')
coastal_greenJ2019 <- wq_calc(raster_stack = file_J2019, alg = 'TurbChip09NIROverGreen', sat = 'sentinel2')
coastal_greenA2019 <- wq_calc(raster_stack = file_A2019, alg = 'TurbChip09NIROverGreen', sat = 'sentinel2')
coastal_greenS2019 <- wq_calc(raster_stack = file_S2019, alg = 'TurbChip09NIROverGreen', sat = 'sentinel2')

#turbidity NIR over red
coastal_red2015 <- wq_calc(raster_stack = file_2015, alg = 'TurbDox02NIRoverRed', sat = 'sentinel2')
coastal_redJ2016 <- wq_calc(raster_stack = file_J2016, alg = 'TurbDox02NIRoverRed', sat = 'sentinel2')
coastal_redA2016 <- wq_calc(raster_stack = file_A2016, alg = 'TurbDox02NIRoverRed', sat = 'sentinel2')
coastal_red2017 <- wq_calc(raster_stack = file_2017, alg = 'TurbDox02NIRoverRed', sat = 'sentinel2')
coastal_redJ2019 <- wq_calc(raster_stack = file_J2019, alg = 'TurbDox02NIRoverRed', sat = 'sentinel2')
coastal_redA2019 <- wq_calc(raster_stack = file_A2019, alg = 'TurbDox02NIRoverRed', sat = 'sentinel2')
coastal_redS2019 <- wq_calc(raster_stack = file_S2019, alg = 'TurbDox02NIRoverRed', sat = 'sentinel2')

#work with Green and NIR turbidity data
crop1_2015 <- crop(coastal_green2015, file_c)
crop1_J2016 <- crop(coastal_greenJ2016, file_c)
crop1_A2016 <- crop(coastal_greenA2016, file_c)
crop1_2017 <- crop(coastal_green2017, file_c)
crop1_J2019 <- crop(coastal_greenJ2019, file_c)
crop1_A2019 <- crop(coastal_greenA2019, file_c)
crop1_S2019 <- crop(coastal_greenS2019, file_c)

#create graphical representation of turbidity
png('turbiditygraphs/G_turb_2015_NRSA.png')
plot(crop1_2015)
dev.off()

png('turbiditygraphs/G_turb_J2016_NRSA.png')
plot(crop1_J2016)
dev.off()

png('turbiditygraphs/G_turb_A2016_NRSA.png')
plot(crop1_A2016)
dev.off()

png('turbiditygraphs/G_turb_2017_NRSA.png')
plot(crop1_2017)
dev.off()

png('turbiditygraphs/G_turb_J2019_NRSA.png')
plot(crop1_J2019)
dev.off()

png('turbiditygraphs/G_turb_A2019_NRSA.png')
plot(crop1_A2019)
dev.off()

png('turbiditygraphs/G_turb_S2019_NRSA.png')
plot(crop1_S2019)
dev.off()

#work with Red and NIR turbidity data
crop2_2015 <- crop(coastal_red2015, file_c)
crop2_J2016 <- crop(coastal_redJ2016, file_c)
crop2_A2016 <- crop(coastal_redA2016, file_c)
crop2_2017 <- crop(coastal_red2017, file_c)
crop2_J2019 <- crop(coastal_redJ2019, file_c)
crop2_A2019 <- crop(coastal_redA2019, file_c)
crop2_S2019 <- crop(coastal_redS2019, file_c)

#create graphical representation of turbidity
png('turbiditygraphs/R_turb_2015_NRSA.png')
plot(crop2_2015)
dev.off()

png('turbiditygraphs/R_turb_J2016_NRSA.png')
plot(crop2_J2016)
dev.off()

png('turbiditygraphs/R_turb_A2016_NRSA.png')
plot(crop2_A2016)
dev.off()

png('turbiditygraphs/R_turb_2017_NRSA.png')
plot(crop2_2017)
dev.off()

png('turbiditygraphs/R_turb_J2019_NRSA.png')
plot(crop2_J2019)
dev.off()

png('turbiditygraphs/R_turb_A2019_NRSA.png')
plot(crop2_A2019)
dev.off()

png('turbiditygraphs/R_turb_S2019_NRSA.png')
plot(crop2_S2019)
dev.off()

