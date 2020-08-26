# required packages
necessary.packages<-c("tidyverse", "here", "raster", "rasterVis", "rgdal", "rgeos", "rgdal", "rgeos", "maptools", "dplyr", "ggplot2", "sf", "bnspatial", "foreach", "waterquality", "Metrics")
already.installed<-necessary.packages%in%installed.packages()[, 'Package']    #asks if the necessary packages are already installed in the library?
if (length(necessary.packages[!already.installed])>=1) {    #if any are NOT installed, download them now.
  install.packages(necessary.packages[!already.installed],dep=T)   #are the dependencies really necessary (there are lots!)?
}
sapply(necessary.packages,function(p) {require(p,quietly=T,character.only=T)})

#set working file
setwd("E:/Gemma/DataProcessing")

ID_CSA <- read.csv("./DataSets/ARCSI_28-08-2015_ID_CSA_dataset.csv")
NR_SA <- read.csv("./DataSets/ARCSI_28-08-2015_NR_SA_dataset.csv")

ID_CSA <- ID_CSA[,2:3]
NR_SA <- NR_SA[,2:3]
date <- ''

tile31UCV <- rbind(NR_SA, ID_CSA)
print(tile31UCV)

tile31UCV[,'Depth']=round(tile31UCV[,'Depth'],1)
tile31UCV <- aggregate(tile31UCV$Lg_Ratio, by=list(Depth=tile31UCV$Depth), FUN=mean)

file <- tile31UCV

col_names <- colnames(file)
names(file)[names(file) == col_names[1]] <- "Depth"
names(file)[names(file) == col_names[2]] <- "Lg_Ratio"

title <- paste0('31UCV - 28-08-2015 ARCSI Data')
xmin <- (min(file$Lg_Ratio))
xmax <- (max(file$Lg_Ratio))
ymin <- (min(file$Depth))
ymax <- (max(file$Depth))
ylabel <- 'Depth (m)'
xlabel <- "Lg(Green) / Lg(Blue)"

fit <- lm(tile31UCV$Depth~tile31UCV$x)
ggplot(tile31UCV, aes(x=x, y=Depth)) + geom_point() + xlim(c(xmin, xmax)) +
  scale_y_continuous(ylabel) +  geom_smooth(method=lm, col="red") +
  ggtitle(title) + theme(plot.title = element_text(hjust=.5)) + labs(x="Ln(Green) / Ln(Blue)", y = "Depth (m)") +
  theme(panel.border=element_rect(colour="black", fill=NA), panel.background = element_blank(),
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12)) +
  geom_label(aes(x=xmin+0.05, y=ymin+3), hjust = 0, label = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared,4),
                                                                  "\nIntercept = ",signif(fit$coef[[1]],4),
                                                                  "\nSlope =",signif(fit$coef[[2]],4)))

fit <- lm(tile31UCV$Depth~tile31UCV$x)
r <- (summary(fit)$adj.r.squared)
fit <- fit$coef
equation <- (fit[2]*(tile31UCV$x)) + fit[1]
print('equation')
RMSE <- rmse(tile31UCV$Depth, equation)
MAE <- mae(tile31UCV$Depth, equation)
fit[1]=round(fit[1],2)
fit[2]=round(fit[2],2)
fit_eq <- paste0('y = ', fit[2], 'x', fit[1] )

tile31UCV_surface <- tile31UCV[tile31UCV$Depth > -12,]

col_names <- colnames(tile31UCV_surface)
names(tile31UCV_surface)[names(tile31UCV_surface) == col_names[1]] <- "Depth"
names(tile31UCV_surface)[names(tile31UCV_surface) == col_names[2]] <- "Lg_Ratio"

title <- paste0('31UCV 0 - 12m below sealevel - 28-08-2015 ARCSI Data')
xmin <- (min(tile31UCV_surface$Lg_Ratio))
xmax <- (max(tile31UCV_surface$Lg_Ratio))
ymin <- (min(tile31UCV_surface$Depth))
ymax <- (max(tile31UCV_surface$Depth))

fit=lm(tile31UCV_surface$Depth~tile31UCV_surface$Lg_Ratio)
ggplot(tile31UCV_surface, aes(x=Lg_Ratio, y=Depth)) + geom_point() + xlim(c(xmin, xmax)) +
  scale_y_continuous(ylabel) +  geom_smooth(method=lm, col="red") +
  ggtitle(title) + theme(plot.title = element_text(hjust=.5)) + labs(x=xlabel, y=ylabel) +
  theme(panel.border=element_rect(colour="black", fill=NA), panel.background = element_blank(),
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12)) +
  geom_label(aes(x=xmin+0.05, y=ymin+3), hjust = 0, label = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared,4),
                                                                  "\nIntercept = ",signif(fit$coef[[1]],4),
                                                                  "\nSlope =",signif(fit$coef[[2]],4)))

fit <- lm(tile31UCV_surface$Depth~tile31UCV_surface$Lg_Ratio)
r <- summary(fit)$adj.r.squared
fit <- fit$coef
equation <- (fit[2]*(tile31UCV_surface$Lg_Ratio)) + fit[1]
print('equation')
RMSE <- rmse(tile31UCV_surface$Depth, equation)
MAE <- mae(tile31UCV_surface$Depth, equation)
fit[1]=round(fit[1],2)
fit[2]=round(fit[2],2)
fit_eq <- paste0('y = ', fit[2], 'x', fit[1] )

tile31UCV_mid <- tile31UCV[tile31UCV$Depth > -20,]

col_names <- colnames(tile31UCV_mid)
names(tile31UCV_mid)[names(tile31UCV_mid) == col_names[1]] <- "Depth"
names(tile31UCV_mid)[names(tile31UCV_mid) == col_names[2]] <- "Lg_Ratio"

title <- paste0('31UCV 0 - 20m below sealevel - 28-08-2015 ARCSI Data')
xmin <- (min(tile31UCV_mid$Lg_Ratio))
xmax <- (max(tile31UCV_mid$Lg_Ratio))
ymin <- (min(tile31UCV_mid$Depth))
ymax <- (max(tile31UCV_mid$Depth))

fit=lm(tile31UCV_mid$Depth~tile31UCV_mid$Lg_Ratio)
ggplot(tile31UCV_mid, aes(x=Lg_Ratio, y=Depth)) + geom_point() + xlim(c(xmin, xmax)) +
  scale_y_continuous(ylabel) +  geom_smooth(method=lm, col="red") +
  ggtitle(title) + theme(plot.title = element_text(hjust=.5)) + labs(x=xlabel, y=ylabel) +
  theme(panel.border=element_rect(colour="black", fill=NA), panel.background = element_blank(),
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12)) +
  geom_label(aes(x=xmin+0.05, y=ymin+3), hjust = 0, label = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared,4),
                                                                  "\nIntercept = ",signif(fit$coef[[1]],4),
                                                                  "\nSlope =",signif(fit$coef[[2]],4)))


col_names <- colnames(tile31UCV_mid)
names(tile31UCV_mid)[names(tile31UCV_mid) == col_names[1]] <- "Depth"
names(tile31UCV_mid)[names(tile31UCV_mid) == col_names[2]] <- "Lg_Ratio"

fit <- lm(tile31UCV_mid$Depth~tile31UCV_mid$Lg_Ratio)
r <- summary(fit)$adj.r.squared
fit <- fit$coef
equation <- (fit[2]*(tile31UCV_mid$Lg_Ratio)) + fit[1]
print('equation')
RMSE <- rmse(tile31UCV_mid$Depth, equation)
MAE <- mae(tile31UCV_mid$Depth, equation)
fit[1]=round(fit[1],2)
fit[2]=round(fit[2],2)
fit_eq <- paste0('y = ', fit[2], 'x', fit[1] )
