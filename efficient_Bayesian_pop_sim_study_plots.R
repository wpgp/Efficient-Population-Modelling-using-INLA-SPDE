###---BAYESIAN GEOSTATISTICAL POP MODEL PAPER SIMULATION STUCY CODES
###---NNANATU ET AL.(2024)
### - Simulation Study 

rm(list=ls())
library(INLA); library(raster); library(maptools)
library(gtools); library(sp); library(spdep)
library(fields); library(mvtnorm); library(gtools)
library(geoR);library(actuar);library(viridisLite)
require(grid);require(gridExtra);require(lattice);require(tidyverse)
library(dplyr); library(sf); library(tmap); library(tmaptools)

# Specify key data paths
path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3/Working/CMR/Chris_N/paper1"
shp_path <- paste0(path, "/application/data/Input_Settlement_Boundaries")
out_path <- paste0(path, "/sim_study/output")

cover <- c(1, 0.8, 0.6, 0.4, 0.2, 0.05) #--percentage coverage of the observations
spat.autocor <- c(0.01, 0.1, 1)#spatial autocorrelation variance

#source_re <- c(0.01, 0.1, 1)#variance of data source random effects
source_var <- list(same = rep(0.01, 5), 
                   diff = c(0.01, 0.04, 0.08, 0.2, 0.7))#variance of data source random effects
dat_type <- c("same", "diff")

#-----------------------------------------------------------------------------------
# Extract posterior fit indices and predictions when data sources have same variability
#-----------------------------------------------------------------------------------
fitsb1 <- fitsb2 <- list()
metsb1 <- metsb2 <- list()
predsb1 <- predsb2 <-list()
for(b in 1:length(spat.autocor))
{ 
  
  result_pathb <- paste0(out_path,"/same/marginal_var_", spat.autocor[b])
  
  c=1 # same data source variability
  result_pathc <- paste0(result_pathb,"_source_var_", dat_type[c])
  
  for(j in 1:length(cover))###----Do the following for each coverage percentage 
  {
    
    result_path2 <- paste0(result_pathc,"_at_", cover[j]*100,"%_coverage")
    
    # fit model fit measures
    mod_fit <- read.csv(paste0(result_path2,"/modfit.csv"))
    mod_fit$cover <- cover[j]*100
    mod_fit$dat_type <- dat_type[c]
    mod_fit$spat_cor <- spat.autocor[b]
    fitsb1[[j]] <- mod_fit
    
    
    # Extract posterior predictions
    pred <- read.csv(paste0(result_path2,"/grid_data.csv"))
    pred$cover <- rep(cover[j]*100, nrow(pred))
    pred$dat_type <- rep(dat_type[c], nrow(pred))
    pred$spat_cor <- rep(spat.autocor[b], nrow(pred))
    predsb1[[j]] <- pred
    
    # Extract more model fit measures
    metsb <- read.csv(paste0(result_path2,"/fit_metrics.csv"))
    metsb$cover <- rep(cover[j]*100, nrow(metsb))
    metsb$dat_type <- rep(dat_type[c], nrow(metsb))
    metsb$spat_cor <- rep(spat.autocor[b], nrow(metsb))
    metsb1[[j]] <- metsb
    
  }
  fitsb2[[b]] <- fitsb1 
  predsb2[[b]] <- predsb1 
  metsb2[[b]] <- metsb1 
}

## model fit measures
unnest_fitsb <- unlist(fitsb2, recursive = FALSE)  #--unnest the list
str(unnest_fitsb)
dim(fitsb <- as.data.frame(matrix(unlist(do.call(rbind, unnest_fitsb)),
                                  nrow=4*18, ncol=7)))# 4 models over 18 combinations & 7 variables
head(fitsb)
colnames(fitsb) <- c("model",names(unnest_fitsb[[1]])[2:7])
head(fitsb)
apply(fitsb[,2:4], 2,class) # check the class of the values
fitsb[,2:4] <- apply(fitsb[,2:4], 2,as.numeric)# convert to numeric
fitsb


## model posterior predictions
unnest_pred <- unlist(predsb2, recursive = FALSE)  #--unnest the list
str(unnest_pred)
names(unnest_pred[[1]])
dim(predsb <- as.data.frame(matrix(unlist(do.call(rbind, unnest_pred)),
                                   nrow=5023*18, ncol=24)))# 5023 grid cells over 18 combinations & 25 variables
head(predsb)
colnames(predsb) <- names(unnest_pred[[1]])
head(predsb)

vars <- c("lon", "lat", paste0("x", 1:6), "dtyp_re1", "dens", "source_re",
          "mean_dens_hat", "mean_pop_hat", "lower_pop_hat", "upper_pop_hat",
          "sd_pop_hat", "cv_pop_hat")
apply(predsb[,vars], 2,class) # all character variables!
predsb[,vars] <- apply(predsb[,vars], 2,as.numeric)# convert to numeric
predsb

## fit metrics
unnest_met <- unlist(metsb2, recursive = FALSE)  #--unnest the list
str(unnest_met)
names(unnest_met[[1]])
dim(metsb <- as.data.frame(matrix(unlist(do.call(rbind, unnest_met)),
                                 nrow=4*18, ncol=5)))# 3 metrics over 18 combinations & 5 variables
##
head(metsb)
colnames(metsb) <- c("metric", "value",names(unnest_met[[1]])[3:5])
head(metsb)
class(metsb)
names(metsb)

# ---Combine the metrics
mod.mets <- metsb
mod.mets2 <- mod.mets %>% filter(!cover %in% c(5,100))
##----- CORR
dim(corr <- mod.mets2[mod.mets2$metric=="corr",])
corr$value <- as.numeric(corr$value)
library(ggpubr)

plot_corr <- ggline(corr, x = "spat_cor", y = "value",
                    error.plot = "value",
                    #facet.by = "cover",
                    panel.labs.font.x = list(size=14),
                    color = "cover",palette = "lancet",
                    point.size=1.5,
                    size=1.4)

rcorr <-  ggpar(plot_corr, xlab="Spatial Variance", ylab="CORR",
                legend = "top", legend.title = "Survey \n Coverage (%)",size=22,
                font.legend=c(14),
                font.label = list(size = 14, face = "bold", color ="red"),
                font.x = c(14),
                font.y = c(14),
                font.main=c(14),
                font.xtickslab =c(14),
                font.ytickslab =c(14),
                # orientation = "reverse",
                xtickslab.rt = 45, ytickslab.rt = 45)
rcorr 


##----- RMSE
dim(rmse <- mod.mets2[mod.mets2$metric=="RMSE",])
rmse$value <- as.numeric(rmse$value)

plot_rmse <- ggline(rmse, x = "spat_cor", y = "value",
                    error.plot = "value",
                    panel.labs.font.x = list(size=14),
                    color = "cover",palette = "lancet",
                    point.size=1.5,
                    #linetype = "pop_cover",
                    size=1.4)

rrmse <-  ggpar(plot_rmse, xlab="Spatial Variance", ylab="RMSE",
               legend = "top", legend.title = "Survey \n Coverage (%)",size=22,
               font.legend=c(14),
               font.label = list(size = 14, face = "bold", color ="red"),
               font.x = c(14),
               font.y = c(14),
               font.main=c(14),
               font.xtickslab =c(14),
               font.ytickslab =c(14),
               # orientation = "reverse",
               xtickslab.rt = 45, ytickslab.rt = 45)
rrmse 

##----- MAE
dim(mae <- mod.mets2[mod.mets2$metric=="MAE",])
mae$value <- as.numeric(mae$value)

plot_mae <- ggline(mae, x = "spat_cor", y = "value",
                    error.plot = "value",
                    #facet.by = "cover",
                    # panel.labs= list(dat_type=c("diff", "same")),
                    panel.labs.font.x = list(size=14),
                    color = "cover",palette = "lancet",
                    point.size=1.5,
                    #linetype = "pop_cover",
                    size=1.4)
rmae <-  ggpar(plot_mae, xlab="Spatial Variance", ylab="MAE",
                legend = "none", 
               font.label = list(size = 14, face = "bold", color ="red"),
               font.x = c(14),
               font.y = c(14),
               font.main=c(14),
               font.xtickslab =c(14),
               font.ytickslab =c(14),
                # orientation = "reverse",
                xtickslab.rt = 45, ytickslab.rt = 45)
rmae 


##----- BIAS
dim(bias <- mod.mets2[mod.mets2$metric=="BIAS",])
bias$value <- as.numeric(bias$value)

plot_bias <- ggline(bias, x = "spat_cor", y = "value",
                   error.plot = "value",
                   panel.labs.font.x = list(size=14),
                   color = "cover",palette = "lancet",
                   point.size=1.5,
                   #linetype = "pop_cover",
                   size=1.4)


rbias <-  ggpar(plot_bias, xlab="Spatial Variance", ylab="BIAS",
                legend = "none", 
                font.label = list(size = 14, face = "bold", color ="red"),
                font.x = c(14),
                font.y = c(14),
                font.main=c(14),
                font.xtickslab =c(14),
                font.ytickslab =c(14),
                # orientation = "reverse",
                xtickslab.rt = 45, ytickslab.rt = 45)
rbias 


ggarrange(rcorr,rrmse, rmae, rbias, 
          ncol = 2, nrow = 2)
#------------------------------------------------------------
#### --- Explore predictions
#----------------------------------------------------------

names(predsb)
vars2put <- c("lon", "lat", "pop", "dens",
              "source", "mean_pop_hat", "lower_pop_hat",
              "upper_pop_hat", "cv_pop_hat", "cover",
              "spat_cor")
mod.pred <- predsb[,vars2put]
mod.pred$cover2 <- factor(mod.pred$cover,
                             levels = c("100", "80", "60", "40", "20", "5"),
                             label = c("Survey Coverage: \n 100% ",
                                       "Survey Coverage: \n 80% ",
                                       "Survey Coverage: \n 60% ",
                                       "Survey Coverage: \n 40% ",
                                       "Survey Coverage: \n 20% ",
                                       "Survey Coverage: \n 5% "))
mod.pred[, vars2put[-c(5, 10,11)]] <- apply(mod.pred[, vars2put[-c(5, 10,11)]], 2, as.numeric)
dim(mod.pred100 <- mod.pred %>% filter(cover==80 & spat_cor ==0.01))

###

# Low
dim(mod.predlow <- mod.pred %>% filter(spat_cor ==0.01))
mod.predlow$pop <- round(mod.predlow$pop/1000) # per 1k
mod.predlow$mean_pop_hat <- round(mod.predlow$mean_pop_hat/1000) # per 1k

plot_low <- ggscatter(mod.predlow, x = "pop", y = "mean_pop_hat",
          add = "reg.line",                         # Add regression line
          facet.by = "cover2",
          conf.int = TRUE,                          # Add confidence interval
          color = "source", 
          palette = "lancet",           # Color by groups "cyl"
          shape = "source"                             # Change point shape by groups "cyl"
) 

rlow <-  ggpar(plot_low, xlab="Simulated Counts ('000)", ylab="Predicted Counts ('000)",
                legend = "right", 
                legend.title = "Data \n Source",size=22,
                font.legend=c(12),
                font.label = list(size = 12, face = "bold", color ="red"),
                font.x = c(12),
                font.y = c(12),
                font.main=c(14),
                font.xtickslab =c(10),
                font.ytickslab =c(10),
                # orientation = "reverse",
                xtickslab.rt = 45, ytickslab.rt = 45)
rlow 

# Medium
dim(mod.predmed <- mod.pred %>% filter(spat_cor ==0.1))
mod.predmed$pop <- round(mod.predmed$pop/1000) # per 1k
mod.predmed$mean_pop_hat <- round(mod.predmed$mean_pop_hat/1000) # per 1k

plot_med <- ggscatter(mod.predmed, x = "pop", y = "mean_pop_hat",
                      add = "reg.line",                         # Add regression line
                      facet.by = "cover2",
                      conf.int = TRUE,                          # Add confidence interval
                      color = "source", 
                      palette = "lancet",           # Color by groups "cyl"
                      shape = "source"                             # Change point shape by groups "cyl"
) 

rmed <-  ggpar(plot_med, xlab="Simulated Counts ('000)", ylab="Predicted Counts ('000)",
               legend = "right", 
               legend.title = "Data \n Source",size=22,
               font.legend=c(12),
               font.label = list(size = 12, face = "bold", color ="red"),
               font.x = c(12),
               font.y = c(12),
               font.main=c(14),
               font.xtickslab =c(10),
               font.ytickslab =c(10),
               # orientation = "reverse",
               xtickslab.rt = 45, ytickslab.rt = 45)
rmed 


# Medium
dim(mod.predhgh <- mod.pred %>% filter(spat_cor ==1))
mod.predhgh$pop <- round(mod.predhgh$pop/1000) # per 1k
mod.predhgh$mean_pop_hat <- round(mod.predhgh$mean_pop_hat/1000) # per 1k

plot_hgh <- ggscatter(mod.predhgh, x = "pop", y = "mean_pop_hat",
                      add = "reg.line",                         # Add regression line
                      facet.by = "cover2",
                      conf.int = TRUE,                          # Add confidence interval
                      color = "source", 
                      palette = "lancet",           # Color by groups "cyl"
                      shape = "source"                             # Change point shape by groups "cyl"
) 

rhgh <-  ggpar(plot_hgh, xlab="Simulated Counts ('000)", ylab="Predicted Counts ('000)",
               legend = "right", 
               legend.title = "Data \n Source",size=22,
               font.legend=c(12),
               font.label = list(size = 12, face = "bold", color ="red"),
               font.x = c(12),
               font.y = c(12),
               font.main=c(14),
               font.xtickslab =c(10),
               font.ytickslab =c(10),
               # orientation = "reverse",
               xtickslab.rt = 45, ytickslab.rt = 45)
rhgh 

### Histograms----------------------------------
plot_hist <- gghistogram(mod.pred, x = "lmean",
                   color = "spat_cor", fill = "spat_cor",
                   add = "mean", #rug = TRUE,
                   facet.by = "cover2",
                   palette = "lancet",
                   bins = 50,
                   alpha = 0.5, ggtheme = theme_bw())
plot_hist

rhist <-  ggpar(plot_hist, xlab="log(Predicted Counts)", ylab="Frequency",
               legend = "right", 
               legend.title = "Spatial \n Variance",size=22,
               font.legend=c(12),
               # palette = c("#00AFBB", "#E7B800", "#FC4E07", "#0D0887FF", "#993333"),
               #palette = "jco",
               font.label = list(size = 12, face = "bold", color ="red"),
               font.x = c(12),
               font.y = c(12),
               font.main=c(14),
               font.xtickslab =c(10),
               font.ytickslab =c(10),
               # orientation = "reverse",
               xtickslab.rt = 45, ytickslab.rt = 45)
rhist 

#---------------------------------------
###  Mapping
#---------------------------------------
#  Load data       
shp <- st_read(paste0(shp_path, "/National/CMR_Boundary.shp")) # national shapefile 
dat1 <- st_read(paste0(shp_path, "/EA/CMR_Data3.shp")) # Combined EA shapefile
dim(dat <- as.data.frame(dat1)); names(dat1)

# visualise
shp %>%
  ggplot2::ggplot() +
  ggplot2::geom_sf()

# create grid
cmr_grid <-
  sf::st_make_grid(shp,
                   cellsize = 10000,
                   what = "polygons",
                   square = T)

# Add grid IDs as data.frame
cmr_grid <-
  sf::st_sf(id = 1:length(lengths(cmr_grid)),
            cmr_grid)

# Visualise
cmr_grid %>%
  ggplot2::ggplot() +
  ggplot2::geom_sf()


# Crop to CMR extents
cmr_grid_cropped <-
  sf::st_intersection(cmr_grid, 
                      shp %>% st_make_valid())


dim(cmr_grid)# 
dim(cmr_grid_cropped) #  

# Visualise
cmr_grid_cropped %>%
  ggplot2::ggplot() +
  ggplot2::geom_sf()

# check CRS - coordinates reference system
st_crs(cmr_grid_cropped)

# Reproject to lon - lat
cmr_shp_tf <- st_transform(cmr_grid_cropped, 
                           crs = "+proj=longlat +datum=WGS84")
st_crs(cmr_shp_tf) # lon-lat

# Convert to spatial object
cmr_shp_sp <- as(st_geometry(cmr_shp_tf),"Spatial")
coords <- data.frame(coordinates(cmr_shp_sp)) # extract the centroids
names(coords)
coords <- coords %>% rename(x = X1, y = X2)
plot(coords)

# visulise
plot(cmr_shp_sp, col="white")
points(coords$x, coords$y, 
       col=2, cex=0.5, pch="*")


# select data to visualise 
dim(mean100 <- mod.pred %>% filter(cover==100, spat_cor == 0.1)) # 100% coverage

length(mean100$mean_pop_hat)
cmr_shp_tf$mean <- mean100$mean_pop_hat # predicted population mean
cmr_shp_tf$pop <- mean100$pop # simulated ('observed') counts
cmr_shp_tf$cv <- mean100$cv # coefficient of variation

# visualise 
library(mapview)
mapview(cmr_shp_tf["mean"])
mapview(cmr_shp_tf["pop"])
plot(cmr_shp_tf["pop"])


# Using tmap
# Observed
tmap_mode(mode="plot")
cmr_pop <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="pop", title="Simulated Counts",
              #style="fixed", 
              legend.hist=F, 
             palette=viridis(100), 
              legend.show =F,
              breaks=c(0, 100, 1000, 5000, 10000, 
                20000, 50000, 100000, 150000, 200000))+
  tm_layout(legend.outside = T, legend.text.size=1.5, legend.title.size=2.5)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "", frame=F) 
cmr_pop

#----- Predicted
cmr_pred <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="mean", title="Simulated Counts",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(100), 
              legend.show =F,
              breaks=c(0, 100, 1000, 5000, 10000, 
                       20000, 50000, 100000, 150000, 200000))+
  tm_layout(legend.outside = T, legend.text.size=1.5, legend.title.size=2.5)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "", frame=F) 
cmr_pred


#----- Coefficient of Variation
cmr_cv <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv", title="Simulated Counts",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F#,
              #breaks=c(0, 100, 1000, 5000, 10000, 
                    #   20000, 50000, 100000, 150000, 200000)
              )+
  tm_layout(legend.outside = T, legend.text.size=1.5, legend.title.size=2.5)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "", frame=F) 
cmr_cv


tmap_arrange(cmr_pop,cmr_pred, cmr_cv,
             ncol = 3)
# ------------------------------------------------------------------
# Extract posterior estimates of coefficient of variation for mapping
#--------------------------------------------------------------------
# 100% Coverage
 # low spatial variance 
dim(dat100L <- mod.pred %>% filter(cover==100, spat_cor == 0.01))
cmr_shp_tf$cv100L <- dat100L$cv
cmr_cv100L <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv100L", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
  )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")

# medium spatial variance 
dim(dat100M <- mod.pred %>% filter(cover==100, spat_cor == 0.1))
cmr_shp_tf$cv100M <- dat100M$cv
cmr_cv100M <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv100M", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
  )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")

# high spatial variance
dim(dat100H <- mod.pred %>% filter(cover==100, spat_cor == 1))
cmr_shp_tf$cv100H <- dat100H$cv
cmr_cv100H <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv100H", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
  )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")


# 80% Coverage
dim(dat80L <- mod.pred %>% filter(cover==80, spat_cor == 0.01))
cmr_shp_tf$cv80L <- dat80L$cv
cmr_cv80L <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv80L", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
  )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")


dim(dat80M <- mod.pred %>% filter(cover==80, spat_cor == 0.1))
cmr_shp_tf$cv80M <- dat80M$cv
cmr_cv80M <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv80M", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
  )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")


dim(dat80H <- mod.pred %>% filter(cover==80, spat_cor == 1))
cmr_shp_tf$cv80H <- dat80H$cv
cmr_cv80H <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv80H", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
  )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")

# 60% Coverage 
dim(dat60L <- mod.pred %>% filter(cover==60, spat_cor == 0.01))
cmr_shp_tf$cv60L <- dat60L$cv
cmr_cv60L <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv60L", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
  )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")


dim(dat60M <- mod.pred %>% filter(cover==60, spat_cor == 0.1))
cmr_shp_tf$cv60M <- dat60M$cv
cmr_cv60M <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv60M", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
  )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")


dim(dat60H <- mod.pred %>% filter(cover==60, spat_cor == 1))
cmr_shp_tf$cv60H <- dat60H$cv
cmr_cv60H <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv60H", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
  )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")


# 40% Coverage
dim(dat40L <- mod.pred %>% filter(cover==40, spat_cor == 0.01))
cmr_shp_tf$cv40L <- dat40L$cv
cmr_cv40L <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv40L", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
  )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")


dim(dat40M <- mod.pred %>% filter(cover==40, spat_cor == 0.1))
cmr_shp_tf$cv40M <- dat40M$cv
cmr_cv40M <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv40M", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
  )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")


dim(dat40H <- mod.pred %>% filter(cover==40, spat_cor == 1))
cmr_shp_tf$cv40H <- dat40H$cv
cmr_cv40H<-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv40H", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
  )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")


# 20% Coverage
dim(dat20L <- mod.pred %>% filter(cover==20, spat_cor == 0.01))
cmr_shp_tf$cv20L <- dat20L$cv
cmr_cv20L <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv20L", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
  )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")


dim(dat20M <- mod.pred %>% filter(cover==20, spat_cor == 0.1))
cmr_shp_tf$cv20M <- dat20M$cv
cmr_cv20M <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv20M", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
  )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")


dim(dat20H <- mod.pred %>% filter(cover==20, spat_cor == 1))
cmr_shp_tf$cv20H <- dat20H$cv
cmr_cv20H<-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv20H", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
  )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")

# 5% Coverage
dim(dat5L <- mod.pred %>% filter(cover==5, spat_cor == 0.01))
cmr_shp_tf$cv5L <- dat5L$cv
cmr_cv100M <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv100M", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
  )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")


dim(dat5M <- mod.pred %>% filter(cover==5, spat_cor == 0.1))
cmr_shp_tf$cv5M <- dat5M$cv
cmr_cv5M <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv5M", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
  )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")


dim(dat5H <- mod.pred %>% filter(cover==5, spat_cor == 1))
cmr_shp_tf$cv5H <- dat5H$cv
cmr_cv5H <-  tm_shape(cmr_shp_tf)+ 
  tm_borders(
    col = "black",
    lwd = 1.5,
    lty = "solid"
  ) +
  tm_polygons(col="cv5H", title="Coefficient of \n Variation",
              #style="fixed", 
              legend.hist=F, 
              palette=plasma(100), 
              legend.show =F,
              breaks=c(0, 0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3)
              )+
  tm_layout(legend.outside = T, legend.text.size=2.5, legend.title.size=3)+
  #tm_compass(position = c("left", "top"), text.size=2, size = 2)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1)+
  tm_layout(main.title = "")


#-----
tmap_arrange(cmr_cv100L,cmr_cv100M, cmr_cv100H, 
             cmr_cv80L,cmr_cv80M, cmr_cv80H,
             cmr_cv60L,cmr_cv60M, cmr_cv60H,
             ncol = 3, nrow = 3)
#----------------------------------------------------------------


