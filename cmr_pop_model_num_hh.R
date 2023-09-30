
####--TITLE: R SCRIPTS FOR MODELLED POPULATION ESTIMATES FOR CAMEROON BASED ON R-INLA----#####
####--METHODS: GEOSTATISTICAL BAYESIAN HIERARCHICAL REGRESSION MODEL---------------------#####
####--AUTHOR: DR CHIBUZOR CHRISTOPHER NNANATU & DR ORTIS YANKEY-------------------------#####
####--INSTITUTION: WORLDPOP, UNIVERSITY OF SOUTHAMPTON----------------------------------##### 
####--DATE: DECEMBER 2022---------------------------------------------------------------####
####@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@####

#$$$$$$$$$---------IMPORTANT NOTE!!!------------------------$$$$$$$$$$$$$$$$$$$$$$$
####################################################################################
#----PLEASE RUN WITH R VERSION 4.0.2 TO REPRODUCE EXACTLY SAME RESULTS-----------
#---ESTIMATES FROM OTHER VERSIONS WILL NORMALLY VARY SLIGHTLY----------------
#---DUE TO DIFFERENT OPTIMIZATION METHODS USED/COMPATIBILITY ISSUES-------------
##################################################################################
#

rm(list=ls()) #----Clear the workspace
packages <- c("raster", "haven", "sf","sp", "tmap","tidyverse","rgdal",
              "lattice", "gridExtra", "devtools", "rlang")
if(length(setdiff(packages, rownames(installed.packages())), type="binary") > 0) { 
  install.packages(setdiff(packages, rownames(installed.packages()))) }


#Installing INLA!!
if(length(setdiff("INLA", rownames(installed.packages()))) > 0){
  install.packages("INLA", type="binary", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
}
library(INLA)
lapply(packages, library, character.only = TRUE) ##--access the libraries

library(spdep)

#install.packages("DClusterm", type="binary")
library(DClusterm)
library(viridis)
###--Specify various file paths

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

out_path <- "//worldpop.files.soton.ac.uk/worldpop/Projects/WP517763_GRID3/Working/CMR/Final_Release"

#path <- paste0(drive_path, "/Chris_N/Rdata/inla/")
data_path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3/Working/CMR/Ortis/Cleaned Data/Combined Data/"#---paths for survey data: 
shp_path <- paste0(drive_path, "/Input Data/Input_Settlement Boundaries/EA")#---paths for shapefile: 
cov_path <- paste0(drive_path, "/Input Data/Input_Covariates/")#---paths for covariates: 
results_path <- paste0(out_path , "/Results/")
admin_path <- paste0(drive_path, "/Input Data/Admin files/")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
###----LOAD DATA       
dat <- read.csv(paste0(data_path,"CMR_Complete_Data4.csv"))# combined survey data
shp <- readOGR(dsn=shp_path, "CMR_Data3") #combined shapefile  
plot(shp); names(dat); table(dat$Dat_Typ); dat$Total_Building_Count
dim(shp); dim(dat) 
##------***************************************************************************--------

####----Map of Cameroon


####


##--convert to spdf
#-install.packages("tmap", type="binary")
library(tmap)
library(tmaptools)
library(rgdal)
library(sf)
library(raster)
#install.packages("cartography")
library('cartography') # mapping dedicated package
#install.packages("OpenStreetMap")
library(OpenStreetMap)

#Specify paths
national_path <- paste0(drive_path, "/Input Data/Input_Settlement Boundaries/National")
region_path <- paste0(drive_path, "/Input Data/Input_Settlement Boundaries/Regional")
dept_path <- paste0(drive_path, "/Input Data/Input_Settlement Boundaries/Departmental")
subdiv_path <- paste0(drive_path, "/Input Data/Input_Settlement Boundaries/Subdivisional")
raster_path <- paste0(drive_path, "/Output_Data/Gridded_final")

adm0 <- readOGR(dsn= paste0(national_path), layer = "CMR_Boundary")
adm1 <- readOGR(dsn= paste0(region_path), layer = "Region_SHP")
adm2 <- readOGR(dsn= paste0(dept_path), layer = "Departement_SHP")
adm3 <- readOGR(dsn= paste0(subdiv_path), layer = "arrondissement_shp")


##----Multiple facets maps by datasource 
tmap_options(check.and.fix = TRUE)
tm_shape(shp) +
  tm_polygons("hh_nm_c", palette = "RdYlBu") +
  #tm_polygons() +
  tm_borders() +
  tm_fill("hh_nm_c",
          palette = get_brewer_pal("YlGnBu"),
          legend.show = F,
          title = "Price in US$ per lb",
          style = "order")+
  tm_facets(by = "Dat_Typ")




tmap_mode("view")
tm_shape(adm1a)+ 
  #tm_polygons(col="libelle", title="Regions")+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=2)+
  #tm_compass(position = c("right", "top"), text.size=2)+
  tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 100, 250))+
  ylab("Latitude")+
  xlab("Longitude")+
  tm_layout(main.title = "(c)")

#tmap_save(t_reg3a, paste0(results_path,'/divisional_estimates_maps3.tiff'), dpi = 300, units="in")


names(adm1a)

######---------------

adm1b <- st_read(dsn= paste0(region_path), layer = "Region_SHP")

adm1b$Region <- adm1b$libelle
ggplot(adm1b) +
  geom_sf(aes(fill = Region))+
  geom_sf_label(aes(label = Region), size=6)+
  ylab("Latitude")+
  xlab("Longitude")+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y =element_text(size=20),
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14))
#



#####
tmap_mode("view")
map_cmr <- tm_basemap("OpenStreetMap")+
  tm_shape(adm1)+
 # tm_polygons(title="Regions")+
  #tm_fill(col="gray")+
  tm_borders(col="red", lwd=2)+
  tm_layout(legend.outside = TRUE)+
  # tm_compass(position = c("left", "bottom"))+
  tm_scale_bar(position = c("right", "bottom"), text.size=3, size=3, breaks=c(0, 100,  250))+
  tm_layout(main.title = "CMR")

map_cmr

tmap_save(map_cmr, paste0(out_path ,"/interactive_map_cmr.html"))



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@-----SOME KEY FUNCTIONS-------------@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
###-1)---Covariates scaling (Z-score)
stdize <- function(x)
{
  stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  return(stdz)
}
#stdize(data)

###-2)---Model fit test metrics
model_metrics <- function(obs, pred, upper, lower)
{
  residual = pred - obs
  INACCURACY = mean(abs(residual), na.rm=T)#MAE
  MSE = mean(residual^2, na.rm=T)
  RMSE = sqrt(MSE)
  BIAS = mean(residual, na.rm=T)
  IMPRECISION = sd(residual, na.rm=T)
  In_IC = mean(obs<upper & obs> lower)*100
  
  output <- list(MAE  = INACCURACY ,
                 RMSE = RMSE,
                 BIAS = abs(BIAS),
                 IMPRECISION = IMPRECISION,
                 In_IC = In_IC)
  return(output)
}
#model_metrics(obs, pred, upper, lower)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#@@@@@@@@@@@@@@@@@@---DATA PREPARATION--@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
####----Exclude observations with zero building counts --- Issues with the BFT raster
dim(dat); dim(shp)
zero <- which(dat$Total_Building_Count==0)
dat2 <- dat[-zero,] #---select data with non zero building counts
shp2 <- shp[-zero,]
dim(shp2); dim(dat2)

#dat2a <- dat[-zero,] 
#@@@@@@@@@@@@----VARIABLES RECODING (OPTIONAL)-----@@@@@@@@@@@@@@@@@@@@
#dat2$source <- factor(dat2$Data_Category)
dat2$set_typ <- factor(dat2$Settlement_Type)
dat2$region <- factor(dat2$Regions)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###----MODEL COVARIATES PREPARATION AND SCALING
names(dat2)
#cov_rows <- 17:59 #---model covariates are on columns 17 to 59 (43 covariates)
cov_rows <- 20:62 #---model covariates are on columns 17 to 59 (43 covariates)
colnames(dat2)[cov_rows] <- paste0("x", 1:43, sep="") #--rename model covariates as x1 - x43

#colnames(dat2a)[cov_rows] <- paste0("x", 1:43, sep="") 

head(dat2[, cov_rows])
dat2[,cov_rows] <- apply(dat2[,cov_rows], 2, stdize)   #z-score
head(dat2[, cov_rows])
head(dat2)

names(dat2)
dat2$av_hh <- dat2$Household_Count/dat2$Total_Building_Count #---density variable 
names(dat2)


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#---##############------VARIABLE SELECTION------###################----------#
library(car) ##--For calculating variance inflation factor (vif)
#######-------For Density=================================================================
#--%%%%%%%%%%%%%%%%%%%
covs_av_hh <- dat2[,c(65,cov_rows)] #--subset for variables selection
covs_av_hh$lav_hh <- log(covs_av_hh$av_hh) #--transform density on logscale
hist(covs_av_hh$lav_hh)
covs_av_hh1 <- covs_av_hh[,-1] %>% drop_na() #- model covariates only without NA's


#---@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##----Covariates selection using GLM STEPWISE REGRESSION
fav_hh <- glm(lav_hh ~., data=covs_av_hh1, family = gaussian)
summary(fav_hh)
stepav_hh <- stepAIC(fav_hh, scale = 0,
                        direction = c("both"),
                        trace = 1, keep = NULL, steps = 1000, use.start = FALSE,
                        k = 2)
stepav_hh

###-----------REFIT THE SELECTED MODEL AND CHOOSE only significant variables with VIF < 5
fav_hh <-  glm(formula = lav_hh ~ x1 + x2 + x3 + x4 + x7 + x9 + x10 + x11 + 
                 x12 + x13 + x14 + x16 + x17 + x19 + x20 + x21 + x22 + x32 + 
                 x33 + x37 + x39 + x40 + x41 + x42 + x43, family = gaussian, 
               data = covs_av_hh1)


summary(fav_hh)

##---VIF
vif_av_hh = vif(fav_hh )
vif_av_hh[which(vif_av_hh < 5)]

##----Selected significant covariates with VIF < 5
cov_av_hh <- c("x2","x3","x17","x19", "x20","x21","x32","x33","x37","x39","x41","x42")



fav_hh2 <-  glm(formula = lav_hh ~ x2 + x17 + x19 + x20 + x21 + x32 + 
                 x37 + x42, family = gaussian, 
               data = covs_av_hh1)


summary(fav_hh2)


cov_av_hh <- c("x2","x3","x17","x19", "x20","x21","x32","x37","x42")
cov_names <- names(dat)[cov_rows]
covariates_av_hh <- cov_names[c(2, 3,17, 19, 20, 21, 32, 37, 42)]
av_hh_covs <- data.frame(cov = cov_av_hh, name=covariates_av_hh)

##----Make correlation plot
require(corrplot)

png(paste0(results_path, "/plots/cor_plots for hh_num model.png"))
corrplot(
  cor(covs_av_hh[,cov_av_hh]),
  method = 'square',
  type = 'upper',
  tl.col = 'black',
  tl.cex = 1,
  col = colorRampPalette(c('purple', 'dark green'))(200)
)
dev.off()

#-************************************************************************
#-********--PREPARE DATA FOR INLA MODEL FITTING--*************************
#-*************************************************************************
##----Extract the coordiates
coords = coordinates(shp2)
dat2$longitude = coords[,1]
dat2$latitude = coords[,2]

#-------Extract boundary for mesh
bdry <- inla.sp2segment(shp2)
bdry$loc <- inla.mesh.map(bdry$loc)

###----Build non-hull mesh
non_convex_bdry <- inla.nonconvex.hull(coords, -0.035, -0.05, resolution = c(100, 100))
mesh <- inla.mesh.2d(boundary = non_convex_bdry, max.edge=c(0.5,1), 
                     offset = c(0.5, 1),
                     cutoff = 0.3)
par(mfrow=c(1,1))

png(paste0(results_path, "/plots/mesh.png"))
plot(mesh) #--plot to 
plot(shp, add=T, col= "grey")
points(coords, cex=0.6, col="red", pch=16)
dev.off()
mesh$n

#--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###---Build projector matrix A
A<-inla.spde.make.A(mesh=mesh,loc=as.matrix(coords));dim(A)

##---Create the SPDE
spde <- inla.spde2.matern(mesh, alpha=2)

##----specify the observation indices for estimation 
iset <- inla.spde.make.index(name = "spatial.field", spde$n.spde)


#-------Fitting the models----------------------------
#------testing various density models: Gamma and LogNormal#

#--Create nesting structure: settlement type within regions 
Zsr <- as(model.matrix( ~ 0 + region:set_typ, data = dat2), "Matrix") 
dat2$IDsr <- 1:nrow(dat2)
dat2$set_reg <- as.factor(apply(Zsr, 1, function(x){names(x)[x == 1]}))#




######-----Population density----------------------------------------------------------------
covars_est_pd <- dat2[,c("x2", "x3", "x17", "x21", "x32", "x34", "x40",
                      "x42", "set_reg", "set_typ", "region", "IDsr")]; dim(covars_est_pd) ##---Population density


names(dat2)
dat2$dens_bldg <- dat2$Imputed_LHHSIZE/dat2$Total_Building_Count #---density variable 
#---Build the stack
stk_est_pd <- inla.stack(data=list(y=dat2$dens_bldg), #the response
                      
                      A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                      
                      effects=list(c(list(Intercept=1), #the Intercept
                                     iset),  #the spatial index
                                   #the covariates
                                   list(covars_est_pd)
                      ), 
                      #this is a quick name so you can call upon easily
                      tag='est')


####-------Population density model-----------------------------------------------------
form9a <- y ~  -1 + Intercept  + x2 + x3 + x17 + x21 + x32 + x34 + x40 + x42 + 
  f(set_reg, model='iid') + f(spatial.field, model=spde)+ f(set_typ, model='iid')+ f(IDsr, model='iid')

mod9a <-inla(form9a, #the formula
             data=inla.stack.data(stk_est_pd,spde=spde),  #the data stack
             family= 'gamma',   #which family the data comes from
             control.predictor=list(A=inla.stack.A(stk_est_pd),compute=TRUE),  #compute gives you the marginals of the linear predictor
             control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
summary(mod9a)

mod9a$summary.fix
mod9a$summary.random
ind9a <-inla.stack.index(stk_est_pd, "est")$data
fit9a <- exp(mod9a$summary.linear.predictor[ind9a,"mean"])
fit9aL <- exp(mod9a$summary.linear.predictor[ind9a,"0.025quant"])
fit9aU <- exp(mod9a$summary.linear.predictor[ind9a,"0.975quant"])
#pred5a <-rpois(nrow(dat2), fit5a*dat2$Total_Building_Count)
sum(pred9a <- round(fit9a*dat2$Total_Building_Count))
sum(pred9aL <- round(fit9aL*dat2$Total_Building_Count))
sum(pred9aU <- round(fit9aU*dat2$Total_Building_Count))
(obs_pred9a <- data.frame(pop=dat2$Imputed_LHHSIZE,pred9a))
apply(obs_pred9a, 2, sum, na.rm=T)

##
dat2$dens_hat <- fit9a
#dat2$dens_hat1 <- stdize(dat2$dens_hat)
dat2$dens_hat1 <- stdize(dat2$dens_hat)

###------Number of households----------------------------------
#----Declare covariates for stack
covars_est <- dat2[, c("x2","x3","x17","x19", "x20","x21","x32","x37","x42", "dens_hat1",
                       "set_reg", "set_typ", "region", "IDsr")]; dim(covars_est)

names(dat2)
#---Build the stack
stk_est <- inla.stack(data=list(y=dat2$av_hh), #the response
                      
                      A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                      
                      effects=list(c(list(Intercept=1), #the Intercept
                                     iset),  #the spatial index
                                   #the covariates
                                   list(covars_est)
                      ), 
                      #this is a quick name so you can call upon easily
                      tag='est')

#--*******************SPECIFY AND FIT VARIOUS MODELS-***************-------------------
#*###---Model 6a-- gamma spatial (with settlement type random effect with covarites only) -----------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
form1a <- y ~ -1 + Intercept + x2 + x3 + x17 + x19 + x20 + x21 + x32 + x37 + x42 + 
  f(IDsr, model='iid') + f(set_typ, model='iid') + f(spatial.field, model=spde)

mod1a <-inla(form1a, #the formula
             data=inla.stack.data(stk_est,spde=spde),  #the data stack
             family= 'gamma',   #which family the data comes from
             control.predictor=list(A=inla.stack.A(stk_est),compute=TRUE),  #compute gives you the marginals of the linear predictor
             control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
summary(mod1a) #----model summary

mod1a$summary.fix #--extract fixed effects 
mod1a$summary.random #--extract random effects 
ind1 <-inla.stack.index(stk_est, "est")$data #--estimation indices 
fit1a <- exp(mod1a$summary.linear.predictor[ind1,"mean"]) #--extract the backtransformed density_hat
pred1a <-rpois(nrow(dat2), fit1a*dat2$Total_Building_Count)
sum(pred1a <- fit1a*dat2$Total_Building_Count) #--obtain the population_hat as the product of density and building count
(obs_pred1a<- data.frame(pop=dat2$Household_Count,pred1a))
apply(obs_pred1a, 2, sum, na.rm=T) ##--compare with observed total




#####------------------------------
form1b <- y ~ -1 + Intercept + x2 + x3 + x17 + x19 + x20 + x21 + x32 + x37 + x42 + f(dens_hat1)+
  f(IDsr, model='iid') + f(set_typ, model='iid') + f(spatial.field, model=spde)

mod1b <-inla(form1b, #the formula
             data=inla.stack.data(stk_est,spde=spde),  #the data stack
             family= 'gamma',   #which family the data comes from
             control.predictor=list(A=inla.stack.A(stk_est),compute=TRUE),  #compute gives you the marginals of the linear predictor
             control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
summary(mod1b) #----model summary

mod1b$summary.fix #--extract fixed effects 
mod1b$summary.random #--extract random effects 
ind1 <-inla.stack.index(stk_est, "est")$data #--estimation indices 
fit1b <- exp(mod1b$summary.linear.predictor[ind1,"mean"]) #--extract the backtransformed density_hat
pred1b <-rpois(nrow(dat2), fit1b*dat2$Total_Building_Count)
sum(pred1b <- fit1b*dat2$Total_Building_Count) #--obtain the population_hat as the product of density and building count
(obs_pred1b<- data.frame(pop=dat2$Household_Count,pred1b))
apply(obs_pred1b, 2, sum, na.rm=T) ##--compare with observed total




####----View the spatial fields of the best fit model
#looking at the spatial field and what it looks like
gproj <- inla.mesh.projector(mesh,  dims = c(300, 300))

#adm0 <- readOGR(dsn= paste0(admin_path, "National"), layer = "CMR_Boundary")

bb <- non_convex_bdry$loc 
library(splancs)
table(xy.in <- inout(gproj$lattice$loc,bb))


g.mean <- inla.mesh.project(gproj, mod1$summary.random$spatial.field$mean)
g.sd <- inla.mesh.project(gproj, mod1$summary.random$spatial.field$sd)

##---Remove points not in the study domain
g.mean[!xy.in] <- g.sd[!xy.in] <- NA

library(viridisLite)
col <- viridis(100)

png(paste0(results_path, "/plots/spatial_random_effects_surface.png"))
grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', cex.lab=2, main='Mean',col.regions = col),
             levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='SD' ,col.regions = col), nrow=1)


#levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='sd' ,col.regions = col), nrow=1)
dev.off()




par(mfrow=c(2,2), mar=c(3,3.5,0,0), mgp=c(1.5, .5, 0), las=0)
plot.default(inla.tmarginal(function(x) exp(x), mod9a$marginals.hy[[4]]),
             type='l', xlab=expression(kappa), ylab='Density', lwd=3)
plot.default(resf$marginals.variance.nominal[[1]], type='l',
             xlab=expression(sigma[x]^2), ylab='Density', lwd=3)
plot.default(resf$marginals.range.nominal[[1]], type='l',
             xlab='Practical range', ylab='Density', lwd=3)



#rm(covariates_dens_bldg, covariates_dens_bldg, covs_dens_bldg1, cov_dens_bldg, dat, pred9b)
#-***********************************************************************************
#-***********************************************************************************
#####------GRID CELL (PIXEL) PREDICTION-------------------###
#---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#---import and prepare the prediction covariates
#-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pred_data <- readRDS(paste0(cov_path, "prediction_data.RDS"))
data <- pred_data
head(pred_data)
names(data)[1:43] <- paste0("x",1:43, sep="") #--rename columns 
head(data)

#--Standardize covs
data[,1:43] <- apply(data[,1:43], 2, stdize)  ##----standardize all model covariates 


###----Prediction Projection matrix 
Aprediction <- inla.spde.make.A(mesh = mesh, loc = cbind(data$Lon, data$Lat))
dim(Aprediction)

#-******-----Obtain national estimates and credible intervals---------**********
#names(data)
#dim(datt2 <- data[, c("Lon", "Lat","x2", "x3", "x17", "x21", "x32", "x34", "x40", "x42", "CMR_buildings_count","Grid_ID",
#                     "CMR_Arrondissement","CMR_Department", "CMR_Regions", 
#                      "CMR_Settlement_Classification")])

#x2 + x3 + x17 + x19 + x20 + x21 + x32 + x37 + x42



#######################################################################################
###--------Population density-------------------------------------------------------
####-------Run the uncertainty estimates 
simDens <- function(model, dat, Aprediction, run)
{
  fixedeff <- dens_hat <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
  #inla.seed = as.integer(runif(1)*.Machine$integer.max)
  inla.seed = 2111236018
  #900099790
  #2103816957
  set.seed(inla.seed)
  print(inla.seed)
  m1.samp <- inla.posterior.sample(run, model, seed = inla.seed, selection=list(x2=1, x3=1, x17=1,
                                                                                x21=1, x32=1, x34=1,
                                                                                x40=1, x42=1),num.threads="1:1")
  
  sfield_nodes_mean <- model$summary.random$spatial.field['mean']
  field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
  for(i in 1:run)
  {
    fixedeff[,i] <- 
      model$summary.fixed['Intercept', 'mean'] +
      m1.samp[[i]]$latent[1,] * dat[,'x2'] +
      m1.samp[[i]]$latent[2,] * dat[,'x3'] +
      m1.samp[[i]]$latent[3,] * dat[,'x17'] +
      m1.samp[[i]]$latent[4,] * dat[,'x21'] +
      m1.samp[[i]]$latent[5,] * dat[,'x32'] + 
      m1.samp[[i]]$latent[6,] * dat[,'x34']  +
      m1.samp[[i]]$latent[7,] * dat[,'x40'] +
      m1.samp[[i]]$latent[8,] * dat[,'x42'] +
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[2]) + #---settlement type and region nested effects
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[5]) + #---settlement type random effect
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[6]) + #---cluster level random effect
      field_mean[,1]
    
    
    ##
    dens_hat[,i] <- exp(fixedeff[,i])
    pop_hat[,i] <- dens_hat[,i]*dat$CMR_buildings_count
    
  }
  #mean_pop_hat <- dat$pop_hat1 #
  mean_dens_hat <- apply(dens_hat, 1, mean, na.rm=T) 
  mean_pop_hat <- apply(pop_hat, 1, mean, na.rm=T) 
  median_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.5), na.rm=T) #
  lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) #
  upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) #
  sd_pop_hat <- apply(pop_hat, 1, sd, na.rm=T) #
  uncert_pop_hat <- (upper_pop_hat - lower_pop_hat)/mean_pop_hat#
  
  
  dat$mean_dens_hat <- mean_dens_hat
  dat$mean_pop_hat <- mean_pop_hat
  dat$median_pop_hat <- median_pop_hat
  dat$lower_pop_hat <- lower_pop_hat
  dat$upper_pop_hat <- upper_pop_hat
  dat$uncert_pop_hat <- uncert_pop_hat
  dat$sd_pop_hat <- sd_pop_hat
  
  
  
  output <- list(pop_hat = pop_hat,
                 est_data = dat)
  
}
#rm(dat.sim, sim.pops_nb, sim.pops, llg.est, prov.est2)
run=100
#rm(vif_dens_bldg,mod7a,mod7b,obs_pred9a,obs_pred9b,obs_pred6b,obs_pred6a,covars_est,adm1,adm1a)
system.time(str(sim.dens1 <- simDens(mod9a,data,Aprediction, run)))








###-------------Number of households-----------------------------------------------

####-------Run the uncertainty estimates 
sim_ave_hh_num1 <- function(model, dat, Aprediction, run)
{
  fixedeff <- av_hh_hat <- hh_hat <- matrix(0, nrow=nrow(dat), ncol = run)
  #inla.seed = as.integer(runif(1)*.Machine$integer.max)
  inla.seed = 2111236018
  #900099790
  #2103816957
  set.seed(inla.seed)
  print(inla.seed)
  m1.samp <- inla.posterior.sample(run, model, seed = inla.seed, selection=list(x2=1, x3=1, x17=1,
                                                                                x19=1, x20=1, x21=1, x32=1,
                                                                                x37=1, x42=1),num.threads="1:1")
  
  sfield_nodes_mean <- model$summary.random$spatial.field['mean']
  field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
  for(i in 1:run)
  {
    fixedeff[,i] <- 
      model$summary.fixed['Intercept', 'mean'] +
      m1.samp[[i]]$latent[1,] * dat[,'x2'] +
      m1.samp[[i]]$latent[2,] * dat[,'x3'] +
      m1.samp[[i]]$latent[3,] * dat[,'x17'] +
      m1.samp[[i]]$latent[4,] * dat[,'x19'] +
      m1.samp[[i]]$latent[5,] * dat[,'x20'] +
      m1.samp[[i]]$latent[6,] * dat[,'x21'] +
      m1.samp[[i]]$latent[7,] * dat[,'x32'] + 
      m1.samp[[i]]$latent[8,] * dat[,'x37']  +
      m1.samp[[i]]$latent[9,] * dat[,'x42'] +
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[3]) + #---settlement type random effect
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[2]) + #---cluster level random effect
      field_mean[,1]
    
    
    ##
    av_hh_hat[,i] <- exp(fixedeff[,i])
    hh_hat[,i] <- av_hh_hat[,i]*dat$CMR_buildings_count
    
  }
  #mean_pop_hat <- dat$pop_hat1 #
  mean_av_hh_hat <- apply(av_hh_hat, 1, mean, na.rm=T) 
  mean_hh_hat <- apply(hh_hat, 1, mean, na.rm=T) 
  median_hh_hat <- apply(hh_hat, 1, quantile, probs=c(0.5), na.rm=T) #
  lower_hh_hat <- apply(hh_hat, 1, quantile, probs=c(0.025), na.rm=T) #
  upper_hh_hat <- apply(hh_hat, 1, quantile, probs=c(0.975), na.rm=T) #
  sd_hh_hat <- apply(hh_hat, 1, sd, na.rm=T) #
  uncert_hh_hat <- (upper_hh_hat - lower_hh_hat)/mean_hh_hat#
  
  
  dat$mean_av_hh_hat <- mean_av_hh_hat
  dat$mean_hh_hat <- mean_hh_hat
  dat$median_hh_hat <- median_hh_hat
  dat$lower_hh_hat <- lower_hh_hat
  dat$upper_hh_hat <- upper_hh_hat
  dat$uncert_hh_hat <- uncert_hh_hat
  dat$sd_hh_hat <- sd_hh_hat
  
  
  
  output <- list(hh_hat = hh_hat,
                 est_data = dat)
  
}
#rm(dat.sim, sim.pops_nb, sim.pops, llg.est, prov.est2)
run=100
#rm(vif_dens_bldg,mod7a,mod7b,obs_pred9a,obs_pred9b,obs_pred6b,obs_pred6a,covars_est,adm1,adm1a)
system.time(str(sim.hh1 <- sim_ave_hh_num1(mod1a,data,Aprediction, run)))



###-------------------------------------------------------
####-------Run the uncertainty estimates 
sim_ave_hh_num2 <- function(model, dat, Aprediction, run)
{
  fixedeff <- av_hh_hat <- hh_hat <- matrix(0, nrow=nrow(dat), ncol = run)
  #inla.seed = as.integer(runif(1)*.Machine$integer.max)
  inla.seed = 2111236018
  #900099790
  #2103816957
  set.seed(inla.seed)
  print(inla.seed)
  m1.samp <- inla.posterior.sample(run, model, seed = inla.seed, selection=list(x2=1, x3=1, x17=1,
                                                                                x19=1, x20=1, x21=1, x32=1,
                                                                                x37=1, x42=1),num.threads="1:1")
  
  sfield_nodes_mean <- model$summary.random$spatial.field['mean']
  field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
  for(i in 1:run)
  {
    fixedeff[,i] <- 
      model$summary.fixed['Intercept', 'mean'] +
      m1.samp[[i]]$latent[1,] * dat[,'x2'] +
      m1.samp[[i]]$latent[2,] * dat[,'x3'] +
      m1.samp[[i]]$latent[3,] * dat[,'x17'] +
      m1.samp[[i]]$latent[4,] * dat[,'x19'] +
      m1.samp[[i]]$latent[5,] * dat[,'x20'] +
      m1.samp[[i]]$latent[6,] * dat[,'x21'] +
      m1.samp[[i]]$latent[7,] * dat[,'x32'] + 
      m1.samp[[i]]$latent[8,] * dat[,'x37']  +
      m1.samp[[i]]$latent[9,] * dat[,'x42'] +
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[4]) + #---settlement type random effect
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[3]) + #---cluster level random effect
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[2]) + #---Population density
      field_mean[,1]
    
    
    ##
    av_hh_hat[,i] <- exp(fixedeff[,i])
    hh_hat[,i] <- av_hh_hat[,i]*dat$CMR_buildings_count
    
  }
  #mean_pop_hat <- dat$pop_hat1 #
  mean_av_hh_hat <- apply(av_hh_hat, 1, mean, na.rm=T) 
  mean_hh_hat <- apply(hh_hat, 1, mean, na.rm=T) 
  median_hh_hat <- apply(hh_hat, 1, quantile, probs=c(0.5), na.rm=T) #
  lower_hh_hat <- apply(hh_hat, 1, quantile, probs=c(0.025), na.rm=T) #
  upper_hh_hat <- apply(hh_hat, 1, quantile, probs=c(0.975), na.rm=T) #
  sd_hh_hat <- apply(hh_hat, 1, sd, na.rm=T) #
  uncert_hh_hat <- (upper_hh_hat - lower_hh_hat)/mean_hh_hat#
  
  
  dat$mean_av_hh_hat <- mean_av_hh_hat
  dat$mean_hh_hat <- mean_hh_hat
  dat$median_hh_hat <- median_hh_hat
  dat$lower_hh_hat <- lower_hh_hat
  dat$upper_hh_hat <- upper_hh_hat
  dat$uncert_hh_hat <- uncert_hh_hat
  dat$sd_hh_hat <- sd_hh_hat
  
  
  
  output <- list(hh_hat = hh_hat,
                 est_data = dat)
  
}
#rm(dat.sim, sim.pops_nb, sim.pops, llg.est, prov.est2)
run=100
#rm(vif_dens_bldg,mod7a,mod7b,obs_pred9a,obs_pred9b,obs_pred6b,obs_pred6a,covars_est,adm1,adm1a)
system.time(str(sim.hh2 <- sim_ave_hh_num2(mod1b,data,Aprediction, run)))
####---Join the posterior sample to the prediction data


data.sim<- data.frame(cbind(data[,c("CMR_Regions", "CMR_Department","CMR_Settlement_Classification",
                                     "CMR_Arrondissement")], sim.dens1$pop_hat))

names(data.sim)



data.sim1<- data.frame(cbind(data[,c("CMR_Regions", "CMR_Department","CMR_Settlement_Classification",
                                     "CMR_Arrondissement")], sim.hh1$hh_hat))

names(data.sim1)


#################################
data.sim2 <- data.frame(cbind(data[,c("CMR_Regions", "CMR_Department","CMR_Settlement_Classification",
                                     "CMR_Arrondissement")], sim.hh2$hh_hat))

names(data.sim2)



###------National pop total
nat_total <- function(dat, run)
{
  h_hat <- dat[,5:(run+4)]
  tots <- apply(h_hat,2, sum, na.rm=T) #Col sums
  
  tot_sd  <- sd(tots, na.rm=T)
  
  tot_mean  <- mean(tots, na.rm=T)
  
  tot_lower <- quantile(tots, probs=c(0.025))
  tot_median <- quantile(tots, probs=c(0.5))
  tot_upper <- quantile(tots, probs=c(0.975))
  
  return(estimates <- data.frame(estimates=unlist(list(total=tot_mean, lower=tot_lower, median=tot_median, upper=tot_upper))))
}

(national <- nat_total(data.sim, run))
(nationalp <- data.frame(total= national[1,],
                         lower = national[2,],
                         median=national[3,],
                         upper=national[4,]))

write.csv(national2, file=paste0(results_path, "/estimates/National_estimates_hh_num.csv"))






##--------Calculate and save National total with uncertainties
nat_total_hh <- function(dat, run)
{
  h_hat <- dat[,5:(run+4)]
  tots <- apply(h_hat,2, sum, na.rm=T) #Col sums
  
  tot_sd  <- sd(tots, na.rm=T)
  
  tot_mean  <- mean(tots, na.rm=T)
  
  tot_lower <- quantile(tots, probs=c(0.025))
  tot_median <- quantile(tots, probs=c(0.5))
  tot_upper <- quantile(tots, probs=c(0.975))
  
  return(estimates <- data.frame(estimates=unlist(list(total=tot_mean, lower=tot_lower, median=tot_median, upper=tot_upper))))
}
(national <- nat_total_hh(data.sim1, run))

(national1 <- data.frame(total= national[1,],
                        lower = national[2,],
                        median=national[3,],
                        upper=national[4,]))


(national <- nat_total_hh(data.sim2, run))
(national2 <- data.frame(total= national[1,],
                        lower = national[2,],
                        median=national[3,],
                        upper=national[4,]))

rbind(national1, national2)
write.csv(national2, file=paste0(results_path, "/estimates/National_estimates_hh_num.csv"))



##---Regional estimates
##----
reg_names <- data.frame(read.csv(paste0(admin_path, "Regions.csv"))) #---region names and codes
reg_names <- reg_names[order(reg_names$id),]

regional_est_hh <- function(datr, run)
{
  #uniR <- unique(datr$CMR_Regions)
  uniR <- unique(reg_names$id)
  regnames <- unique(reg_names$libelle)
  outR <- matrix(0, nrow=length(uniR), ncol=5)
  for(j in uniR)
  {
    reg <- datr[datr$CMR_Regions==j,]
    rtots <- apply(reg[,5:(4+run)], 2, sum, na.rm=T)
    #rtot_mean  <- sum(reg$pop_hat, na.rm=T)
    rtot_mean  <- mean(rtots, na.rm=T)
    rtot_sd <- sd(rtots, na.rm=T)
    #rtot_lower <- rtot_mean - 2*rtot_sd
    #rtot_upper <- rtot_mean  + 2*rtot_sd
    #rtot_median <- quantile(rtots, probs=c(0.5))
    rtot_lower <- quantile(rtots, probs=c(0.025))
    rtot_median <- quantile(rtots, probs=c(0.5))
    rtot_upper <- quantile(rtots, probs=c(0.975))
    rtot_uncert <- (rtot_upper - rtot_lower)/rtot_mean
    
    restimates <- round(c(rtot_mean, rtot_lower, rtot_median,rtot_upper, rtot_uncert),2)
    outR[j,] <- restimates
  }
  outR <- data.frame(outR)
  return(reg_est <- data.frame(id = uniR,
                               names = regnames,
                               total = outR[,1],
                               lower = outR[,2],
                               median = outR[,3],
                               upper = outR[,4],
                               uncertainty = outR[,5]))
}
(regional.est1 <- regional_est_hh(data.sim1, 100))
sum(regional.est1$total)

(regional.est2 <- regional_est_hh(data.sim2, 100))
sum(regional.est2$total)

data.frame(ID=regional.est1 , name=regional.est1 ,
           est1 = regional.est1$total, est2=regional.est1$total)

write.csv(regional.est2, file=paste0(results_path, "/estimates/regional_hh_num.csv"))



##---Divisional estimates
##---
div_names <- data.frame(read.csv(paste0(admin_path, "Department.csv"))) #---division/department names and codes
div_names <- div_names[order(div_names$id),]

divisional_est_hh <- function(datd, run)
{
  #datd <-data.all
  uniD <- unique(div_names$id)
  outD <- matrix(0, nrow=length(uniD), ncol=5)
  divnames <- unique(div_names$libelle)
  for(k in uniD)
  {
    div <- datd[datd$CMR_Department==k,]
    dtots <- apply(div[,5:(4+run)], 2, sum, na.rm=T)
    
    dtot_mean  <- mean(dtots, na.rm=T)
    dtot_sd <- sd(dtots, na.rm=T)
    
    dtot_lower <- quantile(dtots, probs=c(0.025))
    dtot_median <- quantile(dtots, probs=c(0.5))
    dtot_upper <- quantile(dtots, probs=c(0.975))
    dtot_uncert <- (dtot_upper-dtot_lower)/dtot_mean
    destimates <- round(c(dtot_mean, dtot_lower, dtot_median,dtot_upper,  dtot_uncert),2)
    outD[k,] <- destimates
  }
  outD <- data.frame(outD)
  return(div_est <- data.frame(id = uniD,
                               names1 = divnames,
                               total = outD[,1],
                               lower = outD[,2],
                               median = outD[,3],
                               upper = outD[,4],
                               uncertainty = outD[,5]))
}

(divisional.est1 <- divisional_est_hh(data.sim1, 100))
sum(divisional.est1$total)


(divisional.est2<- divisional_est_hh(data.sim2, 100))
sum(divisional.est2$total)
write.csv(divisional.est2, file=paste0(results_path, "/estimates/department_final_hh_num.csv"))




######-----SUBDIVISIONAL
subd_names <- data.frame(read.csv(paste0(admin_path, "Arrondisement.csv"))) #--subdivision names and codes
subd_names <- subd_names[order(subd_names$id),]

sdivisional_est_hh <- function(datsd, run)
{
  #datsd <-data.all
  #uniSD <- unique(datsd$CMR_Arrondissement)
  uniSD <- unique(subd_names$id)
  outSD <- matrix(0, nrow=length(uniSD), ncol=5)
  sdivnames <- unique(subd_names$libelle)
  for(l in uniSD)
  {
    sdiv <- datsd[datsd$CMR_Arrondissement==l,]
    sdtots <- apply(sdiv[,5:(4+run)], 2, sum, na.rm=T)
    
    sdtot_mean  <- mean(sdtots, na.rm=T)
    sdtot_sd <- sd(sdtots, na.rm=T)
    
    sdtot_lower <- quantile(sdtots, probs=c(0.025))
    sdtot_median <- quantile(sdtots, probs=c(0.5))
    sdtot_upper <- quantile(sdtots, probs=c(0.975))
    sdtot_uncert <- (sdtot_upper-sdtot_lower)/sdtot_mean
    
    sdestimates <- round(c(sdtot_mean, sdtot_lower, sdtot_median,sdtot_upper, sdtot_uncert),2)
    outSD[l,] <- sdestimates
  }
  outSD <- data.frame(outSD)
  return(sdiv_est <- data.frame(id =  uniSD,
                                names = sdivnames,
                                total = outSD[,1],
                                lower = outSD[,2],
                                median = outSD[,3],
                                upper = outSD[,4],
                                uncertainty = outSD[,5]))
}
(sdivisional.est1 <- sdivisional_est_hh(data.sim1, 100))
sum(sdivisional.est1$total)


(sdivisional.est2 <- sdivisional_est_hh(data.sim1, 100))
sum(sdivisional.est2$total)
write.csv(sdivisional.est2, file=paste0(results_path, "/estimates/subdivision_final_hh_num.csv"))

#
###########################################################################################

############################

#################################################################################
#####-------------------Cross validation
####################################################################################
#(metrics <- model_metrics(dat2$Imputed_LHHSIZE,  
#                         pred9a, pred9aU,  pred9aL))


cross_val <- function(dat, mod, mod_pred, mesh, spde,
                      A, shp, formula, k_folds)
{
  set.seed(13235)                                     
  N <- nrow(dat)
  ######
  ind_train <- factor(sample(x = rep(1:k_folds, each = floor(N / k_folds)),  # Sample IDs for training data
                             size = N))
  
  table(as.numeric(ind_train)) 
  dat$k_fold <- as.numeric(ind_train)
  coords <- cbind(dat$longitude, dat$latitude)
  
  
  k_uniq <-sort(unique(dat$k_fold))
  metrics_cv <- matrix(0, nrow=length(k_uniq), ncol=5)#5 metrics
  corr_vec <- rep(0, length(k_uniq))
  par(mfrow=c(3,2), mar=c(2,2,1,1))
  plot(dat2$Household_Count, mod_pred, xlab = "Observed", 
       ylab = "Predicted", col=c('dark green','orange'),
       pch=c(16,16), cex.axis=1.5)
  abline(0,1)
  legend("topleft", c("Observed", "Predicted"), col=c("dark green", "orange"), pch=c(16,16),
         bty="n", cex=1.5)
  
  cvs <- list()
  for(i in 1:length(k_uniq))
  {
    #i=1
    print(paste0("fold_", i, sep=""))
    train_ind <- which(dat$k_fold!=k_uniq[i])
    dim(train <- dat[train_ind, ])#---train set for fold i
    dim(test <- dat[-train_ind, ]) #---test set for fold i
    
    train_coords <- coords[train_ind,]
    test_coords <- coords[-train_ind,]
    
    
    ###---Create projection matrices for training and testing datasets
    Ae<-inla.spde.make.A(mesh=mesh,loc=as.matrix(train_coords));dim(Ae) #training
    
    png(paste0(results_path,"/fold_",i,"_sample points for cross-validation.png"))
    plot(shp)
    points(train_coords, pch=16,  col="blue", cex=0.6)
    points(test_coords, pch=2,  col="red", cex=0.6)
    dev.off()
    
    
   # x2 + x3 + x17 + x19 + x20 + x21 + x32 + x37 + x42
    #####------------------------
    covars_train <- train[, c("x2","x3","x17","x19", "x20","x21","x32","x37","x42", "dens_hat1",
                              "set_reg", "set_typ", "region", "IDsr")]; dim(covars_train)
    stk_train <- inla.stack(data=list(y=train$dens_bldg), #the response
                            
                            A=list(Ae,1),  #the A matrix; the 1 is included to make the list(covariates)
                            
                            effects=list(c(list(Intercept=1), #the Intercept
                                           iset),  #the spatial index
                                         #the covariates
                                         list(covars_train)
                            ), 
                            tag='train')
    
    ####-----
    covars_test <- test[, c("x2","x3","x17","x19", "x20","x21","x32","x37","x42", "dens_hat1",
                            "set_reg", "set_typ", "region", "IDsr")]; dim(covars_test)
    
    
    
    ###---Rerun INLA for model test prediction
    model <-inla(formula, #the formula
                 data=inla.stack.data(stk_train,spde=spde),  #the data stack
                 family= 'gamma',   #which family the data comes from
                 control.predictor=list(A=inla.stack.A(stk_train),compute=TRUE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
    summary(model)
    
    
    
    sfield_nodes_mean <- mod$summary.random$spatial.field['mean']
    field_mean <- (A%*% as.data.frame(sfield_nodes_mean)[, 1])
    
    sfield_nodesL <- mod$summary.random$spatial.field['0.025quant']
    fieldL <- (A%*% as.data.frame(sfield_nodesL)[, 1])
    
    sfield_nodesU<- mod$summary.random$spatial.field['0.975quant']
    fieldU <- (A%*% as.data.frame(sfield_nodesU)[, 1])
    
    
   
    ##--------
    fixed <-  
      #m1.samp[[i]]$latent[1,] +
      model$summary.fixed['Intercept', 'mean'] +
      model$summary.fixed['x2', 'mean'] * test[,'x2'] +
      model$summary.fixed['x3', 'mean'] * test[,'x3'] +
      model$summary.fixed['x17', 'mean'] * test[,'x17'] +
      model$summary.fixed['x19', 'mean'] * test[,'x19'] +
      model$summary.fixed['x20', 'mean'] * test[,'x20'] +
      model$summary.fixed['x21', 'mean'] * test[,'x21'] +
      model$summary.fixed['x32', 'mean'] * test[,'x32'] +
      model$summary.fixed['x37', 'mean'] * test[,'x37'] +
      model$summary.fixed['x42', 'mean'] * test[,'x42'] +
      
                 
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$mean[2]) + #---density
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$mean[4]) + #---settlement type random effect
      mod$summary.random$IDsr['mean'][-train_ind,1] + #--uncorrelated spatial random effects
      
      field_mean[-train_ind,1]
    
    dens_ht <- exp(fixed)
    sum(hh_ht <- dens_ht*test$Total_Building_Count)
    
    
    ###
    
    
    ######----Lower
  
    fixedL <-  
      model$summary.fixed['Intercept', '0.025quant'] +
      model$summary.fixed['x2', '0.025quant'] * test[,'x2'] +
      model$summary.fixed['x3', '0.025quant'] * test[,'x3'] +
      model$summary.fixed['x17', '0.025quant'] * test[,'x17'] +
      model$summary.fixed['x19', '0.025quant'] * test[,'x19'] +
      model$summary.fixed['x20', '0.025quant'] * test[,'x20'] +
      model$summary.fixed['x21', '0.025quant'] * test[,'x21'] +
      model$summary.fixed['x32', '0.025quant'] * test[,'x32'] +
      model$summary.fixed['x37', '0.025quant'] * test[,'x37'] +
      model$summary.fixed['x42', '0.025quant'] * test[,'x42'] +
      
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$`0.025quant`[2]) + #---density
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$`0.025quant`[4]) + #---settlement type random effect
      mod$summary.random$IDsr['0.025quant'][-train_ind,1] + #--uncorrelated spatial random effects
      
      fieldL[-train_ind,1]
    
    dens_htL <- exp(fixedL)
    sum(hh_htL <- dens_htL*test$Total_Building_Count)
    
    
    #=========Upper
    fixedU <-   
      model$summary.fixed['Intercept', '0.975quant'] +
      model$summary.fixed['x2', '0.975quant'] * test[,'x2'] +
      model$summary.fixed['x3', '0.975quant'] * test[,'x3'] +
      model$summary.fixed['x17', '0.975quant'] * test[,'x17'] +
      model$summary.fixed['x19', '0.975quant'] * test[,'x19'] +
      model$summary.fixed['x20', '0.975quant'] * test[,'x20'] +
      model$summary.fixed['x21', '0.975quant'] * test[,'x21'] +
      model$summary.fixed['x32', '0.975quant'] * test[,'x32'] +
      model$summary.fixed['x37', '0.975quant'] * test[,'x37'] +
      model$summary.fixed['x42', '0.975quant'] * test[,'x42'] +
      
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$`0.975quant`[2]) + #---density
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$`0.975quant`[4]) + #---settlement type random effect
      mod$summary.random$IDsr['0.975quant'][-train_ind,1] + #--uncorrelated spatial random effects
      
      fieldU[-train_ind,1]
    
    dens_htU <- exp(fixedU)
    sum(hh_htU <- dens_htU*test$Total_Building_Count)
    
    
    
    plot(test$Household_Count, hh_ht, xlab = "Observed", 
         ylab = "Predicted", col=c('dark green','orange'),
         pch=c(16,16), cex.axis=1.5)
    abline(a=0,b=1)
    legend("topleft", c("Observed", "Predicted"), col=c("dark green", "orange"), pch=c(16,16),
           bty="n", cex=1.5) 
    
    
    
    corr_vec[i] <- round(cor(test$Household_Count, hh_ht),4)
    (met <- model_metrics(test$Household_Count,  
                          hh_ht, hh_htU,  hh_htL))
    metrics_cv[i, ] <- as.vector(unlist(met)) 
    
    cvs[[i]] <- data.frame(test=test$Household_Count, pred=hh_ht, lower=hh_htL, upper = hh_htU,
                           fold=rep(k_uniq[i], length(as.vector(hh_ht))))
  }
  
  metrics_cv <- data.frame(Metrics = c("MAE", "RMSE", "Abias","IMPRECISION","%Accuracy","corr"),
                           Fold_1=c(metrics_cv[1,],corr_vec[1]),
                           Fold_2 =c(metrics_cv[2,],corr_vec[2]),
                           Fold_3 =c(metrics_cv[3,],corr_vec[3]),
                           Fold_4 =c(metrics_cv[4,],corr_vec[4]),
                           Fold_5=c(metrics_cv[5,],corr_vec[5]))
  metrics_cv$mean <- c(apply(metrics_cv[,2:6], 1, mean))
  return(list(metrics=metrics_cv, data = cvs))
}

#print(nm)
k_folds=5
(cv <- cross_val(dat2, mod1b, pred1b, mesh, spde,
                 A, shp2, form1b, k_folds))


write.csv(cv, paste0(results_path, "/cross_validation.csv"))


ddt <- do.call("rbind",cv$data)
# Basic scatter plot.

###################


cv$metrics




####################

ddt$fold1 <- factor(ddt$fold)
levels(ddt$fold1) <- paste0("fold", 0:5, sep="")
table(ddt$fold1)

# Library
library(hrbrthemes)

p1 <- ggplot(ddt, aes(x=test, y=pred)) + 
  geom_point( color=c("69b3a2"), size=3) +
  #geom_smooth(method=lm , color="orange", fill=c("#69b3a2"), se=TRUE) +
  geom_smooth(method=lm , color="red", se=FALSE) +
  labs(
    x = "Observed Count (Test)",
    y = "Predicted Count"
  ) +
  theme(#panel.background = element_rect(fill = "#006994", colour = "black"),
    legend.position = "top", legend.spacing.x = unit(1.0, "cm"),
    legend.text = element_text(size = 18, color = "black"),
    legend.title = element_text(size = 20, face = "bold.italic"),
    text=element_text(size=24),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )+
  ggtitle("Out of Sample Cross-Validation")+
  facet_wrap(ddt$fold1)

p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




####-----Prediction------------------------------------
###-----Create the raster files
ref_image <- raster(paste0(cov_path,'CMR_buildings_count.tif'))
#plot(ref_image)
#spplot(ref_image)
#mapview(ref_image)
ref_coords <- cbind(pred_data$Lon, pred_data$Lat)
x <- as.matrix(ref_coords)




#--**********************************************************--#
sum(result <-sim.hh2$est_data$mean_hh_hat, na.rm=T)
resultM <- sim.hh2$est_data$median_hh_hat
resultL <- sim.hh2$est_data$lower_hh_hat
resultU <- sim.hh2$est_data$upper_hh_hat
resultUN <- sim.hh2$est_data$uncert
resultSD <- sim.hh2$est_data$sd_hh_hat

#----Add to grid data
datt2 <- sim.hh2$est_data

datt2$mean <- result
datt2$lower <- resultL
datt2$median <- resultM
datt2$upper <- resultU
datt2$uncert <- resultUN
datt2$sd <- resultSD
##--*****---MEAN--************
z <- as.matrix(result)
cmr_mean = rasterFromXYZ(cbind(x, z))

writeRaster(cmr_mean, filename=paste0(results_path, "gridded_raster/mean_total_final_hh_num.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))


##--*****---MEDIAN--************
zM <- as.matrix(resultM)
cmr_median = rasterFromXYZ(cbind(x, zM))

writeRaster(cmr_median, filename=paste0(results_path, "gridded_raster/median_total_final_hh_num.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))


##--*****---LOWER--************
zL <- as.matrix(resultL)
cmr_lower = rasterFromXYZ(cbind(x, zL))

writeRaster(cmr_lower, filename=paste0(results_path, "gridded_raster/lower_total_final_hh_num.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))

##--*****---UPPER--************
zU <- as.matrix(resultU)
cmr_upper = rasterFromXYZ(cbind(x, zU))

writeRaster(cmr_upper, filename=paste0(results_path, "gridded_raster/upper_total_final_hh_num.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))



##--*****---SD--************
zSD <- as.matrix(resultSD)
cmr_sd = rasterFromXYZ(cbind(x, zSD))

writeRaster(cmr_sd, filename=paste0(results_path, "gridded_raster/sd_total_final_hh_num.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))


##--*****---UNCERTAINTY--************
zUN <- as.matrix(resultUN)
cmr_uncertainty = rasterFromXYZ(cbind(x, zUN))

writeRaster(cmr_uncertainty, filename=paste0(results_path, "gridded_raster/uncertainty_final_hh_num.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))



####
CM1<-rasterize(x, ref_image, field=z, fun='last', background=0)
par(mfrow=c(1,1))
plot(CM1, main = 'Predicted Count')
###------More posterior plots



tmap_mode("view")
#################################
#--Observed regional total
reg_tot_obs <- dat2 %>% dplyr::group_by(Regions) %>%
  dplyr::summarise(count=sum(Household_Count, na.rm=T))

#--Predicted regional total
reg_tot_pred <- datt2 %>% dplyr::group_by(CMR_Regions) %>%
  dplyr::summarise(count=sum(round(mean), na.rm=T))

adm1a <- adm1
adm1a$observed_total <- reg_tot_obs$count[as.numeric(adm1a$id)]
#adm1a$predicted_total <-reg_tot_pred$count[as.numeric(adm1a$id)]
adm1a$predicted_total2 <-regional.est2$total[as.numeric(adm1a$id)]
adm1a$uncertainty <-regional.est2$uncertainty[as.numeric(adm1a$id)]


#library(mapview)
#mapview(adm1a, zcol="observed_total")
#mapview(adm1a, zcol="predicted_total")
tmap_mode("plot")
obs <- tm_shape(adm1a)+ 
  tm_polygons(col="observed_total", title="Observed Count",
              style="quantile", palette="viridis")+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=2)+
  tm_compass(position = c("right", "top"), text.size=2)+
  tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 100, 250))+
  tm_layout(main.title = "(a)") 
obs

##
pred <- tm_shape(adm1a)+ 
  tm_polygons(col="predicted_total2", title="Predicted Count",
              style="quantile", palette="viridis")+
  tm_layout(legend.outside = F, main.title = "(b)",legend.text.size=1, legend.title.size=2)+
  tm_compass(position = c("right", "top"), text.size=2)+
  tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 100, 250))

pred



uncert <- tm_shape(adm1a)+ 
  tm_polygons(col="uncertainty", title="Uncertainty",
              style="quantile",  palette="viridis")+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=2)+
  tm_compass(position = c("right", "top"), text.size=2)+
  tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 100, 250))+
  tm_layout(main.title = "(c)") 

uncert

t_reg1a <-tmap_arrange(obs, pred, uncert, nrow=1)
tmap_save(t_reg1a, paste0(results_path,'/regional_estimates_maps.tiff'), dpi = 300, units="in")
#
##
#############################################################################################
pred1 <- tm_shape(adm1a)+ 
  tm_polygons(col="predicted_total2", title="Predicted Count",
              style="quantile", palette=get_brewer_pal(palette="OrRd", n=9, plot=F))+
  tm_layout(legend.outside = F, legend.text.size=2, legend.title.size=3)+
  tm_compass(position = c("right", "top"), text.size=2)+
  tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 100, 250))+
  tm_layout(main.title = "(a)")

pred1

uncert <- tm_shape(adm1a)+ 
  tm_polygons(col="uncertainty", title="Uncertainty",
              style="quantile",  palette=get_brewer_pal(palette="OrRd", n=9, plot=F))+
  tm_layout(legend.outside = F, legend.text.size=2, legend.title.size=3)+
  tm_compass(position = c("right", "top"), text.size=2)+
  tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 100, 250))+
  tm_layout(main.title = "(b)") 

uncert

t_reg1b <- tmap_arrange(pred1, uncert, nrow=1)
# save map
# save image
tmap_save(t_reg1b, paste0(results_path,'/regional_estimates_maps2.tiff'), dpi = 300, units="in")


########################################
# Error bar plots-----------------------------------------
q <- ggplot(regional.est2) +
  geom_bar( aes(x=names, y=total), stat="identity", fill="#006994", alpha=0.8) +
  geom_errorbar( aes(x=names, ymin=lower, ymax=upper), width=0.4, colour="orange", alpha=0.9, size=1.5) +
  labs(
    x = "Region",
    y = "Population Count"
  ) +
  theme(#panel.background = element_rect(fill = "#006994", colour = "black"),
    legend.position = "top", legend.spacing.x = unit(1.0, "cm"),
    legend.text = element_text(size = 18, color = "black"),
    legend.title = element_text(size = 20, face = "bold.italic"),
    text=element_text(size=24),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )+ggtitle("Regional estimates of  total population and uncertainty")

q + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

q



##----Compare with projected values

reg_est <- read.csv(paste0(drive_path, "/Output_Data/Admin_estimates_final/regional_final_with_projections.csv"))


colors <- c("WorldPop Model" = "black", "NIS Projection" = "red")
q0 <- ggplot(reg_est) +
  geom_bar( aes(x=names, y=total), stat="identity", fill="#006994", alpha=0.8) +
  geom_errorbar( aes(x=names, ymin=lower, ymax=upper), width=0.4, colour="orange", alpha=0.9, size=1.5) +
  geom_point(aes(x=names, y=total, color="WorldPop Model"), shape=18, size = 5, alpha=0.9) +
  geom_point(aes(x=names, y= NIS.Projected, color="NIS Projection"), shape=18, size = 5, alpha=0.9) +
  labs(
    x = "Region",
    y = "Population Count",
    color= "Population Count"
  ) +
  scale_color_manual(values = colors) +
  theme(#panel.background = element_rect(fill = "#006994", colour = "black"),
    legend.position = "top", legend.spacing.x = unit(1.0, "cm"),
    legend.text = element_text(size = 18, color = "black"),
    legend.title = element_text(size = 20, face = "bold.italic"),
    text=element_text(size=24),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )+ggtitle("Regional estimates of  total population and uncertainty")

q0 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



#####################################################################################################

#ggsave(paste0(results_path,'/regional_estimates_maps.tiff'), t_reg, device = "tiff", dpi = 300)
###########################################################################
###------Departmental Estimates----------------------------------------------------------------------------
#########################################################################################################
dept_tot_obs <- dat %>% group_by(Department) %>%
  summarise(count=sum(Household_Count, na.rm=T))

#--Predicted regional total
#dept_tot_pred <- datt2 %>% group_by(CMR_Department) %>%
#  summarise(count=sum(round(mean), na.rm=T))

adm2a <- adm2
adm2a$observed_total <- dept_tot_obs$count[as.numeric(adm2a$id)]
#adm1a$predicted_total <-reg_tot_pred$count[as.numeric(adm1a$id)]
adm2a$predicted_total <-divisional.est$total[as.numeric(adm2a$id)]
adm2a$uncertainty <-divisional.est$uncertainty[as.numeric(adm2a$id)]


library(mapview)
mapview(adm2a, zcol="observed_total")
mapview(adm2a, zcol="predicted_total")
tmap_mode("plot")
obs2a <- tm_shape(adm2a)+ 
  tm_polygons(col="observed_total", title="Observed Count",
              style="quantile", palette="magma", style="cont", n=8)+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=2)+
  tm_compass(position = c("right", "top"), text.size=2)+
  tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 100, 250))+
  tm_layout(main.title = "(a)") 

obs2a

##
pred2a <- tm_shape(adm2a)+ 
  tm_polygons(col="predicted_total", title="Predicted Count",
              style="quantile", palette="magma", style="cont", n=8)+
  tm_layout(legend.outside = F, main.title = "(b)",legend.text.size=1, legend.title.size=2)+
  tm_compass(position = c("right", "top"), text.size=2)+
  tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 100, 250))

pred2a


uncert2 <- tm_shape(adm2a)+ 
  tm_polygons(col="uncertainty", title="Uncertainty",
              style="quantile",  palette="magma", style="cont", n=8)+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=2)+
  tm_compass(position = c("right", "top"), text.size=2)+
  tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 100, 250))+
  tm_layout(main.title = "(c)") 

uncert2

t_reg2a <-tmap_arrange(obs2a, pred2a, uncert2,nrow=1)
tmap_save(t_reg2a, paste0(results_path,'/divisional_estimates_maps2.tiff'), dpi = 300, units="in")
#
##
#############################################################################################
pred2b <- tm_shape(adm2a)+ 
  tm_polygons(col="predicted_total", title="Predicted Count",
              style="quantile", palette="magma", style="cont", n=8)+
  tm_layout(legend.outside = F, legend.text.size=2, legend.title.size=3)+
  tm_compass(position = c("right", "top"), text.size=2)+
  tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 100, 250))+
  tm_layout(main.title = "(a)")

pred2b

uncert2 <- tm_shape(adm2a)+ 
  tm_polygons(col="uncertainty", title="Uncertainty",
              style="quantile",  palette="magma", style="cont", n=8)+
  tm_layout(legend.outside = F, legend.text.size=2, legend.title.size=3)+
  tm_compass(position = c("right", "top"), text.size=2)+
  tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 100, 250))+
  tm_layout(main.title = "(b)") 

uncert2

t_reg2b <- tmap_arrange(pred2b, uncert2, nrow=1)
# save map
# save image
tmap_save(t_reg2b, paste0(results_path,'/divisional_estimates_maps2b.tiff'), dpi = 300, units="in")


########################################
# Error bar plots-----------------------------------------
q2 <- ggplot(divisional.est) +
  geom_bar( aes(x=names1, y=total), stat="identity", fill="#006994", alpha=0.8) +
  geom_errorbar( aes(x=names1, ymin=lower, ymax=upper), width=0.4, colour="orange", alpha=0.9, size=1.5) +
  labs(
    x = "Department",
    y = "Population Count"
  ) +
  theme(#panel.background = element_rect(fill = "#006994", colour = "black"),
    legend.position = "top", legend.spacing.x = unit(1.0, "cm"),
    legend.text = element_text(size = 18, color = "black"),
    legend.title = element_text(size = 20, face = "bold.italic"),
    text=element_text(size=24),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )+ggtitle("Divisional estimates of  total population and uncertainty")

q2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



##########################################################################################
###-----SUBDIVISIONAL------------------------------------------------------------
#########################################################################################
sdiv_tot_obs <- dat %>% group_by(Arrondissement) %>%
  summarise(count=sum(Household_Count, na.rm=T))

#--Predicted regional total
#sdiv_tot_pred <- datt2 %>% group_by(datt2$CMR_Arrondissement) %>%
#summarise(count=sum(round(mean), na.rm=T))

adm3a <- adm3
adm3a$observed_total <- sdiv_tot_obs$count[as.numeric(adm3a$id)]
#adm1a$predicted_total <-reg_tot_pred$count[as.numeric(adm1a$id)]
adm3a$predicted_total <-sdivisional.est2$total[as.numeric(adm3a$id)]
adm3a$uncertainty <-sdivisional.est$uncertainty[as.numeric(adm3a$id)]


library(mapview)
mapview(adm3a, zcol="observed_total")
mapview(adm3a, zcol="predicted_total")
tmap_mode("plot")
obs3a <- tm_shape(adm3a)+ 
  tm_polygons(col="observed_total", title="Observed Count",
              style="quantile", palette="plasma", style="cont", n=8)+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=2)+
  tm_compass(position = c("right", "top"), text.size=2)+
  tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 100, 250))+
  tm_layout(main.title = "(a)") 

obs3a 
##
pred3a <- tm_shape(adm3a)+ 
  tm_polygons(col="predicted_total", title="Predicted Count",
              style="quantile", palette="plasma", style="cont", n=8)+
  tm_layout(legend.outside = F, main.title = "(b)",legend.text.size=1, legend.title.size=2)+
  tm_compass(position = c("right", "top"), text.size=2)+
  tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 100, 250))

pred3a

uncert3 <- tm_shape(adm3a)+ 
  tm_polygons(col="uncertainty", title="Uncertainty",
              style="quantile",  palette="plasma", style="cont", n=8)+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=2)+
  tm_compass(position = c("right", "top"), text.size=2)+
  tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 100, 250))+
  tm_layout(main.title = "(c)")

uncert3

t_reg3a <-tmap_arrange(obs3a, pred3a, uncert3, nrow=1)
tmap_save(t_reg3a, paste0(results_path,'/divisional_estimates_maps3.tiff'), dpi = 300, units="in")
#


#############################################################################################
pred3b <- tm_shape(adm3a)+ 
  tm_polygons(col="predicted_total", title="Predicted Count",
              style="quantile", palette=get_brewer_pal(palette="OrRd", n=9, plot=F))+
  tm_layout(legend.outside = F, legend.text.size=2, legend.title.size=3)+
  tm_compass(position = c("right", "top"), text.size=2)+
  tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 100, 250))+
  tm_layout(main.title = "(a)")

pred3b

uncert3 <- tm_shape(adm3a)+ 
  tm_polygons(col="uncertainty", title="Uncertainty",
              style="quantile",  palette=get_brewer_pal(palette="OrRd", n=9, plot=F))+
  tm_layout(legend.outside = F, legend.text.size=2, legend.title.size=3)+
  tm_compass(position = c("right", "top"), text.size=2)+
  tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 100, 250))+
  tm_layout(main.title = "(b)")

uncert3

t_reg3b <- tmap_arrange(pred3b, uncert3, nrow=1)
# save map
# save image
tmap_save(t_reg3b, paste0(results_path,'/divisional_estimates_maps2b.tiff'), dpi = 300, units="in")
####------------------------------------------------------------------------

########################################
xy <- cbind(datt2$Lon, datt2$Lat)
spdf <- SpatialPointsDataFrame(coords = xy, data = datt2,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

spplot(spdf, "mean")

class(spdf)
gridded(spdf) <- T
spplot(spdf,"mean")

library(lattice)
library(sp)
spplot(spdf, zcol = names(spdf), scales = list(draw = FALSE))


#tm_shape(spdf)+ 
# tm_polygons(col="mean", title="Predicted Count",
#              style="quantile", palette=get_brewer_pal(palette="OrRd", n=9, plot=F)) #--ran out of memory!




ggplot(data = ddt, aes(x = pred, y =test)) +
  geom_point(aes(shape = factor(fold))) +
  geom_point(aes(color = factor(fold))) +
  geom_smooth(method = "lm", 
              se = FALSE, 
              aes(color = factor(fold))) +
  labs(title = "Survival Time from Malignant Melanoma",
       x = "Age (in years)",
       y = "Survival Time (in days)")


######################################
#--Model level
list_marginals <- mod9a$marginals.fitted.values

marginals <- data.frame(do.call(rbind, list_marginals))
marginals$hospital <- rep(names(list_marginals),
                          times = sapply(list_marginals, nrow))

plot(mod9a, plot.fixed.effects=TRUE, plot.lincomb=FALSE, plot.random.effects=FALSE,
     plot.hyperparameters=FALSE, plot.predictor=FALSE, plot.q=FALSE, plot.cpo=FALSE,
     single=FALSE)


plot(mod9a, plot.fixed.effects=FALSE, plot.lincomb=FALSE, plot.random.effects=FALSE,
     plot.hyperparameters=TRUE, plot.predictor=FALSE, plot.q=FALSE, plot.cpo=FALSE,
     single=FALSE)




######

par(mar=c(5,5,1,0.5))
toplot = c("0.025quant", "mean", "0.975quant")
matplot(mod9a$summary.random$set_reg[, toplot],
        lty =c(2,1,2), type="l", col=1, ylab="")



par(mar=c(5,5,1.5,0.5))
matplot(mod9a$summary.linear.predictor,
        lty =c(2,1,2), type="l", col=1, ylim=c(0,1), ylab="")







##############################################


#c("x2", "x3", "x17", "x21", "x32", "x34", "x40","x42")



##########################
##################################
######################################

require(ggplot2)
# Convert the variable dose from a numeric to a factor variable

# Change color by groups
#dp <- ggplot(ddt, aes(x=fold, y=pred, fill=fold)) + 
#geom_violin(trim=FALSE)+
#geom_boxplot(width=0.1, fill="white")+
#labs(title="Plot of cross-validation predictions by fold",x="k-Fold", y = "Predicted")
#dp + theme_classic()


#----Make further posterior plots

###---hISTOGRAMS

dp1 <- ggplot(ddt, aes(y=log(pred), fill=fold)) + 
  geom_density()+
  geom_histogram()+
  facet_wrap(~fold) +
  scale_fill_brewer(palette="Dark2") + theme_minimal()  
dp1



#----Make Violing aNd boxplots 
ddt$Fold <- factor(ddt$fold)
ggplot(ddt, aes(x=Fold,y=log(pred), fill=Fold)) + 
  geom_violin(trim=TRUE)+
  geom_boxplot(width=0.1, fill="white")+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.title.x =element_text(size=16),
        axis.title.y =element_text(size=16),
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=16)) +#change legend text font size
  coord_flip()




ggplot(ddt, aes(y=log(pred), fill=Fold)) + 
  geom_density(adjust=1.5, alpha=.4)+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x =element_text(size=15),
        axis.title.y =element_text(size=15),
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) +#change legend text font size
  coord_flip()


###----violin plot for the gridded density
datt2$Region <- factor(datt2$CMR_Regions)
table(datt2$Region)
levels(datt2$Region) <- reg_names$libelle


sum(datt2$mean_pop_hat[datt2$Region=="Sud"], na.rm=T)
mean(datt2$mean_pop_hat[datt2$Region=="Sud"], na.rm=T)



datt2$Scaled_Pop <- stdize(datt2$mean_pop_hat) 
sum(datt2$mean_pop_hat[datt2$Region=="Ouest"], na.rm=T)
mean(datt2$mean_pop_hat[datt2$Region=="Ouest"], na.rm=T)



table(datt2$Region)
# Use single color
ggplot(datt2, aes(y=log(mean_dens_hat), fill=Region)) +
  geom_histogram(trim=FALSE, fill='white', color="black")+
  #geom_boxplot(width=0.1) + 
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x =element_text(size=15),
        axis.title.y =element_text(size=15),
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) #change legend text font size


##---Boxplot for poplation  
ggplot(datt2, aes(x=Region,y=mean_pop_hat, fill=Region)) +
  # geom_violin(trim=FALSE, fill='white', color="black")+
  geom_boxplot(width=0.1) + 
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x =element_text(size=15),
        axis.title.y =element_text(size=15),
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) #change legend text font size





# Change colors
#p+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))





dat2$dens_pred <- fit9a
dens <- ggplot(dat2, aes(y=log(dens_pred))) + 
  geom_histogram(trim=FALSE, color="#E69F00") + 
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x =element_text(size=15),
        axis.title.y =element_text(size=15),
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14))+ #change legend text font size
  
  
  dat2$pop_pred <- pred9a
pop <- ggplot(dat2, aes(y=log(pop_pred))) + 
  geom_histogram(trim=FALSE, color="#E69F00") + 
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x =element_text(size=15),
        axis.title.y =element_text(size=15),
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) +#change legend text font size
  coord_flip()

library(ggpubr)
ggarrange(dens, pop, nrow=1)



###-----------GRID CELL
densg <- ggplot(datt2, aes(y=log(mean_dens_hat))) + 
  geom_histogram(trim=FALSE, color="#E69F00") + 
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x =element_text(size=15),
        axis.title.y =element_text(size=15),
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14))+ #change legend text font size
  coord_flip()


popg <- ggplot(datt2, aes(y=log(mean_pop_hat))) + 
  geom_histogram(trim=FALSE, color="#E69F00") + 
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x =element_text(size=15),
        axis.title.y =element_text(size=15),
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) +#change legend text font size
  coord_flip()

library(ggpubr)
ggarrange(densg, popg, nrow=1)
##------------


dp2 <- ggplot(dts, aes(y=mean_dens_hat)) + 
  geom_histogram(trim=FALSE)
dp2 + scale_fill_brewer(palette="Dark2") + theme_minimal()  




geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of cross-validation predictions by fold",x="k-Fold", y = "Predicted")
dp + theme_classic()


#####
# Continusous colors
dp + scale_fill_brewer(palette="Blues") + theme_classic()
# Discrete colors
dp1 + scale_fill_brewer(palette="Dark2") + theme_minimal()
# Gradient colors
dp + scale_fill_brewer(palette="RdBu") + theme_minimal()


#######

n_folds=5
par(mfrow=c(2,3), mar=c(4,4,2,2))
ddt$fold1 <- as.numeric(ddt$fold)
for(i in 1:n_folds)
{
  dts <- ddt[ddt$fold1==i,]
  with(dts, plot(test, pred, xlab = "Observed", type="p",
                 ylab = "Predicted", col=c("red", "#56B4E9"), cex=1.5,
                 pch=c(16,16), cex.axis=2, cex.lab=2))
  abline(0,1)
  # hack: we draw arrows but with very special "arrowheads"
  #with(dts, arrows(test, lower, test, upper, length=0.05, angle=90, code=3))
  legend("topleft", c("Observed", "Predicted"), col=c("red", "#56B4E9"), pch=c(16,16),
         bty="n", cex=1.5, text.width=2) 
}


set_reg<-data.frame(mod9a$summary.random$set_reg)
set_reg$ID2 <- as.numeric(as.factor(set_reg$ID))

###---region-settlement type interaction plots
par(mfrow=c(1,1), mar=c(12,5,2,2))
with(set_reg, plot(ID2, mean, xlab = "", type="p",
                   ylab = "Effect size", col="#56B4E9",
                   ylim=c(min(X0.025quant), max(X0.975quant)),
                   pch=c(16,16), cex.axis=2, cex.lab=2, cex=2, xaxt="n"))
axis(las=2,1, labels=set_reg$ID, at =1:40, cex.axis=1.4)
# hack: we draw arrows but with very special "arrowheads"
with(set_reg, arrows(ID2, X0.025quant, ID2, X0.975quant, 
                     length=0.05, angle=90, code=3, lwd=2))






###########--Settlement type only
set_typ<-data.frame(mod9a$summary.random$set_typ)
set_typ$ID2 <- as.numeric(as.factor(set_typ$ID))


par(mfrow=c(1,1), mar=c(12,5,2,2))
with(set_typ, plot(ID2, mean, xlab = "", type="p",
                   ylab = "Effect size", col="#56B4E9",#col="#56B4E9",
                   ylim=c(min(X0.025quant), max(X0.975quant)),
                   pch=c(16,16), cex.axis=2, cex.lab=2, cex=2, xaxt="n"))
axis(las=1,1, labels=set_typ$ID, at =1:4, cex.axis=1.8)
# hack: we draw arrows but with very special "arrowheads"
with(set_typ, arrows(ID2, X0.025quant, ID2, X0.975quant, 
                     length=0.05, angle=90, code=3, lwd=2))


###----------------------------

normalize <- function(x)
{
  xx <- x/sum(x)
  return(xx)
}


###
par(mfrow=c(1,1), mar=c(11,3,2,2))
regional <- regional.est
regional[, c(4,3,6)] <- apply(regional[, c(4,3,6)], 2, normalize)


with(regional, plot(id, total, xlab = "", type="p",
                    ylab = "Effect size", col="red",#col="#56B4E9",
                    ylim=c(min(lower), max(upper)),
                    pch=c(16,16), cex.axis=2, cex=1.5, xaxt="n"))
#abline(h=0, lwd=2, lty=2)
axis(las=2,1, labels=regional$names, at =1:10, cex.axis=1.8)
# hack: we draw arrows but with very special "arrowheads"
with(regional, arrows(id, lower, id, upper, 
                      length=0.05, angle=90, code=3, lwd=2))



####
###
par(mfrow=c(1,1), mar=c(11,3,2,2))
div <- divisional.est
div[, c(4,3,6)] <- apply(div[, c(4,3,6)], 2, normalize)


with(div, plot(id, total, xlab = "", type="p",
               ylab = "Effect size", col="red",#col="#56B4E9",
               ylim=c(min(lower), max(upper)),
               pch=c(16,16), cex.axis=2, cex=1.5, xaxt="n"))
#abline(h=0, lwd=2, lty=2)
axis(las=2,1, labels=div$names1, at =1:58, cex.axis=1.8)
# hack: we draw arrows but with very special "arrowheads"
with(div, arrows(id, lower, id, upper, 
                 length=0.05, angle=90, code=3, lwd=2))



## Draw the x-axis labels.
ggplot(data, aes(x, y)) +        # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper))



qplot(x = test, y = pred, color = fold, data = ddt) +
  geom_smooth(method = "lm", se=TRUE)



qplot(x = test, y = pred, data = ddt) +
  geom_smooth(method = "lm", se=TRUE, level = 0.90)+
  facet_wrap(~fold1)




ddt %>% 
  ggplot(aes(x = test, y = pred)) +
  geom_point(colour = "midnightblue", alpha = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper))+
  geom_smooth(method = "lm", se = F, colour = "brown") +
  facet_wrap(~fold1)
labs(x = "Speed", y = "Distance", color = NULL) +
  theme_test() +
  theme(legend.position = c(0.08, 0.85))


#########################################

myfit <-mod9a

###Obtain the posterior marginal distribution
res.zzy4.f <- inla.spde2.result(myfit, 'spatial.field', spde, do.transf=TRUE)

#####--Posterior marginal estimates of model parameters and 95% credible interval (Table 2)
c(mean=inla.emarginal(function(x) x, res.zzy4.f$marginals.variance.nominal[[1]]), 
  q=inla.hpdmarginal(0.95, res.zzy4.f$marginals.variance.nominal[[1]]))[c(1,2,3)]#marginal variance
c(mean=inla.emarginal(function(x) x, res.zzy4.f$marginals.range[[1]]), 
  q=inla.hpdmarginal(0.95, res.zzy4.f$marginals.range[[1]]))[c(1,2,3)] # nominal range

#####----------------------------------------------------------------------------------------

####---Posterior marginal plots of model parameters for Occurrence model diagnoistics (Figure S3, supp mat)
par(mfrow=c(2,2), mar=c(4,3.5,2,2), mgp=c(2.5, 1, 0), las=0)
plot.default(inla.tmarginal(function(x) exp(x), myfit$marginals.hy[[4]]),
             type='l', xlab=expression(kappa), ylab='Density', lwd=3, cex.lab=1.5, cex.axis=1.5)
plot.default(res.zzy4.f$marginals.variance.nominal[[1]], type='l',
             xlab=expression(sigma[c]^2), ylab='Density', lwd=3, cex.lab=1.5, cex.axis=1.5)
plot.default(res.zzy4.f$marginals.range.nominal[[1]], type='l',
             xlab='Practical range', ylab='Density', lwd=3, cex.lab=1.5, cex.axis=1.5)
#----------------------------------------------------------


##----Variance explained

sigma1 <- 1/ c(mean=inla.emarginal(function(x) x, res.zzy4.f$marginals.variance.nominal[[1]]), 
               q=inla.hpdmarginal(0.95, res.zzy4.f$marginals.variance.nominal[[1]]))[1]#---SPATILA
sigma2 <- 1/mod9a$summary.hyperpar$mean[6]  ----IID


sigma1/(sigma1+sigma2)# 94.8%




###################################################################################

save.image(paste0(results_path, "/posterior_samples/CMR_MODEL_main_workspace_final.Rdata"))
