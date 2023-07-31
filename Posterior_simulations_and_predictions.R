
####--TITLE: R SCRIPTS FOR MODELLED POPULATION ESTIMATES FOR CAMEROON BASED ON R-INLA----#####
####--METHODS: GEOSTATISTICAL BAYESIAN HIERARCHICAL REGRESSION MODEL---------------------#####
####--AUTHOR: DR CHIBUZOR CHRISTOPHER NNANATU-------------------------#####
####--INSTITUTION: WORLDPOP, UNIVERSITY OF SOUTHAMPTON----------------------------------##### 
####--DATE: DECEMBER 2022---------------------------------------------------------------####
####@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@####


rm(list=ls()) #----Clear the workspace
packages <- c("readxl", "xlsx","raster", "haven", "sf","sp", "tmap","tidyverse","rgdal",
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

drive_path <- "Z:/Projects/WP517763_GRID3/Working/CMR"
path <- paste0(drive_path, "/Chris_N/Rdata/inla/")
data_path <- paste0(drive_path, "/Ortis/Cleaned Data/Combined Data/")#---paths for survey data: 
shp_path <- paste0(drive_path, "/Ortis/Cleaned Data/Combined Data")#---paths for shapefile: 
cov_path <- paste0(drive_path, "/Chris_N/covariates/")#---paths for covariates: 
results_path <- paste0(drive_path, "/Chris_N/results")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
###----LOAD DATA       
dat <- read.csv(paste0(data_path,"CMR_Complete_Data3.csv"))# combined survey data
shp <- readOGR(dsn=shp_path, "CMR_Data3") #combined shapefile  
plot(shp); names(dat); table(dat$Dat_Typ); dat$Total_Building_Count
dim(shp); dim(dat) 
##------***************************************************************************--------


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

#@@@@@@@@@@@@@@@@@@---DATA PREPARATION--@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
####----Exclude observations with zero building counts --- Issues with the BFT raster
dim(dat); dim(shp)
zero <- which(dat$Total_Building_Count==0)
dat2 <- dat[-zero,] #---select data with non zero building counts
shp2 <- shp[-zero,]
dim(shp2); dim(dat2)
#@@@@@@@@@@@@----VARIABLES RECODING (OPTIONAL)-----@@@@@@@@@@@@@@@@@@@@
#dat2$source <- factor(dat2$Data_Category)
dat2$set_typ <- factor(dat2$Settlement_Type)
dat2$region <- factor(dat2$Regions)


#-************************************************************************
#-********--PREPARE DATA FOR INLA MODEL FITTING USING THE BEST COVARIATES**********************
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


#----Declare covariates for stack
covars_est <- dat2[,c("x2", "x3", "x17", "x21", "x32", "x34", "x40",
                      "x42", "set_reg", "set_typ", "region", "IDsr")]; dim(covars_est)

names(dat2)
#---Build the stack
stk_est <- inla.stack(data=list(y=dat2$dens_bldg), #the response
                      
                      A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                      
                      effects=list(c(list(Intercept=1), #the Intercept
                                     iset),  #the spatial index
                                   #the covariates
                                   list(covars_est)
                      ), 
                      #this is a quick name so you can call upon easily
                      tag='est')



##---Refit best fit model
form4 <- y ~  -1 + Intercept  + x2 + x3 + x17 + x21 + x32 + x34 + x40 + x42 + 
  f(set_reg, model='iid') + f(spatial.field, model=spde)+ f(set_typ, model='iid')+ f(IDsr, model='iid')

mod4 <-inla(form4, #the formula
            data=inla.stack.data(stk_est,spde=spde),  #the data stack
            family= 'gamma',   #which family the data comes from
            control.predictor=list(A=inla.stack.A(stk_est),compute=TRUE),  #compute gives you the marginals of the linear predictor
            control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
            verbose = FALSE) #can include verbose=TRUE to see the log of the model runs


#---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#---import and prepare the prediction covariates
#-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pred_data <- readRDS("Z:/Projects/WP517763_GRID3/Working/CMR/Chris_N/data/Rdata/inla/prediction_data.RDS")#Load gridded prediction covariates
data <- pred_data
head(pred_data)
names(data)[1:43] <- paste0("x",1:43, sep="") #--rename columns 
head(data)

#--Standardize covs
data[,1:43] <- apply(data[,1:43], 2, stdize)  ##----standardize all model covariates 


###----Prediction Projection matrix 
Aprediction <- inla.spde.make.A(mesh = mesh, loc = cbind(data$Lon, data$Lat))
dim(Aprediction)



####-------Run postrior simulation to obtain stable estimates and sampling distribution of totals
simDens <- function(model, dat, Aprediction, run)
{
  fixedeff <- dens_hat <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
  #inla.seed = as.integer(runif(1)*.Machine$integer.max)
  inla.seed = 2111236018

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
  
    dens_hat[,i] <- exp(fixedeff[,i]) #--predcited density
    pop_hat[,i] <- dens_hat[,i]*dat$CMR_buildings_count #--predicted count
    
  }
  mean_dens_hat <- apply(pop_hat, 1, mean, na.rm=T) 
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
run=100
system.time(str(sim.dens1 <- simDens(mod4,data,Aprediction, run)))




####---Join the posterior sample to the prediction data
data.sim <- data.frame(cbind(data[,c("CMR_Regions", "CMR_Department","CMR_Settlement_Classification",
                                     "CMR_Arrondissement")], sim.dens1$pop_hat))

names(data.sim)


##--------Calculate and save National total with uncertainties
nat_total <- function(dat, run)
{
  p_hat <- dat[,5:(run+4)]
  tots <- apply(p_hat,2, sum, na.rm=T) #Col sums
  
  tot_sd  <- sd(tots, na.rm=T)
  
  tot_mean  <- mean(tots, na.rm=T)
  
  tot_lower <- quantile(tots, probs=c(0.025))
  tot_median <- quantile(tots, probs=c(0.5))
  tot_upper <- quantile(tots, probs=c(0.975))
  
  return(estimates <- data.frame(estimates=unlist(list(total=tot_mean, lower=tot_lower, median=tot_median, upper=tot_upper))))
}
(national <- nat_total(data.sim, run))
(national <- data.frame(total= national[1,],
                        lower = national[2,],
                        median=national[3,],
                        upper=national[4,]))
#write.csv(national, file=paste0(results_path, "/estimates/National_estimates_final.csv"))



##---Regional estimates
##----
#reg_names <- data.frame(read.csv("Z:/Projects/WP517763_GRID3/Working/CMR/Data_release/Regions.csv")) #---region names and codes
#reg_names <- reg_names[order(reg_names$id),]

regional_est <- function(datr, run)
{
  #uniR <- unique(datr$CMR_Regions)
  uniR <- unique(reg_names$id)
  regnames <- unique(reg_names$libelle)
  outR <- matrix(0, nrow=length(uniR), ncol=5)
  for(j in uniR)
  {
    reg <- datr[datr$CMR_Regions==j,]
    rtots <- apply(reg[,5:(4+run)], 2, sum, na.rm=T)

    rtot_mean  <- mean(rtots, na.rm=T)
    rtot_sd <- sd(rtots, na.rm=T)

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
(regional.est <- regional_est(data.sim, 100))
sum(regional.est$total)
#write.csv(regional.est, file=paste0(results_path, "/estimates/regional_final.csv"))



##---Divisional estimates
##---
#div_names <- data.frame(read.csv("Z:/Projects/WP517763_GRID3/Working/CMR/Data_release/Department.csv")) #---division/department names and codes
#div_names <- div_names[order(div_names$id),]

divisional_est <- function(datd, run)
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

(divisional.est <- divisional_est(data.sim, 100))
sum(divisional.est$total)
#write.csv(divisional.est, file=paste0(results_path, "/estimates/department_final.csv"))





######-----SUBDIVISIONAL
#subd_names <- data.frame(read.csv("Z:/Projects/WP517763_GRID3/Working/CMR/Data_release/Arrondisement.csv")) #--subdivision names and codes
#subd_names <- subd_names[order(subd_names$id),]

sdivisional_est <- function(datsd, run)
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
(sdivisional.est <- sdivisional_est(data.sim, 100))
sum(sdivisional.est$total)
#write.csv(sdivisional.est, file=paste0(results_path, "/estimates/subdivision_final.csv"))


###########################################################
###-----Create the output raster files
ref_coords <- cbind(data$Lon, data$Lat) #--GPS of the grid cells
x <- as.matrix(ref_coords)


##---Extract stats of interest from the gridded posterior samples
sum(result <-sim.dens1$est_data$mean_pop_hat, na.rm=T)# --- Mean
resultM <- sim.dens1$est_data$median_pop_hat #----Median
resultL <- sim.dens1$est_data$lower_pop_hat #--Lower Bound
resultU <- sim.dens1$est_data$upper_pop_hat  #--Upper Bound
resultUN <- sim.dens1$est_data$uncert  #---Uncertainty
resultSD <- sim.dens1$est_data$sd_pop_hat #--Standard Deviation

##----Add estimates to grid data
datt2 <- sim.dens1$est_data
datt2$mean <- result
datt2$lower <- resultL
datt2$median <- resultM
datt2$upper <- resultU
datt2$uncert <- resultUN
datt2$sd <- resultSD


##--*****---MEAN--************
z <- as.matrix(result)
cmr_mean = rasterFromXYZ(cbind(x, z))
writeRaster(cmr_mean, filename=paste0(results_path, "/gridded_raster/mean_total_final.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))


##--*****---MEDIAN--************
zM <- as.matrix(resultM)
cmr_median = rasterFromXYZ(cbind(x, zM))
writeRaster(cmr_median, filename=paste0(results_path, "/gridded_raster/median_total_final.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))


##--*****---LOWER--************
zL <- as.matrix(resultL)
cmr_lower = rasterFromXYZ(cbind(x, zL))
writeRaster(cmr_lower, filename=paste0(results_path, "/gridded_raster/lower_total_final.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))

##--*****---UPPER--************
zU <- as.matrix(resultU)
cmr_upper = rasterFromXYZ(cbind(x, zU))
writeRaster(cmr_upper, filename=paste0(results_path, "/gridded_raster/upper_total_final.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))



##--*****---SD--************
zSD <- as.matrix(resultSD)
cmr_sd = rasterFromXYZ(cbind(x, zSD))
writeRaster(cmr_sd, filename=paste0(results_path, "/gridded_raster/sd_total_final.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))


##--*****---UNCERTAINTY--************
zUN <- as.matrix(resultUN)
cmr_uncertainty = rasterFromXYZ(cbind(x, zUN))
writeRaster(cmr_uncertainty, filename=paste0(results_path, "/gridded_raster/uncertainty_final.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))


