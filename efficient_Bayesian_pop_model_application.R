####--TITLE: Efficient Bayesian Hierarchical Population Estimation Modelling using INLA-SPDE: 
#            Integrating Multiple Data Sources and Spatial Autocorrelation
####--AUTHOR: DR chris Nnanatu
####--INSTITUTION: WORLDPOP, UNIVERSITY OF SOUTHAMPTON------------------
####--DATE: DECEMBER 2022----------------------------------------------


#--------IMPORTANT NOTE!!!--------------------------------------------
#---------------------------------------------------------------------
#   PLEASE RUN WITH R VERSION 4.0.2 TO REPRODUCE EXACTLY SAME RESULTS
#   ESTIMATES FROM OTHER VERSIONS WILL NORMALLY VARY SLIGHTLY
#   DUE TO DIFFERENT OPTIMIZATION METHODS USED/COMPATIBILITY ISSUES
#---------------------------------------------------------------------

rm(list=ls()) #----Clear the workspace

# Load key libraries
packages <- c("raster", "haven", "sf","sp", "tmap","tidyverse",
              "lattice", "gridExtra", "devtools", "rlang")
if(length(setdiff(packages, rownames(installed.packages())), type="binary") > 0) { 
  install.packages(setdiff(packages, rownames(installed.packages()))) }


#Instal INLA
if(length(setdiff("INLA", rownames(installed.packages()))) > 0){
  install.packages("INLA", type="binary", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
}
library(INLA)
lapply(packages, library, character.only = TRUE) ##--access the libraries

library(spdep)

library(DClusterm)
library(viridis)
library(tmap)
library(tmaptools)
#---------------------------------------------------------------------------------
###--Specify key file paths
drive_path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3/Working/CMR/Chris_N/paper1/application/"

data_path <- paste0(drive_path, "data/")# all datasets 
output_path <- paste0(drive_path, "output/")# all results
#surv_path <- paste0(data_path, "Input_Survey_Data/")# combined household listing data
#shp_path <- paste0(data_path, "Input_Settlement_Boundaries/EA/")# shapefile - boundaries  
#cov_path <- paste0(data_path, "Input_Covariates/")# Model geospatial covariates
#admin_path <- paste0(data_path, "Admin_files/") # Other administrative files 
#---------------------------------------------------------------------------------

# Load and explore data       
dat <- read.csv(paste0(data_path,"Input_Survey_Data/CMR_Complete_Data3.csv"))# combined survey data
shp <- st_read(paste0(data_path, "Input_Settlement_Boundaries/EA/CMR_Data3.shp")) #combined shapefile  
names(shp)

# visualise the imputed household size
plot(shp["I_LHHSI"]) 
names(dat); table(dat$Dat_Typ); dat$Total_Building_Count
dim(shp); dim(dat) 


#  visualise data sources using tmap 
tmap_options(check.and.fix = TRUE)
tm_shape(shp) +
  tm_polygons("O_LHHSI", palette = "RdYlBu") +
  tm_borders() +
  tm_fill("O_LHHSI",
          palette = get_brewer_pal("YlGnBu"),
          legend.show = F,
          title = "Price in US$ per lb",
          style = "order")+
  tm_facets(by = "Dat_Typ")

# -----------------------------
#  Covariates scaling (Z-score)
# ---------------------------
stdize <- function(x)
{
  stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  return(stdz)
}
#stdize(data)

#--------------------------------------------------------------------------------
#   1)  INPUT DATA PREPARATION--
#-------------------------------------------------------------------------------
dim(dat); dim(shp)

# Some variable recoding
dat$bldg <- dat$Total_Building_Count # observed building count
(zerob <- which(dat$bldg==0)) # check if there are samples without building counts

dat$pop <- dat$Imputed_LHHSIZE # Observed population counts
(zerop <- which(dat$pop==0)) # check if there are samples without population count

dat$dens <- dat$pop/dat$bldg #---density variable 
dat$dens[is.infinite(dat$dens)] = NA # check if there are infinite terms and set them to NA
                                 # Infinite terms arise when zero serves as a denominator

dat$set_typ <- factor(dat$Settlement_Type) # settlement type
dat$region <- factor(dat$Regions) # regions 


### Extract standardized continuous geospatial covariates 
names(dat)
cov_rows <- 20:62 # model covariates are on columns 20 to 62 (43 covariates)
colnames(dat)[cov_rows] <- paste0("x", 1:43, sep="") # rename model covariates as x1 - x43

# standardize/scale 
dat[,cov_rows] <- apply(dat[,cov_rows], 2, stdize) 

# subset only the renamed geospatial covariates
dim(mcovs <- dat %>% dplyr::select(starts_with("x")) %>%
      dplyr::select(-X)) 
 
#----------------------------------------------------------------
# Below GLM-based covariates selection will be performed.
# - GLM model is fitted
# - Step-wise selection is performed 
# - selected covariates are further assessed
# - only significant covariates with 
# - variance inflation factor (vif) less than 5
# - are selected for Bayesian modelling 
#----------------------------------------------------------------

# Add response variable to the covariates only data
mcovs$resp <- log(dat$dens)


dim(mcovs <- mcovs %>% drop_na()) # Drop all NA's 

##----Covariates selection using GLM STEPWISE REGRESSION
fdens <- glm(resp ~., data=mcovs, family = gaussian)
summary(fdens)


step_dens <- stepAIC(fdens, scale = 0,
                        direction = c("both"),
                        trace = 1, keep = NULL, steps = 1000, use.start = FALSE,
                        k = 2)
step_dens

###-----------REFIT THE SELECTED MODEL AND CHOOSE only significant variables with VIF < 5
fdens2 <-  glm(formula = resp ~ x1 + x2 + x3 + x4 + x7 + x9 + x10 + x12 + 
                 x13 + x14 + x15 + x16 + x17 + x20 + x21 + x22 + x32 + x33 + 
                 x34 + x39 + x40 + x42 + x43, family = gaussian, data = mcovs)
summary(fdens2)

# Calculate the vif 
library(car) # required for the vif
vif1 = vif(fdens2)
vif1[(vif1 < 5)]

#
#x2       x3      x17      x20      x21      x32      x33      x34      x39      x40      x42 
#2.536825 4.517798 1.626241 4.387028 4.019219 1.670609 2.018067 2.391184 1.537882 3.953380 3.675827 

# Refit with the selected covariates and retain only the statistically significant covs
fdens3 <-  glm(formula = resp ~ x2 + x3 + x17 + x21 + x32 + # x20, x33 & x39 dropped (NS)
                 x34 + x40 + x42, family = gaussian, data = mcovs)
summary(fdens3)

##----Selected significant covariates with VIF < 5
covS <- c("x2","x3","x17","x21","x32","x34","x40","x42")

# -------------------------------------------------------------
# The names of the  selected covariates are given below
# --------------------------------------------------------------
# x2:  Conflicts data (ACLED)
# x3:  Explosions data (ACLED)
# x17: Distance to water bodies
# x21: Distance to herbaceous areas
# x32: Distance to local roads
# x34: Distance to market places
# x40: Slope
# x42: Night-time lights intensity
# --------------------------------------------------------------


# Visualise the matrix of correlation coefficients of the 
# selected covs
require(corrplot) # for correlation plot

#png(paste0(results_path, "/plots/cor_plots.png"))
corrplot(
  cor(mcovs[,covS]),
  method = 'square',
  type = 'upper',
  tl.col = 'black',
  tl.cex = 1,
  col = colorRampPalette(c('purple', 'dark green'))(200)
)
#dev.off()

#---------------------------------------------------------
### Test for Spatial Autoccorrelation using Moran's I stat
# Define the neighbors (use queen case)
names(shp)
#library(sf)
library(spdep)
#library(tmap)
# check for invalid geometries, and fix if any
table(st_is_valid(shp))
dim(valid <- st_make_valid(shp))
shp <- valid
table(st_is_valid(shp))

# Check for empty polygons, and remove if any
table(st_is_empty(shp))

# build neighbourhood structures 
nb <- poly2nb(shp, queen=TRUE)
lw <- nb2listw(nb, style="W", zero.policy=TRUE)

# Carry out the Moran's I test
moran.test(shp$O_LHHSI,lw, alternative="greater", zero.policy=TRUE)
              # O_LHHSI: is the population size variable
              # Moran's I stat = 0.149; p-value = 0.0002
              # there is spatial clsutering in the data
#-------------------------------------------------------------
# Below, the selected model covariates will be used to fit
# Bayesian Geostatistical Model using INLA-SPDE
# The codes of the key steps are provided
#-------------------------------------------------------------

# Select the coordinates of the EA centroids
shp2 <- as(st_geometry(shp), "Spatial")
dat$lon <- coordinates(shp2)[,1] # longitude
dat$lat <- coordinates(shp2)[,2] # latitude

coords = cbind(dat$lon, dat$lat)

# Build non-conver hull mesh
non_convex_bdry <- inla.nonconvex.hull(coords, -0.035, -0.05, resolution = c(100, 100))
mesh <- inla.mesh.2d(boundary = non_convex_bdry, max.edge=c(0.5,1), 
                     offset = c(0.5, 1),
                     cutoff = 0.3)

#png(paste0(results_path, "/plots/mesh.png"))
plot(mesh) #--plot to 
plot(shp, add=T, col= "grey")
points(coords, cex=0.6, col="red", pch=16)
#dev.off()
mesh$n # Number of mesh nodes

###---Build projector matrix A
A<-inla.spde.make.A(mesh=mesh,loc=as.matrix(coords));dim(A)

##---Create the SPDE
spde <- inla.spde2.matern(mesh, alpha=2)

##----specify the observation indices for estimation 
iset <- inla.spde.make.index(name = "spatial.field", spde$n.spde)

#------------------------------------------------------------------
# 2) Fitting Bayesian Hierarchical models 
#------------------------------------------------------------------
# create settlement type - regions interaction term
Zsr <- as(model.matrix( ~ 0 + region:set_typ, data = dat), "Matrix") 
dat$IDsr <- 1:nrow(dat)
dat$set_reg <- as.factor(apply(Zsr, 1, function(x){names(x)[x == 1]}))#interaction term


# Extract model covariates and variables for for data stack 
covars_est <- dat[,c("x2", "x3", "x17", "x21", "x32", "x34", "x40",
                      "x42", "set_reg", "set_typ", "region", "IDsr")]; dim(covars_est)

names(dat)


#---Build data stack
stk_est <- inla.stack(data=list(y=dat$dens), # the response variable 
                      
                      A=list(A,1),  # the A matrix; the 1 is included to make the list(covariates)
                      
                      effects=list(c(list(Intercept=1), # the Intercept
                                     iset),  # the spatial index
                                   # the covariates
                                   list(covars_est)
                      ), 
                      #this is a quick name so you can call upon easily
                      tag='est')


# Model 1: covs + Spatial autocor + settlement type
form1 <- y ~ -1 + Intercept + x2 + x3 + x17 + x21 + x32 + x34 + x40 + x42 + 
  f(IDsr, model='iid') + f(set_typ, model='iid')

mod1 <-inla(form1, #the formula
             data=inla.stack.data(stk_est,spde=spde),  #the data stack
             family= 'gamma',   #which family the data comes from
             control.predictor=list(A=inla.stack.A(stk_est),compute=TRUE),  #compute gives you the marginals of the linear predictor
             control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs


# Model 2: covs + Spatial autocor + settlement type  + region random effects
form2 <- y ~ -1 + Intercept + x2 + x3 + x17 + x21 + x32 + x34 + x40 + x42 + 
  f(IDsr, model='iid')  + f(region, model='iid') + f(spatial.field, model=spde)

mod2<-inla(form2, #the formula
             data=inla.stack.data(stk_est,spde=spde),  #the data stack
             family= 'gamma',   #which family the data comes from
             control.predictor=list(A=inla.stack.A(stk_est),compute=TRUE),  #compute gives you the marginals of the linear predictor
             control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs



# Model 3 : covs + Spatial autocor + settlement type - region interactions effects
form3 <- y ~  -1 + Intercept + x2 + x3 + x17 + x21 + x32 + x34 + x40 + x42 + 
  f(IDsr, model='iid') + f(set_reg, model='iid') + f(spatial.field, model=spde)

mod3 <-inla(form3, #the formula
             data=inla.stack.data(stk_est,spde=spde),  #the data stack
             family= 'gamma',   #which family the data comes from
             control.predictor=list(A=inla.stack.A(stk_est),compute=TRUE),  #compute gives you the marginals of the linear predictor
             control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs


# Model 4 : covs + Spatial autocor + settlement type - region interactions effects + region
form4 <- y ~  -1 + Intercept  + x2 + x3 + x17 + x21 + x32 + x34 + x40 + x42 + 
  f(set_reg, model='iid') + f(spatial.field, model=spde)+ f(set_typ, model='iid')+ f(IDsr, model='iid')

mod4 <-inla(form4, #the formula
             data=inla.stack.data(stk_est,spde=spde),  #the data stack
             family= 'gamma',   #which family the data comes from
             control.predictor=list(A=inla.stack.A(stk_est),compute=TRUE),  #compute gives you the marginals of the linear predictor
             control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs


## Run Model fit checks Based on DIC
DIC <- t(c(mod1 = mod1$dic$dic, mod2 = mod2$dic$dic,
           mod3 = mod3$dic$dic, mod4 = mod4$dic$dic))

WAIC <- t(c(mod1 = mod1$waic$waic, mod2 = mod2$waic$waic,
           mod3 = mod3$waic$waic, mod4 = mod4$waic$waic))

CPO <- t(c(mod1 = -sum(log(mod1$cpo$cpo), na.rm=T), mod2 = -sum(log(mod2$cpo$cpo), na.rm=T),
            mod3 = -sum(log(mod3$cpo$cpo), na.rm=T), mod4 = -sum(log(mod4$cpo$cpo), na.rm=T)))

(mod.select <- data.frame(DIC = t(DIC),
                         WAIC = t(WAIC),
                         CPO = t(CPO)))# Model 4 provided the best fit, thus, selected

# extract the fixed effects 
fixed.effects <- round(mod4$summary.fixed,4) # to 4 dps
fixed.effects

# Extract the back transformed linear predictor based on the best fit model (model 4)
ind4 <-inla.stack.index(stk_est, "est")$data # indices of the estimation mesh nodes 
fit4 <- exp(mod4$summary.linear.predictor[ind4,"mean"]) # back transformed linear predictor

# Obtain posterior estimates at EA level 
sum(pred4 <- round(fit4*dat$bldg))

# visualise the spatial fields of the best fit model
gproj <- inla.mesh.projector(mesh,  dims = c(300, 300))
bb <- non_convex_bdry$loc # the coordiates 

library(splancs) # for inout function for mapping
table(xy.in <- inout(gproj$lattice$loc,bb))
# FALSE  TRUE 
# 59771 30229


##---Remove points not in the study domain
g.mean <- inla.mesh.project(gproj, mod4$summary.random$spatial.field$mean)
g.sd <- inla.mesh.project(gproj, mod4$summary.random$spatial.field$sd)

g.mean[!xy.in] <- g.sd[!xy.in] <- NA

library(viridisLite)
col <- viridis(100)
grid.arrange(levelplot(g.mean, scales=list(draw=F), 
                       xlab='', ylab='', cex.lab=2, 
                       main='Mean',col.regions = col),
             levelplot(g.sd, scal=list(draw=F), xla='', 
                       yla='', main='SD',
                       col.regions = col), nrow=1)


#--------------------------------------------------------------------------
# The codes below provide the guides on how to:
#  1) run posterior simulations based on the parameters of the best fit model
#  2) Ensure there is a large enough (at least 30 iteration) sample size 
#  3) Predict population density parameters simulated at each iterations
#  4) Obtain the predicted population counts and save the output 
#  5) Save the full density and population counts matrices
#  6) Obtain the posterior statistics as well as uncertainty quantification 
#--------------------------------------------------------------------------
# Load the stack of prediction covariates
pred_data <- readRDS(paste0(data_path, "Input_Covariates/prediction_data.RDS"))
data <- pred_data
head(pred_data)
names(data)[1:43] <- paste0("x",1:43, sep="") #--rename covariates columns 
head(data)

#      Standardize covs
data[,1:43] <- apply(data[,1:43], 2, stdize)  ##----standardize all model covariates 

#    Prediction covs matrix
Apred <- inla.spde.make.A(mesh = mesh, loc = cbind(data$Lon, data$Lat))
dim(Apred)


####-------Run the uncertainty estimates 
simDens <- function(model, dat, Aprediction, run)
{
  fixedeff <- dens_hat <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
  #inla.seed = as.integer(runif(1)*.Machine$integer.max)
  inla.seed = 2111236018
  set.seed(inla.seed)
  print(inla.seed)
  m1.samp <- inla.posterior.sample(run, model, 
                                   seed = inla.seed, 
                                   selection=list(x2 = 1, 
                                                  x3 = 1, 
                                                  x17 = 1, 
                                                  x21 = 1,
                                                  x32 = 1,
                                                  x34 = 1,
                                                  x40 = 1, 
                                                  x42 = 1),
                                   num.threads="1:1")

  sfield_nodes_mean <- model$summary.random$spatial.field['mean']
  field_mean <- (Apred%*% as.data.frame(sfield_nodes_mean)[, 1])
  for(i in 1:run)
  {
    fixedeff[,i] <- 
      model$summary.fixed['Intercept', 'mean'] +
      m1.samp[[i]]$latent[1,] * dat[,'x2']  +
      m1.samp[[i]]$latent[2,] * dat[,'x3']  +
      m1.samp[[i]]$latent[3,] * dat[,'x17'] +
      m1.samp[[i]]$latent[4,] * dat[,'x21'] +
      m1.samp[[i]]$latent[5,] * dat[,'x32'] + 
      m1.samp[[i]]$latent[6,] * dat[,'x34'] +
      m1.samp[[i]]$latent[7,] * dat[,'x40'] +
      m1.samp[[i]]$latent[8,] * dat[,'x42'] +
      
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[2]) + #---settlement type and region nested effects
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[5]) + #---settlement type random effect
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[6]) + #---cluster level random effect
      field_mean[,1]
#
    dens_hat[,i] <- exp(fixedeff[,i])
    pop_hat[,i] <- dens_hat[,i]*dat$bld
    
  }
  dat$mean_dens_hat <- apply(dens_hat, 1, mean, na.rm=T) 
  dat$mean_pop_hat  <- apply(pop_hat, 1, mean, na.rm=T) #
  dat$lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) #
  dat$upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) #
  dat$sd_pop_hat <- apply(pop_hat, 1, sd, na.rm=T) #
  dat$cv <- dat$sd_pop_hat/dat$mean_pop_hat#
  
  output <- list(pop_hat = pop_hat,
                 est_data = dat)
  
}

run=100

data$bld <- data$CMR_buildings_count
system.time(str(sim.dens1 <- simDens(mod4,data,Apred, run)))


#  Join the posterior sample to the prediction data
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
  
  return(estimates <- data.frame(estimates=unlist(list(total=tot_mean, 
                                                       lower=tot_lower, 
                                                       median=tot_median, 
                                                       upper=tot_upper))))
}
(national <- nat_total(data.sim, run))
(national <- data.frame(total= national[1,],
                        lower = national[2,],
                        median=national[3,],
                        upper=national[4,]))
#write.csv(national, file=paste0(results_path, "/estimates/National_estimates_final.csv"))

##---Regional estimates
reg_names <- data.frame(read.csv(paste0(data_path, "Admin_Data/Regions.csv"))) #---region names and codes
reg_names <- reg_names[order(reg_names$id),]

regional_est <- function(datr, run)
{
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
div_names <- data.frame(read.csv(paste0(data_path, "Admin_Data/Department.csv"))) #---division/department names and codes
div_names <- div_names[order(div_names$id),]

divisional_est <- function(datd, run)
{
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

##------------------------------------------------------------------
#    Model Cross Validation (K-fold)
#-------------------------------------------------------------------
# Function: Model fit metrics 
mod_metrics2 <- function(obs, pred)
{
  residual = pred - obs
  MAE = mean(abs(residual), na.rm=T)# Mean Absolute Error
  MSE = mean(residual^2, na.rm=T) # Mean Square Error
  RMSE = sqrt(MSE) # Root Mean Square Error
  BIAS = mean(residual, na.rm=T) # Bias
  CORR = cor(obs[!is.na(obs)], pred[!is.na(obs)]) # Correlation Coefficient
  output <- list(MAE  = MAE,
                 RMSE = RMSE,
                 BIAS = abs(BIAS),
                 CC=CORR)
  return(output)
}

#----Extract settement type effects
set_t <- function(dat, st)
{
  uniq <- unique(dat$set_typ)
  uniq[1]
  for(i in  1:nrow(dat))
  {
    
    for(j in 1:3)
    {
      if(dat$set_typ[i]==uniq[j]) dat$set_typ2[i] = st[j]
    }
    
  }
  dat$set_typ2
}


cross_validate <- function(dat, n.folds, mod, form, A, seed)
{
#--------------------------------------------
# dat: the input survey data containing the all the variables
# n.folds: number of test (k) folds to use
# mod: the best model of the full or reference data
# A: the projection  matrix used in training the full data model
# seed: a random sample seed to make results reproducible 
#--------------------------------------------
seed = 13235
set.seed(seed) 
#dat <- dat2 # survey data
N <- nrow(dat)
#n.folds <- 5

######
table(ind_train <- factor(sample(x = rep(1:n.folds, each = floor(N / n.folds)),  # Sample IDs for training data
                                 size = N)))

table(as.numeric(ind_train)) 
dat$k_fold <- as.numeric(ind_train)
coords = cbind(dat$lon, dat$lat)


k_uniq <-sort(unique(dat$k_fold))



#---------------------------------------------------------------
#                   in-sample
#---------------------------------------------------------------

met_list_in <- list()
pred_list_in <- list()
for(i in 1:length(k_uniq))
{
  
  print(paste0("in-sample cross-validation using fold ", i, sep=""))
  test_ind <- which(dat$k_fold==k_uniq[i])
  dim(test <- dat[test_ind, ]) #---test set for fold i
  
  
  train_coords <- coords
  test_coords <- coords[test_ind,]
  
  
  # spatial random effects based on the full data best model
  sfield_nodes_mean <- mod$summary.random$spatial.field['mean']
  field_mean <- (A%*% as.data.frame(sfield_nodes_mean)[, 1])
  
  
  ##--------
  fixed <-  
    mod$summary.fixed['Intercept', 'mean'] +
    mod$summary.fixed['x2', 'mean'] * test[,'x2'] +
    mod$summary.fixed['x3', 'mean'] * test[,'x3'] +
    mod$summary.fixed['x17', 'mean'] * test[,'x17'] +
    mod$summary.fixed['x21', 'mean'] * test[,'x21'] +
    mod$summary.fixed['x32', 'mean'] * test[,'x32'] +
    mod$summary.fixed['x34', 'mean'] * test[,'x34'] + 
    mod$summary.fixed['x40', 'mean'] * test[,'x40']  +
    mod$summary.fixed['x42', 'mean'] * test[,'x42'] +
    
    
    rnorm(nrow(test), 0, 1/mod$summary.hyperpar$mean[2]) + #---settlement type and region nested effects
    rnorm(nrow(test), 0, 1/mod$summary.hyperpar$mean[5]) + #---settlement type random effect
    
    
    mod$summary.random$IDsr['mean'][test_ind,1] + #--uncorrelated spatial random effects
    
    field_mean[test_ind,1]
  
  dens_ht <- exp(fixed)
  sum(pop_ht <- dens_ht*test$bld)
  

  # visualise samples
  #par(mfrow =c(1,2))
  #plot(shp)
  #points(train_coords, col="blue", pch =15, cex=0.6)
  #points(test_coords, col="orange", pch=15, cex=0.6)
  
  
  par(mfrow =c(1,1))
  plot(test$obs, pop_ht, xlab = "Observed", 
       ylab = "Predicted", col=c('blue','orange'),
       pch=c(16,16), cex.axis=1.5)
  abline(0,1)
  legend("topleft", c("Observed", "Predicted"), col=c("blue", "orange"), pch=c(16,16),
         bty="n", cex=1.5) 
  
  
  
  # calculate fit metrics
  met_in <- mod_metrics2(test$resp,  
                          pop_ht)
  
  met_list_in[[i]]<- unlist(met_in)
  pred_list_in[[i]] <- data.frame(obs = test$obs, pred = pop_ht,
                                  fold = rep(i, length(test$obs)),
                                  data = rep("insample", length(test$obs)))
}
met_list_in_dat <- do.call(rbind,met_list_in)
metrics_in <- apply(met_list_in_dat, 2, mean)
pred_list_in_dat <- do.call(rbind,pred_list_in)

#-----------------------------------------------------------
#               out - of -sample
#-----------------------------------------------------------
met_list_out <- list()
pred_list_out <- list()
for(i in 1:length(k_uniq))
{
  
  print(paste0("out-of-sample cross-validation using fold ", i, sep=""))
  train_ind <- which(dat$k_fold!=k_uniq[i])
  test_ind <- which(dat$k_fold==k_uniq[i])
  dim(train <- dat[train_ind, ])#---train set for fold i
  dim(test <- dat[test_ind, ]) #---test set for fold i
  
  
  train_coords <- coords[train_ind,]
  test_coords <- coords[test_ind,]
  
  
  ###---Create projection matrices for training and testing datasets
  Ae<-inla.spde.make.A(mesh=mesh,loc=as.matrix(train_coords));dim(Ae) #training
  
  
  #####------------------------
  covars_train <- train[,c("x2", "x3", "x17", "x21", "x32", "x34", "x40",
                           "x42", "set_reg", "set_typ", "region", "IDsr")]; dim(covars_train)
  stk_train <- inla.stack(data=list(y=train$dens), #the response
                          
                          A=list(Ae,1),  #the A matrix; the 1 is included to make the list(covariates)
                          
                          effects=list(c(list(Intercept=1), #the Intercept
                                         iset),  #the spatial index
                                       #the covariates
                                       list(covars_train)
                          ), 
                          tag='train')
  
  
  ###---Rerun INLA for model test prediction
  model <-inla(form, #the formula
               data=inla.stack.data(stk_train,spde=spde),  #the data stack
               family= 'gamma',   #which family the data comes from
               control.predictor=list(A=inla.stack.A(stk_train),compute=TRUE),  #compute gives you the marginals of the linear predictor
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
               verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
  summary(model)
  
  
  # Extract Spatial random effects from the full data best model
  sfield_nodes_mean <- mod$summary.random$spatial.field['mean']
  field_mean <- (A%*% as.data.frame(sfield_nodes_mean)[, 1])
  
  
  ##--------
  fixed <-  
    model$summary.fixed['Intercept', 'mean'] +
    model$summary.fixed['x2', 'mean'] * test[,'x2'] +
    model$summary.fixed['x3', 'mean'] * test[,'x3'] +
    model$summary.fixed['x17', 'mean'] * test[,'x17'] +
    model$summary.fixed['x21', 'mean'] * test[,'x21'] +
    model$summary.fixed['x32', 'mean'] * test[,'x32'] +
    model$summary.fixed['x34', 'mean'] * test[,'x34'] + 
    model$summary.fixed['x40', 'mean'] * test[,'x40']  +
    model$summary.fixed['x42', 'mean'] * test[,'x42'] +
    
   
    rnorm(nrow(test), 0, 1/model$summary.hyperpar$mean[2]) + #---settlement type and region nested effects
    rnorm(nrow(test), 0, 1/model$summary.hyperpar$mean[5]) + #---settlement type random effect
    
    
    mod$summary.random$IDsr['mean'][test_ind,1] + #--uncorrelated spatial random effects
    
    field_mean[test_ind,1]
  
  dens_ht <- exp(fixed)
  sum(pop_ht <- dens_ht*test$bld)
  

  # visualise samples
 # par(mfrow =c(1,2))
  #plot(shp)
  #points(train_coords, col="blue", pch =15, cex=0.6)
  #points(test_coords, col="orange", pch=15, cex=0.6)
  
  par(mfrow =c(1,1))
  plot(test$obs, pop_ht, xlab = "Observed", 
       ylab = "Predicted", col=c('blue','orange'),
       pch=c(16,16), cex.axis=1.5)
  abline(0,1)
  legend("topleft", c("Observed", "Predicted"), col=c("blue", "orange"), pch=c(16,16),
         bty="n", cex=1.5) 
  
  
  
  # calculate fit metrics
  met_out <- mod_metrics2(test$resp,  
                           pop_ht)
  
  met_list_out[[i]]<- unlist(met_out)
  pred_list_out[[i]] <- data.frame(obs = test$obs, pred = pop_ht,
                                  fold = rep(i, length(test$obs)),
                                  data = rep("outsample", length(test$obs)))
}
met_list_out_dat <- do.call(rbind,met_list_out)
metrics_out <- apply(met_list_out_dat, 2, mean) # fit metrics

pred_list_out_dat <- do.call(rbind,pred_list_out)# predictions

cv_mets <- rbind(metrics_in, metrics_out)
output <- list( met_list_in_dat = met_list_in_dat,
                met_list_out_dat = met_list_out_dat,
               pred_dat = rbind(pred_list_in_dat, pred_list_out_dat),
               cv_metrics = rbind(metrics_in, metrics_out))
}

dat <- dat
dat$obs <- dat$Imputed_LHHSIZE # observed household size

(cross_val <- cross_validate(dat, n.folds = 5, 
                             mod = mod4, 
                             form = form4,
                             A, 
                             seed = 13235))
  
cross_val$met_list_in_dat  # in-sample metrics per fold 
cross_val$met_list_out_dat  # out-of-sample metrics per fold
cross_val$cv_metrics    # combined averaged metrics
cross_val$pred_dat  # combined prediction data

#-------------------------------------------------------------
# Scatter plots of model cross-validation results
#------------------------------------------------------------
pred.data <- cross_val$pred_dat
pred.data$Fold <- factor(pred.data$fold)
pred.data$data <- factor(pred.data$data,
                         levels = c("insample", "outsample"),
                         labels = c("In-Sample", "Out-of-Sample"))

table(pred.data$data)
names(pred.data)

library(ggpubr)
plot_cval <- ggscatter(pred.data, x = "obs", y = "pred",
                      add = "reg.line",                         # Add regression line
                      facet.by = "data",
                      conf.int = TRUE,                          # Add confidence interval
                      color = "Fold", 
                      palette = "lancet",           # Color by groups "cyl"
                      shape = "Fold"                             # Change point shape by groups "cyl"
) 

rcval <-  ggpar(plot_cval, xlab="Observed counts", ylab="Predicted Counts",
               legend = "top", 
               legend.title = "Fold (k{=5}-fold)",size=20,
               font.legend=c(12),
               font.label = list(size = 12, face = "bold", color ="red"),
               font.x = c(12),
               font.y = c(12),
               font.main=c(14),
               font.xtickslab =c(12),
               font.ytickslab =c(12),
               # orientation = "reverse",
               xtickslab.rt = 45, ytickslab.rt = 45)
rcval


# -----------------------------------------------------------
# Prepare and save the grid vector predictions as raster files
# -----------------------------------------------------------

# specify the reference coordinates for the predictions
ref_coords <- cbind(pred_data$Lon, pred_data$Lat)
x <- as.matrix(ref_coords)

# Extract the grid cell predictions
sum(result <-sim.dens1$est_data$mean_pop_hat, na.rm=T) # mean
resultL <- sim.dens1$est_data$lower_pop_hat # lower bound
resultU <- sim.dens1$est_data$upper_pop_hat # upper bound
resultCV <- sim.dens1$est_data$cv # coefficient of variation

# add grid cell values and write the raster files

#  Mean
z <- as.matrix(result)
cmr_mean = rasterFromXYZ(cbind(x, z))

writeRaster(cmr_mean, filename=paste0(output_path, "gridded_mean_total.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))


# Lower 
zL <- as.matrix(resultL)
cmr_lower = rasterFromXYZ(cbind(x, zL))

writeRaster(cmr_lower, filename=paste0(output_path, "gridded_lower_total.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))


# Upper
zU <- as.matrix(resultU)
cmr_upper = rasterFromXYZ(cbind(x, zU))

writeRaster(cmr_upper, filename=paste0(output_path, "gridded_upper_total.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))

# Coefficient of Variation 
zcv <- as.matrix(resultCV)
cmr_cv = rasterFromXYZ(cbind(x, zcv))

writeRaster(cmr_cv, filename=paste0(output_path, "/gridded_cv_total.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))





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

# Specify specific admin shapefiles paths
national_path <- paste0(data_path, "Input_Settlement_Boundaries/National")
region_path <- paste0(data_path, "Input_Settlement_Boundaries/Regional")
dept_path <- paste0(data_path, "Input_Settlement_Boundaries/Departmental")
sdept_path <- paste0(data_path, "Input_Settlement_Boundaries/Subdivisional")

# load shapefiles
adm0 <- st_read(paste0(national_path,"/CMR_Boundary.shp"))
adm1 <- st_read(paste0(region_path, "/Region_SHP.shp"))
adm2 <- st_read(paste0(dept_path, "/Departement_SHP.shp"))
adm3 <- st_read(paste0(sdept_path, "/arrondissement_shp.shp"))


# read in the saved raster files
mean_tot <- raster(paste0(output_path, "/gridded_mean_total.tif"))
mean_low <- raster(paste0(output_path, "/gridded_lower_total.tif"))
mean_upp <- raster(paste0(output_path, "/gridded_upper_total.tif"))
mean_cv <- raster(paste0(output_path, "/gridded_cv_total.tif"))


# Visualise the raster data
adm0tr = st_transform(adm0, coords = c('lon', 'lat'), crs=4326) # transform to lonlat
adm0s <- as(st_geometry(adm0tr), "Spatial") # convert to spatial object
plot(adm0s, col="white")
plot(mean_tot, add=T) # add the raster of mean population counts


# Zoom in to Yaounde and Doula
adm3tr = st_transform(adm3, coords = c('lon', 'lat'), crs=4326)

# Yaounde
#library(raster)
dim(shp_ynde <- adm3tr %>% filter(libelle %in% paste0("Yaoundé ", 5)))
ext_ynde <- extent(11.52, 11.56, 3.84, 3.88)

# mean
plot(mean_tot, ext = extent(ext_ynde),
     legend=T, col = viridis(100),  asp = NA, cex.axis=1.4)

# cv
plot(mean_cv, ext = extent(ext_ynde),
     legend=T, col = magma(100),  asp = NA, cex.axis=1.4)



## Doula
dim(shp_dou <- adm3tr %>% filter(libelle %in% paste0("Douala ", 3)))
ext_dou <- extent(9.72, 9.76, 4.02, 4.06)

# mean
plot(mean_tot, ext = extent(ext_dou),
     legend=T, col = viridis(100),  asp = NA, cex.axis=1.4)
# cv
plot(mean_cv, ext = extent(ext_dou),
     legend=T, col = magma(100),  asp = NA, cex.axis=1.4)

##

#---------------------------------------------------------------
## Mapping Administrative totals
#----------------------------------------------------------------

# Regional

# Extract Observed count and add to regional shapefile


#tmap_mode("view")
#interactive1 <- tm_basemap("OpenStreetMap")+
  ##tm_basemap("OpenStreetMap")+
 # tm_shape(adm1)+
  #tm_fill(col="gray")+
 # tm_borders(col="red", lwd=2)+
 # tm_layout(legend.outside = TRUE)+
  # tm_compass(position = c("left", "bottom"))+
  #tm_scale_bar(position = c("left", "bottom"))+
  #tm_layout(main.title = "CMR")

#interactive1

#tmap_save(interactive1, "interactive1.html")
#tmap_mode("plot")


datt <- sim.dens1$est_data
#--Observed regional total
names(dat)
reg_obs <- dat %>% dplyr::group_by(region) %>%
  dplyr::summarise(count=sum(pop, na.rm=T))

#--Predicted regional total
reg_pred <- datt %>% dplyr::group_by(CMR_Regions) %>%
  dplyr::summarise(count=sum(round(mean_pop_hat), na.rm=T))

adm1a <- adm1
adm1a$obs <- reg_obs$count[as.numeric(adm1a$id)]
adm1a$pred <-regional.est$total[as.numeric(adm1a$id)]
adm1a$predL <-regional.est$lower[as.numeric(adm1a$id)]
adm1a$predU <-regional.est$upper[as.numeric(adm1a$id)]
adm1a$uncert <-regional.est$uncertainty[as.numeric(adm1a$id)]


# Map regional values 
# Mean
tmap_mode("plot")
rpred <- tm_shape(adm1a)+ 
  tm_polygons(col="pred", title="Predicted Count",
              style="cont", n=8, palette="Blues",
              legend.show = F,
              breaks = c(700000, 1000000, 2000000,
                         3000000, 5000000, 7000000))+
  tm_borders(col="black", lwd=2)+
  tm_layout(frame=F)
rpred


# lower
rpred.lower <- tm_shape(adm1a)+ 
  tm_polygons(col="predL", title="Predicted Count",
              style="cont", n=8, palette="Blues",
              legend.show = F,
              breaks = c(700000, 1000000, 2000000,
                         3000000, 5000000, 7000000))+
  tm_borders(col="black", lwd=2)+
  tm_layout(frame=F)
rpred.lower

# upper
rpred.upper <- tm_shape(adm1a)+ 
  tm_polygons(col="predU", title="Predicted Count",
              style="cont", n=8, palette="Blues",
              legend.show = F,
              breaks = c(700000, 1000000, 2000000,
                         3000000, 5000000, 7000000))+
  tm_borders(col="black", lwd=2)+
  tm_layout(frame=F)
rpred.upper

t_reg <-tmap_arrange(rpred, rpred.lower, 
                     rpred.upper,nrow=1)
#------------------------------------------
# -- Divisional maps
#-----------------------------------------
dept_obs <- dat %>% group_by(Department) %>%
  summarise(count=sum(pop, na.rm=T))

adm2a <- adm2
adm2a$obs <- dept_obs$count[as.numeric(adm2a$id)]

# ---------------------------------------------------
adm2a$pred <-divisional.est$total[as.numeric(adm2a$id)]
adm2a$predL <-divisional.est$lower[as.numeric(adm2a$id)]
adm2a$predU <-divisional.est$upper[as.numeric(adm2a$id)]
adm2a$uncert <-divisional.est$uncertainty[as.numeric(adm2a$id)]

# Mean
tmap_mode("plot")
dpred <- tm_shape(adm2a)+ 
  tm_polygons(col="pred", title="Predicted Count",
              style="cont", n=6, palette=plasma(100),
              legend.show =F,
              breaks = c(500000, 1000000, 1500000,
                         3000000, 4000000, 5000000)
              )+
  tm_borders(col="black", lwd=2)+
  tm_layout(frame=F)
dpred


# lower
dpred.lower <- tm_shape(adm2a)+ 
  tm_polygons(col="predL", title="Predicted Count",
              style="cont", n=6, palette=plasma(100),
              legend.show = F,
              breaks = c(500000, 1000000, 1500000,
                         3000000, 4000000, 5000000)
              )+
  tm_borders(col="black", lwd=2)+
  tm_layout(frame=F)
dpred.lower

# upper
dpred.upper <- tm_shape(adm2a)+ 
  tm_polygons(col="predU", title="Predicted Count",
              style="cont", n=6, palette=plasma(100),
              legend.show = F,
              breaks = c(500000, 1000000, 1500000,
                         3000000, 4000000, 5000000)
              )+
  tm_borders(col="black", lwd=2)+
  tm_layout(frame=F)
dpred.upper

t_dep <-tmap_arrange(dpred, dpred.lower, 
                     dpred.upper,nrow=1)



# Compare the projected regional estimates with the modelled
library(ggplot2)
library(scales)
theme_set(theme_classic())

# prep data
df <- read.csv(paste0(output_path, "regional_final_with_projections_07_06_23.csv"))
df <- df[, c("names", "NIS.Projected", "total")]
names(df)
df$NIS.Projected<- df$NIS.Projected/1000
df$total <- df$total/1000


# Projected 
df1 <- df[,c("names", "NIS.Projected")]
df1$Method <- rep("Projected", nrow(df1))
colnames(df1) <- c("Region", "Estimate", "Method")


# Modelled 
df2 <- df[,c("names", "total")]
df2$Method <- rep("Modelled", nrow(df2))
colnames(df2) <- c("Region", "Estimate", "Method")



# combined
df3 <- rbind(df1, df2)
df3$Method <- factor(df3$Method, levels = c("Projected", "Modelled"))
# Plot the individuals
# Use a consistent y range
ymax <- max(df3$Estimate)
ymin <- min(df3$Estimate)
plot_slp <- ggplot(df3, aes(x=Method, y=Estimate, colour=Region, group=Region)) +
  geom_line(size = 1) + geom_point(shape=21, fill="white", size=2) + 
  geom_vline(xintercept=c("Projected", "Modelled"), 
             linetype='dashed', color=c('grey', 'grey'))+
  scale_y_continuous(breaks=seq(500,7000,500))
#ylim(ymin,ymax)

rslp <-  ggpar(plot_slp, xlab="Method", ylab="Estimated Counts ('000)",
               legend = "right", 
               legend.title = "Region",size=20,
               font.legend=c(12),
               font.label = list(size = 12, face = "bold", color ="red"),
               font.x = c(12),
               font.y = c(12),
               font.main=c(14),
               font.xtickslab =c(12),
               font.ytickslab =c(12),
               # orientation = "reverse",
               #xtickslab.rt = 45, 
               ytickslab.rt = 45)
rslp

save.image(paste0(output_path, "cmr_application.Rdata"))
