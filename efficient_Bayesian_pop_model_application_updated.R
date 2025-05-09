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
if(length(setdiff(packages, rownames(installed.packages()))) > 0) { 
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

# Load and explore data       
githublink <- "https://raw.github.com/wpgp/Efficient-Population-Modelling-using-INLA-SPDE/main/CMR_input_data.RData"
load(url(githublink))
out_path <- tempdir()# please specify your output folder here

ls() # view the loaded files: Only need - dat, shp, df, div_names & reg_names 
# dat: the input demographic data with observed population counts and covariates
str(dat)
#  shp: the associated shapefiles fo the combined household listing datasets
str(shp)
# df: contains NIS 2022 regional projections and predicted counts - for post analysis
str(df)
# div_names & reg_names are required for calculating agggregated totals
str(div_names); str(reg_names)



# visualise the imputed household size
plot(shp["I_LHHSI"]) 
names(dat); table(dat$Dat_Typ); dat$Total_Building_Count
dim(shp); dim(dat) 


#  visualise data sources using tmap 
library(tmap)
library(RColorBrewer)
library(tmaptools)
tmap_options(check.and.fix = TRUE)
tm_shape(shp) +
  tm_polygons("O_LHHSI", palette = "RdYlBu",
              legend.show = F) +
  tm_borders() +
  tm_fill("O_LHHSI",
          palette = get_brewer_pal("YlGnBu"),
          
          title = "Price in US$ per lb",
          style = "order")+
  tm_scale_bar(position = c("right","bottom"), 
               text.size=2)+
  tm_compass(position = c("left","top"), 
               text.size=1)+
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


names(dat)

#-----------------------------------------------
# add data source random effect
#------------------------------------------
dat$source <- factor(dat$Dat_Typ)

# Extract model covariates and variables for for data stack 
covars_est <- dat[,c("x2", "x3", "x17", "x21", "x32", "x34", "x40",
                     "x42", "set_reg", "set_typ", "region", "source","IDsr")]; dim(covars_est)

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
  f(IDsr, model='iid') + f(set_typ, model='iid') + f(source, model='iid')

mod1 <-inla(form1, #the formula
            data=inla.stack.data(stk_est,spde=spde),  #the data stack
            family= 'gamma',   #which family the data comes from
            control.predictor=list(A=inla.stack.A(stk_est),compute=TRUE),  #compute gives you the marginals of the linear predictor
            control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
            verbose = FALSE) #can include verbose=TRUE to see the log of the model runs


# Model 2: covs + Spatial autocor + settlement type  + region random effects
form2 <- y ~ -1 + Intercept + x2 + x3 + x17 + x21 + x32 + x34 + x40 + x42 + 
  f(IDsr, model='iid')  + f(region, model='iid') + f(spatial.field, model=spde) + f(source, model='iid')

mod2<-inla(form2, #the formula
           data=inla.stack.data(stk_est,spde=spde),  #the data stack
           family= 'gamma',   #which family the data comes from
           control.predictor=list(A=inla.stack.A(stk_est),compute=TRUE),  #compute gives you the marginals of the linear predictor
           control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
           verbose = FALSE) #can include verbose=TRUE to see the log of the model runs



# Model 3 : covs + Spatial autocor + settlement type - region interactions effects
form3 <- y ~  -1 + Intercept + x2 + x3 + x17 + x21 + x32 + x34 + x40 + x42 + 
  f(IDsr, model='iid') + f(set_reg, model='iid') + f(spatial.field, model=spde) + f(source, model='iid')

mod3 <-inla(form3, #the formula
            data=inla.stack.data(stk_est,spde=spde),  #the data stack
            family= 'gamma',   #which family the data comes from
            control.predictor=list(A=inla.stack.A(stk_est),compute=TRUE),  #compute gives you the marginals of the linear predictor
            control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
            verbose = FALSE) #can include verbose=TRUE to see the log of the model runs


# Model 4 : covs + Spatial autocor + settlement type - region interactions effects + region
form4 <- y ~  -1 + Intercept  + x2 + x3 + x17 + x21 + x32 + x34 + x40 + x42 + 
  f(set_reg, model='iid') + f(spatial.field, model=spde)+ f(set_typ, model='iid')+ f(IDsr, model='iid') + f(source, model='iid')

mod4 <-inla(form4, #the formula
            data=inla.stack.data(stk_est,spde=spde),  #the data stack
            family= 'gamma',   #which family the data comes from
            control.predictor=list(A=inla.stack.A(stk_est),compute=TRUE),  #compute gives you the marginals of the linear predictor
            control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
            verbose = FALSE) #can include verbose=TRUE to see the log of the model runs


## Run Model fit checks Based on DIC
(DIC <- t(c(mod1 = mod1$dic$dic, mod2 = mod2$dic$dic,
           mod3 = mod3$dic$dic, mod4 = mod4$dic$dic)))
#    mod1     mod2     mod3     mod4
# 1979.103 1952.221 1955.517 1928.952


# non- spatial

# fit the best fit model without spatial random effect
form4b <- dens ~  1 + x2 + x3 + x17 + x21 + x32 + x34 + x40 + x42 + 
  f(set_reg, model='iid')+ f(set_typ, model='iid')+ f(IDsr, model='iid') #+ f(source, model='iid')

mod4b <-inla(form4b, #the formula
             data=dat,  #the data stack
             family= 'gamma',   #which family the data comes from
             control.predictor=list(compute=TRUE),  #compute gives you the marginals of the linear predictor
             control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs

# extract the fixed effects  for model 4
fixed.effects <- round(mod4$summary.fixed,4) # to 4 dps
fixed.effects

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


#------------------------------------------------------------
# Cross-validation for spatial model
#-----------------------------------------------------------
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
    
    

    
    
    ##--------
    # mean --------------------------------------------------------
    # spatial random effects based on the full data best model
    sfield_nodes_mean <- mod$summary.random$spatial.field['mean']
    field_mean <- (A%*% as.data.frame(sfield_nodes_mean)[, 1])
    
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
      rnorm(nrow(test), 0, 1/mod$summary.hyperpar$mean[7]) + #---data source random effect
      
      mod$summary.random$IDsr['mean'][test_ind,1] + #--uncorrelated spatial random effect
    
      field_mean[test_ind,1]
    
       dens_ht <- exp(fixed)
       sum(pop_ht <- dens_ht*test$bld)
    
    
    
    #### LOWER---------------------------------
      # spatial random effects based on the full data best model
       sfield_nodes_lower <- mod$summary.random$spatial.field['0.025quant']
       field_lower <- (A%*% as.data.frame(sfield_nodes_lower)[, 1])
       
    fixedL <-  
      mod$summary.fixed['Intercept', '0.025quant'] +
      mod$summary.fixed['x2', '0.025quant'] * test[,'x2'] +
      mod$summary.fixed['x3', '0.025quant'] * test[,'x3'] +
      mod$summary.fixed['x17', '0.025quant'] * test[,'x17'] +
      mod$summary.fixed['x21', '0.025quant'] * test[,'x21'] +
      mod$summary.fixed['x32', '0.025quant'] * test[,'x32'] +
      mod$summary.fixed['x34', '0.025quant'] * test[,'x34'] + 
      mod$summary.fixed['x40', '0.025quant'] * test[,'x40']  +
      mod$summary.fixed['x42', '0.025quant'] * test[,'x42'] +
      
      
      rnorm(nrow(test), 0, 1/mod$summary.hyperpar$'0.025quant'[2]) + #---settlement type and region nested effects
      rnorm(nrow(test), 0, 1/mod$summary.hyperpar$'0.025quant'[5]) + #---settlement type random effect
      rnorm(nrow(test), 0, 1/mod$summary.hyperpar$'0.025quant'[7]) + #---data source random effect
      
      mod$summary.random$IDsr['0.025quant'][test_ind,1] + #--uncorrelated spatial random effects
    
      field_lower[test_ind,1]
       
      dens_htL <- exp(fixedL)
      pop_htL <- dens_htL*test$bld
    
    
    
    ### upper ---------------------------------------------
      
      # spatial random effects based on the full data best model
      sfield_nodes_upper <- mod$summary.random$spatial.field['0.975quant']
      field_upper <- (A%*% as.data.frame(sfield_nodes_upper)[, 1])
      
      fixedU <-  
        mod$summary.fixed['Intercept', '0.975quant'] +
        mod$summary.fixed['x2', '0.975quant'] * test[,'x2'] +
        mod$summary.fixed['x3', '0.975quant'] * test[,'x3'] +
        mod$summary.fixed['x17', '0.975quant'] * test[,'x17'] +
        mod$summary.fixed['x21', '0.975quant'] * test[,'x21'] +
        mod$summary.fixed['x32', '0.975quant'] * test[,'x32'] +
        mod$summary.fixed['x34', '0.975quant'] * test[,'x34'] + 
        mod$summary.fixed['x40', '0.975quant'] * test[,'x40']  +
        mod$summary.fixed['x42', '0.975quant'] * test[,'x42'] +
        
        
        rnorm(nrow(test), 0, 1/mod$summary.hyperpar$'0.975quant'[2]) + #---settlement type and region nested effects
        rnorm(nrow(test), 0, 1/mod$summary.hyperpar$'0.975quant'[5]) + #---settlement type random effect
        rnorm(nrow(test), 0, 1/mod$summary.hyperpar$'0.975quant'[7]) + #---data source random effect
        
        mod$summary.random$IDsr['0.975quant'][test_ind,1] + #--uncorrelated spatial random effects

      field_upper[test_ind,1]
      
      dens_htU <- exp(fixedU)
      pop_htU <- dens_htU*test$bld
    
    
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
    met_in <- mod_metrics2(test$pop,  
                           pop_ht)
    
    met_list_in[[i]]<- unlist(met_in)
    pred_list_in[[i]] <- data.frame(obs = test$obs, pred = pop_ht,
                                    lower = pop_htL,
                                    upper = pop_htU,
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
                             "x42", "set_reg", "set_typ", "region", "source","IDsr")]; dim(covars_train)
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
    
    
    ##--------
    # mean --------------------------------------------------------
    # spatial random effects based on the full data best model
    sfield_nodes_mean <- model$summary.random$spatial.field['mean']
    field_mean <- (A%*% as.data.frame(sfield_nodes_mean)[, 1])
    
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
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$mean[7]) + #---data source random effect
      
      mod$summary.random$IDsr['mean'][test_ind,1] + #--uncorrelated spatial random effects

    field_mean[test_ind,1]
    
    dens_ht <- exp(fixed)
    sum(pop_ht <- dens_ht*test$bld)
    
    
    
    #### LOWER---------------------------------
    
    # spatial random effects based on the full data best model
    sfield_nodes_lower <- model$summary.random$spatial.field['0.025quant']
    field_lower <- (A%*% as.data.frame(sfield_nodes_lower)[, 1])
    
    fixedL <-  
      model$summary.fixed['Intercept', '0.025quant'] +
      model$summary.fixed['x2', '0.025quant'] * test[,'x2'] +
      model$summary.fixed['x3', '0.025quant'] * test[,'x3'] +
      model$summary.fixed['x17', '0.025quant'] * test[,'x17'] +
      model$summary.fixed['x21', '0.025quant'] * test[,'x21'] +
      model$summary.fixed['x32', '0.025quant'] * test[,'x32'] +
      model$summary.fixed['x34', '0.025quant'] * test[,'x34'] + 
      model$summary.fixed['x40', '0.025quant'] * test[,'x40']  +
      model$summary.fixed['x42', '0.025quant'] * test[,'x42'] +
      
      
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$'0.025quant'[2]) + #---settlement type and region nested effects
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$'0.025quant'[5]) + #---settlement type random effect
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$'0.025quant'[7]) + #---data source random effect
      
      model$summary.random$IDsr['0.025quant'][test_ind,1] + #--uncorrelated spatial random effects
 
    field_lower[test_ind,1]
    
    dens_htL <- exp(fixedL)
    pop_htL <- dens_htL*test$bld
    
    
    
    ### upper ---------------------------------------------
    
    # spatial random effects based on the full data best model
    sfield_nodes_upper <- model$summary.random$spatial.field['0.975quant']
    field_upper <- (A%*% as.data.frame(sfield_nodes_upper)[, 1])
    
    
    fixedU <-  
      model$summary.fixed['Intercept', '0.975quant'] +
      model$summary.fixed['x2', '0.975quant'] * test[,'x2'] +
      model$summary.fixed['x3', '0.975quant'] * test[,'x3'] +
      model$summary.fixed['x17', '0.975quant'] * test[,'x17'] +
      model$summary.fixed['x21', '0.975quant'] * test[,'x21'] +
      model$summary.fixed['x32', '0.975quant'] * test[,'x32'] +
      model$summary.fixed['x34', '0.975quant'] * test[,'x34'] + 
      model$summary.fixed['x40', '0.975quant'] * test[,'x40']  +
      model$summary.fixed['x42', '0.975quant'] * test[,'x42'] +
      
      
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$'0.975quant'[2]) + #---settlement type and region nested effects
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$'0.975quant'[5]) + #---settlement type random effect
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$'0.975quant'[7]) + #---data source random effect
      
      model$summary.random$IDsr['0.975quant'][test_ind,1] + #--uncorrelated spatial random effects

    field_upper[test_ind,1]
    
    dens_htU <- exp(fixedU)
    pop_htU <- dens_htU*test$bld
    
    
    
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
    met_out <- mod_metrics2(test$pop,  
                            pop_ht)
    
    met_list_out[[i]]<- unlist(met_out)
    pred_list_out[[i]] <- data.frame(obs = test$obs, pred = pop_ht,
                                     lower = pop_htL,
                                     upper = pop_htU,
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


#------------------------------------------------------------
# Cross-validation for non-spatial model
#-----------------------------------------------------------

cross_validate_ns <- function(dat, n.folds, mod, form, seed)
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
  
  
  k_uniq <-sort(unique(dat$k_fold))
  
  # mod = mod4b
  
  #---------------------------------------------------------------
  #                   in-sample
  #---------------------------------------------------------------
  
  met_list_in <- list()
  pred_list_in <- list()
  for(i in 1:length(k_uniq))
  {
    # i=2
    print(paste0("in-sample cross-validation using fold ", i, sep=""))
    test_ind <- which(dat$k_fold==k_uniq[i])
    dim(test <- dat[test_ind, ]) #---test set for fold i
    
    
    train_coords <- coords
    test_coords <- coords[test_ind,]
    
    
    ##--------
    # mean -----------------------------------------------------
    fixed <-  
      mod$summary.fixed['(Intercept)', 'mean'] * 1 +
      mod$summary.fixed['x2', 'mean'] * test[,'x2'] +
      mod$summary.fixed['x3', 'mean'] * test[,'x3'] +
      mod$summary.fixed['x17', 'mean'] * test[,'x17'] +
      mod$summary.fixed['x21', 'mean'] * test[,'x21'] +
      mod$summary.fixed['x32', 'mean'] * test[,'x32'] +
      mod$summary.fixed['x34', 'mean'] * test[,'x34'] + 
      mod$summary.fixed['x40', 'mean'] * test[,'x40']  +
      mod$summary.fixed['x42', 'mean'] * test[,'x42'] +
      
      
      rnorm(nrow(test), 0, 1/mod$summary.hyperpar$mean[2]) + #---settlement type and region nested effects
      rnorm(nrow(test), 0, 1/mod$summary.hyperpar$mean[3]) + #---settlement type random effect
      #rnorm(nrow(test), 0, 1/mod$summary.hyperpar$mean[5]) + #---data source random effect
      
      mod$summary.random$IDsr$mean[test_ind] 
    
    dens_ht <- exp(fixed)
    sum(pop_ht <- dens_ht*test$bld)
    
    
    
    # lower------------------------------------------------
    fixedL <-  
      mod$summary.fixed['(Intercept)', '0.025quant'] * 1 +
      mod$summary.fixed['x2', '0.025quant'] * test[,'x2'] +
      mod$summary.fixed['x3', '0.025quant'] * test[,'x3'] +
      mod$summary.fixed['x17', '0.025quant'] * test[,'x17'] +
      mod$summary.fixed['x21', '0.025quant'] * test[,'x21'] +
      mod$summary.fixed['x32', '0.025quant'] * test[,'x32'] +
      mod$summary.fixed['x34', '0.025quant'] * test[,'x34'] + 
      mod$summary.fixed['x40', '0.025quant'] * test[,'x40']  +
      mod$summary.fixed['x42', '0.025quant'] * test[,'x42'] +
      
      
      rnorm(nrow(test), 0, 1/mod$summary.hyperpar$'0.025quant'[2]) + #---settlement type and region nested effects
      rnorm(nrow(test), 0, 1/mod$summary.hyperpar$'0.025quant'[3]) + #---settlement type random effect
  
      mod$summary.random$IDsr$'0.025quant'[test_ind] 
    
    dens_htL <- exp(fixedL)
    sum(pop_htL <- dens_htL*test$bld)
    
    
    # upper -----------------------------------------------------------
    fixedU <-  
      mod$summary.fixed['(Intercept)', '0.975quant'] * 1 +
      mod$summary.fixed['x2', '0.975quant'] * test[,'x2'] +
      mod$summary.fixed['x3', '0.975quant'] * test[,'x3'] +
      mod$summary.fixed['x17', '0.975quant'] * test[,'x17'] +
      mod$summary.fixed['x21', '0.975quant'] * test[,'x21'] +
      mod$summary.fixed['x32', '0.975quant'] * test[,'x32'] +
      mod$summary.fixed['x34', '0.975quant'] * test[,'x34'] + 
      mod$summary.fixed['x40', '0.975quant'] * test[,'x40']  +
      mod$summary.fixed['x42', '0.975quant'] * test[,'x42'] +
      
      
      rnorm(nrow(test), 0, 1/mod$summary.hyperpar$'0.975quant'[2]) + #---settlement type and region nested effects
      rnorm(nrow(test), 0, 1/mod$summary.hyperpar$'0.975quant'[3]) + #---settlement type random effect
      mod$summary.random$IDsr$'0.975quant'[test_ind] 
    
    dens_htU <- exp(fixedU)
    sum(pop_htU <- dens_htU*test$bld)
    
    # visualise samples

    
    par(mfrow =c(1,1))
    plot(test$obs, pop_ht, xlab = "Observed", 
         ylab = "Predicted", col=c('blue','orange'),
         pch=c(16,16), cex.axis=1.5)
    abline(0,1)
    legend("topleft", c("Observed", "Predicted"), col=c("blue", "orange"), pch=c(16,16),
           bty="n", cex=1.5) 
    
    
    
    # calculate fit metrics
    met_in <- mod_metrics2(test$pop,  
                           pop_ht)
    
    met_list_in[[i]]<- unlist(met_in)
    pred_list_in[[i]] <- data.frame(obs = test$obs, pred = pop_ht,
                                    lower = pop_htL,
                                    upper = pop_htU,
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
    
    # i =2
    print(paste0("out-of-sample cross-validation using fold ", i, sep=""))
    train_ind <- which(dat$k_fold!=k_uniq[i])
    test_ind <- which(dat$k_fold==k_uniq[i])
    dim(train <- dat[train_ind, ])#---train set for fold i
    dim(test <- dat[test_ind, ]) #---test set for fold i
    
    # form = form4b
    ###---Rerun INLA for model test prediction
    model <-inla(form, #the formula
                 data=dat,  #the data stack
                 family= 'gamma',   #which family the data comes from
                 control.predictor=list(compute=TRUE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
    summary(model)
    
    
    # mean -----------------------------------------------------
    fixed <-  
      model$summary.fixed['(Intercept)', 'mean'] * 1 +
      model$summary.fixed['x2', 'mean'] * test[,'x2'] +
      model$summary.fixed['x3', 'mean'] * test[,'x3'] +
      model$summary.fixed['x17', 'mean'] * test[,'x17'] +
      model$summary.fixed['x21', 'mean'] * test[,'x21'] +
      model$summary.fixed['x32', 'mean'] * test[,'x32'] +
      model$summary.fixed['x34', 'mean'] * test[,'x34'] + 
      model$summary.fixed['x40', 'mean'] * test[,'x40']  +
      model$summary.fixed['x42', 'mean'] * test[,'x42'] +
      
      
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$mean[2]) + #---settlement type and region nested effects
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$mean[3]) + #---settlement type random effect
      #rnorm(nrow(test), 0, 1/model$summary.hyperpar$mean[5]) + #---data source random effect
      
      model$summary.random$IDsr$mean[test_ind] 
    
    dens_ht <- exp(fixed)
    sum(pop_ht <- dens_ht*test$bld)
    
    
    
    # lower------------------------------------------------
    fixedL <-  
      model$summary.fixed['(Intercept)', '0.025quant'] * 1 +
      model$summary.fixed['x2', '0.025quant'] * test[,'x2'] +
      model$summary.fixed['x3', '0.025quant'] * test[,'x3'] +
      model$summary.fixed['x17', '0.025quant'] * test[,'x17'] +
      model$summary.fixed['x21', '0.025quant'] * test[,'x21'] +
      model$summary.fixed['x32', '0.025quant'] * test[,'x32'] +
      model$summary.fixed['x34', '0.025quant'] * test[,'x34'] + 
      model$summary.fixed['x40', '0.025quant'] * test[,'x40']  +
      model$summary.fixed['x42', '0.025quant'] * test[,'x42'] +
      
      
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$'0.025quant'[2]) + #---settlement type and region nested effects
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$'0.025quant'[3]) + #---settlement type random effect
      
      model$summary.random$IDsr$'0.025quant'[test_ind] 
    
    dens_htL <- exp(fixedL)
    pop_htL <- dens_htL*test$bld
    
    
    # upper -----------------------------------------------------------
    fixedU <-  
      model$summary.fixed['(Intercept)', '0.975quant'] * 1 +
      model$summary.fixed['x2', '0.975quant'] * test[,'x2'] +
      model$summary.fixed['x3', '0.975quant'] * test[,'x3'] +
      model$summary.fixed['x17', '0.975quant'] * test[,'x17'] +
      model$summary.fixed['x21', '0.975quant'] * test[,'x21'] +
      model$summary.fixed['x32', '0.975quant'] * test[,'x32'] +
      model$summary.fixed['x34', '0.975quant'] * test[,'x34'] + 
      model$summary.fixed['x40', '0.975quant'] * test[,'x40']  +
      model$summary.fixed['x42', '0.975quant'] * test[,'x42'] +
      
      
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$'0.975quant'[2]) + #---settlement type and region nested effects
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$'0.975quant'[3]) + #---settlement type random effect
      model$summary.random$IDsr$'0.975quant'[test_ind] 
    
    dens_htU <- exp(fixedU)
    pop_htU <- dens_htU*test$bld
    
    
    
    par(mfrow =c(1,1))
    plot(test$obs, pop_ht, xlab = "Observed", 
         ylab = "Predicted", col=c('blue','orange'),
         pch=c(16,16), cex.axis=1.5)
    abline(0,1)
    legend("topleft", c("Observed", "Predicted"), col=c("blue", "orange"), pch=c(16,16),
           bty="n", cex=1.5) 
    
    
    
    # calculate fit metrics
    met_out <- mod_metrics2(test$pop,  
                            pop_ht)
    
    met_list_out[[i]]<- unlist(met_out)
    pred_list_out[[i]] <- data.frame(obs = test$pop, pred = pop_ht,
                                     lower = pop_htL,
                                     upper = pop_htU,
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

(cross_valb <- cross_validate_ns(dat, n.folds = 5, 
                                 mod = mod4b, 
                                 form = form4b, 
                                 seed = 13235))

cross_valb$met_list_in_dat  # in-sample metrics per fold 
cross_valb$met_list_out_dat  # out-of-sample metrics per fold
cross_valb$cv_metrics    # combined averaged metrics
cross_valb$pred_dat  # combined prediction data-


# combined metrics
(com_metrics <- rbind(cross_val$cv_metrics,cross_valb$cv_metrics))

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
                       #facet.by = "data",
                       conf.int = TRUE,                          # Add confidence interval
                       color = "data", 
                       palette = "lancet"#,           # Color by groups "cyl"
                       #shape = "Fold"                             # Change point shape by groups "cyl"
) 
pred.data <- pred.data %>% filter(pred<5000 & obs < 5000) %>% 
  filter(obs!=0 & pred!=0)# remove outliers to get nice plots

plot_cval <-ggplot(pred.data, aes(x=obs, y=pred))+
  geom_point() +
  geom_smooth(method="lm") +
  theme_bw()+
  theme(strip.text = element_text(size=20))+
  facet_wrap(~data)

rcval <-  ggpar(plot_cval, xlab="Observed counts", ylab="Predicted Counts",
                legend = "top", 
                legend.title = "Fold (k{=5}-fold)",size=18,
                font.legend=c(18),
                font.label = list(size = 18, face = "bold", color ="red"),
                font.x = c(18),
                font.y = c(18),
                font.main=c(18),
                font.xtickslab =c(18),
                font.ytickslab =c(18),
                # orientation = "reverse",
                xtickslab.rt = 45, ytickslab.rt = 45)
rcval



##
plot_cval4 <- pred.data %>%
    ggplot( aes(x=Fold, y=pred, fill=Fold)) +
    geom_violin(width=1.4) +
    geom_boxplot(width=0.4, color="grey", alpha=0.2,
                 # custom boxes
                 color="blue",
                 fill="magenta",
                 alpha=0.2,
                 
                 # Notch?
                 notch=TRUE,
                 notchwidth = 0.8) +
    scale_fill_viridis(discrete = TRUE) +
    facet_wrap(~data) +
  theme_bw()+
  theme(strip.text = element_text(size=20))

rcval4 <-  ggpar(plot_cval4, xlab="Data Fold", ylab="Predicted Counts",
                 legend = "none", 
                 #legend.title = "Dataset",size=20,
                 legend.title = list(color = "Dataset", 
                                     linetype = "data", shape = "Fold"),
                 font.legend=c(18),
                 font.label = list(size = 18, face = "bold", color ="red"),
                 font.x = c(18),
                 palette = "lancet",
                 font.y = c(18),
                 font.main=c(18),
                 font.xtickslab =c(18),
                 font.ytickslab =c(18),
                 # orientation = "reverse",
                 xtickslab.rt = 45, ytickslab.rt = 45)
rcval4

ggarrange(rcval4, rcval,
          nrow=2, ncol=1)



#---------------------------------------------------------
## For the base method (non-spatial)
#-------------------------------------------------------
##
# remove outliers to get nice plots
pred.datab <- cross_valb$pred_dat
pred.datab$Fold <- factor(pred.datab$fold)
pred.datab$data <- factor(pred.datab$data,
                         levels = c("insample", "outsample"),
                         labels = c("In-Sample", "Out-of-Sample"))


pred.data2b <- pred.datab %>% dplyr::filter(pred<5000)
plot_cval4b <- ggboxplot(pred.data2b, y = "pred",
                        #add = "reg.line",                         # Add regression line
                        facet.by = "data",
                        #conf.int = TRUE,                          # Add confidence interval
                        color = "Fold", 
                        palette = "lancet",        # Color by groups "cyl"
                        shape = "Fold"                             # Change point shape by groups "cyl"
) 

plot_cval4b <- pred.data2b %>% dplyr::filter(pred<5000) %>%
  ggplot( aes(x=Fold, y=pred, fill=Fold)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.4, color="grey", alpha=0.2,
               # custom boxes
               color="blue",
               fill="magenta",
               alpha=0.2,
               
               # Notch?
               notch=TRUE,
               notchwidth = 0.8) +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~data) +
  theme_bw()+
  theme(strip.text = element_text(size=20))


rcval4b <-  ggpar(plot_cval4b, xlab="Data Fold", ylab="Predicted Counts",
                 legend = "none", 
                 #legend.title = "Dataset",size=20,
                 legend.title = list(color = "Dataset", 
                                     linetype = "data", shape = "Fold"),
                 font.legend=c(18),
                 font.label = list(size = 18, face = "bold", color ="red"),
                 font.x = c(18),
                 palette = "jco",
                 font.y = c(18),
                 font.main=c(18),
                 font.xtickslab =c(18),
                 font.ytickslab =c(18),
                 # orientation = "reverse",
                 xtickslab.rt = 45, ytickslab.rt = 45)
rcval4b
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
# download the prediction stack "prediction_data.RDS" from Google drive here: https://drive.google.com/file/d/1290hqUnBHhQS0I3iijddj34solTG-S_S/view?usp=sharing
# ....data_path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3/Working/CMR/Chris_N/paper1/submission"
pred_data <- readRDS(paste0(data_path, "/prediction_data.RDS"))
data <- pred_data
head(data)


#      Standardize covs
covss <- c("x2", "x3", "x17", "x21", "x32", "x34", "x40", "x42")
data[,covss] <- apply(data[,covss], 2, stdize)  ##----standardize all model covariates 

#    Prediction covs matrix
Apred <- inla.spde.make.A(mesh = mesh, loc = cbind(data$Lon, data$Lat))
dim(Apred)


#
summary(mod4)
####----mod4####-------Run the uncertainty estimates - spatial
sim_spatial <- function(model, dat, Aprediction, run)
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
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[7]) + #---data type random effects
      field_mean[,1]
    
    
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



# for non-spatial model
sim_nonspatial <- function(model, dat, run)
{
  fixedeff <- dens_hat <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
  #inla.seed = as.integer(runif(1)*.Machine$integer.max)
  inla.seed = 211123018
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

  for(i in 1:run)
  {
    fixedeff[,i] <- 
      model$summary.fixed['(Intercept)', 'mean'] +
      m1.samp[[i]]$latent[1,] * dat[,'x2']  +
      m1.samp[[i]]$latent[2,] * dat[,'x3']  +
      m1.samp[[i]]$latent[3,] * dat[,'x17'] +
      m1.samp[[i]]$latent[4,] * dat[,'x21'] +
      m1.samp[[i]]$latent[5,] * dat[,'x32'] + 
      m1.samp[[i]]$latent[6,] * dat[,'x34'] +
      m1.samp[[i]]$latent[7,] * dat[,'x40'] +
      m1.samp[[i]]$latent[8,] * dat[,'x42'] +
      
      
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[2]) + #---settlement type and region nested effects
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[3]) + #---settlement type random effect
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[4])  #---cluster level random effect
    
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
system.time(str(sim.spat <- sim_spatial(mod4,data,Apred, run))) # spatial
system.time(str(sim.nspat <- sim_nonspatial(mod4b,data,run))) # non-spatial

sum(sim.spat$est_dat$mean_pop_hat,na.rm=T)
sum(sim.nspat$est_dat$mean_pop_hat,na.rm=T)

#  Join the posterior sample to the prediction data
names(data)
data.sim <- data.frame(cbind(data[,c("CMR_Regions", "CMR_Department","CMR_Settlement_Classification")], sim.spat$pop_hat))
data.sim.ns <- data.frame(cbind(data[,c("CMR_Regions", "CMR_Department","CMR_Settlement_Classification")], sim.nspat$pop_hat))

##--------Calculate and save National total with uncertainties
nat_total <- function(dat, run)
{
  p_hat <- dat[,4:(run+3)]
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

(national2 <- nat_total(data.sim.ns, run))
(national2 <- data.frame(total= national2[1,],
                        lower = national2[2,],
                        median=national2[3,],
                        upper=national2[4,]))

rbind(national, national2)

#write.csv(national, file=paste0(results_path, "/estimates/National_estimates_final.csv"))

##---Regional estimates
regional_est <- function(datr, run)
{
  uniR <- unique(reg_names$id)
  regnames <- unique(reg_names$libelle)
  outR <- matrix(0, nrow=length(uniR), ncol=5)
  for(j in uniR)
  {
    reg <- datr[datr$CMR_Regions==j,]
    rtots <- apply(reg[,4:(3+run)], 2, sum, na.rm=T)
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
regional.est$method <- rep("Proposed", nrow(regional.est))

(regional.est2 <- regional_est(data.sim.ns, 100))
sum(regional.est2$total)
regional.est2$method <- rep("Base", nrow(regional.est2))

reg_dat <- rbind(regional.est, regional.est2) # combined regional data

#write.csv(regional.est, file=paste0(results_path, "/estimates/regional_final.csv"))

##---Divisional estimates
divisional_est <- function(datd, run)
{
  uniD <- unique(div_names$id)
  uniR <- unique(reg_names$id)
  outD <- matrix(0, nrow=length(uniD), ncol=5)
  divnames <- unique(div_names$libelle)
  regnames <- unique(reg_names$libelle)
  reg <- rep(0, length(uniD))
  for(k in uniD)
  {
    div <- datd[datd$CMR_Department==k,]
    dtots <- apply(div[,4:(3+run)], 2, sum, na.rm=T)
    regn <- unique(datd$CMR_Regions[datd$CMR_Department==k])
    reg[k] = regnames[regn]
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
                               names2 = reg,
                               total = outD[,1],
                               lower = outD[,2],
                               median = outD[,3],
                               upper = outD[,4],
                               uncertainty = outD[,5]))
}

(divisional.est <- divisional_est(data.sim, 100))# Spatial
sum(divisional.est$total)
divisional.est$method <- rep("Proposed", nrow(divisional.est))

(divisional.est2 <- divisional_est(data.sim.ns, 100)) # Non spatial
sum(divisional.est2$total)
divisional.est2$method <- rep("Proposed", nrow(divisional.est2))

dep_dat <- rbind(divisional.est,divisional.est2) 

#write.csv(divisional.est, file=paste0(results_path, "/estimates/department_final.csv"))


# Error bar plots
regional.est$total2 <- round(regional.est$total/1000)
regional.est$lower2 <- round(regional.est$lower/1000)
regional.est$upper2 <- round(regional.est$upper/1000)

pd <- ggplot(regional.est, aes(x=reorder(names, -total), 
                      y=total2)) +
  geom_errorbar(aes(ymin = lower2, ymax = upper2),
                position = position_dodge(width = 0.4), width = 0.4,
                size=0.9)+
  geom_point(size = 2, color = "red", alpha=0.9) +
  theme_bw()#+
  #coord_flip() 
  
  
pdh3 <- ggpar(pd, ylab="Estimated population count ('000)", xlab="Region",
              legend = "none", legend.title = "Region:",
              font.legend=c(14),
              # font.label = list(size = 12, face = "bold", color ="red"),
              #palette = "pnj",
              font.x = c(14),
              font.y = c(14),
              font.main=c(14),
              font.xtickslab =c(12),
               font.ytickslab =c(14),
              xtickslab.rt = 45)

pdh3
# -----------------------------------------------------------
# Prepare and save the grid vector predictions as raster files
# -----------------------------------------------------------

# specify the reference coordinates for the predictions
ref_coords <- cbind(pred_data$Lon, pred_data$Lat)
x <- as.matrix(ref_coords)

# Extract the grid cell predictions
sum(result <-sim.spat$est_data$mean_pop_hat, na.rm=T) # mean
resultL <- sim.spat$est_data$lower_pop_hat # lower bound
resultU <- sim.spat$est_data$upper_pop_hat # upper bound
resultCV <- sim.spat$est_data$cv # coefficient of variation

# add grid cell values and write the raster files

# specify your output_path for the raster file 
#  Mean

output_path <- data_path
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


# Compare the projected regional estimates with the modelled
library(ggplot2)
library(scales)
theme_set(theme_classic())

# prep data
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
  geom_line(size = 1.5) + geom_point(shape=21, fill="grey", size=3) + 
  geom_vline(xintercept=c("Projected", "Modelled"), 
             linetype='dashed', color=c('black', 'black'))+
  scale_y_continuous(breaks=seq(500,7000,500))+
  theme_minimal()
#ylim(ymin,ymax)

rslp <-  ggpar(plot_slp, xlab="Method", ylab="Estimated Counts ('000)",
               legend = "top", 
               legend.title = "Region",size=20,
               font.legend=c(12),
               palette = "jco",
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



##### Compare with barplots
 
barp <- df3 %>%
  ggplot( aes(x=Region, y=Estimate, fill=Method)) +
   #ggplot(data = df_continent, aes(x = reorder(continent, -mean_temp), y = mean_temp, fill = as.factor(Year))) +
   geom_bar(stat = "identity", color = "black", 
            position = position_dodge()) +
   scale_y_continuous(breaks=seq(500,7000,1000))+
   theme_minimal()+
   coord_flip()


rcval5 <-  ggpar(barp, xlab="Region", ylab="Population Estimate('000)",
                 legend = "right", 
                 #legend.title = "Dataset",size=20,
                 legend.title = list(color = "Method", 
                                     linetype = "data", shape = "Fold"),
                 font.legend=c(12),
                 font.label = list(size = 12, face = "bold", color ="red"),
                 font.x = c(12),
                 palette = "lancet",
                 font.y = c(12),
                 font.main=c(14),
                 font.xtickslab =c(12),
                 font.ytickslab =c(12),
                 # orientation = "reverse",
                 xtickslab.rt = 45, ytickslab.rt = 45)
rcval5

#----------------------------------------------------------------
# Visualise the gridded data
library(terra)
#r <- raster(paste0(output_path, "gridded_mean_total.tif"))
r <- rast(paste0(output_path, "/gridded_mean_total.tif"))
rcv <- rast(paste0(output_path, "/gridded_cv_total.tif"))

####
cmr_grd.adj <-r

library(exactextractr)
r_df <- data.frame(r)
head(r_df)
summary(r_df)
#gridded_mean_total
#Min.   :  0.7773  
#1st Qu.:  4.0386  
#Median :  9.0014  
#Mean   : 21.6166  
#3rd Qu.: 21.4526  
#Max.   :798.6568 



cmr_shp <- st_read(paste0(data_path, "/CMR_Boundary.shp"))
library(basemaps)
library(ggplot2)
library(tidyterra)
library(raster)
library(ggmap)
library(mapview)
set_defaults(cmr_grd.adj, map_service = "esri", map_type = "world_imagery")
# mapping the predicted means
crs(cmr_grd.adj) = "EPSG:3857"
#plot(cmr_grd.adj)
x <- basemap_raster(cmr_grd.adj, map_service = "esri", map_type = "world_imagery")
x_terr <- rast(x) 
crs(x_terr)
x_lonlat <- project(x_terr, "EPSG:4326")

crs(x_lonlat)# check
plot(x_lonlat)

shpp <- st_transform(cmr_shp, crs=4326) # Transform
shpp <- as(st_geometry(shpp), "Spatial") # Change to spatial object
plot(shpp) # checks
crs(shpp)#
#"dark green", "orange", "white"
library(viridis)
#cmr_grd.adj <- round(r, 2)
brks <- seq(1,820.66, by=89)
arg <- list(brks, labels=brks)
zlim=c(0.78 ,798.66)
plot(cmr_grd.adj,
     legend=T, col = colorRampPalette(c("dark green", "grey", "orange"))(255),
     asp = NA, cex.axis=1.4, axes=T,#zlim=zlim,
     breaks=brks)


plot(x_lonlat, add=T, legend=F)
#plot(shpp,add=T, col="black", lwd=2)
plot(shpp,add=T, lwd=2, bg = "transparent")
plot(cmr_grd.adj,
     legend=T, col = colorRampPalette(c("dark green", "grey", "orange"))(255),
     asp = NA, 
     cex.axis=1.4, add=T, breaks = brks, zlim=zlim)

# add Douala
ext_dou <- extent(9.6, 9.9, 3.9, 4.2)
plot(ext_dou , add=TRUE, col='yellow', lwd=2)

# add Yaounde
ext_ynde <- extent(11.4, 11.7, 3.7, 4)
plot(ext_ynde , add=TRUE, col='yellow', lwd=2)

# plot Douala and Yaounde counts separately

plot(cmr_grd.adj, ext = extent(ext_dou),
     legend=T, col = colorRampPalette(c("dark green", "white", "orange"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)

#
plot(cmr_grd.adj, ext = extent(ext_ynde),
     legend=F, col = colorRampPalette(c("dark green", "white", "orange"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)


#### Zoom in a bit further
plot(cmr_grd.adj,
     legend=T, col = colorRampPalette(c("green", "grey", "orange"))(255),
     asp = NA, cex.axis=1.4, axes=F,
     gridded=T, breaks = brks, zlim=zlim)
plot(x_lonlat, add=T, legend=F)
plot(shpp,add=T, col="black", lwd=2)
plot(cmr_grd.adj,
     legend=F, col = colorRampPalette(c("green", "grey", "orange"))(255),
     asp = NA, 
     cex.axis=1.4, add=T, breaks = brks, zlim=zlim)

# add Douala
ext_dou <- extent(9.6, 9.9, 3.9, 4.2)/12
ext_dou <- extent(9.69, 9.72, 4.02, 4.05)
plot(ext_dou , add=TRUE, col='yellow', lwd=3)

# add Yaounde
ext_ynde <- extent(11.4, 11.7, 3.7, 4)/6
ext_ynde <- extent(11.52, 11.55, 3.85, 3.88)
plot(ext_ynde , add=TRUE, col='yellow', lwd=3)

# plot Douala and Yaounde counts separately

plot(cmr_grd.adj, ext = extent(ext_dou),
     legend=F, col = colorRampPalette(c("dark green", "grey", "orange"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)

#
plot(cmr_grd.adj, ext = extent(ext_ynde),
     legend=F, col = colorRampPalette(c("dark green", "white", "orange"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)


#### Mapping the coefficient of variation
# mapping the predicted means
#utm_crs <- "EPSG:3857"
#cmr_utm <- project(rcv, utm_crs)
#res(cmr_utm)
xcv <- basemap_raster(rcv, map_service = "esri", map_type = "world_imagery")
x_cv <- rast(xcv) 
#res(xcv) <- c(50,50)
x_lonlat2 <- project(x_cv, "EPSG:4326")

#crs(x_lonlat2)# check
plot(x_lonlat2)

# reproject shapefile
shpp <- st_transform(cmr_shp, crs=4326) # Transform
shpp <- as(st_geometry(shpp), "Spatial") # Change to spatial object
plot(shpp) # checks
#crs(shpp)#

# set breaks and plot limits
brks2 <- seq(0.20,0.8, by=0.1)
zlim2=c(0.20,0.8)
plot(rcv,
     legend=T, col=colorRampPalette(c("white", "orange", "red"))(255),
     asp = NA, cex.axis=1.4, axes=F, breaks = brks2, zlim=zlim2)
plot(x_lonlat2, add=T, legend=F)
plot(shpp,add=T, col="black", lwd=2)
plot(rcv,
     legend=F, col=colorRampPalette(c("white", "orange", "red"))(255),
     asp = NA,cex.axis=1.4, add=T, breaks = brks2, zlim=zlim2)


# add Douala
ext_dou <- extent(9.6, 9.9, 3.9, 4.2)
ext_dou <- extent(9.69, 9.72, 4.02, 4.05)
#ext_dou <- extent(10.7E+5, 11.0E+5, 4.35E+5, 4.65E+5)
plot(ext_dou , add=TRUE, col='yellow', lwd=2)

# add Yaounde
#ext_ynde <- extent(12.69E+5, 12.99E+5, 4.1E+5, 4.4E+5)
#plot(ext_ynde , add=TRUE, col='yellow', lwd=2)

ext_ynde <- extent(11.4, 11.7, 3.7, 4)
ext_ynde <- extent(11.52, 11.55, 3.85, 3.88)
plot(ext_ynde , add=TRUE, col='yellow', lwd=2)

# plot Douala and Yaounde counts separately

plot(rcv, ext = extent(ext_dou),
     legend=F,col=colorRampPalette(c("white", "orange", "red"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black", breaks = brks2, zlim=zlim2,
     axes=F)

#
plot(rcv, ext = extent(ext_ynde),
     legend=F, col=colorRampPalette(c("white", "orange", "red"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks2, zlim=zlim2)


#save.image(paste0(output_path, "cmr_application.Rdata"))
