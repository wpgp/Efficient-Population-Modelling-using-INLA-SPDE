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
              "lattice", "gridExtra", "devtools", "rlang", "DClusterm",
              "viridis", "tmaptools", "spdep", "ggplot2", "ggpubr",
              "psych", "knitr", "car", "MASS", "terra")


# install.packages("psych", type="binary")
#if(length(setdiff(packages, rownames(installed.packages()))) > 0) { 
  #install.packages(setdiff(packages, rownames(installed.packages(type="binary")))) }


#Install INLA
#if(length(setdiff("INLA", rownames(installed.packages()))) > 0){
#  install.packages("INLA", type="binary", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#}
library(INLA)
lapply(packages, library, character.only = TRUE) ##--access the libraries

#---------------------------------------------------------------------------------
###--Specify key file paths

# Load and explore data       
githublink <- "https://raw.github.com/wpgp/Efficient-Population-Modelling-using-INLA-SPDE/main/CMR_input_data.RData"
load(url(githublink))

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

#--------------------------------------------------------------------------------
#   1)  INPUT DATA PREPARATION--
#-------------------------------------------------------------------------------
dim(dat); dim(shp)

# Some variable recoding
dat$bldg <- dat$Total_Building_Count # observed building count
(zerob <- which(dat$bldg==0)) # check if there are samples without building counts

dat$pop <- dat$Imputed_LHHSIZE # Observed population counts
(zerop <- which(dat$pop==0)) # check if there are samples without population count

dat$dens <- dat$pop/dat$bldg #---people per building
dat$dens_hh <- dat$hh_nm_c/dat$bldg #---households per building


dat$dens[is.infinite(dat$dens)] = NA # check if there are infinite terms and set them to NA
dat$dens_hh[is.infinite(dat$dens_hh)] = NA # check if there are infinite terms and set them to NA

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

names(dat)


## some descriptive stats
# household number
sum(dat$hh_nm_c, na.rm=T)# 509,628 hhlds
summary(dat$hh_nm_c, na.rm=T)# min = 2, max = 1459; median 208, mean = 222.5 hhlds


# population count
sum(dat$Imputed_LHHSIZE, na.rm=T)# 2,587,459
summary(dat$Imputed_LHHSIZE, na.rm=T)# min = 8, max = 26,689; median 979.5, mean = 1129.9 people

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

mcovs$resp <- log(dat$dens) # # population density model
mcovs$resp2 <- log(dat$dens_hh) ## household number model


dim(mcovs <- mcovs %>% drop_na()) # Drop all NA's 

##----Covariates selection using GLM STEPWISE REGRESSION
dim(mcovs)

# select covariates for population density
fdens <- glm(resp ~., data=mcovs[,-45], family = gaussian)# exclude 'resp2'
summary(fdens)

# run stepwise regression
step_dens <- stepAIC(fdens, scale = 0,
                     direction = c("both"),
                     trace = 1, keep = NULL, steps = 1000, use.start = FALSE,
                     k = 2)
step_dens

###-----------REFIT THE SELECTED MODEL AND CHOOSE only significant variables with VIF < 5
fdens2 <-  glm(formula = resp ~ x2 + x3 + x17 + x21 + x32 + 
                 x34 + x40 + x42, family = gaussian, data = mcovs[,-45])
summary(fdens2)

# Calculate the vif 
#library(car) # required for the vif
vif1 = vif(fdens2)
vif1[(vif1 < 5)]

covs1 <- c("x2","x3","x17","x21","x32","x34","x40","x42") # for population density

### select covariates for the number of households
fdensb <- glm(resp2 ~., data=mcovs[,-44], family = gaussian)# exclude 'resp'
summary(fdensb)

# run stepwise regression
step_densb <- stepAIC(fdensb, scale = 0,
                     direction = c("both"),
                     trace = 1, keep = NULL, steps = 1000, use.start = FALSE,
                     k = 2)
step_densb

###-----------REFIT THE SELECTED MODEL AND CHOOSE only significant variables with VIF < 5
fdensb2<-  glm(formula = resp2 ~  x2 + x17 + x19 + x20 + x21  + x32 + 
                 x37 + x42, family = gaussian, 
               data = mcovs[, -44])

summary(fdensb2)

# Calculate the vif 
#library(car) # required for the vif
vifb1 = vif(fdensb2)
vifb1[(vifb1 < 5)]
##----Selected significant covariates with VIF < 5
covs2 <- c("x2","x17","x19","x20","x21","x32","x37","x42")# for number of households

# -------------------------------------------------------------
# The names of the  selected covariates are given below
# --------------------------------------------------------------
# x2:  Conflicts data (ACLED)
# x3:  Explosions data (ACLED)
# x17: Distance to water bodies
# x19: esaccilc_dst040_100m_2015
# x20: esaccilc_dst130_100m_2015
# x21: Distance to herbaceous areas
# x32: Distance to local roads
# x34: Distance to market places
# x37: osm_dst_railways
# x40: Slope
# x42: Night-time lights intensity
# --------------------------------------------------------------


# Visualise the matrix of correlation coefficients of the 
# selected covs
require(corrplot) # for correlation plot

#png(paste0(results_path, "/plots/cor_plots.png"))
# for population density
corrplot(
  cor(mcovs[,covs1]),
  method = 'circle',
  type = 'lower',
  tl.col = 'black',
  tl.cex = 1,
  col = colorRampPalette(c('orange', 'brown'))(200)
)
#dev.off()


# for number of households density
corrplot(
  cor(mcovs[,covs2]),
  method = 'circle',
  type = 'lower',
  tl.col = 'black',
  tl.cex = 1,
  col = colorRampPalette(c('orange', 'brown'))(200)
)
#dev.off()
#---------------------------------------------------------
### Test for Spatial Autoccorrelation using Moran's I stat
# Define the neighbors (use queen case)
names(shp)
#library(sf)
#library(spdep)
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
covars_est <- dat[,c("x2", "x3", "x17", "x19","x20", "x21", "x32", "x34","x37", "x40",
                     "x42", "set_reg", "set_typ", "region", "IDsr")]; dim(covars_est)

names(dat)


#---Build data stacks

# population density stack
stk_est1 <- inla.stack(data=list(y=dat$dens), # the response variable 
                      
                      A=list(A,1),  # the A matrix; the 1 is included to make the list(covariates)
                      
                      effects=list(c(list(Intercept=1), # the Intercept
                                     iset),  # the spatial index
                                   # the covariates
                                   list(covars_est)
                      ), 
                      #this is a quick name so you can call upon easily
                      tag='est')



# households per building stack
stk_est2 <- inla.stack(data=list(y=dat$dens_hh), # the response variable 
                      
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
            data=inla.stack.data(stk_est1,spde=spde),  #the data stack
            family= 'gamma',   #which family the data comes from
            control.predictor=list(A=inla.stack.A(stk_est1),compute=TRUE),  #compute gives you the marginals of the linear predictor
            control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
            verbose = FALSE) #can include verbose=TRUE to see the log of the model runs


## households per bldg
form1b <- y ~ -1 + Intercept + x2 + x17 + x19 + x20 + x21  + x32 + x37 + x42 + 
  f(IDsr, model='iid') + f(set_typ, model='iid')

mod1b <-inla(form1b, #the formula
            data=inla.stack.data(stk_est2,spde=spde),  #the data stack
            family= 'gamma',   #which family the data comes from
            control.predictor=list(A=inla.stack.A(stk_est2),compute=TRUE),  #compute gives you the marginals of the linear predictor
            control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
            verbose = FALSE) #can include verbose=TRUE to see the log of the model runs


# Model 2 : covs + Spatial autocor + settlement type - region interactions effects
form2 <- y ~  -1 + Intercept + x2 + x3 + x17 + x21 + x32 + x34 + x40 + x42 + 
  f(IDsr, model='iid') + f(set_reg, model='iid') + f(spatial.field, model=spde)

mod2 <-inla(form2, #the formula
            data=inla.stack.data(stk_est1,spde=spde),  #the data stack
            family= 'gamma',   #which family the data comes from
            control.predictor=list(A=inla.stack.A(stk_est1),compute=TRUE),  #compute gives you the marginals of the linear predictor
            control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
            verbose = FALSE) #can include verbose=TRUE to see the log of the model runs


# households per bldg
form2b <- y ~  -1 + Intercept + x2 + x17 + x19 + x20 + x21  + x32 + x37 + x42 + 
  f(IDsr, model='iid') + f(set_reg, model='iid') + f(spatial.field, model=spde)

mod2b <-inla(form2b, #the formula
            data=inla.stack.data(stk_est2,spde=spde),  #the data stack
            family= 'gamma',   #which family the data comes from
            control.predictor=list(A=inla.stack.A(stk_est2),compute=TRUE),  #compute gives you the marginals of the linear predictor
            control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
            verbose = FALSE) #can include verbose=TRUE to see the log of the model runs



# Model 3 : covs + Spatial autocor + settlement type - region interactions effects + region
form3 <- y ~  -1 + Intercept  + x2 + x3 + x17 + x21 + x32 + x34 + x40 + x42 + 
  f(set_reg, model='iid') + f(spatial.field, model=spde)+ f(set_typ, model='iid')+ f(IDsr, model='iid')

mod3 <-inla(form3, #the formula
            data=inla.stack.data(stk_est1,spde=spde),  #the data stack
            family= 'gamma',   #which family the data comes from
            control.predictor=list(A=inla.stack.A(stk_est1),compute=TRUE),  #compute gives you the marginals of the linear predictor
            control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
            verbose = FALSE) #can include verbose=TRUE to see the log of the model runs

# households per bldg
form3b <- y ~  -1 + Intercept  + x2 + x17 + x19 + x20 + x21  + x32 + x37 + x42 + 
  f(set_reg, model='iid') + f(spatial.field, model=spde)+ f(set_typ, model='iid')+ f(IDsr, model='iid')

mod3b <-inla(form3b, #the formula
            data=inla.stack.data(stk_est2,spde=spde),  #the data stack
            family= 'gamma',   #which family the data comes from
            control.predictor=list(A=inla.stack.A(stk_est2),compute=TRUE),  #compute gives you the marginals of the linear predictor
            control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
            verbose = FALSE) #can include verbose=TRUE to see the log of the model runs

## Run Model fit checks Based on DIC
DIC <- t(c(mod1 = mod1$dic$dic, mod2 = mod2$dic$dic,# pop density
           mod3 = mod3$dic$dic)) 
DIC2 <- t(c(mod1 = mod1b$dic$dic, mod2 = mod2b$dic$dic,# households per building
           mod3 = mod3b$dic$dic))


WAIC <- t(c(mod1 = mod1$waic$waic, mod2 = mod2$waic$waic,# pop density
            mod3 = mod3$waic$waic))
WAIC2 <- t(c(mod1 = mod1b$waic$waic, mod2 = mod2b$waic$waic,# households per building
            mod3 = mod3b$waic$waic))


CPO <- t(c(mod1 = -sum(log(mod1$cpo$cpo), na.rm=T), mod2 = -sum(log(mod2$cpo$cpo), na.rm=T),
           mod3 = -sum(log(mod3$cpo$cpo), na.rm=T)))
CPO2 <- t(c(mod1 = -sum(log(mod1b$cpo$cpo), na.rm=T), mod2 = -sum(log(mod2b$cpo$cpo), na.rm=T),
           mod3 = -sum(log(mod3b$cpo$cpo), na.rm=T)))


out_path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3/Working/CMR/Chris_N/paper2/output"
# population density model fit checks
(mod.select <- data.frame(DIC = t(DIC),
                          WAIC = t(WAIC),
                          CPO = t(CPO)))# Model 4 provided the best fit, thus, selected
#        DIC     WAIC      CPO
# mod1 1953.595 1453.742 6143.408
# mod2 1344.970 1825.464 6779.899
# mod3 1946.788 1324.180 5728.780 *** best model

# number of households model fit checks
(mod.select2 <- data.frame(DIC = t(DIC2),
                          WAIC = t(WAIC2),
                          CPO = t(CPO2)))# Model 4 provided the best fit, thus, selected
#         DIC      WAIC       CPO
# mod1 -5243.679 -5729.389 2457.2033
# mod2 -5266.214 -5585.049 2567.2805
# mod3 -5184.568 -6352.544  968.1896 *** best model

#write.csv(mod.select, paste0(out_path, "/pop_density_fit.csv"))
#write.csv(mod.select2, paste0(out_path, "/number_of_hholds_fit.csv"))
# Model 3 provided the best fit for both population density and number of hholds models
# extract the fixed effects 
fixed.effects <- round(mod3$summary.fixed,4) # to 4 dps
fixed.effects


# households per bldg
fixed.effects2 <- round(mod3b$summary.fixed,4) # to 4 dps
fixed.effects2

write.csv(fixed.effects, paste0(out_path, "/pop_density_betas.csv"))
write.csv(fixed.effects2, paste0(out_path, "/number_of_hholds_betas.csv"))
# Extract the back transformed linear predictor based on the best fit model (model 4)
ind3 <-inla.stack.index(stk_est1, "est")$data # indices of the estimation mesh nodes 
fit3 <- exp(mod3$summary.linear.predictor[ind3,"mean"]) # back transformed linear predictor
# Obtain posterior estimates at EA level 
sum(pred3 <- round(fit3*dat$bldg))


# households per bldg
ind3b <-inla.stack.index(stk_est2, "est")$data # indices of the estimation mesh nodes 
fit3b <- exp(mod3b$summary.linear.predictor[ind3b,"mean"]) # back transformed linear predictor
# Obtain posterior estimates at EA level 
sum(pred3b <- round(fit3b*dat$bldg))

# visualise the spatial fields of the best fit model
gproj <- inla.mesh.projector(mesh,  dims = c(300, 300))
bb <- non_convex_bdry$loc # the coordiates 

library(splancs) # for inout function for mapping
table(xy.in <- inout(gproj$lattice$loc,bb))
# FALSE  TRUE 
# 59771 30229


##---Remove points not in the study domain
g.mean <- inla.mesh.project(gproj, mod3$summary.random$spatial.field$mean)
g.sd <- inla.mesh.project(gproj, mod3$summary.random$spatial.field$sd)

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
#  2) Ideally, ensure there is a large enough (at least 30 iteration) sample size 
#  3) Predict population density parameters simulated at each iterations
#  4) Obtain the predicted population counts and save the output 
#  5) Save the full density and population counts matrices
#  6) Obtain the posterior statistics as well as uncertainty quantification 
#--------------------------------------------------------------------------
# Load the stack of prediction covariates
# download the prediction stack "prediction_data.RDS" from Google drive here: https://drive.google.com/file/d/1290hqUnBHhQS0I3iijddj34solTG-S_S/view?usp=sharing
# ....data_path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3/Working/CMR/Chris_N/paper1/submission"
pred_data <- readRDS(paste0(out_path, "/prediction_data.RDS"))#already standardized
data <- pred_data
head(data)
#    Prediction covs matrix
Apred <- inla.spde.make.A(mesh = mesh, loc = cbind(data$Lon, data$Lat))
dim(Apred)


####-------Run the posterior stimulation and uncertainty quantification 
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

##  number of households
simhhb <- function(model, dat, Aprediction, run)
{
  fixedeff <- dens2_hat <- hh_hat <- matrix(0, nrow=nrow(dat), ncol = run)
  #inla.seed = as.integer(runif(1)*.Machine$integer.max)
  inla.seed = 2111236018
  set.seed(inla.seed)
  print(inla.seed)
  m1.samp <- inla.posterior.sample(run, model, 
                                   seed = inla.seed, 
                                   selection=list(x2 = 1, 
                                                  x17 = 1, 
                                                  x19 = 1,
                                                  x20 = 1,
                                                  x21 = 1,
                                                  x32 = 1,
                                                  x37 = 1,
                                                  x42 = 1),
                                   num.threads="1:1")
  
  sfield_nodes_mean <- model$summary.random$spatial.field['mean']
  field_mean <- (Apred%*% as.data.frame(sfield_nodes_mean)[, 1])
  for(i in 1:run)
  {
    fixedeff[,i] <- 
      model$summary.fixed['Intercept', 'mean'] +
      m1.samp[[i]]$latent[1,] * dat[,'x2']  +
      m1.samp[[i]]$latent[2,] * dat[,'x17'] +
      m1.samp[[i]]$latent[3,] * dat[,'x19'] +
      m1.samp[[i]]$latent[4,] * dat[,'x20'] +
      m1.samp[[i]]$latent[5,] * dat[,'x21'] +
      m1.samp[[i]]$latent[6,] * dat[,'x32'] + 
      m1.samp[[i]]$latent[7,] * dat[,'x37'] +
      m1.samp[[i]]$latent[8,] * dat[,'x42'] +
      
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[2]) + #---settlement type and region nested effects
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[5]) + #---settlement type random effect
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[6]) + #---cluster level random effect
      field_mean[,1]
    #
    dens2_hat[,i] <- exp(fixedeff[,i])
    hh_hat[,i] <- dens2_hat[,i]*dat$bld
    
  }
  dat$mean_dens2_hat <- apply(dens2_hat, 1, mean, na.rm=T) 
  dat$mean_hh_hat  <- apply(hh_hat, 1, mean, na.rm=T) #
  dat$lower_hh_hat <- apply(hh_hat, 1, quantile, probs=c(0.025), na.rm=T) #
  dat$upper_hh_hat <- apply(hh_hat, 1, quantile, probs=c(0.975), na.rm=T) #
  dat$sd_hh_hat <- apply(hh_hat, 1, sd, na.rm=T) #
  dat$cv <- dat$sd_hh_hat/dat$mean_hh_hat#
  
  output <- list(hh_hat = hh_hat,
                 est_data = dat)
  
}
run=100

data$bld <- data$CMR_buildings_count
system.time(str(sim.dens1 <- simDens(mod3,data,Apred, run)))
#  Join the posterior sample to the prediction data
data.sim <- data.frame(cbind(data[,c("CMR_Regions", "CMR_Department","CMR_Settlement_Classification",
                                     "CMR_Arrondissement")], sim.dens1$pop_hat))


system.time(str(sim.hh <- simhhb(mod3b,data,Apred, run)))
data.sim2 <- data.frame(cbind(data[,c("CMR_Regions", "CMR_Department","CMR_Settlement_Classification",
                                     "CMR_Arrondissement")], sim.hh$hh_hat))

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

##----------
(national2 <- nat_total(data.sim2, run))
(national2 <- data.frame(total= national2[1,],
                        lower = national2[2,],
                        median=national2[3,],
                        upper=national2[4,]))

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

(regional.est2 <- regional_est(data.sim2, 100))
sum(regional.est2$total)
#write.csv(regional.est, file=paste0(results_path, "/estimates/regional_final.csv"))

##---Divisional estimates
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

##
(divisional.est2 <- divisional_est(data.sim2, 100))
sum(divisional.est2$total)
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
    met_in <- mod_metrics2(test$pop,  
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
    met_out <- mod_metrics2(test$obs,  
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
                             mod = mod3, 
                             form = form3,
                             A, 
                             seed = 13235))

cross_val$met_list_in_dat  # in-sample metrics per fold 
cross_val$met_list_out_dat  # out-of-sample metrics per fold
cross_val$cv_metrics    # combined averaged metrics
cross_val$pred_dat  # combined prediction data



#####------------------------------------------------------------
# cross - validation for household number model
#----------------------------------------------------------------
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
      mod$summary.fixed['x17', 'mean'] * test[,'x17'] +
      mod$summary.fixed['x19', 'mean'] * test[,'x19'] +
      mod$summary.fixed['x20', 'mean'] * test[,'x20'] +
      mod$summary.fixed['x21', 'mean'] * test[,'x21'] +
      mod$summary.fixed['x32', 'mean'] * test[,'x32'] +
      mod$summary.fixed['x37', 'mean'] * test[,'x37'] + 
      mod$summary.fixed['x42', 'mean'] * test[,'x42'] +
      
      
      rnorm(nrow(test), 0, 1/mod$summary.hyperpar$mean[2]) + #---settlement type and region nested effects
      rnorm(nrow(test), 0, 1/mod$summary.hyperpar$mean[5]) + #---settlement type random effect
      
      
      mod$summary.random$IDsr['mean'][test_ind,1] + #--uncorrelated spatial random effects
      
      field_mean[test_ind,1]
    
    dens_hh <- exp(fixed)
    sum(pop_hh <- dens_hh*test$bld)
    
    
    # visualise samples
    #par(mfrow =c(1,2))
    #plot(shp)
    #points(train_coords, col="blue", pch =15, cex=0.6)
    #points(test_coords, col="orange", pch=15, cex=0.6)
    
    
    par(mfrow =c(1,1))
    plot(test$obs2, pop_hh, xlab = "Observed", 
         ylab = "Predicted", col=c('blue','orange'),
         pch=c(16,16), cex.axis=1.5)
    abline(0,1)
    legend("topleft", c("Observed", "Predicted"), col=c("blue", "orange"), pch=c(16,16),
           bty="n", cex=1.5) 
    
    
    
    # calculate fit metrics
    met_in <- mod_metrics2(test$obs2,  
                           pop_hh)
    
    met_list_in[[i]]<- unlist(met_in)
    pred_list_in[[i]] <- data.frame(obs = test$obs2, pred = pop_hh,
                                    fold = rep(i, length(test$obs2)),
                                    data = rep("insample", length(test$obs2)))
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
    covars_train <- train[,c("x2", "x17", "x19", "x20", "x21", "x32", "x37",
                             "x42", "set_reg", "set_typ", "region", "IDsr")]; dim(covars_train)
    
    stk_train <- inla.stack(data=list(y=train$dens_hh), #the response
                            
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
      mod$summary.fixed['Intercept', 'mean'] +
      mod$summary.fixed['x2', 'mean'] * test[,'x2'] +
      mod$summary.fixed['x17', 'mean'] * test[,'x17'] +
      mod$summary.fixed['x19', 'mean'] * test[,'x19'] +
      mod$summary.fixed['x20', 'mean'] * test[,'x20'] +
      mod$summary.fixed['x21', 'mean'] * test[,'x21'] +
      mod$summary.fixed['x32', 'mean'] * test[,'x32'] +
      mod$summary.fixed['x37', 'mean'] * test[,'x37'] + 
      mod$summary.fixed['x42', 'mean'] * test[,'x42'] +
      
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$mean[2]) + #---settlement type and region nested effects
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$mean[5]) + #---settlement type random effect
      
      
      mod$summary.random$IDsr['mean'][test_ind,1] + #--uncorrelated spatial random effects
      
      field_mean[test_ind,1]
    
    dens_hh <- exp(fixed)
    sum(pop_hh <- dens_hh*test$bld)
    
    
    # visualise samples
    # par(mfrow =c(1,2))
    #plot(shp)
    #points(train_coords, col="blue", pch =15, cex=0.6)
    #points(test_coords, col="orange", pch=15, cex=0.6)
    
    par(mfrow =c(1,1))
    plot(test$obs2, pop_hh, xlab = "Observed", 
         ylab = "Predicted", col=c('blue','orange'),
         pch=c(16,16), cex.axis=1.5)
    abline(0,1)
    legend("topleft", c("Observed", "Predicted"), col=c("blue", "orange"), pch=c(16,16),
           bty="n", cex=1.5) 
    
    
    
    # calculate fit metrics
    met_out <- mod_metrics2(test$obs2,  
                            pop_hh)
    
    met_list_out[[i]]<- unlist(met_out)
    pred_list_out[[i]] <- data.frame(obs = test$obs2, pred = pop_hh,
                                     fold = rep(i, length(test$obs2)),
                                     data = rep("outsample", length(test$obs2)))
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
dat$obs2 <- dat$hh_nm_c # observed household size

(cross_val2 <- cross_validate(dat, n.folds = 5, 
                             mod = mod3b, 
                             form = form3b,
                             A, 
                             seed = 13235))

cross_val2$met_list_in_dat  # in-sample metrics per fold 
cross_val2$met_list_out_dat  # out-of-sample metrics per fold
cross_val2$cv_metrics    # combined averaged metrics
cross_val2$pred_dat  # combined prediction data

####
write.csv(cross_val$cv_metrics, paste0(out_path, "/pop_density_cv_metrics.csv"))
write.csv(cross_val2$cv_metrics, paste0(out_path, "/number_of_hholds_cv_metrics.csv"))
#-------------------------------------------------------------
# Scatter plots of model cross-validation results
#------------------------------------------------------------
pred.data <- cross_val$pred_dat
pred.data<- pred.data %>% dplyr::filter(pred<2000)

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
                       palette = "lancet"#,           # Color by groups "cyl"
                       #shape = "Fold"                             # Change point shape by groups "cyl"
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


# Individually
plot_cval2 <- ggscatter(pred.data, x = "obs", y = "pred",
                        add = "reg.line",                         # Add regression line
                        facet.by = "Fold",
                        #conf.int = TRUE,                          # Add confidence interval
                        color = "data", 
                        palette = "jco"#,           # Color by groups "cyl"
                        #shape = "Fold"                             # Change point shape by groups "cyl"
) 
rcval2 <-  ggpar(plot_cval2, xlab="Observed counts", ylab="Predicted Counts",
                 legend = "right", 
                 #legend.title = "Dataset",size=20,
                 legend.title = list(color = "Dataset", 
                                     linetype = "data", shape = "Fold"),
                 font.legend=c(12),
                 font.label = list(size = 12, face = "bold", color ="red"),
                 font.x = c(12),
                 font.y = c(12),
                 font.main=c(14),
                 font.xtickslab =c(12),
                 font.ytickslab =c(12),
                 # orientation = "reverse",
                 xtickslab.rt = 45, ytickslab.rt = 45)
rcval2


### Altogether
plot_cval3 <- ggscatter(pred.data, x = "obs", y = "pred",
                        add = "reg.line",                         # Add regression line
                        #facet.by = "Fold",
                        #conf.int = TRUE,                          # Add confidence interval
                        color = "data", 
                        palette = c("orange", "dark green")#,
                        #palette = "pnj",           # Color by groups "cyl"
                        #shape = "Fold"                             # Change point shape by groups "cyl"
) 

library(ggpubr)
rcval3 <-  ggpar(plot_cval3, xlab="Observed counts", ylab="Predicted Counts",
                 legend = "top", 
                 #legend.title = "Dataset",size=20,
                 legend.title = list(color = "Dataset", 
                                     linetype = "data", shape = "Fold"),
                 font.legend=c(12),
                 font.label = list(size = 12, face = "bold", color ="red"),
                 font.x = c(12),
                 font.y = c(12),
                 font.main=c(14),
                 font.xtickslab =c(12),
                 font.ytickslab =c(12),
                 # orientation = "reverse",
                 xtickslab.rt = 45, ytickslab.rt = 45)
rcval3


### - Density plots for cros validation for pop density
dens_plt <-ggplot(data=pred.data,aes(obs,pred)) + # swapped 
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black',
                 contour_var = "ndensity") + 
  scale_fill_continuous(type="viridis") +
  geom_smooth(method=lm,linetype=1,colour="black",se=F) + 
  scale_x_continuous(breaks=seq(0, 2500,500))+
  scale_y_continuous(breaks=seq(0, 2500,500))+
  guides(alpha="none") +
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))+
 # geom_point(col="red", size = 0.05) +
  facet_wrap(~data)

rdens_pop <-  ggpar(dens_plt, xlab="Observed population count", 
                    ylab="Predicted population count",
                   legend = "right", 
                   legend.title = "Density",
                   font.legend=c(15),
                   font.label = list(size = 15, face = "bold", color ="red"),
                   #font.x = c(12),
                  # palette = "lancet",
                   font.y = c(15),
                   font.x = c(15),
                   font.main=c(14),
                   font.xtickslab =c(14),
                   font.ytickslab =c(14),
                   xtickslab.rt = 45, ytickslab.rt = 45)
rdens_pop

#-------------------------------
#### fOr household number models
#------------------------------------------

dim(pred.datab <- cross_val2$pred_dat)
pred.datab<- pred.datab %>% dplyr::filter(pred<500)

pred.datab$Fold <- factor(pred.datab$fold)
pred.datab$data <- factor(pred.datab$data,
                         levels = c("insample", "outsample"),
                         labels = c("In-Sample", "Out-of-Sample"))
dens_hh <-ggplot(data=pred.datab,aes(obs,pred)) + # swapped 
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black',
                 contour_var = "ndensity") + 
  scale_fill_continuous(type="viridis") +
  geom_smooth(method=lm,linetype=1,colour="black",se=F) + 
  scale_x_continuous(breaks=seq(0, 600,150))+
  scale_y_continuous(breaks=seq(0, 600,150))+
  guides(alpha="none") +
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))+
 # geom_point(col="red", size = 0.05) +
  facet_wrap(~data)

rdens_hh<-  ggpar(dens_hh, xlab="Observed household count", 
                  ylab="Predicted household count",
                    legend = "right", 
                    legend.title = "Density",
                    font.legend=c(15),
                  font.label = list(size = 15, face = "bold", color ="red"),
                  #font.x = c(12),
                 # palette = "lancet",
                  font.y = c(15),
                  font.x = c(15),
                  font.main=c(14),
                  font.xtickslab =c(14),
                  font.ytickslab =c(14),
                    xtickslab.rt = 45, ytickslab.rt = 45)
rdens_hh
ggarrange(rdens_pop, rdens_hh+ rremove("x.text"), 
          labels = c("(A)", "(B)"),
          nrow=2)


# --------------------------------------------------------
# Notched box plots
#-----------------------------------------------------

# for population count
plot_cval4 <- pred.data %>% dplyr::filter(pred<2000) %>%
    ggplot( aes(x=Fold, y=pred, fill=Fold)) +
    #geom_violin(width=1.4) +
    geom_boxplot(width=0.8, #color="grey", 
                 # custom boxes
                # color="blue",
                # fill="magenta",
                 alpha=0.8,
                 
                 # Notch?
                 notch=TRUE,
                 notchwidth = 0.8) +
    #scale_fill_viridis(discrete = TRUE) +
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))+
  facet_wrap(~data)
  

rcval4 <-  ggpar(plot_cval4, xlab="Data Fold", ylab="Predicted population count",
                 legend = "none", 
                 #legend.title = "Dataset",size=20,
                 legend.title = list(color = "Dataset", 
                                     linetype = "data", shape = "Fold"),
                 font.legend=c(15),
                 font.label = list(size = 15, face = "bold", color ="red"),
                 #font.x = c(12),
                 palette = "lancet",
                 font.y = c(15),
                 font.x = c(15),
                 font.main=c(14),
                 font.xtickslab =c(14),
                 font.ytickslab =c(14),
                 # orientation = "reverse",
                 xtickslab.rt = 45, ytickslab.rt = 45)
rcval4

#---------------------------------------------
# for number of households

plot_cval4b <- pred.datab %>% dplyr::filter(pred<2000) %>%
  ggplot( aes(x=Fold, y=pred, fill=Fold)) +
  #geom_violin(width=1.4) +
  geom_boxplot(width=0.8, #color="grey", 
               # custom boxes
               # color="blue",
               # fill="magenta",
               alpha=0.8,
               
               # Notch?
               notch=TRUE,
               notchwidth = 0.8) +
  #scale_fill_viridis(discrete = TRUE) +
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))+
  facet_wrap(~data)


rcval4b <-  ggpar(plot_cval4b, xlab="Data Fold", ylab="Predicted household count",
                 legend = "none", 
                 #legend.title = "Dataset",size=20,
                 legend.title = list(color = "Dataset", 
                                     linetype = "data", shape = "Fold"),
                 font.legend=c(15),
                 font.label = list(size = 15, face = "bold", color ="red"),
                 #font.x = c(12),
                 palette = "lancet",
                 font.y = c(15),
                 font.x = c(15),
                 font.main=c(14),
                 font.xtickslab =c(14),
                 font.ytickslab =c(14),
                 # orientation = "reverse",
                 xtickslab.rt = 45, ytickslab.rt = 45)
rcval4b


ggarrange(rcval4, rcval4b+ rremove("x.text"), 
          labels = c("(A)", "(B)"),
          nrow=2)


# p
#
ggarrange(rdens_hh, rcval4b+ rremove("x.text"), 
          labels = c("(A)", "(B)"),
          nrow=2)


# -----------------------------------------------------------
# Prepare and save the grid vector predictions as raster files
# -----------------------------------------------------------

# specify the reference coordinates for the predictions
ref_coords <- cbind(pred_data$Lon, pred_data$Lat)
x <- as.matrix(ref_coords)

# Extract the grid cell population predictions
sum(result <-sim.dens1$est_data$mean_pop_hat, na.rm=T) # mean
resultL <- sim.dens1$est_data$lower_pop_hat # lower bound
resultU <- sim.dens1$est_data$upper_pop_hat # upper bound
resultCV <- sim.dens1$est_data$cv # coefficient of variation


# Extract the grid cell household count predictions
sum(resultb <-sim.hh$est_data$mean_hh_hat, na.rm=T) # mean
resultLb <- sim.hh$est_data$lower_hh_hat # lower bound
resultUb <- sim.hh$est_data$upper_hh_hat # upper bound
resultCVb <- sim.hh$est_data$cv # coefficient of variation


# add grid cell values and write the raster files

# specify your output_path for the raster file 
#  Mean

z <- as.matrix(result)
cmr_mean = rasterFromXYZ(cbind(x, z))

zb <- as.matrix(resultb)
cmr_meanb = rasterFromXYZ(cbind(x, zb))

# predicted population count
writeRaster(cmr_mean, filename=paste0(out_path, "/gridded_mean_pop_total.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))

# predicted household count
writeRaster(cmr_meanb, filename=paste0(out_path, "/gridded_mean_hh_total.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))


# Lower 
zL <- as.matrix(resultL)
cmr_lower = rasterFromXYZ(cbind(x, zL))

zLb <- as.matrix(resultLb)
cmr_lowerb = rasterFromXYZ(cbind(x, zLb))


writeRaster(cmr_lower, filename=paste0(out_path, "/gridded_lower_pop_total.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))
writeRaster(cmr_lowerb, filename=paste0(out_path, "/gridded_lower_hh_total.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))


# Upper
zU <- as.matrix(resultU)
cmr_upper = rasterFromXYZ(cbind(x, zU))

zUb <- as.matrix(resultUb)
cmr_upperb = rasterFromXYZ(cbind(x, zUb))

writeRaster(cmr_upper, filename=paste0(out_path, "/gridded_upper_pop_total.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))
writeRaster(cmr_upperb, filename=paste0(out_path, "/gridded_upper_hh_total.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))

# Coefficient of Variation 
zcv <- as.matrix(resultCV)
cmr_cv = rasterFromXYZ(cbind(x, zcv))

##
zcvb <- as.matrix(resultCVb)
cmr_cvb = rasterFromXYZ(cbind(x, zcvb))

writeRaster(cmr_cv, filename=paste0(out_path, "/gridded_cv_pop_total.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))
writeRaster(cmr_cvb, filename=paste0(out_path, "/gridded_cv_hh_total.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))


# Compare the projected regional estimates with the modelled
library(ggplot2)
library(scales)
theme_set(theme_classic())

# prep data
names(df)
df 
df$NIS.Projected<- df$NIS.Projected/1000
df$total <- df$total/1000
plot(df$total, df$NIS.Projected)
abline(lm(df$NIS.Projected ~ df$total))

# add the percentage difference between modelled and projected

df$diff <- round(((df$NIS.Projected - df$total)/df$NIS.Projected)*100)
# Create the bar plot for the percentage difference
df$Change <- as.factor(ifelse(df$diff > 0, "Higher", "Lower"))
ggplot(df, aes(x = names, y = diff, fill = Change)) +
  geom_bar(stat = "identity") +
  #scale_fill_manual(values = c("Higher" = "blue", "Lower" = "red")) +
  labs(title = "Bar Plot of the percentage changes \n between the modelled and projected population estimates",
       x = "Region",
       y = "Percentage difference") +
  #scale_x_continuous(breaks=seq(0, 100,10))+
  scale_y_continuous(breaks=seq(-100, 100,10))+
  ylim(-100, 100)+
  theme_minimal()



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
               font.legend=c(15),
               palette = "jco",
               font.label = list(size = 15, face = "bold", color ="red"),
               font.x = c(15),
               font.y = c(15),
               font.main=c(14),
               font.xtickslab =c(14),
               font.ytickslab =c(14),
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


#####   Bar plots with Error bar
main_drive <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3/Working"
data_drive <- paste0(main_drive,"/CMR/Chris_N/paper2")
adm_csv <- paste0(data_drive, "/CMR_Population_and_household_Estimates/Admin_estimates_final")

#
# Load in regional and departmental shapefiles
adm1 <- st_read(paste0(data_drive,"/Regional/Region_SHP.shp"))#regional
plot(adm1["libelle"])

adm2 <- st_read(paste0(data_drive,"/Departmental/Departement_SHP.shp"))# departmental
plot(adm2["libelle"])
###-----------------
# regional data with projections
 dff <- read.csv(paste0(adm_csv, "/regional_final_with_projections_07_06_23.csv"))

dff$names1 <- as.factor(dff$names)
dff$total <- dff$total/1000
dff$lower <- dff$lower/1000
dff$upper <- dff$upper/1000
dff$NIS.Projected <- dff$NIS.Projected/1000

levels(dff$names1)[1] <- "Extrême Nord"
#myPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")))

dff1 <- diff
pcp <- ggplot(dff, aes(x=reorder(names1, -total), 
                       y=total, fill = total)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 1, 
           #fill= "grey",
           colour = "blue") +
  scale_fill_viridis_c() +  # Use a beautiful color palette
  geom_point(size = 3, color = "blue", alpha=0.8) +  
  geom_point(aes(x = names1, y = NIS.Projected),
             colour = "red",size = 3) + # Add points
   geom_errorbar(aes(ymin = lower, ymax = upper, colour="95% Uncertainty"),
                position = position_dodge(width = 0.4), width = 0.35,
                size=0.7, colour = "magenta")+
  #scale_fill_viridis_c(option = "viridis")+
  #coord_flip() + 
  theme_minimal()
#scale_fill_manual(values=c("dark blue", "light blue"))

pdcp <- ggpar(pcp, ylab="Estimated population ('000)", xlab="Region",
              legend = "none", legend.title = "Population \n Estimate ('000)",
              font.legend=c(12),
              font.label = list(size = 12, face = "bold", color ="red"),
              #palette = "pnj",
              font.x = c(14),
              font.y = c(14),
              font.main=c(14),
              font.xtickslab =c(14),
              font.ytickslab =c(14),
              xtickslab.rt = 45, ytickslab.rt = 45)

pdcp

#----------------------------------------------------------------
# Visualise the gridded data
library(terra)
#r <- raster(paste0(output_path, "gridded_mean_total.tif"))
grid_path <- paste0(data_drive, "/CMR_Population_and_household_Estimates/Gridded_final")

# population count grid data
r <- rast(paste0(grid_path, "/CMR_population_v1_0_gridded.tif"))
r.un <- rast(paste0(grid_path, "/CMR_population_v1_0_uncertainty.tif"))


# number of households grid data
rh <- rast(paste0(grid_path, "/CMR_population_v1_0_hh_num.tif"))
r.unh <- rast(paste0(grid_path, "/CMR_population_v1_0_hh_num_uncertainty.tif"))
####
cmr_grd.adj <-rh

# install.packages("exactextractr", type="binary")
#library(exactextractr)
r_df <- data.frame(rh)
head(r_df)
summary(r_df)
#hist(log(r_df$CMR_population_v1_0_hh_num))
#> summary(r_df)
#CMR_population_v1_0_hh_num
#Min.   :  0.1678          
#1st Qu.:  0.7390          
#Median :  1.6393          
#Mean   :  3.9273          
#3rd Qu.:  3.8753          
#Max.   :168.8653

rdf <- data.frame(r)
head(rdf)
summary(rdf)
#gridded_mean_total
#Min.   :  0.7773  
#1st Qu.:  4.0386  
#Median :  9.0014  
#Mean   : 21.6166  
#3rd Qu.: 21.4526  
#Max.   :798.6568 


cmr_shp <- st_read(paste0(data_path, "/CMR_Boundary.shp"))
# install.packages("basemaps")
library(basemaps)
library(ggplot2)
library(tidyterra)
library(raster)
#install.packages("ggmap")
library(ggmap)
library(mapview)
set_defaults(rh, map_service = "esri", map_type = "world_street_map")
# mapping the predicted means of number of households
#crs(r) = "EPSG:3857"
#plot(cmr_grd.adj)
# get_maptypes()
#ext(rh) <- c(6.5, 17, 1, 14)# EXTEND THE EXTENT
#crs(rh) = "EPSG:4326"
x <- basemap_raster(rh, map_service = "esri", map_type = "natgeo_world_map")
ext(x)
x_terr <- rast(x) 

plot(x_terr)
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
brks <- seq(1,180.5, by=21)
arg <- list(brks, labels=brks)
zlim=c(1 ,169)
plot(rh,
     legend=F, col = colorRampPalette(c("dark blue", "red", "yellow"))(255),
     asp = NA, cex.axis=1.4, axes=T,#zlim=zlim,
     breaks=brks)
plot(x_lonlat, add=T, legend=F)
plot(shpp,add=T, col="black", lwd=2)#
#plot(shpp,add=T, lwd=2,bg = "transparent")
plot(rh,
     legend=F, col = colorRampPalette(c("dark blue", "red", "yellow"))(255),
     asp = NA, 
     cex.axis=1.4, add=T, breaks = brks, zlim=zlim)

# add Douala
ext_dou <- extent(9.6, 9.9, 3.9, 4.2)
plot(ext_dou , add=TRUE, col='red', lwd=2)

# add Yaounde
ext_ynde <- extent(11.4, 11.7, 3.7, 4)
plot(ext_ynde , add=TRUE, col='red', lwd=2)


# plot Douala and Yaounde counts separately
plot(rh, ext = extent(ext_dou),
     legend=F, col = colorRampPalette(c("dark blue", "red", "yellow"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)

#
plot(rh, ext = extent(ext_ynde),
     legend=F, col = colorRampPalette(c("dark blue", "red", "yellow"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)


#### Zoom in a bit further
# add Douala
plot(rh, ext = extent(ext_dou),
     legend=F, col = colorRampPalette(c("dark blue", "red", "yellow"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)

#ext_dou <- extent(9.6, 9.9, 3.9, 4.2)/12
ext_dou2 <- extent(9.69, 9.72, 4.02, 4.05)
plot(ext_dou2 , add=TRUE, col='green', lwd=3)

# add Yaounde
plot(rh, ext = extent(ext_ynde),
     legend=F, col = colorRampPalette(c("dark blue", "red", "yellow"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)
#ext_ynde <- extent(11.4, 11.7, 3.7, 4)/6
ext_ynde2 <- extent(11.52, 11.55, 3.85, 3.88)
plot(ext_ynde2 , add=TRUE, col='green', lwd=3)

# plot Douala and Yaounde ZOOMED IN counts separately

plot(rh, ext = extent(ext_dou2),
     legend=F, col = colorRampPalette(c("dark blue", "red", "yellow"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)

#
plot(rh, ext = extent(ext_ynde2),
     legend=F, col = colorRampPalette(c("dark blue", "red", "yellow"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)


#### Mapping the uncertainty for number of households
# mapping the predicted means
#utm_crs <- "EPSG:3857"
#cmr_utm <- project(rcv, utm_crs)
#res(cmr_utm)
xbh <- basemap_raster(r.unh, map_service = "esri", map_type = "natgeo_world_map")
ext(xbh)
x_terrbh <- rast(xbh) 

plot(x_terrbh)
crs(x_terrbh)
x_lonlat <- project(x_terrbh, "EPSG:4326")

crs(x_lonlat)# check
plot(x_lonlat)

shpp <- st_transform(cmr_shp, crs=4326) # Transform
shpp <- as(st_geometry(shpp), "Spatial") # Change to spatial object
plot(shpp) # checks
crs(shpp)#
#"dark green", "orange", "white"
library(viridis)
#cmr_grd.adj <- round(r, 2)
brks <- seq(0.5,3.7, by=0.4)
arg <- list(brks, labels=brks)
zlim=c(0.5,3.7)
plot(r.unh,
     legend=F, col = colorRampPalette(c("green", "dark blue", "red"))(255),
     asp = NA, cex.axis=1.4, axes=T,#zlim=zlim,
     breaks=brks)
plot(x_lonlat, add=T, legend=F)
plot(shpp,add=T, col="black", lwd=2)#
#plot(shpp,add=T, lwd=2,bg = "transparent")
plot(r.unh,
     legend=F, col = colorRampPalette(c("green", "dark blue", "red"))(255),
     asp = NA, 
     cex.axis=1.4, add=T, breaks = brks, zlim=zlim)

# add Douala
ext_dou <- extent(9.6, 9.9, 3.9, 4.2)
plot(ext_dou , add=TRUE, col='red', lwd=2)

# add Yaounde
ext_ynde <- extent(11.4, 11.7, 3.7, 4)
plot(ext_ynde , add=TRUE, col='red', lwd=2)


# plot Douala and Yaounde counts separately
plot(r.unh, ext = extent(ext_dou),
     legend=F, col = colorRampPalette(c("green", "dark blue", "red"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)

#
plot(r.unh, ext = extent(ext_ynde),
     legend=F, col = colorRampPalette(c("green", "dark blue", "red"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)


#### Zoom in a bit further
# add Douala
plot(r.unh, ext = extent(ext_dou),
     legend=F, col = colorRampPalette(c("green", "dark blue", "red"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)

#ext_dou <- extent(9.6, 9.9, 3.9, 4.2)/12
ext_dou2 <- extent(9.69, 9.72, 4.02, 4.05)
plot(ext_dou2 , add=TRUE, col='magenta', lwd=3)

# add Yaounde
plot(r.unh, ext = extent(ext_ynde),
     legend=F, col = colorRampPalette(c("green", "dark blue", "red"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)
#ext_ynde <- extent(11.4, 11.7, 3.7, 4)/6
ext_ynde2 <- extent(11.52, 11.55, 3.85, 3.88)
plot(ext_ynde2 , add=TRUE, col='magenta', lwd=3)

# plot Douala and Yaounde ZOOMED IN counts separately

plot(r.unh, ext = extent(ext_dou2),
     legend=F, col = colorRampPalette(c("green", "dark blue", "red"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)

#
plot(r.unh, ext = extent(ext_ynde2),
     legend=F, col = colorRampPalette(c("green", "dark blue", "red"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)


### mapping the predicted means of population count
#crs(r) = "EPSG:3857"
#plot(cmr_grd.adj)
# get_maptypes()
#ext(r) <- c(6.5, 17, 1, 14)# EXTEND THE EXTENT
#crs(r) = "EPSG:4326"
x <- basemap_raster(r, map_service = "esri", map_type = "natgeo_world_map")
ext(x)
x_terr <- rast(x) 

plot(x_terr)
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
brks <- seq(1,801, by=80)
arg <- list(brks, labels=brks)
zlim=c(1 ,799)
plot(r,
     legend=F, col = colorRampPalette(c("dark blue", "orange", "yellow"))(255),
     asp = NA, cex.axis=1.4, axes=T,#zlim=zlim,
     breaks=brks)
plot(x_lonlat, add=T, legend=F)
plot(shpp,add=T, col="black", lwd=2)#

#plot(shpp,add=T, lwd=2,bg = "transparent")
plot(r,
     legend=F, col = colorRampPalette(c("dark blue", "orange", "yellow"))(255),
     asp = NA, 
     cex.axis=1.4, add=T, breaks = brks, zlim=zlim)

# add Douala
ext_dou <- extent(9.6, 9.9, 3.9, 4.2)
plot(ext_dou , add=TRUE, col='red', lwd=2)

# add Yaounde
ext_ynde <- extent(11.4, 11.7, 3.7, 4)
plot(ext_ynde , add=TRUE, col='red', lwd=2)


# plot Douala and Yaounde counts separately
plot(r, ext = extent(ext_dou),
     legend=F, col = colorRampPalette(c("dark blue", "orange", "yellow"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)

#
plot(r, ext = extent(ext_ynde),
     legend=F, col = colorRampPalette(c("dark blue", "orange", "yellow"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)


#### Zoom in a bit further
# add Douala
plot(r, ext = extent(ext_dou),
     legend=F, col = colorRampPalette(c("dark blue", "orange", "yellow"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)

#ext_dou <- extent(9.6, 9.9, 3.9, 4.2)/12
ext_dou2 <- extent(9.69, 9.72, 4.02, 4.05)
plot(ext_dou2 , add=TRUE, col='green', lwd=3)

# add Yaounde
plot(r, ext = extent(ext_ynde),
     legend=F, col = colorRampPalette(c("dark blue", "orange", "yellow"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)
#ext_ynde <- extent(11.4, 11.7, 3.7, 4)/6
ext_ynde2 <- extent(11.52, 11.55, 3.85, 3.88)
plot(ext_ynde2 , add=TRUE, col='green', lwd=3)

# plot Douala and Yaounde ZOOMED IN counts separately

plot(r, ext = extent(ext_dou2),
     legend=F, col = colorRampPalette(c("dark blue", "orange", "yellow"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)

#
plot(r, ext = extent(ext_ynde2),
     legend=F, col = colorRampPalette(c("dark blue", "orange", "yellow"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)


#### Mapping the uncertainty for predicted population
# mapping the predicted means
#utm_crs <- "EPSG:3857"
#cmr_utm <- project(rcv, utm_crs)
#res(cmr_utm)
xb <- basemap_raster(r.un, map_service = "esri", map_type = "natgeo_world_map")
ext(xb)
x_terrb <- rast(xb) 

plot(x_terrb)
crs(x_terrb)
x_lonlat <- project(x_terrb, "EPSG:4326")

crs(x_lonlat)# check
plot(x_lonlat)

shpp <- st_transform(cmr_shp, crs=4326) # Transform
shpp <- as(st_geometry(shpp), "Spatial") # Change to spatial object
plot(shpp) # checks
crs(shpp)#
#"dark green", "orange", "white"
library(viridis)
#cmr_grd.adj <- round(r, 2)
brks <- seq(0.8,4.8, by=0.4)
arg <- list(brks, labels=brks)
zlim=c(0.8,4.8)
plot(r.un,
     legend=F, col = colorRampPalette(c("green", "dark blue", "red"))(255),
     asp = NA, cex.axis=1.4, axes=T,#zlim=zlim,
     breaks=brks)
plot(x_lonlat, add=T, legend=F)
plot(shpp,add=T, col="black", lwd=2)#
#plot(shpp,add=T, lwd=2,bg = "transparent")
plot(r.un,
     legend=F, col = colorRampPalette(c("green", "dark blue", "red"))(255),
     asp = NA, 
     cex.axis=1.4, add=T, breaks = brks, zlim=zlim)

# add Douala
ext_dou <- extent(9.6, 9.9, 3.9, 4.2)
plot(ext_dou , add=TRUE, col='red', lwd=2)

# add Yaounde
ext_ynde <- extent(11.4, 11.7, 3.7, 4)
plot(ext_ynde , add=TRUE, col='red', lwd=2)


# plot Douala and Yaounde counts separately
plot(r.un, ext = extent(ext_dou),
     legend=F, col = colorRampPalette(c("green", "dark blue", "red"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)

#
plot(r.un, ext = extent(ext_ynde),
     legend=F, col = colorRampPalette(c("green", "dark blue", "red"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)


#### Zoom in a bit further
# add Douala
plot(r.un, ext = extent(ext_dou),
     legend=F, col = colorRampPalette(c("green", "dark blue", "red"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)

#ext_dou <- extent(9.6, 9.9, 3.9, 4.2)/12
ext_dou2 <- extent(9.69, 9.72, 4.02, 4.05)
plot(ext_dou2 , add=TRUE, col='magenta', lwd=3)

# add Yaounde
plot(r.un, ext = extent(ext_ynde),
     legend=F, col = colorRampPalette(c("green", "dark blue", "red"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)
#ext_ynde <- extent(11.4, 11.7, 3.7, 4)/6
ext_ynde2 <- extent(11.52, 11.55, 3.85, 3.88)
plot(ext_ynde2 , add=TRUE, col='magenta', lwd=3)

# plot Douala and Yaounde ZOOMED IN counts separately

plot(r.un, ext = extent(ext_dou2),
     legend=F, col = colorRampPalette(c("green", "dark blue", "red"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)

#
plot(r.un, ext = extent(ext_ynde2),
     legend=F, col = colorRampPalette(c("green", "dark blue", "red"))(255),
     asp = NA, cex.axis=1.4, 
     colNA="black",
     axes=F, breaks = brks, zlim=zlim)



##  aggregated map plots
#----------------------------------------------------------
### Map the regional household counts
#-----------------------------------------------------------------
# install.packages("prettymapr")
library(ggmap)
library(prettymapr)
library(raster)
library(ggplot2)
library(ggspatial)

# Households
# regions
dim(reg_hh <- read.csv(paste0(adm_csv, "/regional_hh_num_07_06_23.csv")))

reg_hh$names1 <- as.factor(reg_hh$names)
reg_hh$total <- reg_hh$total/1000
reg_hh$lower <- reg_hh$lower/1000
reg_hh$upper <- reg_hh$upper/1000

levels(reg_hh$names1)[1] <- "Extrême Nord"
#myPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")))


pch <- ggplot(reg_hh, aes(x=reorder(names1, -total), 
                       y=total, fill = total)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 1, 
           #fill= "grey",
           colour = "blue") +
  scale_fill_viridis_c(option="viridis") +  # Use a beautiful color palette
  geom_point(size = 3, color = "red", alpha=0.8) +  
  geom_errorbar(aes(ymin = lower, ymax = upper, colour="95% Uncertainty"),
                position = position_dodge(width = 0.4), width = 0.35,
                size=0.9, colour = "grey")+
  #scale_fill_viridis_c(option = "viridis")+
  #coord_flip() + 
  theme_minimal()
#scale_fill_manual(values=c("dark blue", "light blue"))

pdch <- ggpar(pch, ylab="Estimated number of households", xlab="Region",
              legend = "none", legend.title = "Household \n Count('000)",
              font.legend=c(12),
              font.label = list(size = 12, face = "bold", color ="red"),
              #palette = "pnj",
              font.x = c(14),
              font.y = c(14),
              font.main=c(14),
              font.xtickslab =c(14),
              font.ytickslab =c(14),
              xtickslab.rt = 45, ytickslab.rt = 45)

pdch

# departments
dim(dep_hh <- read.csv(paste0(adm_csv, "/department_final_hh_num_07_06_23.csv"))) 

dep_hh$names2 <- as.factor(dep_hh$names1)
dep_hh$total <- dep_hh$total/1000
dep_hh$lower <- dep_hh$lower/1000
dep_hh$upper <- dep_hh$upper/1000
pdh <- ggplot(dep_hh, aes(x=reorder(names2, -total), 
                          y=total, fill = total)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, colour="95% Uncertainty"),
                position = position_dodge(width = 0.4), width = 0.75,
                size=0.7, colour = "darkblue")+
  geom_point(size = 1, color = "red", alpha=0.8) +
  coord_flip() + 
  theme_minimal()
pdh2 <- ggpar(pdh, ylab="Estimated number of households", xlab="Department",
              legend = "none", legend.title = "Household \n Count('000)",
              font.legend=c(12),
             # font.label = list(size = 12, face = "bold", color ="red"),
              #palette = "pnj",
              font.x = c(12),
              font.y = c(12),
              font.main=c(12),
              font.xtickslab =c(10),
             # font.ytickslab =c(14),
              xtickslab.rt = 45)

pdh2



# popultaion counts 
# regions
#dim(reg <- read.csv(paste0(adm_csv, "/regional_final_with_projections_07_06_23.csv")))

dim(reg <- read.csv(paste0(adm_csv, "/regional_07_06_23.csv")))
reg$names1 <- as.factor(reg$names)
reg$total <- reg$total/1000
reg$lower <- reg$lower/1000
reg$upper <- reg$upper/1000

levels(reg$names1)[1] <- "Extrême Nord"

pc <- ggplot(reg, aes(x=reorder(names1, -total), 
                          y=total, fill = total)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 1, 
           #fill= "grey",
           colour = "blue") +
  scale_fill_viridis_c(option="inferno") +  # Use a beautiful color palette
  geom_point(size = 3, color = "red", alpha=0.8) +  
  geom_errorbar(aes(ymin = lower, ymax = upper, colour="95% Uncertainty"),
                position = position_dodge(width = 0.4), width = 0.35,
                size=0.9, colour = "grey")+
  #scale_fill_viridis_c(option = "viridis")+
  #coord_flip() + 
  theme_minimal()
#scale_fill_manual(values=c("dark blue", "light blue"))

pdc <- ggpar(pc, ylab="Estimated Population", xlab="Region",
              legend = "none", legend.title = "Population \n Count('000)",
              font.legend=c(12),
              font.label = list(size = 12, face = "bold", color ="red"),
              #palette = "pnj",
              font.x = c(14),
              font.y = c(14),
              font.main=c(14),
              font.xtickslab =c(14),
              font.ytickslab =c(14),
              xtickslab.rt = 45, ytickslab.rt = 45)

pdc



# departments
dim(dep <- read.csv(paste0(adm_csv, "/department_07_06_23.csv")))
dep$names2 <- as.factor(dep$names1)
dep$total <- dep$total/1000
dep$lower <- dep$lower/1000
dep$upper <- dep$upper/1000
pd <- ggplot(dep, aes(x=reorder(names2, -total), 
                          y=total, fill = total)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, colour="95% Uncertainty"),
                position = position_dodge(width = 0.4), width = 0.75,
                size=0.7, colour = "dark green")+
  geom_point(size = 1, color = "red", alpha=0.8) +
  coord_flip() + 
  theme_minimal()
pdh3 <- ggpar(pd, ylab="Estimated population count", xlab="Department",
              legend = "none", legend.title = "Population \n Count('000)",
              font.legend=c(12),
              # font.label = list(size = 12, face = "bold", color ="red"),
              #palette = "pnj",
              font.x = c(12),
              font.y = c(12),
              font.main=c(12),
              font.xtickslab =c(10),
              # font.ytickslab =c(14),
              xtickslab.rt = 45)

pdh3

# Join household number data with pop count data
# region
regg <- merge(reg, reg_hh, by = "names")
regg$libelle <- factor(regg$names)
levels(regg$libelle)[1] <- "Extrême Nord"

# depatment
depp <- merge(dep, dep_hh, by="names1")
plot(depp$total.x, depp$total.y)


# Join estimates to the shapefile

# region
regg_shp <- merge(adm1, regg, by = "libelle")
plot(regg_shp["total.y"])
plot(regg_shp["total.x"])

# department
depp$libelle <- depp$names1
depp_shp <- merge(adm2, depp, by = "libelle")
plot(depp_shp["total.y"])
plot(depp_shp["total.x"])


## Total number of households per region 
regg_shp$hhold_cnt <- round(regg_shp$total.y/1000)
min(regg_shp$hhold_cnt); max(regg_shp$hhold_cnt)
regg_shp$brks1 <- cut(regg_shp$hhold_cnt, 
                      breaks=c(200, 400, 600, 800,1000,1200), 
                      labels=c("200 - 400", "400- 600","600 - 800", 
                               "800 - 1000", "1000 - 1200"))
  hh_reg <- ggplot() + 
  annotation_map_tile(type = "osm") +  # OpenStreetMap basemap
  geom_sf(data = regg_shp, aes(fill = brks1)) +
  scale_fill_viridis_d(option="plasma","Household \n count('000)")+
  #ggspatial::annotation_scale(location="br")+
  #ggspatial::annotation_north_arrow(location="topleft")+
  theme_minimal() +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right")+
  theme(plot.title=element_text(face="bold",size=12))+
  theme(plot.subtitle=element_text(size=11))+ 
  theme(legend.text = element_text(size = 14))+  # Increase legend text size
  theme(legend.title = element_text(size = 16))  # Increase legend title size

hh_reg


## Uncertainty in the total number of households per region
regg_shp$hhold_un <- regg_shp$uncertainty.y
min(regg_shp$hhold_un); max(regg_shp$hhold_un)
regg_shp$brks2 <- cut(regg_shp$hhold_un, 
                      breaks=c(0.2, 0.35, 0.5, 0.65,0.8,0.95), 
                      labels=c("0.2 - 0.35", "0.35- 0.5","0.5 - 0.65", 
                               "0.65 - 0.8", "0.8 - 0.95"))
hh_reg_un <- ggplot() + 
  annotation_map_tile(type = "osm") +  # OpenStreetMap basemap
  geom_sf(data = regg_shp, aes(fill = brks2)) +
  scale_fill_brewer(palette = "RdYlGn","Uncertainty",
                    direction = -1) +
  #ggspatial::annotation_scale(location="br")+
  #ggspatial::annotation_north_arrow(location="topleft")+
  theme_minimal() +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right")+
  theme(plot.title=element_text(face="bold",size=12))+
  theme(plot.subtitle=element_text(size=11))+
  theme(legend.text = element_text(size = 14))+  # Increase legend text size
  theme(legend.title = element_text(size = 16))  # Increase legend title size
hh_reg_un


#RdYlGn
# Household counts: Departmental 
depp_shp$hhold_cnt <- round(depp_shp$total.y/1000)
min(depp_shp$hhold_cnt); max(depp_shp$hhold_cnt)
depp_shp$brks <- cut(depp_shp$hhold_cnt, 
                     breaks=c(10, 30, 50, 100,200,400,700), 
                     labels=c("10 - 30", "30 - 50", "50 - 100", "100 - 200", "200 - 400", "400 - 700"))
hh_dep <- ggplot() + 
  annotation_map_tile(type = "osm") +  # OpenStreetMap basemap
  geom_sf(data = depp_shp, aes(fill = brks)) +
  scale_fill_viridis_d(option="plasma","Household \n count ('000)")+
  theme_minimal() +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title=element_text(face="bold",size=12))+
  theme(legend.text = element_text(size = 14))+  # Increase legend text size
  theme(legend.title = element_text(size = 16))  # Increase legend title size
hh_dep



# Uncertainty in the number of households estimation
depp_shp$hhold_un <- depp_shp$uncertainty.y
min(depp_shp$hhold_un); max(depp_shp$hhold_un)
depp_shp$brks2 <- cut(depp_shp$hhold_un, 
                      breaks=c(0.2, 0.35, 0.5, 0.75, 0.9, 1.05, 1.3), 
                      labels=c("0.2 - 0.35", "0.35 - 0.5", "0.5 - 0.75", 
                               "0.75 - 0.9", "0.9 - 1.05", "1.05 - 1.3"))
  hh_dep_un <- ggplot() + 
  annotation_map_tile(type = "osm") +  # OpenStreetMap basemap
  geom_sf(data = depp_shp, aes(fill = brks2)) +
  scale_fill_brewer(palette = "RdYlGn","Uncertainty",
                    direction = -1) +
  theme_minimal() +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right")+
  theme(plot.title=element_text(face="bold",size=12))+
  theme(legend.text = element_text(size = 14))+  # Increase legend text size
  theme(legend.title = element_text(size = 16))  # Increase legend title size
hh_dep_un


## Population Count
## Total population per region 
regg_shp$pop_cnt <- round(regg_shp$total.x/1000)
min(regg_shp$pop_cnt); max(regg_shp$pop_cnt)
regg_shp$brks3 <- cut(regg_shp$pop_cnt, 
                      breaks=c(1000, 1500, 2000, 3000,4000,5500), 
                      labels=c("1000 - 1500", "1500 - 2000","2000 - 3000", 
                               "3000 - 4000", "4000 - 5500"))
pop_reg <- ggplot() + 
  annotation_map_tile(type = "osm") +  # OpenStreetMap basemap
  geom_sf(data = regg_shp, aes(fill = pop_cnt)) +
  scale_fill_viridis_c(option="viridis","Population \n count('000)")+
  theme_minimal() +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right")+
  theme(plot.title=element_text(face="bold",size=12))+
  theme(plot.subtitle=element_text(size=11))+ 
  theme(legend.text = element_text(size = 14))+  # Increase legend text size
  theme(legend.title = element_text(size = 16))  # Increase legend title size

pop_reg


## Uncertainty in the total population per region
regg_shp$pop_un <- regg_shp$uncertainty.x
min(regg_shp$pop_un); max(regg_shp$pop_un)

pop_reg_un <- ggplot() + 
  annotation_map_tile(type = "osm") +  # OpenStreetMap basemap
  geom_sf(data = regg_shp, aes(fill = pop_un)) +
   scale_fill_viridis_c(option="inferno","Uncertainty")+
  theme_minimal() +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right")+
  theme(plot.title=element_text(face="bold",size=12))+
  theme(legend.text = element_text(size = 14))+  # Increase legend text size
  theme(legend.title = element_text(size = 16))  # Increase legend title size
pop_reg_un


#RdYlGn
# Population count: Departmental 
depp_shp$pop_cnt <- round(depp_shp$total.x/1000)
min(depp_shp$pop_cnt); max(depp_shp$pop_cnt)

pop_dep <- ggplot() + 
  annotation_map_tile(type = "osm") +  # OpenStreetMap basemap
  geom_sf(data = depp_shp, aes(fill = pop_cnt)) +
  scale_fill_viridis_c(option="viridis","Population \n count('000)")+
 # scale_fill_viridis_d(option="viridis","Population \n count ('000)")+
  theme_minimal() +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title=element_text(face="bold",size=12))+
  theme(legend.text = element_text(size = 14))+  # Increase legend text size
  theme(legend.title = element_text(size = 16))  # Increase legend title size
pop_dep



# Uncertainty in population estimation 

depp_shp <- st_make_valid(depp_shp)
depp_shp$pop_un <- depp_shp$uncertainty.x
min(depp_shp$pop_un); max(depp_shp$pop_un)

pop_dep_un <- ggplot() + 
  annotation_map_tile(type = "osm") +  # OpenStreetMap basemap
  geom_sf(data = depp_shp, aes(fill = pop_un)) +
  #scale_fill_brewer(palette = "OrRd","Uncertainty") +
  scale_fill_viridis_c(option="inferno","Uncertainty")+
  theme_minimal() +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(plot.title=element_text(face="bold",size=12))+
  theme(legend.text = element_text(size = 14))+  # Increase legend text size
  theme(legend.title = element_text(size = 16)) # Increase legend title size
pop_dep_un

#
ggarrange(pop_reg, pop_reg_un+ rremove("x.text"), 
          labels = c("(A)", "(B)"),
          nrow=2)

depp_shp$uncertainty.x[is.na(depp_shp$uncertainty.x)]
##
# install.packages("patchwork")
library(patchwork)

# Combine the maps
combined_maps <- pop_reg + pop_dep + hh_reg + hh_dep

# Print the combined maps
print(combined_maps)

#save.image(paste0(output_path, "cmr_application.Rdata"))