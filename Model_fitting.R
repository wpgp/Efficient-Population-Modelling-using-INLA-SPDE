####--TITLE: R SCRIPTS FOR MODELLED POPULATION ESTIMATES FOR CAMEROON BASED ON R-INLA----#####
####--METHODS: GEOSTATISTICAL BAYESIAN HIERARCHICAL REGRESSION MODEL---------------------#####
####--AUTHOR: DR CHIBUZOR CHRISTOPHER NNANATU -------------------------#####
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%visualize observed household totals
#names(shp2)
png(paste0(results_path, "/plots/observed_totals_by_EA.png"))
spplot(shp, "O_LHHSI", 
       col.regions = rev(magma(16))) #gray.colors(16, 0.9, 0.4))
dev.off()


library(tmap)
library(tmaptools)
##----Multiple facets maps by datasource 
tmap_options(check.and.fix = TRUE)
tm_shape(shp) +
  tm_polygons("O_LHHSI", palette = "RdYlBu") +
  #tm_polygons() +
  tm_borders() +
  tm_fill("O_LHHSI",
          palette = get_brewer_pal("YlGnBu"),
          legend.show = F,
          title = "Price in US$ per lb",
          style = "order")+
  tm_facets(by = "Dat_Typ")


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

#--*******************SPECIFY AND FIT VARIOUS MODELS-***************-------------------
#*###---Model 1-- gamma spatial (with settlement type random effect with covariates only) -----------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
form1 <- y ~ -1 + Intercept + x2 + x3 + x17 + x21 + x32 + x34 + x40 + x42 + 
  f(IDsr, model='iid') + f(set_typ, model='iid') + f(spatial.field, model=spde)

mod1<-inla(form1, #the formula
             data=inla.stack.data(stk_est,spde=spde),  #the data stack
             family= 'gamma',   #which family the data comes from
             control.predictor=list(A=inla.stack.A(stk_est),compute=TRUE),  #compute gives you the marginals of the linear predictor
             control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs




#--@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
###---Model 2-- gamma spatial (with settlement type and region random effects no nesting) -----------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
form72 <- y ~ -1 + Intercept + x2 + x3 + x17 + x21 + x32 + x34 + x40 + x42 + 
  f(IDsr, model='iid') + f(set_typ, model='iid') + f(region, model='iid') + f(spatial.field, model=spde)

mod2 <-inla(form2, #the formula
             data=inla.stack.data(stk_est,spde=spde),  #the data stack
             family= 'gamma',   #which family the data comes from
             control.predictor=list(A=inla.stack.A(stk_est),compute=TRUE),  #compute gives you the marginals of the linear predictor
             control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs




#-----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###---Model 3-- gamma spatial (with settlement type and region nesting) -----------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
form3 <- y ~  -1 + Intercept + x2 + x3 + x17 + x21 + x32 + x34 + x40 + x42 + 
  f(IDsr, model='iid') + f(set_reg, model='iid') + f(spatial.field, model=spde)

mod3 <-inla(form3, #the formula
             data=inla.stack.data(stk_est,spde=spde),  #the data stack
             family= 'gamma',   #which family the data comes from
             control.predictor=list(A=inla.stack.A(stk_est),compute=TRUE),  #compute gives you the marginals of the linear predictor
             control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs



#----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###---Model 4-- Lognormal spatial (with settlement type and region nesting and settlement type random effect) -----------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
form4 <- y ~  -1 + Intercept  + x2 + x3 + x17 + x21 + x32 + x34 + x40 + x42 + 
  f(set_reg, model='iid') + f(spatial.field, model=spde)+ f(set_typ, model='iid')+ f(IDsr, model='iid')

mod4 <-inla(form4, #the formula
             data=inla.stack.data(stk_est,spde=spde),  #the data stack
             family= 'gamma',   #which family the data comes from
             control.predictor=list(A=inla.stack.A(stk_est),compute=TRUE),  #compute gives you the marginals of the linear predictor
             control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
summary(mod4)

mod4$summary.fix
mod4$summary.random
ind4 <-inla.stack.index(stk_est, "est")$data

#-Predicted density
fit4 <- exp(mod4$summary.linear.predictor[ind4,"mean"])

#--predicted count 
sum(pred4 <- round(fit4*dat2$Total_Building_Count))

##---Model selection based on WAIC
##-----Model fit: WAIC
(waic=c(mod1$waic$waic, mod2$waic$waic,mod3$waic$waic, mod4$waic$waic))

#-- model 4 provided the best fit with the lowest WAIC. 
##---Scatter plot for model checks
plot(dat2$Imputed_LHHSIZE, pred4, xlab = "Observed", ylab = "Predicted", col='red')
abline(a=0, b=1)
cor(dat2$Imputed_LHHSIZE, pred4)  


####----View the spatial fields of the best fit model
gproj <- inla.mesh.projector(mesh,  dims = c(300, 300))
g.mean <- inla.mesh.project(gproj, mod9a$summary.random$spatial.field$mean)
g.sd <- inla.mesh.project(gproj, mod9a$summary.random$spatial.field$sd)

grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', main='mean',col.regions = heat.colors(16, rev=T)),
             levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='sd' ,col.regions = heat.colors(16, rev=T)), nrow=1)
