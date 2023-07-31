####--TITLE: R SCRIPTS FOR MODELLED POPULATION ESTIMATES FOR CAMEROON BASED ON R-INLA----#####
####--METHODS: GEOSTATISTICAL BAYESIAN HIERARCHICAL REGRESSION MODEL---------------------#####
####--AUTHOR: DR CHIBUZOR CHRISTOPHER NNANATU-------------------------#####
####--INSTITUTION: WORLDPOP, UNIVERSITY OF SOUTHAMPTON----------------------------------##### 
####--DATE: March 2023---------------------------------------------------------------####
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###----MODEL COVARIATES PREPARATION AND SCALING
names(dat2)

cov_rows <- 20:62 #---model covariates are on columns 20 to 62 (43 covariates)
colnames(dat2)[cov_rows] <- paste0("x", 1:43, sep="") #--rename model covariates as x1 - x43

head(dat2[, cov_rows])
dat2[,cov_rows] <- apply(dat2[,cov_rows], 2, stdize)   #z-score
head(dat2[, cov_rows])
head(dat2)

names(dat2)
dat2$dens_bldg <- dat2$Imputed_LHHSIZE/dat2$Total_Building_Count #---density variable 
names(dat2)


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#---##############------VARIABLE SELECTION------###################----------#
library(car) ##--For calculating variance inflation factor (vif)
#######-------For Density=================================================================
#--%%%%%%%%%%%%%%%%%%%
covs_dens_bldg <- dat2[,c(65,cov_rows)] #--subset for variables selection
covs_dens_bldg$ldens_bldg <- log(covs_dens_bldg$dens_bldg) #--transform density on logscale
hist(covs_dens_bldg$ldens_bldg)
covs_dens_bldg1 <- covs_dens_bldg[,-1] %>% drop_na() #- model covariates only without NA's


#---@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##----Covariates selection using GLM STEPWISE REGRESSION
fdens_bldg <- glm(ldens_bldg ~., data=covs_dens_bldg1, family = gaussian)
summary(fdens_bldg)
stepden_bldg <- stepAIC(fdens_bldg, scale = 0,
                        direction = c("both"),
                        trace = 1, keep = NULL, steps = 1000, use.start = FALSE,
                        k = 2)
stepden_bldg

###-----------REFIT THE SELECTED MODEL AND CHOOSE only significant variables with VIF < 5
fdens_bldg <-  glm(formula = ldens_bldg ~ x1 + x2 + x3 + x4 + x7 + x9 + x10 + 
                     x12 + x13 + x14 + x15 + x16 + x17 + x20 + x21 + x22 + x32 + 
                     x33 + x34 + x39 + x40 + x42 + x43, family = gaussian, data = covs_dens_bldg1)


summary(fdens_bldg)

##---VIF
vif_dens_bldg = vif(fdens_bldg)
vif_dens_bldg[which(vif_dens_bldg < 5)]

##----Selected significant covariates with VIF < 5
cov_dens_bldg <- c("x2","x3","x17","x21","x32","x34","x40","x42")

cov_names <- names(dat)[cov_rows]
covariates_dens_bldg <- cov_names[c(2, 3, 17, 21, 32, 34, 40, 42)]
dens_bldg_covs <- data.frame(cov = cov_dens_bldg, name=covariates_dens_bldg )

##----Make correlation plot
require(corrplot)

png(paste0(results_path, "/plots/cor_plots.png"))
corrplot(
  cor(covs_dens_bldg[,cov_dens_bldg]),
  method = 'square',
  type = 'upper',
  tl.col = 'black',
  tl.cex = 1,
  col = colorRampPalette(c('purple', 'dark green'))(200)
)
dev.off()
