
###---BAYESIAN GEOSTATISTICAL POP MODEL PAPER SIMULATION STUDY CODES
###---NNANATU ET AL.(2024)
### - Simulation Study 

rm(list=ls())

library(INLA); library(raster); 
library(gtools); library(sp); library(spdep)
library(fields); library(mvtnorm); library(gtools)
library(geoR);library(actuar);library(viridisLite)
require(grid);require(gridExtra);require(lattice);require(tidyverse)
library(dplyr); library(sf); library(tmap); library(tmaptools)


out_path <- tempdir()
githublink <- "https://raw.github.com/wpgp/Efficient-Population-Modelling-using-INLA-SPDE/main/simulation_shapefile.RData"
load(url(githublink))


ls()      

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


dim(cmr_grid)# 44032 grid cells at 5km by 5km
dim(cmr_grid_cropped) # 19480 grid cells left after croping to CMR extents

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


####---Model fit metrics function
model_metrics <- function(obs, pred)
{
  residual = pred - obs
  MAE = mean(abs(residual), na.rm=T)#MAE
  MAPE = (1/length(obs))*sum(abs((obs-pred)/obs))*100#MAPE
  MSE = mean(residual^2, na.rm=T)
  RMSE = sqrt(MSE)
  BIAS = mean(residual, na.rm=T)
  corr = cor(obs[!is.na(obs)],pred[!is.na(obs)])
  
  output <- list(#MAPE = MAPE,
    MAE  = MAE ,
    BIAS = abs(BIAS),
    RMSE = RMSE,
    corr = corr)
  return(output)
}

#model_metrics(obs, pred, upper, lower)


#---------------------------------------------------------------------------------------
cover <- seq(0.2,1, by=0.2)
spat.autocor <- c(0.75, 1.25, 1.75, 2.75, 3.75) 
source_var <- list(same = rep(0.001, 5), 
                   diff = c(0.0015, 0.002, 0.0017, 0.001, 0.0032))#variance of data source random effects

#source_re[[2]]
#---Create the data frame 
dim(dat.all <- as.data.frame(coords))

n.dat_typ = 5 # number of different data sources 
prop.dat_typ <- c(0.19, 0.16, 0.16, 0.17, 0.32) # proportion of the different data sources
n.sample <- nrow(dat.all)
table(dat_typ1 <- sample(1:n.dat_typ,n.sample, prob=prop.dat_typ, rep=T))


###---Add to dataset
table(dat.all$dat_typ_ID <- as.factor(dat_typ1))


####------Build the mesh
# dim(coord)
coord <- cbind(coords$x, coords$y) # must be a matrix object
bnd <- inla.nonconvex.hull(coord, -0.035, -0.04, resolution = c(100, 100))
mesh <- inla.mesh.2d(boundary = bnd, max.edge=c(0.2,1), 
                     offset = c(0.2, 0.7),
                     cutoff = 0.2)
par(mfrow=c(1,1))
plot(mesh)
mesh$n # 1133 mesh nodes
#--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##---Specify SPDE parameters
rho <- 0.2 #--range
nu <- 1 #--smooth parameter

for(b in 1:length(spat.autocor))
{ 
  
  # b=1
  result_pathb <- paste0(out_path,"/marginal_var_", spat.autocor[b])
  if (file.exists(result_pathb)){
    #setwd(file.path(result_pathb))
  } else {
    dir.create(file.path(result_pathb))
    #setwd(file.path(result_pathb))
  }
  
  sigma0 <- spat.autocor[b]#--marginal variance  
  kappa0 <- sqrt(8*nu)/rho #--scale parameters
  (tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)) #--precision
  
  #---SPDE
  spde <- inla.spde2.matern(mesh, 
                            B.tau = matrix(c(log(tau0), -1, +1), nrow=1, ncol=3),
                            B.kappa = matrix(c(log(kappa0), 0, -1), nrow=1, ncol=3),
                            theta.prior.mean = c(0,0),
                            theta.prior.prec = c(0.1, 0.1))
  
  
  #--Precision Matrix
  Q <- inla.spde2.precision(spde=spde, theta=c(0,0))
  
  #---Simulate the GMRF
  sam <- as.vector(inla.qsample(
    n = 1, Q = Q, seed=100))
  #length(sam)
  
  ###---Build projector matrix A
  A <- inla.spde.make.A(mesh=mesh, loc=coord);dim(A)
  
  # Spatial Random effects
  sp.raneff <- as.vector(A %*% sam)
  length(sp.raneff)
  hist(sp.raneff)
  
  ##----specify the observation indices for estimation 
  iset <- inla.spde.make.index(name = "spatial.field", spde$n.spde)
  
  
  #-- generate the Covariates 
  nn <- nrow(dat.all)
  covs = cbind(1,runif(nn,1,10),rnorm(nn, 2, 2),rpois(nn,4),rnorm(nn, 4, 3), runif(nn))
  apply(covs, 2, mean) #---check averages
  dim(dat.all)
  
  # Standardization 
  stdize <- function(x)
  {
    #return(stdze <- (x-min(x, na.rm=T))/(max(x, na.rm =T)-min(x, na.rm =T)))
    return(stdze <- (x-mean(x, na.rm=T))/(sd(x, na.rm =T)))
  }
  
  covs_std <- cbind(covs[,1],apply(covs[,-1], 2, stdize)) 
  
  # Add covariates to data
  head(dat.all)
  rownames(dat.all) <- NULL
  dim(ddat <- cbind(dat.all, covs_std))
  ddat <- data.frame(ddat)
  colnames(ddat) <- c("lon", "lat", "source", "x1", "x2", "x3", "x4", "x5", "x6")
  names(ddat)
  head(ddat)
  ###############
  #Parameters 
  #
  #dim(zpred <- covs_std)
  dim(zpred <- covs)
  
  require(MASS)
  #--generate data type random effects 
  
  
  dtyp_re1 <- as.vector(t(mvrnorm(1, rep(0, n.dat_typ), 
                                  diag(source_var[[2]]))))
  
  
  ##-----Extract data types effects
  data_type_re <- function(da2, da1, ID)
  {
    typ <- rep(1, nrow(da2))
    typ_ID <- as.numeric(ID)
    uniq <- unique(typ_ID)
    for(k in  1:nrow(da2))
    {
      
      for(j in uniq)
      {
        if(typ_ID[k]==uniq[j]) typ[k] = da1[j]
      }
      
    }
    typ
  }
  
  ##--Add to dataset
  ddat$dtyp_re1 <- data_type_re(ddat, dtyp_re1,ddat$source)#--data type re with varying effects
  
  
  dty_re1 <- ddat$dtyp_re1
  # print(AdmUnit[i]^2)
  #####----------ALTERNATIVELY---------########
  #---Simulate Pop Count and Building Count separately and calculate pop density
  #---see the alternative codes below
  ###---Simulate building count
  sigma_b <- 0.55
  epsb <- rnorm(nrow(coord), 0, sigma_b)#--iid random effect for building count
  #betaB <- c(3.51, 0.46, 0.85, 0.21, 0.58, 0.235) #---betas- fixed effects
  betaB <- c(2.21, 0.06, 0.15, 0.21, 0.18, 0.27) 
  bld <- lambdaB <- numeric(nn) #
  for (i in 1:nn)
  {
    lambdaB[i] <- exp(zpred[i,1]*betaB[1] + zpred[i,2]*betaB[2] + 
                        zpred[i,3]*betaB[3] + zpred[i,4]*betaB[4] + 
                        zpred[i,5]*betaB[5] + zpred[i,6]*betaB[6] + sp.raneff[i])
    bld[i] <- rpois(1, lambdaB[i])
  }
  bld
  min(bld); max(bld)
  hist(bld); hist(log(bld))
  mean(bld); var(bld); 
  summary(bld)
  
  
  # set any 0 buildings to 1
  bld[bld==0]=1
  
  ###---Simulate Population count
  sigma_p <- 0.015
  epsp <- rnorm(nrow(coord), 0, sigma_p) 
  betaP <- c(3.50, 0.41, 0.08, 0.04, 0.15, 0.22) #--betas - fixed effects
  
  pop <- lambdaP <- numeric(nn)

  
  
  for (i in 1:nn)
  {
    lambdaP[i] <- exp(zpred[i,1]*betaP[1] + zpred[i,2]*betaP[2] + 
                        zpred[i,3]*betaP[3] + zpred[i,4]*betaP[4] + 
                        zpred[i,5]*betaP[5] + zpred[i,6]*betaP[6]  + sp.raneff[i]+ dty_re1[i]
                      )
    pop[i] <- rpois(1, lambdaP[i])
  }
  sum(pop)
  hist(pop); hist(log(pop))
  mean(pop); var(pop)
  summary(pop)
  # set any 0 population to 1
  pop[pop==0]=1
  
  #--------Add to dataset
  ddat$bld <- bld
  ddat$pop <- pop
  
  # define population density
  ddat$dens <- ddat$pop/ddat$bld #----population density
  hist(ddat$dens); hist(log(ddat$dens))
  
  
  datam <- ddat
  head(ddat)
  names(datam); head(datam)
  datam[,paste0("x",1:6)] <- covs
  head(datam)
  ##---------------------------------------
  ##---Repeat under various proportions of missing observations
  ###--------------------------------------------------------
  for(j in 1:length(cover))###----Do the following for each coverage percentage 
  {
    #set.seed(5060)
    #j = 2
    print(c(spat.autocor[b],paste0(cover[j]*100, "%")))
    
    
    result_path2 <- paste0(result_pathb,"_at_", cover[j]*100,"%_coverage")
    if (file.exists(result_path2)){
      # setwd(file.path(result_path2))
    } else {
      dir.create(file.path(result_path2))
      #setwd(file.path(result_path2))
    } 
    
    
    ##---Select the subset of the data as observed at random.
    pp <- cover[j]
    opts <- sample(nrow(datam), pp*nrow(datam))
    opts <- sort(opts)
    dim(dat.est  <-  datam) #---Locations with observations
    #dim(dat.est  <-  data.frame(datam[opts,])) #---Locations with observations
    
    dat.est$dens[-opts] = NA
    dat.est$eps <- 1:nrow(dat.est)
    names(dat.est)
    #dat.est$dens <- dat.est$pop/dat.est$bld
    dat.est$dens[is.infinite(dat.est$dens)] = NA
    # non Spatial 
    dat.est[,c("x2", "x3", "x4", "x5", "x6")] <- apply(dat.est[,c("x2", "x3", "x4", "x5", "x6")],2, stdize)
    
    table(dat.est$source1 <- as.factor(dat.est$source))
    
    form1 <- dens ~ 1 + x2 + x3 + x4 + x5 + x6 #+ 
      #f(eps, model='iid') #+
      #f(source1, model='iid')
    
    mod1<-inla(form1, #the formula
               data=dat.est,  #the data stack
               family= 'gamma',   #which family the data comes from
               control.predictor=list(compute=TRUE),  #compute gives you the marginals of the linear predictor
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
               verbose = FALSE) 
    
    
    
    ##--------------
    coord8 <- cbind(dat.est$lon, dat.est$lat)
    dim(coord8)
    
    
    ##---Build mesh at grid cell level at the coverage level
    bnd <- inla.nonconvex.hull(coord8, -0.035, -0.04, resolution = c(100, 100))
    mesh <- inla.mesh.2d(boundary = bnd, max.edge=c(0.2,1), 
                         offset = c(0.2, 0.7),
                         cutoff = 0.2)
    par(mfrow=c(1,1))
    plot(mesh)
    mesh$n # 1133 mesh nodes
    
    ###---Build projector matrix A 
    A<-inla.spde.make.A(mesh=mesh,loc=as.matrix(coord8));dim(A)
    
    ##---Create the SPDE
    spde <- inla.spde2.matern(mesh, alpha=2)
    
    ##----specify the observation indices for estimation 
    iset<- inla.spde.make.index(name = "spatial.field", spde$n.spde)
    
   covars.est <- dat.est[,c("x1","x2", "x3", "x4", "x5", "x6", "source1", "eps")]; dim(covars.est) ##---Population density
    head(covars.est)
  
    head(dat.est)
    ##---Build data stack 
    stk.est <- inla.stack(data=list(y=dat.est$dens), #the response
                          
                          A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                          
                          effects=list(c(list(Intercept=1), #the Intercept
                                         iset),  #the spatial index
                                       #the covariates
                                       list(covars.est)
                          ), 
                          #this is a quick name so you can call upon easily
                          tag='est')
    
    # Spatial 
    
    form2 <- y ~ -1 + Intercept + x2 + x3 + x4 + x5 + x6 + 
      #f(eps,  model = 'iid') +
      f(source1, model='iid') +
      f(spatial.field, model=spde)
    
    mod2 <-inla(form2, #the formula
                data=inla.stack.data(stk.est,spde=spde),  #the data stack
                family= 'gamma',   #which family the data comes from
                control.predictor=list(A=inla.stack.A(stk.est),compute=TRUE),  #compute gives you the marginals of the linear predictor
                control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                verbose = FALSE) 
    
    
    # Initial model fit checks
    (DIC <- t(c(mod1=mod1$dic$dic, mod2=mod2$dic$dic)))
    
    
    
    # extract posterior values from the non spatial model      
    index <-inla.stack.index(stk.est, "est")$data #---extract the data location indices
    pred2<- exp(mod2$summary.linear.predictor[index ,"mean"]) #--predicted mean count
    sum(dat.est$fit2 <- pred2*dat.est$bld)
    sum(dat.est$pop)
    plot(dat.est$pop, dat.est$fit2)
    cor(dat.est$pop, dat.est$fit2)
    
    
    
    # extract posterior values from the non spatial model
    pred1<- exp(mod1$summary.linear.predictor$mean) #--predicted mean count
    sum(dat.est$fit1 <- pred1*dat.est$bld)
    sum(dat.est$pop)
    plot(dat.est$pop, dat.est$fit1)
    cor(dat.est$pop, dat.est$fit1)
    
    met1a <- unlist(model_metrics(dat.est$pop, 
                                  dat.est$fit1))
    met2a <- unlist(model_metrics(dat.est$pop, 
                                  dat.est$fit2))
    
    
    print(metrics1a <- rbind(met1a,met2a))
    #----------------------------------------------------------------------
    #  Run posterior simulation and predictions simulataneously
    #-----------------------------------------------------------------
    
    # Prediction Projection matrix 
    Aprd <- inla.spde.make.A(mesh = mesh, loc = as.matrix(coord))
    dim(Aprd)#--prediction projectionnmatrix
    
    dat.pred <- datam
    

    #mean(dat.pred$source_re); var(dat.pred$source_re)
    dat.pred[,c("x2", "x3", "x4", "x5", "x6")] <- apply(dat.pred[,c("x2", "x3", "x4", "x5", "x6")],2, stdize)
    
    head(datam)
    head(dat.pred)
    ####-------Run posterior and grid cell prediction (GRID2GRID)
    sim_spatial <- function(model, dat, Apred, run) # spatial model
    {
      fixedeff <- dens_hat <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
    #inla.seed = 1663336
      inla.seed = as.integer(runif(1)*.Machine$integer.max)
      set.seed(inla.seed)
      print(inla.seed)
      m1.samp <- inla.posterior.sample(run, model, 
                                       seed = inla.seed, 
                                       selection=list(x2=1, x3=1, x4=1,x5=1, x6=1),
                                       num.threads="1:1")
      
      sfield_nodes_mean <- model$summary.random$spatial.field['mean']
      field_mean <- (Apred%*% as.data.frame(sfield_nodes_mean)[, 1])
      for(i in 1:run)
      {
        fixedeff[,i] <- 
          model$summary.fixed['Intercept', 'mean'] +
          m1.samp[[i]]$latent[1,] * dat[,'x2'] +
          m1.samp[[i]]$latent[2,] * dat[,'x3'] +
          m1.samp[[i]]$latent[3,] * dat[,'x4'] +
          m1.samp[[i]]$latent[4,] * dat[,'x5'] +
          m1.samp[[i]]$latent[5,] * dat[,'x6'] + 
          
          #dat$source_re2 + #Data source random effect
          #rnorm(nrow(dat), 0, 1/model$summary.hyperpar$mean[2]) + #IID
          rnorm(nrow(dat), 0, 1/model$summary.hyperpar$mean[2]) + # data source
          field_mean[,1]
        
        ##
        dens_hat[,i] <- exp(fixedeff[,i])
        pop_hat[,i] <- dens_hat[,i]*dat$bld
        
      }
      
      dat$mean_dens_hat <- apply(dens_hat, 1, mean, na.rm=T) # mean population density
      dat$mean_pop_hat <- apply(pop_hat, 1, mean, na.rm=T) # mean population count 
      dat$lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) # lower bound
      dat$upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) # upper bound
      dat$sd_pop_hat <- apply(pop_hat, 1, sd, na.rm=T) # standard deviation
      dat$cv_pop_hat <- dat$sd_pop_hat/dat$mean_pop_hat # coefficient of variation
      
      
      output <- list(est_data = dat)
      
    }
    
    
    
    sim_nonspatial <- function(model, dat, run) # non spatial model
    {
      fixedeff <- dens_hat <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
      #inla.seed = 16646
      inla.seed = as.integer(runif(1)*.Machine$integer.max)
      set.seed(inla.seed)
      print(inla.seed)
      m1.samp <- inla.posterior.sample(run, model, 
                                       seed = inla.seed, 
                                       selection=list(x2=1, x3=1, x4=1,x5=1, x6=1),
                                       num.threads="1:1")
      
      for(i in 1:run)
      {
        fixedeff[,i] <- 
          model$summary.fixed['(Intercept)', 'mean'] +
          m1.samp[[i]]$latent[1,] * dat[,'x2'] +
          m1.samp[[i]]$latent[2,] * dat[,'x3'] +
          m1.samp[[i]]$latent[3,] * dat[,'x4'] +
          m1.samp[[i]]$latent[4,] * dat[,'x5'] +
          m1.samp[[i]]$latent[5,] * dat[,'x6'] #+ 
          
          #dat$source_re1 + #Data source random effect
         # rnorm(nrow(dat), 0, 1/model$summary.hyperpar$mean[2]) #+# IID
          #rnorm(nrow(dat), 0, 1/model$summary.hyperpar$mean[3]) # data source
        
        ##
        dens_hat[,i] <- exp(fixedeff[,i])
        pop_hat[,i] <- dens_hat[,i]*dat$bld
        
      }
      
      dat$mean_dens_hat <- apply(dens_hat, 1, mean, na.rm=T) # mean population density
      dat$mean_pop_hat <- apply(pop_hat, 1, mean, na.rm=T) # mean population count 
      dat$lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) # lower bound
      dat$upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) # upper bound
      dat$sd_pop_hat <- apply(pop_hat, 1, sd, na.rm=T) # standard deviation
      dat$cv_pop_hat <- dat$sd_pop_hat/dat$mean_pop_hat # coefficient of variation
      
      
      output <- list(est_data = dat)
      
    }
    run=100
    
    
    system.time(str(sim.nspat <- sim_nonspatial(mod1,dat.pred, run))) # non-spatial
    system.time(str(sim.spat <- sim_spatial(mod2,dat.pred,Aprd, run))) # spatial
    
    dat1 <- sim.nspat$est_data
    dat1$method <- rep("Base", nrow(dat1))
    dat2 <- sim.spat$est_data
    dat2$method <- rep("Proposed", nrow(dat2))
    
    # combine
    ddtt <- rbind(dat1,dat2)
    write.csv(ddtt, file=paste0(result_path2,"/post_dat.csv"))
    
    ##---Calculate model fit and model cross-validation metrics  
    met1 <- unlist(model_metrics(dat1$pop, 
                                 dat1$mean_pop_hat))
    met2 <- unlist(model_metrics(dat2$pop, 
                                 dat2$mean_pop_hat))
    
    
  metrics <- rbind(met1,met2)
    
    par(mfrow=c(1,2))
    plot(dat1$pop, dat1$mean_pop_hat)
    plot(dat2$pop, dat2$mean_pop_hat)
    par(mfrow=c(1,1))
    write.csv(metrics1a, file=paste0(result_path2,"/fit_metrics.csv"))
    write.csv(metrics, file=paste0(result_path2,"/fit_metrics2.csv"))
  }
  
}

 
#Run posterior plots for different data source variance models
#-----------------------------------------------------------------------------------
# Extract posterior fit indices and predictions when data sources have same variability
#-----------------------------------------------------------------------------------
metsb1 <- metsb2 <- list()
for(b in 1:length(spat.autocor))
{ 
  
  result_pathb <- paste0(out_path,"/marginal_var_", spat.autocor[b])
  
  for(j in 1:length(cover))###----Do the following for each coverage percentage 
  {
    
    result_path2 <- paste0(result_pathb,"_at_", cover[j]*100,"%_coverage")
    
    # fit model fit measures
    mod_fit <- read.csv(paste0(result_path2,"/fit_metrics.csv"))
    mod_fit$cover <- rep(cover[j]*100,2)
    mod_fit$spat_cor <- rep(spat.autocor[b],2)
    mod_fit$method <- c("Base", "Proposed")
    metsb1[[j]] <- mod_fit
  }
  metsb2[[b]] <- metsb1 
}


## fit metrics
unnest_met <- unlist(metsb2, recursive = FALSE)  #--unnest the list
mod.mets <- do.call(rbind, unnest_met)

# exclude 5% and 100% missingned for more meaningful comparison
metrics<- mod.mets %>% dplyr::select(-X) %>% filter(!cover %in% c(5,100))

write.csv(metrics, paste0(out_path,"/combined_fit_metrics.csv"))#save the combined data

#---Covert to long format
# install.packages("reshape2")
library(reshape2)
dim(met_long <- melt(metrics, id.vars=c("method","cover", "spat_cor"),
                     value.name="estimate", variable.name = "metric"))

met_long$method = factor(met_long$method)
met_long$cover = factor(met_long$cover)
met_long$spat_cor = factor(met_long$spat_cor)


met_long <- met_long %>% filter(spat_cor !="0.75")
levels(met_long$method) <- c("Base","Proposed")
write.csv(met_long, "combined_fit_metrics_long.csv", row.names=FALSE)

library(ggpubr)



##----- CORR
dim(corr <- met_long[met_long$metric=="corr",])

plot_corr <- ggline(corr, x = "spat_cor", y = "estimate",
                    error.plot = "estimate",
                    facet.by = "method",
                    panel.labs.font.x = list(size=25),
                    color = "cover",
                    palette = "aaas",
                    point.size=2,
                    #linetype = "method",
                    size=2)+
  theme_bw()+
  theme(strip.text = element_text(size=25))

rcorr <-  ggpar(plot_corr, xlab="Spatial Variance", ylab="Correlation Coefficient (CC)",
                legend = "top", legend.title = "Survey Coverage (%)",size=25,
                font.legend=c(25),
                font.label = list(size = 20, face = "bold", color ="red"),
                font.x = c(25),
                font.y = c(25),
                font.main=c(20),
                font.xtickslab =c(20),
                font.ytickslab =c(20),
                # orientation = "reverse",
                xtickslab.rt = 45, ytickslab.rt = 45)
rcorr 



##----- RMSE
dim(rmse <- met_long[met_long$metric=="RMSE",])

plot_rmse <- ggline(rmse, x = "spat_cor", y = "estimate",
                    error.plot = "estimate",
                    facet.by = "method",
                    panel.labs.font.x = list(size=25),
                    color = "cover",
                    palette = "aaas",
                    point.size=2,
                    #linetype = "method",
                    size=2)+
  theme_bw()+
  theme(strip.text = element_text(size=25))

rrmse <-  ggpar(plot_rmse, xlab="Spatial Variance", ylab="Root Mean Square Error (RMSE)",
                legend = "top", legend.title = "Survey Coverage (%)",size=25,
                font.legend=c(25),
                font.label = list(size = 20, face = "bold", color ="red"),
                font.x = c(25),
                font.y = c(25),
                font.main=c(20),
                font.xtickslab =c(20),
                font.ytickslab =c(20),
                # orientation = "reverse",
                xtickslab.rt = 45, ytickslab.rt = 45)
rrmse 

##----- MAE
dim(mae <- met_long[met_long$metric=="MAE",])

plot_mae <- ggline(mae, x = "spat_cor", y = "estimate",
                   error.plot = "estimate",
                   facet.by = "method",
                   # panel.labs= list(dat_type=c("diff", "same")),
                   panel.labs.font.x = list(size=25),
                   color = "cover",palette = "aaas",
                   point.size=2,
                   #linetype = "pop_cover",
                   size=2)+
  theme_bw()+
  theme(strip.text = element_text(size=25))
rmae <-  ggpar(plot_mae, xlab="Spatial Variance", ylab="Mean Absolute Error (MAE)",
               legend = "none", 
               font.legend=c(25),
               font.label = list(size = 20, face = "bold", color ="red"),
               font.x = c(25),
               font.y = c(25),
               font.main=c(22),
               font.xtickslab =c(22),
               font.ytickslab =c(22),
               # orientation = "reverse",
               xtickslab.rt = 45, ytickslab.rt = 45)
rmae 



##----- BIAS
dim(bias <- met_long[met_long$metric=="BIAS",])

plot_bias <- ggline(bias, x = "spat_cor", y = "estimate",
                   error.plot = "estimate",
                   facet.by = "method",
                   # panel.labs= list(dat_type=c("diff", "same")),
                   panel.labs.font.x = list(size=25),
                   color = "cover",palette = "aaas",
                   point.size=2,
                   #linetype = "pop_cover",
                   size=2)+
  theme_bw()+
  theme(strip.text = element_text(size=25))
rbias <-  ggpar(plot_bias, xlab="Spatial Variance", ylab="Absolute Bias (Abias)",
               legend = "none", 
               font.legend=c(25),
               font.label = list(size = 20, face = "bold", color ="red"),
               font.x = c(25),
               font.y = c(25),
               font.main=c(22),
               font.xtickslab =c(22),
               font.ytickslab =c(22),
               # orientation = "reverse",
               xtickslab.rt = 45, ytickslab.rt = 45)
rbias 

ggarrange(rrmse,rmae, 
          ncol = 1, nrow = 2)
#------------------------------------------------------------
#### --- Perecentage reductions in RMSE and MAE
#----------------------------------------------------------
### Reduction in relative error measures
rrem <- as.data.frame(read.csv(paste0(out_path, "/error_reduction.csv"))) # calculate this from the excel fivle saved in the outpute folder
rr_dat <- rrem %>% drop_na(rrmae)%>% dplyr::select(c(cover,
                                                     spat_cor,
                                                     rrmae, rrrmse))


rr_dat <- rr_dat %>% filter(spat_cor !="0.75")
rr_dat$rrmae <- round(rr_dat$rrmae*100, 2)
rr_dat$rrrmse <- round(rr_dat$rrrmse*100, 2)
rr_dat$spat_cor <- factor(rr_dat$spat_cor)
rr_dat$cover <- factor(rr_dat$cover)
rr_dat$cover <- factor(rr_dat$cover,
                       labels=c("20%", "40%",
                                "60%", "80%"))


# rrmae
bar_mae <- rr_dat %>%
  ggplot( aes(x=cover, y=rrmae, fill=spat_cor)) +
  geom_bar(stat = "identity", 
           position = position_dodge()) +
  #scale_x_continuous(breaks = seq(20,90, by=10))+
  #scale_y_continuous(breaks=seq(500,7000,1000))+
  #coord_flip()+
  theme_bw()+
  theme(strip.text = element_text(size = 16),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))

rbar_mae <-  ggpar(bar_mae, xlab="Proportion of population observed(%)", 
                   ylab="Reduction in relative MAE(%)",
                   legend = "none", legend.title=element_text("Spatial Variance:"),
                   font.legend=c(25),
                   palette = c("aaas"),
                   font.label = list(size = 25, face = "bold", color ="red"),
                   font.x = c(25),
                   font.y = c(25),
                   font.main=c(25),
                   font.xtickslab =c(22),
                   font.ytickslab =c(22),
                   xtickslab.rt = 45, ytickslab.rt = 45)
rbar_mae



# rrrmse
bar_rmse <- rr_dat %>%
  ggplot( aes(x=cover, y=rrrmse, fill=spat_cor)) +
  geom_bar(stat = "identity", 
           position = position_dodge()) +
  #scale_y_continuous(breaks=seq(500,7000,1000))+
  #coord_flip()+
  theme_bw()+
  theme(strip.text = element_text(size = 16),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))

rbar_rmse <-  ggpar(bar_rmse, xlab="Proportion of population observed(%)", 
                    ylab="Reduction in relative RMSE(%)",
                    legend = "top", legend.title=element_text("Spatial Variance:"),
                    font.legend=c(25),
                    palette = c("aaas"),
                    font.label = list(size = 25, face = "bold", color ="red"),
                    font.x = c(25),
                    font.y = c(25),
                    font.main=c(25),
                    font.xtickslab =c(22),
                    font.ytickslab =c(22),
                    xtickslab.rt = 45, ytickslab.rt = 45)
rbar_rmse



ggarrange(rbar_rmse, rbar_mae, 
          ncol = 1, nrow = 2)


min(rr_dat$rrmae); max(rr_dat$rrmae)
min(rr_dat$rrrmse); max(rr_dat$rrrmse)

# for models where all the data sources have the same variance
# reduction in relative mae ranged from 17.65% to 25.35%
# while reduction in relative rmse ranged from 13.60% to 27.02%

# for models where all the data sources have the different variances
# reduction in relative mae ranged from 21.39% to 31.21%
# while reduction in relative rmse ranged from 18.71% to 28.08%


#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
### Same data source variance
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3/Working/CMR/Chris_N/paper1"
shp_path <- paste0(path, "/application/data/Input_Settlement_Boundaries")
out_path <- paste0(path, "/sim_study/output_revised2")
#  Load data       

for(b in 1:length(spat.autocor))
{ 
  
  # b=1
  result_pathb <- paste0(out_path,"/marginal_var_", spat.autocor[b])
  if (file.exists(result_pathb)){
    #setwd(file.path(result_pathb))
  } else {
    dir.create(file.path(result_pathb))
    #setwd(file.path(result_pathb))
  }
  
  sigma0 <- spat.autocor[b]#--marginal variance  
  kappa0 <- sqrt(8*nu)/rho #--scale parameters
  (tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)) #--precision
  
  #---SPDE
  spde <- inla.spde2.matern(mesh, 
                            B.tau = matrix(c(log(tau0), -1, +1), nrow=1, ncol=3),
                            B.kappa = matrix(c(log(kappa0), 0, -1), nrow=1, ncol=3),
                            theta.prior.mean = c(0,0),
                            theta.prior.prec = c(0.1, 0.1))
  
  
  #--Precision Matrix
  Q <- inla.spde2.precision(spde=spde, theta=c(0,0))
  
  #---Simulate the GMRF
  sam <- as.vector(inla.qsample(
    n = 1, Q = Q, seed=100))
  #length(sam)
  
  ###---Build projector matrix A
  A <- inla.spde.make.A(mesh=mesh, loc=coord);dim(A)
  
  # Spatial Random effects
  sp.raneff <- as.vector(A %*% sam)
  length(sp.raneff)
  hist(sp.raneff)
  
  ##----specify the observation indices for estimation 
  iset <- inla.spde.make.index(name = "spatial.field", spde$n.spde)
  
  
  #-- generate the Covariates 
  nn <- nrow(dat.all)
  covs = cbind(1,runif(nn,1,10),rnorm(nn, 2, 2),rpois(nn,4),rnorm(nn, 4, 3), runif(nn))
  apply(covs, 2, mean) #---check averages
  dim(dat.all)
  
  # Standardization 
  stdize <- function(x)
  {
    #return(stdze <- (x-min(x, na.rm=T))/(max(x, na.rm =T)-min(x, na.rm =T)))
    return(stdze <- (x-mean(x, na.rm=T))/(sd(x, na.rm =T)))
  }
  
  covs_std <- cbind(covs[,1],apply(covs[,-1], 2, stdize)) 
  
  # Add covariates to data
  head(dat.all)
  rownames(dat.all) <- NULL
  dim(ddat <- cbind(dat.all, covs_std))
  ddat <- data.frame(ddat)
  colnames(ddat) <- c("lon", "lat", "source", "x1", "x2", "x3", "x4", "x5", "x6")
  names(ddat)
  head(ddat)
  ###############
  #Parameters 
  #
  #dim(zpred <- covs_std)
  dim(zpred <- covs)
  
  require(MASS)
  #--generate data type random effects 
  
  
  dtyp_re1 <- as.vector(t(mvrnorm(1, rep(0, n.dat_typ), 
                                  diag(source_var[[1]]))))# same data source variance
  
  
  ##-----Extract data types effects
  data_type_re <- function(da2, da1, ID)
  {
    typ <- rep(1, nrow(da2))
    typ_ID <- as.numeric(ID)
    uniq <- unique(typ_ID)
    for(k in  1:nrow(da2))
    {
      
      for(j in uniq)
      {
        if(typ_ID[k]==uniq[j]) typ[k] = da1[j]
      }
      
    }
    typ
  }
  
  ##--Add to dataset
  ddat$dtyp_re1 <- data_type_re(ddat, dtyp_re1,ddat$source)#--data type re with varying effects
  
  
  dty_re1 <- ddat$dtyp_re1
  # print(AdmUnit[i]^2)
  #####----------ALTERNATIVELY---------########
  #---Simulate Pop Count and Building Count separately and calculate pop density
  #---see the alternative codes below
  ###---Simulate building count
  sigma_b <- 0.55
  epsb <- rnorm(nrow(coord), 0, sigma_b)#--iid random effect for building count
  #betaB <- c(3.51, 0.46, 0.85, 0.21, 0.58, 0.235) #---betas- fixed effects
  betaB <- c(2.21, 0.06, 0.15, 0.21, 0.18, 0.27) 
  bld <- lambdaB <- numeric(nn) #
  for (i in 1:nn)
  {
    lambdaB[i] <- exp(zpred[i,1]*betaB[1] + zpred[i,2]*betaB[2] + 
                        zpred[i,3]*betaB[3] + zpred[i,4]*betaB[4] + 
                        zpred[i,5]*betaB[5] + zpred[i,6]*betaB[6] + sp.raneff[i])
    bld[i] <- rpois(1, lambdaB[i])
  }
  bld
  min(bld); max(bld)
  hist(bld); hist(log(bld))
  mean(bld); var(bld); 
  summary(bld)
  
  
  # set any 0 buildings to 1
  bld[bld==0]=1
  
  ###---Simulate Population count
  sigma_p <- 0.015
  epsp <- rnorm(nrow(coord), 0, sigma_p) 
  betaP <- c(3.50, 0.41, 0.08, 0.04, 0.15, 0.22) #--betas - fixed effects
  
  pop <- lambdaP <- numeric(nn)
  
  
  
  for (i in 1:nn)
  {
    lambdaP[i] <- exp(zpred[i,1]*betaP[1] + zpred[i,2]*betaP[2] + 
                        zpred[i,3]*betaP[3] + zpred[i,4]*betaP[4] + 
                        zpred[i,5]*betaP[5] + zpred[i,6]*betaP[6]  + sp.raneff[i]+ dty_re1[i]
    )
    pop[i] <- rpois(1, lambdaP[i])
  }
  sum(pop)
  hist(pop); hist(log(pop))
  mean(pop); var(pop)
  summary(pop)
  # set any 0 population to 1
  pop[pop==0]=1
  
  #--------Add to dataset
  ddat$bld <- bld
  ddat$pop <- pop
  
  # define population density
  ddat$dens <- ddat$pop/ddat$bld #----population density
  hist(ddat$dens); hist(log(ddat$dens))
  
  
  datam <- ddat
  head(ddat)
  names(datam); head(datam)
  datam[,paste0("x",1:6)] <- covs
  head(datam)
  ##---------------------------------------
  ##---Repeat under various proportions of missing observations
  ###--------------------------------------------------------
  for(j in 1:length(cover))###----Do the following for each coverage percentage 
  {
    #set.seed(5060)
    #j = 2
    print(c(spat.autocor[b],paste0(cover[j]*100, "%")))
    
    
    result_path2 <- paste0(result_pathb,"_at_", cover[j]*100,"%_coverage")
    if (file.exists(result_path2)){
      # setwd(file.path(result_path2))
    } else {
      dir.create(file.path(result_path2))
      #setwd(file.path(result_path2))
    } 
    
    
    ##---Select the subset of the data as observed at random.
    pp <- cover[j]
    opts <- sample(nrow(datam), pp*nrow(datam))
    opts <- sort(opts)
    dim(dat.est  <-  datam) #---Locations with observations
    #dim(dat.est  <-  data.frame(datam[opts,])) #---Locations with observations
    
    dat.est$dens[-opts] = NA
    dat.est$eps <- 1:nrow(dat.est)
    names(dat.est)
    #dat.est$dens <- dat.est$pop/dat.est$bld
    dat.est$dens[is.infinite(dat.est$dens)] = NA
    # non Spatial 
    dat.est[,c("x2", "x3", "x4", "x5", "x6")] <- apply(dat.est[,c("x2", "x3", "x4", "x5", "x6")],2, stdize)
    
    table(dat.est$source1 <- as.factor(dat.est$source))
    
    form1 <- dens ~ 1 + x2 + x3 + x4 + x5 + x6 #+ 
    #f(eps, model='iid') #+
    #f(source1, model='iid')
    
    mod1<-inla(form1, #the formula
               data=dat.est,  #the data stack
               family= 'gamma',   #which family the data comes from
               control.predictor=list(compute=TRUE),  #compute gives you the marginals of the linear predictor
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
               verbose = FALSE) 
    
    
    
    ##--------------
    coord8 <- cbind(dat.est$lon, dat.est$lat)
    dim(coord8)
    
    
    ##---Build mesh at grid cell level at the coverage level
    bnd <- inla.nonconvex.hull(coord8, -0.035, -0.04, resolution = c(100, 100))
    mesh <- inla.mesh.2d(boundary = bnd, max.edge=c(0.2,1), 
                         offset = c(0.2, 0.7),
                         cutoff = 0.2)
    par(mfrow=c(1,1))
    plot(mesh)
    mesh$n # 1133 mesh nodes
    
    ###---Build projector matrix A 
    A<-inla.spde.make.A(mesh=mesh,loc=as.matrix(coord8));dim(A)
    
    ##---Create the SPDE
    spde <- inla.spde2.matern(mesh, alpha=2)
    
    ##----specify the observation indices for estimation 
    iset<- inla.spde.make.index(name = "spatial.field", spde$n.spde)
    
    covars.est <- dat.est[,c("x1","x2", "x3", "x4", "x5", "x6", "source1", "eps")]; dim(covars.est) ##---Population density
    head(covars.est)
    
    head(dat.est)
    ##---Build data stack 
    stk.est <- inla.stack(data=list(y=dat.est$dens), #the response
                          
                          A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                          
                          effects=list(c(list(Intercept=1), #the Intercept
                                         iset),  #the spatial index
                                       #the covariates
                                       list(covars.est)
                          ), 
                          #this is a quick name so you can call upon easily
                          tag='est')
    
    # Spatial 
    
    form2 <- y ~ -1 + Intercept + x2 + x3 + x4 + x5 + x6 + 
      #f(eps,  model = 'iid') +
      f(source1, model='iid') +
      f(spatial.field, model=spde)
    
    mod2 <-inla(form2, #the formula
                data=inla.stack.data(stk.est,spde=spde),  #the data stack
                family= 'gamma',   #which family the data comes from
                control.predictor=list(A=inla.stack.A(stk.est),compute=TRUE),  #compute gives you the marginals of the linear predictor
                control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                verbose = FALSE) 
    
    
    # Initial model fit checks
    (DIC <- t(c(mod1=mod1$dic$dic, mod2=mod2$dic$dic)))
    
    
    
    # extract posterior values from the non spatial model      
    index <-inla.stack.index(stk.est, "est")$data #---extract the data location indices
    pred2<- exp(mod2$summary.linear.predictor[index ,"mean"]) #--predicted mean count
    sum(dat.est$fit2 <- pred2*dat.est$bld)
    sum(dat.est$pop)
    plot(dat.est$pop, dat.est$fit2)
    cor(dat.est$pop, dat.est$fit2)
    
    
    
    # extract posterior values from the non spatial model
    pred1<- exp(mod1$summary.linear.predictor$mean) #--predicted mean count
    sum(dat.est$fit1 <- pred1*dat.est$bld)
    sum(dat.est$pop)
    plot(dat.est$pop, dat.est$fit1)
    cor(dat.est$pop, dat.est$fit1)
    
    met1a <- unlist(model_metrics(dat.est$pop, 
                                  dat.est$fit1))
    met2a <- unlist(model_metrics(dat.est$pop, 
                                  dat.est$fit2))
    
    
    print(metrics1a <- rbind(met1a,met2a))
    #----------------------------------------------------------------------
    #  Run posterior simulation and predictions simulataneously
    #-----------------------------------------------------------------
    
    # Prediction Projection matrix 
    Aprd <- inla.spde.make.A(mesh = mesh, loc = as.matrix(coord))
    dim(Aprd)#--prediction projectionnmatrix
    
    dat.pred <- datam
    
    
    #mean(dat.pred$source_re); var(dat.pred$source_re)
    dat.pred[,c("x2", "x3", "x4", "x5", "x6")] <- apply(dat.pred[,c("x2", "x3", "x4", "x5", "x6")],2, stdize)
    
    head(datam)
    head(dat.pred)
    ####-------Run posterior and grid cell prediction (GRID2GRID)
    sim_spatial <- function(model, dat, Apred, run) # spatial model
    {
      fixedeff <- dens_hat <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
      #inla.seed = 1663336
      inla.seed = as.integer(runif(1)*.Machine$integer.max)
      set.seed(inla.seed)
      print(inla.seed)
      m1.samp <- inla.posterior.sample(run, model, 
                                       seed = inla.seed, 
                                       selection=list(x2=1, x3=1, x4=1,x5=1, x6=1),
                                       num.threads="1:1")
      
      sfield_nodes_mean <- model$summary.random$spatial.field['mean']
      field_mean <- (Apred%*% as.data.frame(sfield_nodes_mean)[, 1])
      for(i in 1:run)
      {
        fixedeff[,i] <- 
          model$summary.fixed['Intercept', 'mean'] +
          m1.samp[[i]]$latent[1,] * dat[,'x2'] +
          m1.samp[[i]]$latent[2,] * dat[,'x3'] +
          m1.samp[[i]]$latent[3,] * dat[,'x4'] +
          m1.samp[[i]]$latent[4,] * dat[,'x5'] +
          m1.samp[[i]]$latent[5,] * dat[,'x6'] + 
          
          #dat$source_re2 + #Data source random effect
          #rnorm(nrow(dat), 0, 1/model$summary.hyperpar$mean[2]) + #IID
          rnorm(nrow(dat), 0, 1/model$summary.hyperpar$mean[2]) + # data source
          field_mean[,1]
        
        ##
        dens_hat[,i] <- exp(fixedeff[,i])
        pop_hat[,i] <- dens_hat[,i]*dat$bld
        
      }
      
      dat$mean_dens_hat <- apply(dens_hat, 1, mean, na.rm=T) # mean population density
      dat$mean_pop_hat <- apply(pop_hat, 1, mean, na.rm=T) # mean population count 
      dat$lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) # lower bound
      dat$upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) # upper bound
      dat$sd_pop_hat <- apply(pop_hat, 1, sd, na.rm=T) # standard deviation
      dat$cv_pop_hat <- dat$sd_pop_hat/dat$mean_pop_hat # coefficient of variation
      
      
      output <- list(est_data = dat)
      
    }
    
    
    
    sim_nonspatial <- function(model, dat, run) # non spatial model
    {
      fixedeff <- dens_hat <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
      #inla.seed = 16646
      inla.seed = as.integer(runif(1)*.Machine$integer.max)
      set.seed(inla.seed)
      print(inla.seed)
      m1.samp <- inla.posterior.sample(run, model, 
                                       seed = inla.seed, 
                                       selection=list(x2=1, x3=1, x4=1,x5=1, x6=1),
                                       num.threads="1:1")
      
      for(i in 1:run)
      {
        fixedeff[,i] <- 
          model$summary.fixed['(Intercept)', 'mean'] +
          m1.samp[[i]]$latent[1,] * dat[,'x2'] +
          m1.samp[[i]]$latent[2,] * dat[,'x3'] +
          m1.samp[[i]]$latent[3,] * dat[,'x4'] +
          m1.samp[[i]]$latent[4,] * dat[,'x5'] +
          m1.samp[[i]]$latent[5,] * dat[,'x6'] #+ 
        
        #dat$source_re1 + #Data source random effect
        # rnorm(nrow(dat), 0, 1/model$summary.hyperpar$mean[2]) #+# IID
        #rnorm(nrow(dat), 0, 1/model$summary.hyperpar$mean[3]) # data source
        
        ##
        dens_hat[,i] <- exp(fixedeff[,i])
        pop_hat[,i] <- dens_hat[,i]*dat$bld
        
      }
      
      dat$mean_dens_hat <- apply(dens_hat, 1, mean, na.rm=T) # mean population density
      dat$mean_pop_hat <- apply(pop_hat, 1, mean, na.rm=T) # mean population count 
      dat$lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) # lower bound
      dat$upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) # upper bound
      dat$sd_pop_hat <- apply(pop_hat, 1, sd, na.rm=T) # standard deviation
      dat$cv_pop_hat <- dat$sd_pop_hat/dat$mean_pop_hat # coefficient of variation
      
      
      output <- list(est_data = dat)
      
    }
    run=100
    
    
    system.time(str(sim.nspat <- sim_nonspatial(mod1,dat.pred, run))) # non-spatial
    system.time(str(sim.spat <- sim_spatial(mod2,dat.pred,Aprd, run))) # spatial
    
    dat1 <- sim.nspat$est_data
    dat1$method <- rep("Base", nrow(dat1))
    dat2 <- sim.spat$est_data
    dat2$method <- rep("Proposed", nrow(dat2))
    
    # combine
    ddtt <- rbind(dat1,dat2)
    write.csv(ddtt, file=paste0(result_path2,"/post_dat.csv"))
    
    ##---Calculate model fit and model cross-validation metrics  
    met1 <- unlist(model_metrics(dat1$pop, 
                                 dat1$mean_pop_hat))
    met2 <- unlist(model_metrics(dat2$pop, 
                                 dat2$mean_pop_hat))
    
    
    metrics <- rbind(met1,met2)
    
    par(mfrow=c(1,2))
    plot(dat1$pop, dat1$mean_pop_hat)
    plot(dat2$pop, dat2$mean_pop_hat)
    par(mfrow=c(1,1))
    write.csv(metrics1a, file=paste0(result_path2,"/fit_metrics.csv"))
    write.csv(metrics, file=paste0(result_path2,"/fit_metrics2.csv"))
  }
  
}


#Run posterior plots for different data source variance models
# ==========================================================
# Specify key data paths
path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3/Working/CMR/Chris_N/paper1"
shp_path <- paste0(path, "/application/data/Input_Settlement_Boundaries")
out_path <- paste0(path, "/sim_study/output_revised2")

#-----------------------------------------------------------------------------------
# Extract posterior fit indices and predictions when data sources have same variability
#-----------------------------------------------------------------------------------
metsb1 <- metsb2 <- list()
for(b in 1:length(spat.autocor))
{ 
  
  result_pathb <- paste0(out_path,"/marginal_var_", spat.autocor[b])
  
  for(j in 1:length(cover))###----Do the following for each coverage percentage 
  {
    
    result_path2 <- paste0(result_pathb,"_at_", cover[j]*100,"%_coverage")
    
    # fit model fit measures
    mod_fit <- read.csv(paste0(result_path2,"/fit_metrics.csv"))
    mod_fit$cover <- rep(cover[j]*100,2)
    mod_fit$spat_cor <- rep(spat.autocor[b],2)
    mod_fit$method <- c("Base", "Proposed")
    metsb1[[j]] <- mod_fit
  }
  metsb2[[b]] <- metsb1 
}


## fit metrics
unnest_met <- unlist(metsb2, recursive = FALSE)  #--unnest the list
mod.mets <- do.call(rbind, unnest_met)

# exclude 5% and 100% missingned for more meaningful comparison
metrics<- mod.mets %>% dplyr::select(-X) %>% filter(!cover %in% c(5,100))

write.csv(metrics, paste0(out_path,"/combined_fit_metrics.csv"))#save the combined data

#---Covert to long format
# install.packages("reshape2")
library(reshape2)
dim(met_long <- melt(metrics, id.vars=c("method","cover", "spat_cor"),
                     value.name="estimate", variable.name = "metric"))

met_long$method = factor(met_long$method)
met_long$cover = factor(met_long$cover)
met_long$spat_cor = factor(met_long$spat_cor)


met_long <- met_long %>% filter(spat_cor !="0.75")
levels(met_long$method) <- c("Base","Proposed")
write.csv(met_long, "combined_fit_metrics_long.csv", row.names=FALSE)

library(ggpubr)



##----- CORR
dim(corr <- met_long[met_long$metric=="corr",])

plot_corr <- ggline(corr, x = "spat_cor", y = "estimate",
                    error.plot = "estimate",
                    facet.by = "method",
                    panel.labs.font.x = list(size=25),
                    color = "cover",
                    palette = "aaas",
                    point.size=2,
                    #linetype = "method",
                    size=2)+
  theme_bw()+
  theme(strip.text = element_text(size=25))

rcorr <-  ggpar(plot_corr, xlab="Spatial Variance", ylab="Correlation Coefficient (CC)",
                legend = "top", legend.title = "Survey Coverage (%)",size=25,
                font.legend=c(25),
                font.label = list(size = 20, face = "bold", color ="red"),
                font.x = c(25),
                font.y = c(25),
                font.main=c(20),
                font.xtickslab =c(20),
                font.ytickslab =c(20),
                # orientation = "reverse",
                xtickslab.rt = 45, ytickslab.rt = 45)
rcorr 



##----- RMSE
dim(rmse <- met_long[met_long$metric=="RMSE",])

plot_rmse <- ggline(rmse, x = "spat_cor", y = "estimate",
                    error.plot = "estimate",
                    facet.by = "method",
                    panel.labs.font.x = list(size=25),
                    color = "cover",
                    palette = "jco",
                    point.size=2,
                    #linetype = "method",
                    size=2)+
  theme_bw()+
  theme(strip.text = element_text(size=25))

rrmse <-  ggpar(plot_rmse, xlab="Spatial Variance", ylab="Root Mean Square Error (RMSE)",
                legend = "top", legend.title = "Survey Coverage (%)",size=25,
                font.legend=c(25),
                font.label = list(size = 20, face = "bold", color ="red"),
                font.x = c(25),
                font.y = c(25),
                font.main=c(20),
                font.xtickslab =c(20),
                font.ytickslab =c(20),
                # orientation = "reverse",
                xtickslab.rt = 45, ytickslab.rt = 45)
rrmse 

##----- MAE
dim(mae <- met_long[met_long$metric=="MAE",])

plot_mae <- ggline(mae, x = "spat_cor", y = "estimate",
                   error.plot = "estimate",
                   facet.by = "method",
                   # panel.labs= list(dat_type=c("diff", "same")),
                   panel.labs.font.x = list(size=25),
                   color = "cover",palette = "jco",
                   point.size=2,
                   #linetype = "pop_cover",
                   size=2)+
  theme_bw()+
  theme(strip.text = element_text(size=25))
rmae <-  ggpar(plot_mae, xlab="Spatial Variance", ylab="Mean Absolute Error (MAE)",
               legend = "none", 
               font.legend=c(25),
               font.label = list(size = 20, face = "bold", color ="red"),
               font.x = c(25),
               font.y = c(25),
               font.main=c(22),
               font.xtickslab =c(22),
               font.ytickslab =c(22),
               # orientation = "reverse",
               xtickslab.rt = 45, ytickslab.rt = 45)
rmae 



##----- BIAS
dim(bias <- met_long[met_long$metric=="BIAS",])

plot_bias <- ggline(bias, x = "spat_cor", y = "estimate",
                    error.plot = "estimate",
                    facet.by = "method",
                    # panel.labs= list(dat_type=c("diff", "same")),
                    panel.labs.font.x = list(size=25),
                    color = "cover",palette = "aaas",
                    point.size=2,
                    #linetype = "pop_cover",
                    size=2)+
  theme_bw()+
  theme(strip.text = element_text(size=25))
rbias <-  ggpar(plot_bias, xlab="Spatial Variance", ylab="Absolute Bias (Abias)",
                legend = "none", 
                font.legend=c(25),
                font.label = list(size = 20, face = "bold", color ="red"),
                font.x = c(25),
                font.y = c(25),
                font.main=c(22),
                font.xtickslab =c(22),
                font.ytickslab =c(22),
                # orientation = "reverse",
                xtickslab.rt = 45, ytickslab.rt = 45)
rbias 

ggarrange(rrmse,rmae, 
          ncol = 1, nrow = 2)
#------------------------------------------------------------
#### --- Perecentage reductions in RMSE and MAE
#----------------------------------------------------------
### Reduction in relative error measures
rrem <- as.data.frame(read.csv(paste0(out_path, "/error_reduction.csv")))
rr_dat <- rrem %>% drop_na(rrmae)%>% dplyr::select(c(cover,
                                                     spat_cor,
                                                     rrmae, rrrmse))


rr_dat <- rr_dat %>% filter(spat_cor !="0.75")
rr_dat$rrmae <- round(rr_dat$rrmae*100, 2)
rr_dat$rrrmse <- round(rr_dat$rrrmse*100, 2)
rr_dat$spat_cor <- factor(rr_dat$spat_cor)
rr_dat$cover <- factor(rr_dat$cover)
rr_dat$cover <- factor(rr_dat$cover,
                       labels=c("20%", "40%",
                                "60%", "80%"))


# rrmae
bar_mae <- rr_dat %>%
  ggplot( aes(x=cover, y=rrmae, fill=spat_cor)) +
  geom_bar(stat = "identity", 
           position = position_dodge()) +
  #scale_x_continuous(breaks = seq(20,90, by=10))+
  #scale_y_continuous(breaks=seq(500,7000,1000))+
  #coord_flip()+
  theme_bw()+
  theme(strip.text = element_text(size = 16),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))

rbar_mae <-  ggpar(bar_mae, xlab="Proportion of population observed(%)", 
                   ylab="Reduction in relative MAE(%)",
                   legend = "none", legend.title=element_text("Spatial Variance:"),
                   font.legend=c(25),
                   palette = c("jco"),
                   font.label = list(size = 25, face = "bold", color ="red"),
                   font.x = c(25),
                   font.y = c(25),
                   font.main=c(25),
                   font.xtickslab =c(22),
                   font.ytickslab =c(22),
                   xtickslab.rt = 45, ytickslab.rt = 45)
rbar_mae



# rrrmse
bar_rmse <- rr_dat %>%
  ggplot( aes(x=cover, y=rrrmse, fill=spat_cor)) +
  geom_bar(stat = "identity", 
           position = position_dodge()) +
  #scale_y_continuous(breaks=seq(500,7000,1000))+
  #coord_flip()+
  theme_bw()+
  theme(strip.text = element_text(size = 16),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))

rbar_rmse <-  ggpar(bar_rmse, xlab="Proportion of population observed(%)", 
                    ylab="Reduction in relative RMSE(%)",
                    legend = "top", legend.title=element_text("Spatial Variance:"),
                    font.legend=c(25),
                    palette = c("jco"),
                    font.label = list(size = 25, face = "bold", color ="red"),
                    font.x = c(25),
                    font.y = c(25),
                    font.main=c(25),
                    font.xtickslab =c(22),
                    font.ytickslab =c(22),
                    xtickslab.rt = 45, ytickslab.rt = 45)
rbar_rmse



ggarrange(rbar_rmse, rbar_mae, 
          ncol = 1, nrow = 2)


min(rr_dat$rrmae); max(rr_dat$rrmae)
min(rr_dat$rrrmse); max(rr_dat$rrrmse)


