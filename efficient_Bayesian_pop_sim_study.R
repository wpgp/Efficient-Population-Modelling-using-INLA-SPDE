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
model_metrics <- function(obs, pred, upper, lower)
{
  residual = pred - obs
  MAE = mean(abs(residual), na.rm=T) # Mean Absolute Error
  MSE = mean(residual^2, na.rm=T) # Mean Square Error 
  RMSE = sqrt(MSE) # Root Mean Square Error 
  BIAS = mean(residual, na.rm=T) # Bias
  
  output <- list(MAE  = MAE ,
                 RMSE = RMSE,
                 BIAS = abs(BIAS))
  return(output)
}

#model_metrics(obs, pred, upper, lower)



cover <- c(1, 0.8, 0.6, 0.4, 0.2, 0.05) #--percentage coverage of the observations
spat.autocor <- c(0.01, 0.1, 1)#spatial autocorrelation variance

#source_re <- c(0.01, 0.1, 1)#variance of data source random effects
source_var <- list(same = rep(0.01, 5), 
                  diff = c(0.01, 0.04, 0.08, 0.2, 0.7))#variance of data source random effects
dat_type <- c("same", "diff")
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
  r0 <- 0.3 #--range
  nu <- 1  #--smooth parameter
  
  for(b in 1:length(spat.autocor))
  { 
  
 # b=1
  result_pathb <- paste0(out_path,"/marginal_var_", spat.autocor[b])
  if (file.exists(result_pathb)){
    setwd(file.path(result_pathb))
  } else {
    dir.create(file.path(result_pathb))
    setwd(file.path(result_pathb))
  }
  
  sigma0 <- spat.autocor[b] #--marginal variance  C(0.01, 0.1, 1)
  kappa0 <- sqrt(8*nu)/r0 #--scale parameters
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
  names(dat)
  covs = cbind(1,runif(nn,1,10),rnorm(nn, 2, 5),rpois(nn,4),rnorm(nn, 4, 5), runif(nn))
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

  for(c in 1:length(dat_type))
  { 
  #c=2
  result_pathc <- paste0(result_pathb,"_source_var_", dat_type[c])
  if (file.exists(result_pathc)){
    setwd(file.path(result_pathc))
  } else {
    dir.create(file.path(result_pathc))
    setwd(file.path(result_pathc))
  }
  
  dtyp_re1 <- as.vector(t(mvrnorm(1, rep(0, n.dat_typ), 
                        diag(source_var[[c]]))))
  
  
  ##-----Extract data types effects
  data_type_re <- function(dat1, dat2, ID)
  {
    typ <- rep(1, nrow(dat1))
    typ_ID <- as.numeric(ID)
    uniq <- unique(typ_ID)
    for(k in  1:nrow(dat1))
    {
      
      for(j in uniq)
      {
        if(typ_ID[k]==uniq[j]) typ[k] = dat2[j]
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
  sigma_b <- 0.2
  epsb <- rnorm(nrow(coord), 0, sigma_b)#--iid random effect for building count
  #betaB <- c(3.51, 0.46, 0.85, 0.21, 0.58, 0.235) #---betas- fixed effects
  betaB <- c(2.21, 0.06, 0.15, 0.21, 0.18, 0.27) 
  bld <- lambdaB <- numeric(nn) #
  for (i in 1:nn)
  {
    lambdaB[i] <- exp(zpred[i,1]*betaB[1] + zpred[i,2]*betaB[2] + 
                        zpred[i,3]*betaB[3] + zpred[i,4]*betaB[4] + 
                        zpred[i,5]*betaB[5] + zpred[i,6]*betaB[6] + sp.raneff[i] + epsb[i])
    bld[i] <- rpois(1, lambdaB[i])
  }
  bld
  min(bld); max(bld)
  hist(bld); hist(log(bld))
  mean(bld); var(bld); 
  summary(bld)
  summary(dat$Bldng_C)
  
  # set any 0 buildings to 1
  bld[bld==0]=1
  
  
  ###---Simulate Population count
  sigma_p <- 0.85
  epsp <- rnorm(nrow(coord), 0, sigma_p) 
  betaP <- c(3.50, 0.41, 0.08, 0.04, 0.15, 0.22) #--betas - fixed effects
  
  pop <- lambdaP <- numeric(nn)
  for (i in 1:nn)
  {
    lambdaP[i] <- exp(zpred[i,1]*betaP[1] + zpred[i,2]*betaP[2] + 
                        zpred[i,3]*betaP[3] + zpred[i,4]*betaP[4] + 
                        zpred[i,5]*betaP[5] + zpred[i,6]*betaP[6] + dty_re1[i] + sp.raneff[i] + epsp[i])
    pop[i] <- rpois(1, lambdaP[i])
  }
  sum(pop)
  hist(pop); hist(log(pop))
  mean(pop); var(pop)
  summary(pop)
  summary(dat$I_LHHSI)
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
    set.seed(566200)
    #j = 1
    print(c(spat.autocor[b], dat_type[c],paste0(cover[j]*100, "%")))
    
  
    result_path2 <- paste0(result_pathc,"_at_", cover[j]*100,"%_coverage")
    if (file.exists(result_path2)){
      setwd(file.path(result_path2))
    } else {
      dir.create(file.path(result_path2))
      setwd(file.path(result_path2))
    } 
    
    
    ##---Select the subset of the data as observed at random.
    pp8 <- cover[j]
    opts8 <- sample(nrow(coord), pp8*nrow(coord))
    opts8 <- sort(opts8)
    dim(dat8 <-  data.frame(datam[opts8,])) #---Locations with observations
    
    ##--------------
    coord8 <- cbind(dat8$lon, dat8$lat)
    dim(coord8)
    
    
    ##---Build mesh at grid cell level at the coverage level
    bnd8 <- inla.nonconvex.hull(coord8, -0.035, -0.04, resolution = c(100, 100))
    mesh8 <- inla.mesh.2d(boundary = bnd8, max.edge=c(0.2,1), 
                         offset = c(0.2, 0.7),
                         cutoff = 0.2)
    par(mfrow=c(1,1))
    plot(mesh8)
    mesh8$n # 1133 mesh nodes
    
    ###---Build projector matrix A 
    A8<-inla.spde.make.A(mesh=mesh8,loc=as.matrix(coord8));dim(A8)
    
    ##---Create the SPDE
    spde8 <- inla.spde2.matern(mesh8, alpha=2)
    
    ##----specify the observation indices for estimation 
    iset8<- inla.spde.make.index(name = "spatial.field", spde8$n.spde)
    
    
    dat.est8 <- dat8
    dat.est8$eps <- 1:nrow(dat.est8)
    head(dat8)
    names(dat.est8)
    
    sum(datam$pop); sum(dat.est8$pop)
    names(dat.est8)
    head(dat.est8)
    table(dat.est8$source1 <- as.factor(dat.est8$source))
    
    dat.est8[,c("x2", "x3", "x4", "x5", "x6")] <- apply(dat.est8[,c("x2", "x3", "x4", "x5", "x6")],2, stdize)
    covars.est8 <- dat.est8[,c("x1","x2", "x3", "x4", "x5", "x6", "source1", "eps")]; dim(covars.est8) ##---Population density
    head(covars.est8)
    dat.est8$dens <- dat.est8$pop/dat.est8$bld
    dat.est8$dens[is.infinite(dat.est8$dens)] = NA
    head(dat.est8)
    ##---Build data stack 
    stk.est8 <- inla.stack(data=list(y=dat.est8$dens), #the response
                           
                           A=list(A8,1),  #the A matrix; the 1 is included to make the list(covariates)
                           
                           effects=list(c(list(Intercept=1), #the Intercept
                                          iset8),  #the spatial index
                                        #the covariates
                                        list(covars.est8)
                           ), 
                           #this is a quick name so you can call upon easily
                           tag='est')
    
  
    ######---------INLA Models
    # Model 1
    form8a <- y ~ -1 + Intercept + x2 + x3 + x4 + x5 + x6 + f(eps, model='iid')
    mod8a <-inla(form8a, #the formula
                    data=inla.stack.data(stk.est8,spde=spde8),  #the data stack
                    family= 'gamma',   #which family the data comes from
                    control.predictor=list(A=inla.stack.A(stk.est8),compute=TRUE),  #compute gives you the marginals of the linear predictor
                    control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                    verbose = FALSE) 
    #summary(mod8a)
    #index <-inla.stack.index(stk.est8, "est")$data #---extract the data location indices
    #pred8a <- exp(mod8a$summary.linear.predictor[index ,"mean"]) #--predicted mean count
    #sum(fit8a <- pred8a*dat.est8$bld)
   # sum(dat.est8$pop)
    
    # Model 2
    form8b <- y ~ -1 + Intercept + x2 + x3 + x4 + x5 + x6 + 
      f(spatial.field, model=spde8) + f(eps, model='iid')
    mod8b <-inla(form8b, #the formula
                     data=inla.stack.data(stk.est8,spde=spde8),  #the data stack
                     family= 'gamma',   #which family the data comes from
                     control.predictor=list(A=inla.stack.A(stk.est8),compute=TRUE),  #compute gives you the marginals of the linear predictor
                     control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                     verbose = FALSE) 
    #summary(mod8b)
    
    
    # Model 3
    form8c <- y ~ -1 + Intercept + x2 + x3 + x4 + x5 + x6 + 
      f(source1, model='iid') + f(eps, model='iid')
    mod8c <-inla(form8c, #the formula
                     data=inla.stack.data(stk.est8,spde=spde8),  #the data stack
                     family= 'gamma',   #which family the data comes from
                     control.predictor=list(A=inla.stack.A(stk.est8),compute=TRUE),  #compute gives you the marginals of the linear predictor
                     control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                     verbose = FALSE) 
    #summary(mod8c)
    
    
    
    # Model 4
    form8d <- y ~ -1 + Intercept + x2 + x3 + x4 + x5 + x6 + 
      f(eps, model='iid') +
      f(source1, model='iid') +
    f(spatial.field, model=spde8)
    
    mod8d<-inla(form8d, #the formula
                     data=inla.stack.data(stk.est8,spde=spde8),  #the data stack
                     family= 'gamma',   #which family the data comes from
                     control.predictor=list(A=inla.stack.A(stk.est8),compute=TRUE),  #compute gives you the marginals of the linear predictor
                     control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                     verbose = FALSE) 
    #summary(mod8d)
    
    
   # Model fit checks -----------------------------------------------
    
    DIC <- t(c(mod1 = mod8a$dic$dic,
               mod2 = mod8b$dic$dic,
               mod3 = mod8c$dic$dic,
               mod4 = mod8d$dic$dic))
    
    WAIC <- t(c(mod1 = mod8a$waic$waic,
               mod2 = mod8b$waic$waic,
               mod3 = mod8c$waic$waic,
               mod4 = mod8d$waic$waic))
    
    CPO <- t(c(mod1 = -sum(log(mod8a$cpo$cpo)),
               mod2 = -sum(log(mod8b$cpo$cpo)),
               mod3 = -sum(log(mod8c$cpo$cpo)),
               mod4 = -sum(log(mod8d$cpo$cpo))))
    
    (modfit <- data.frame(DIC=t(DIC), WAIC=t(WAIC), CPO=t(CPO)))
    
    write.csv(modfit, file=paste0(result_path2,"/modfit.csv"))
    
    ####----View the spatial fields of the best fit model
    #looking at the spatial field and what it looks like
    gproj8 <- inla.mesh.projector(mesh8,  dims = c(300, 300))
    
    col <- viridis(100)
    
    g.mean8 <- inla.mesh.project(gproj8, mod8d$summary.random$spatial.field$mean)
    g.sd8 <- inla.mesh.project(gproj8, mod8d$summary.random$spatial.field$sd)
    
    
    grid.arrange(levelplot(g.mean8, scales=list(draw=F), xlab='', ylab='', cex.lab=2, main='Mean',col.regions = col),
                 levelplot(g.sd8, scal=list(draw=F), xla='', yla='', main='SD' ,col.regions = col), nrow=1)
    
    index <-inla.stack.index(stk.est8, "est")$data #---extract the data location indices
    pred8d <- exp(mod8d$summary.linear.predictor[index ,"mean"]) #--predicted mean count
     sum(fit8d <- pred8d*dat.est8$bld)
     sum(dat.est8$pop)
     plot(dat.est8$pop, fit8d)
     cor(dat.est8$pop, fit8d)
  
#----------------------------------------------------------------------
    #  Run posterior simulation and predictions simulataneously
#-----------------------------------------------------------------
    
    # Prediction Projection matrix 
    Aprd <- inla.spde.make.A(mesh = mesh8, loc = as.matrix(coord))
    dim(Aprd)#--prediction projectionnmatrix
    
    dat.pred <- datam
  
    ##-----Extract data types effects
    #dat.pred$source2 <- as.numeric(dat.pred$source)
    dat.pred$source_re <- data_type_re(dat.pred, 
                                       mod8d$summary.random$source1$mean,
                                       dat.pred$source)
    
    mean(dat.pred$source_re); var(dat.pred$source_re)
    dat.pred[,c("x2", "x3", "x4", "x5", "x6")] <- apply(dat.pred[,c("x2", "x3", "x4", "x5", "x6")],2, stdize)
    
    head(datam)
    head(dat.pred)
    ####-------Run posterior and grid cell prediction (GRID2GRID)
    post_sim <- function(model, dat, Apred, run)
    {
      fixedeff <- dens_hat <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
      inla.seed = 166
      #inla.seed = as.integer(runif(1)*.Machine$integer.max)
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
          
          dat$source_re + #Data source random effect
          #rnorm(nrow(dat), 0, 1/model$summary.hyperpar$mean[3]) + #Data source random effect
          rnorm(nrow(dat), 0, 1/model$summary.hyperpar$mean[2]) + #IID
          #model$summary.random$eps$mean[opts8]+ #IID
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
    run=100
    system.time(str(sim.dens8 <- post_sim(mod8d,dat.pred,Aprd, run)))

#summary(mod8d)
    sum(sim.dens8$est_data$mean_pop_hat)
    sum(sim.dens8$est_data$pop)

    write.csv(sim.dens8$est_data, file=paste0(result_path2,"/grid_data.csv"))
    
##---Calculate model fit and model cross-validation metrics  
   met <- model_metrics(sim.dens8$est_data$pop, 
                  sim.dens8$est_data$mean_pop_hat)
   (metrics <- round(unlist(met),3))
   plot(sim.dens8$est_data$pop, 
        sim.dens8$est_data$mean_pop_hat)
   write.csv(metrics, file=paste0(result_path2,"/fit_metrics.csv"))
  }
  
}
}


#save.image(paste0(results_path, "/posterior_samples/CMR_MODEL_main_workspace_final.Rdata"))
