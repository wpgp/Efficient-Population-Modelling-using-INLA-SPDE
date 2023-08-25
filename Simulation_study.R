####--TITLE: R SCRIPTS FOR MODELLED POPULATION ESTIMATES BASED ON R-INLA:SIMULATION STUDY----#####
####--METHODS: GEOSTATISTICAL BAYESIAN HIERARCHICAL REGRESSION MODEL---------------------#####
####--AUTHOR: DR CHIBUZOR CHRISTOPHER NNANATU ------------------------#####
####--INSTITUTION: WORLDPOP, UNIVERSITY OF SOUTHAMPTON----------------------------------##### 
####--DATE: DECEMBER 2022---------------------------------------------------------------####
####@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@####
###--LOAD KEY PACKAGES AND LIBRARIES
library(INLA); library(raster); library(maptools)
library(gtools); library(sp); library(spdep)
library(fields); library(mvtnorm); library(gtools)
library(geoR);library(actuar);library(viridisLite)
require(grid);require(gridExtra);require(lattice);require(tidyverse)

##---Specify files and outputs paths
path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3/Working/CMR/Chris_N/codes/INLA_models/final_models/updated/sim_study2"
file_path <- paste0(path, "/files")
out_path <- paste0(path, "/outputs")
setwd(path)
#getwd()

AdmUnit <- c(8, 10, 20, 30)#square root of the number of areal units
cover <- c(0.8, 0.6, 0.4, 0.2) #--percentage coverage of the observations
pg <- 120# = 14400 prediction grids

####---Model fit metrics function
model_metrics <- function(obs, pred, upper, lower)
{
  residual = pred - obs
  INACCURACY = mean(abs(residual), na.rm=T)#MAE
  MSE = mean(residual^2, na.rm=T)
  RMSE = sqrt(MSE)
  BIAS = mean(residual, na.rm=T)
  corr = cor(obs[!is.na(obs)],pred[!is.na(obs)])
  
  output <- list(MAE  = INACCURACY ,
                 RMSE = RMSE,
                 BIAS = abs(BIAS),
                 corr = corr)
  return(output)
}
#model_metrics(obs, pred, upper, lower)

###---loop over each number of areal unit (adminunit) level
for(i in AdmUnit)
{
  
i=8
  print(i)
  result_path <- paste0(out_path,"/estimates_for_", i^2,"_area_units")
  if (file.exists(result_path)){
   setwd(file.path(result_path))
  } else {
    dir.create(file.path(result_path))
    setwd(file.path(result_path))
  }
  
  #---Select AdmUnit Size
  Area <- as(raster(nrow=i, ncol=i, xmn=0, xmx=1, ymn=0, ymx=1), "SpatialPolygons")
  proj4string(Area) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  #
  tiff(paste0(result_path,"/shape_plot.tiff"), width = 4, height = 4, units = 'in', res = 300)
  plot(Area, lwd=2)
  dev.off()
  
  #---Create raster
  pgd = raster(nrow=pg, ncol=pg, xmn=0, xmx=1, ymn=0, ymx=1) #--prediction grid
  #res(pgd) = c(0.00083, 0.00083) #--approx 100m by 100m resolution
  pgd1 = pgd
  
  
  ######
  coord.1 = xyFromCell(pgd1, 1:ncell(pgd1), spatial=FALSE)
  coord = coord.1
  
  
  ####----Generate points within the areas
  spol = Area
  c.sp=SpatialPoints(coord)
  proj4string(c.sp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  sp.1 <- rep(NA, nrow(coord))
  for(i in 1:length(spol)){
    sp.1[as.vector(which(!is.na(over(c.sp, spol[i]))))] <- i
  }
  spred <- sp.1
  #table(spred)
 
  
  #plot(coord)
  #plot(Area, add=T, lwd=2)
  
  
  dim(dat.all <- as.data.frame(coord))
  table(dat.all$regionID <- as.factor(spred))
  
  n.dat_typ = 5
  prop.dat_typ <- c(0.19, 0.16, 0.16, 0.17, 0.32)
  n.sample <- nrow(dat.all)
  table(dat_typ1 <- sample(1:n.dat_typ,n.sample, prob=prop.dat_typ, rep=T))
  ###---Add to dataset
  table(dat.all$dat_typ_ID <- as.factor(dat_typ1))
  table(dat.all$clusterID <- as.factor(1:nrow(coord)))
  
  
####------Build the mesh
  dim(coord)
  bnd <- inla.nonconvex.hull(coord, -0.035, -0.04, resolution = c(100, 100))
  meshf <- inla.mesh.2d(boundary = bnd, 
                          offset=c(0.1, 0.15), max.edge=c(0.06472101,0.3)) 
  #plot(meshf)
  #points(coord, col="light blue", pch="*", cex=1.5)
  #plot(Area, add=T, lwd=1.5)
  #plot(meshf, add=T)
  #meshf$n
  
  
  #--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##---Specify SPDE parameters
  r0 <- 0.3 #--range
  nu <- 1  #--smooth parameter
  sigma0 <- 1 #--marginal variance
  kappa0 <- sqrt(8*nu)/r0 #--scale parameters
  tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0) #--precision
  
  #---SPDE
  spde <- inla.spde2.matern(meshf, 
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
  A <- inla.spde.make.A(mesh=meshf, loc=coord);dim(A)
  S.pred <- as.vector(A %*% sam)
  #hist(S.pred)
  
  
  #####---Spatial random effects
  tau2 =1
  dim(Q); class(Q)
  Q <- as.matrix(Q)
  sig <- solve(Q)
  mm <- meshf$n
  phi.a = rmvnorm(1,rep(0,mm), sig)
  phi.a = as.vector(phi.a)  #
  
  #Construct phi for each grid cell
  phi.pred <- numeric()
  for(i in 1:length(spred)) phi.pred[i] <- phi.a[spred[i]]
  
  
  ##----specify the observation indices for estimation 
  iset <- inla.spde.make.index(name = "spatial.field", spde$n.spde)
  #length(iset)
  
  
  ##---Nugget effect/iid term for population count
  sigma_e <- 0.05
  ee.pred <- rnorm(nrow(coord), 0, sigma_e) 
  
  
  #-- generate the Covariates 
  nn <- nrow(dat.all)
  covs = cbind(1,runif(nn),rnorm(nn),rpois(nn,2),rnorm(nn), runif(nn))
  
  apply(covs, 2, mean) #---check averages
  dim(dat.all)
  
  dim(ddat <- cbind(dat.all, covs))
  ddat <- data.frame(ddat)
  #head(ddat)
  #names(ddat)
  
  ###############
  #Parameters 
  dim(zpred <- covs)
  
  #dim(zpred)
  #For prediction grid
  #nrow(ddat)
  #mean(bld) = 369
  
  require(MASS)
  #--generate data type random effects 
  dtyp_re1 <- t(mvrnorm(1, rep(0, n.dat_typ), diag(c(0.06, 0.01, 0.08, 0.05, 0.02))))
 
  ##-----Extract data types effects
  data_type_re <- function(dat1, dat2, ID)
  {
  typ <- rep(1, nrow(dat1))
  typ_ID <- as.numeric(ID)
  uniq <- unique(typ_ID)
  for(i in  1:nrow(dat1))
  {
    
    for(j in uniq)
    {
      if(typ_ID[i]==uniq[j]) typ[i] = dat2[j]
    }
    
  }
  typ
  }
  dat.all$dtyp_re1 <- data_type_re(dat.all, dtyp_re1,dat.all$dat_typ_ID)#--data type re with varying effects
 
  ##--Add to dataset
  
  dty_re1 <- dat.all$dtyp_re1

  beta <- c(1.5, 0.76, 0.2, 0.2, -0.028, 0.5)
  rr=0.01
  resp <- mu <- numeric(nn) #
  for (i in 1:nn)
  {
    mu[i] <- exp(zpred[i,1]*beta[1] + zpred[i,2]*beta[2] + zpred[i,3]*beta[3] + zpred[i,4]*beta[4] + 
        zpred[i,5]*beta[5] + zpred[i,6]*beta[6] + S.pred[i] + phi.pred[i]+ dty_re1[i] + ee.pred[i])
    resp[i] <- rgamma(1, mu[i]^2/rr,mu[i]/rr)
  }
  
  
mean(resp)
  #####----------ALTERNATIVELY---------########
  #---Simulate Pop Count and Building Count separately and calculate pop density
  #---see the alternative codes below
  ###---Simulate building count
  sigma_b <- 0.05
  epsb <- rnorm(nrow(coord), 0, sigma_b)#--iid random effect for building count
  betaB <- c(3.1, 0.16, 0.25, 0.21, 0.18, 0.0935) #---betas- fixed effects
  
  bld <- lambdaB <- numeric(nn) #
  for (i in 1:nn)
  {
    lambdaB[i] <- exp(zpred[i,1]*betaB[1] + zpred[i,2]*betaB[2] + 
                        zpred[i,3]*betaB[3] + zpred[i,4]*betaB[4] + 
                        zpred[i,5]*betaB[5] + zpred[i,6]*betaB[6] + S.pred[i] + phi.pred[i] + dty_re1[i] + epsb[i])
    bld[i] <- rpois(1, lambdaB[i])
  }
  bld
  min(bld)
  hist(bld); hist(log(bld))
  mean(bld); var(bld); 
  
  
  ###---Simulate Population count
  sigma_p <- 0.05
  epsp <- rnorm(nrow(coord), 0, sigma_p) 
  betaP <- c(6.85, 0.01, 0.12, 0.02, 0.003, 0.012) #--betas - fixed effects
  
  pop <- lambdaP <- numeric(nn)
  for (i in 1:nn)
  {
    lambdaP[i] <- exp(zpred[i,1]*betaP[1] + zpred[i,2]*betaP[2] + 
                        zpred[i,3]*betaP[3] + zpred[i,4]*betaP[4] + 
                        zpred[i,5]*betaP[5] + zpred[i,6]*betaP[6] + S.pred[i] + phi.pred[i] + epsp[i])
    pop[i] <- rpois(1, lambdaP[i])
  }
  pop
  hist(pop); hist(log(pop))
  mean(pop); var(pop)
  
  ###--------Add to dataset
  ddat$bld <- bld
  ddat$pop <- pop
  ddat$dens <- ddat$pop/ddat$bld #----population density
  hist(ddat$dens); hist(log(ddat$dens))
  
  
  
  ###---rename x and y variables and select
  names(ddat)
  datam <- ddat %>% 
    mutate(Lon = x, Lat=y) %>%
    dplyr::select(-x, -y)
  
  datam <- as.data.frame(datam)
  
  class(datam)
  names(datam)
  dim(ddat)

#--rename as necessary
  datam$resp <- datam$dens
  hist(datam$resp)
  
  
  
  names(datam)
  #----Scale covariates for stack
  stdize <- function(x)
  {
    stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
    return(stdz)
  }
  
  ##--prepare covariates fro INLA modelling
  names(dat.fit)
  dim(dat.fit <- datam)
  class(dat.fit$regionID)
  class(dat.fit$dat_typ_ID)
  dat.fit[, c("X2", "X3", "X4", "X5", "X6")] <- apply(dat.fit[, c("X2", "X3", "X4", "X5", "X6")], 2, stdize)
  
  covars.fit <- dat.fit[,c("X2", "X3", "X4", "X5", "X6","dat_typ_ID", "regionID", "clusterID")]; dim(covars.fit) ##---Population density
  
  names(dat.fit)
  
  #---Build the stack
  stk.fit <- inla.stack(data=list(y=dat.fit$dens), #the response
                        
                        A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                        
                        effects=list(c(list(Intercept=1), #the Intercept
                                       iset),  #the spatial index
                                     #the covariates
                                     list(covars.fit)
                        ), 
                        #this is a quick name so you can call upon easily
                        tag='est')
  
  
  
  ######--------- fit the initial model at grid cell level with 100% observations  
  form.fit <- y ~ -1 + Intercept + X2 + X3 + X4 + X5 + X6 + f(spatial.field, model=spde) + 
   f(dat_typ_ID, model='iid') + f(regionID, model='iid') 
  mod.fit <-inla(form.fit, #the formula
                 data=inla.stack.data(stk.fit,spde=spde),  #the data stack
                 family= 'gamma',   #which family the data comes from
                 control.predictor=list(A=inla.stack.A(stk.fit),compute=TRUE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 verbose = FALSE) 
  summary(mod.fit)
  
  
  #saveRDS(mod.fit, file=paste0(result_path,"/full_grid_model.rds"))
  ind.fit <-inla.stack.index(stk.fit, "est")$data
  fit1 <- exp(mod.fit$summary.linear.predictor[ind.fit,"mean"])
  sum(pred1 <- round(fit1*dat.fit$bld))

  
  dat.fit$pdens1 <-fit1 
  dat.fit$pred1 <- pred1
  
  #########################################################################################################
  ####----View the spatial fields of the best fit model
  #looking at the spatial field and what it looks like
  gproj <- inla.mesh.projector(meshf,  dims = c(300, 300))
  
  col <- viridis(100)
  
  g.mean <- inla.mesh.project(gproj, mod.fit$summary.random$spatial.field$mean)
  g.sd <- inla.mesh.project(gproj, mod.fit$summary.random$spatial.field$sd)
  
  
  
  grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', cex.lab=2, main='Mean',col.regions = col),
               levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='SD' ,col.regions = col), nrow=1)
  
  #####===============
  ##-----Extract regional effects
  dat.fit2 <- dat.fit
  dat.fit2$regionID <- dat.all$regionID
  dat.fit2$reg <- data_type_re(dat.fit, mod.fit$summary.random$regionID$mean,dat.fit$regionID)

  
  ##-----Extract data types effects
  dim(dat.fit2)
  names(dat.fit2)
  dat.fit2$dtyp <- data_type_re(dat.fit2, mod.fit$summary.random$dat_typ_ID$mean,dat.fit$dat_typ_ID)
  
  
  ####-------Posterior simulation and grid cell prediction (GRID2GRID)
  post_sim <- function(model, dat, Aprediction, run)
  {
    fixedeff <- dens_hat <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
    inla.seed = as.integer(runif(1)*.Machine$integer.max)
    #inla.seed = 
    
    set.seed(inla.seed)
    print(inla.seed)
    m1.samp <- inla.posterior.sample(run, model, seed = inla.seed, selection=list(X2=1, X3=1, X4=1,
                                                                                  X5=1, X6=1),num.threads="1:1")
    
    sfield_nodes_mean <- model$summary.random$spatial.field['mean']
    field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
    for(i in 1:run)
    {
      fixedeff[,i] <- 
        model$summary.fixed['Intercept', 'mean'] +
        m1.samp[[i]]$latent[1,] * dat[,'X2'] +
        m1.samp[[i]]$latent[2,] * dat[,'X3'] +
        m1.samp[[i]]$latent[3,] * dat[,'X4'] +
        m1.samp[[i]]$latent[4,] * dat[,'X5'] +
        m1.samp[[i]]$latent[5,] * dat[,'X6'] + 
        dat$reg + #--Areal level random effects random effect
        dat$dtyp + #--data source random effect

        field_mean[,1] #---spatial random effect
      
      ##
      dens_hat[,i] <- exp(fixedeff[,i]) #--density
      pop_hat[,i] <- dens_hat[,i]*dat$bld #--population count
      
    }
    
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
  run=200#--Note that run=100 or above would suffice but the larger the better
  system.time(str(sim.dens <- post_sim(mod.fit,dat.fit2,A, run)))
  write.csv(sim.dens$est_data, file=paste0(result_path,"/grid_data_gg.csv"))
  
 
  ###----Model cross validation and fit metrics
  cpo <- mod.fit$cpo$cpo #--for leave one out cross validation
  lcpo <- -sum(log(cpo), na.rm=T)
  dic <-  mod.fit$dic$dic #--for model fit 
  mets <- model_metrics(dat.fit$pop, sim.dens$est_data$mean_pop_hat, sim.dens$est_data$upper_pop_hat, 
                        sim.dens$est_data$lower_pop_hat)
  metrics <- list(lcpo=lcpo, dic=dic, mets=mets)
  capture.output(metrics, file=paste0(result_path, "/metrics_gg.txt"))
  
  
  
  ###------
  dat.fit$pdens_gg <- sim.dens$est_data$mean_dens_hat
  dat.fit$ppop_gg <- sim.dens$est_data$mean_pop_hat
  mean.pred <- sim.dens$est_data$mean_pop_hat
  
  
  
  ###---------
  dat.fit2$pred <- mean.pred 
  
  
  
  ####---Join the posterior sample to the prediction data
  data.sim <- data.frame(cbind(dat.fit2[,c("regionID")], sim.dens$pop_hat))
  names(data.sim)[1] <- "regionID"
  
  
  ##----Calculate and save National total with uncertainties
  nat_total <- function(dat, run)
  {
    p_hat <- dat[,2:(run+1)]
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
  
  sum(pred1 <- round(fit1*dat.fit$bld))
  write.csv(national, file=paste0(result_path,"/national_estimates_gg.csv"))
  # write.csv(national, file=paste0(results_path, "/estimates/National_estimates_final.csv"))
  
  
  ##---Areal level estimates with uncertanties
  admin_est <- function(datr, run)
  {
    uniR <-unique(datr$regionID)
    outR <- matrix(0, nrow=length(uniR), ncol=5)
    for(j in uniR)
    {
      reg <- datr[datr$regionID==j,]
      rtots <- apply(reg[,2:(1+run)], 2, sum, na.rm=T)
      #rtot_mean  <- sum(reg$pop_hat, na.rm=T)
      rtot_mean  <- mean(rtots, na.rm=T)
      rtot_sd <- sd(rtots, na.rm=T)
      rtot_lower <- quantile(rtots, probs=c(0.025))
      rtot_median <- quantile(rtots, probs=c(0.5))
      rtot_upper <- quantile(rtots, probs=c(0.975))
      rtot_uncert <- (rtot_upper - rtot_lower)/rtot_mean
      
      restimates <- round(c(rtot_mean, rtot_lower, rtot_median,rtot_upper, rtot_uncert),4)
      outR[j,] <- restimates
    }
    outR <- data.frame(outR)
    return(reg_est <- data.frame(ID = uniR,
                                 total = outR[,1],
                                 lower = outR[,2],
                                 median = outR[,3],
                                 upper = outR[,4],
                                 uncertainty = outR[,5]))
  }
  (admin.est <- admin_est(data.sim, run))
  sum(admin.est$total)
  write.csv(admin.est, file=paste0(result_path,"/area_estimates_gg.csv"))
  
  
##----Aggregate observation to areal units
  ##---Group by region and data type
  dat.agg1 <- data.frame(datam1 <-  datam %>%
                          group_by(regionID, dat_typ_ID) %>%
                          summarise(bld2 = sum(bld), pop2 = sum(pop), X2 = mean(X2),
                                    X3=mean(X3), X4= mean(X4), X5=mean(X5), X6 = mean(X6))) 
                        
  head(dat.agg1)
  
  
  ###---Group by area (or region)
  dat.agg2 <- data.frame(datamm <-  dat.agg1  %>%
                           group_by(regionID) %>%
                           summarise(bld3 = sum(bld2), pop3= sum(pop2), X2 = mean(X2),
                                     X3=mean(X3), X4= mean(X4), X5=mean(X5), X6 = mean(X6))) 
  
  head(dat.agg2)
  dat.agg2$Lon <- coordinates(Area)[,1]
  dat.agg2$Lat <- coordinates(Area)[,2]
  ####

  
  ##---Add Lon-Lat to dat.agg1
  dat.agg1$Lon <- rep(1, nrow(dat.agg1))
  dat.agg1$Lat <- rep(1, nrow(dat.agg1))
  
  for(i in 1:nrow(dat.agg2))
  {
    dat.agg1$Lon[dat.agg1$regionID == i] = dat.agg2$Lon[dat.agg2$regionID == i]
    dat.agg1$Lat[dat.agg1$regionID == i] = dat.agg2$Lat[dat.agg2$regionID == i]
  }

  
  ###--Check-------------
  par(mfrow=c(1,2))
  plot(dat.agg1$Lon, dat.agg1$Lat)
  plot(dat.agg2$Lon, dat.agg2$Lat)
  dim(dat.agg1); dim(dat.agg2)
  dat.agg <- dat.agg1
  names(dat.agg)
 
  ##--plot(dat.agg$Lon, dat.agg$Lat)
  names(dat.agg)
  write.csv(dat.agg, file=paste0(result_path,"/area_data.csv"))
  
  ##########
  #par(mfrow=c(1,1))
  coordb <- cbind(dat.agg$Lon, dat.agg$Lat)
  #plot(coordb)
  #plot(Area, add=T)
  
  
  
  ##---Build mesh for areal-level unit
  #meshb: fine triangulated mesh
  bndb <- inla.nonconvex.hull(coordb, -0.08, -0.7, resolution = c(100, 100))
  meshb<- inla.mesh.2d(boundary = bndb, 
                       offset=c(0.1, 0.15), max.edge=c(0.07,0.3)) 
  
  #plot(meshb)
  #plot(Area, add=T, col="light green", lwd=1.5)
  #plot(meshb, add=T)
  #points(coordb, pch="*", col="red", cex=1)
  #points(coord[-opts,], pch="*", col="blue", cex=1.5)
  meshb$n
  
  
  #--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ###---Build projector matrix A for areal level data
  Ab<-inla.spde.make.A(mesh=meshb,loc=as.matrix(coordb));dim(Ab)
  
  ##---Create the SPDE for areal level data
  spdeb <- inla.spde2.matern(meshb, alpha=2)
  
  ##----specify the observation indices for estimation 
  isetb<- inla.spde.make.index(name = "spatial.field", spdeb$n.spde)
  
  
  
  ###----Define the density variable 
  dat.agg$dens <- dat.agg$pop2/dat.agg$bld2
  
  
  ##---Scale the areal-level covariates as using the same scaling at the grid cell level

  stdize <- function(x)
  {
    stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
    return(stdz)
  }
  
  dat.estb <- dat.agg
  names(dat.agg)
  dat.estb[, c("X2", "X3", "X4", "X5", "X6")] <- apply(dat.estb[, c("X2", "X3", "X4", "X5", "X6")], 2, stdize)
  
  covars.estb <- dat.estb[,c("X2", "X3", "X4", "X5", "X6", "regionID", "dat_typ_ID")]; dim(covars.estb) ##---Population density
  
  
  
  ##---Build the data stack for the areal unit 
  stk.estb <- inla.stack(data=list(y=dat.estb$dens), #the response
                         
                         A=list(Ab,1),  #the A matrix; the 1 is included to make the list(covariates)
                         
                         effects=list(c(list(Intercept=1), #the Intercept
                                        isetb),  #the spatial index
                                      #the covariates
                                      list(covars.estb)
                         ), 
                         #this is a quick name so you can call upon easily
                         tag='est')
  

  ##---------areal-level model 
  form.estb <- y ~ -1 + Intercept + X2 + X3 + X4 + X5 + X6 + f(spatial.field, model=spdeb) +
    f(regionID, model='iid')+ f(dat_typ_ID, model='iid')
  mod.estb <-inla(form.estb, #the formula
                  data=inla.stack.data(stk.estb,spde=spdeb),  #the data stack
                  family= 'gamma',   #which family the data comes from
                  control.predictor=list(A=inla.stack.A(stk.estb),compute=TRUE),  #compute gives you the marginals of the linear predictor
                  control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                  verbose = FALSE) 
  summary(mod.estb)
  
  #saveRDS(mod.estb, file=paste0(result_path,"/full_area_model.rds"))
  
  indb <-inla.stack.index(stk.estb, "est")$data
  fitb <- exp(mod.estb$summary.linear.predictor[indb,"mean"])
  sum(predb <- round(fitb*dat.estb$bld2))
      
  #(obs_predb <- data.frame(pop=dat.estb$pop2,predb))
  #apply(obs_predb, 2, sum, na.rm=T)
  #plot(dat.estb$pop2,predb)
  #cor(dat.estb$pop2,predb)
  #abline(a=0, b=1)
  
###---add predictions to the areal dataset 
  dat.agg$pred1 <- predb
  dat.agg$dens1 <- fitb
  dat.agg$ppop_ga <- admin.est$total
  
  #plot(dat.agg$pop2, dat.agg$ppop_ga)
  #abline(a=0, b=1)
  #cor(dat.agg$pop2,dat.agg$ppop_ga)
  ####-----Make prediction
  
###----Prediction Projection matrix 
  #dim(dat.pred <- datam[-opts,])
  dim(dat.predb <- dat.fit)
  dim(coord.predb <- cbind(dat.predb$Lon, dat.predb$Lat))
  Apredictionb<- inla.spde.make.A(mesh = meshb, loc = coord.predb); dim(Apredictionb)
  dat.fit[, c("X2", "X3", "X4", "X5", "X6")]
  covars.estb <- dat.predb[,c("X1", "X2", "X3", "X4", "X5", "X6", "regionID", "dat_typ_ID")]; dim(covars.estb) ##---Population density
  
  
##----Extract areal level random effects----
  dat.predb$reg <- data_type_re(dat.predb, mod.estb$summary.random$regionID$mean,dat.predb$regionID)
  
##-----Extract data source random effects-------
  dat.predb$dat_typ <- data_type_re(dat.predb, mod.estb$summary.random$dat_typ_ID$mean,dat.predb$dat_typ_ID)
  
####-------Run posterior simulation and grid cell prediction (AREA2GRID)
  post_simb <- function(model, dat, Aprediction, run)
  {
    fixedeff <- dens_hat <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
    inla.seed = as.integer(runif(1)*.Machine$integer.max)
    #inla.seed = 
    
    set.seed(inla.seed)
    print(inla.seed)
    m1.samp <- inla.posterior.sample(run, model, seed = inla.seed, selection=list(X2=1, X3=1, X4=1,
                                                                                  X5=1, X6=1),num.threads="1:1")
    
    sfield_nodes_mean <- model$summary.random$spatial.field['mean']
    field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
    for(i in 1:run)
    {
      fixedeff[,i] <- 
        model$summary.fixed['Intercept', 'mean'] +
        m1.samp[[i]]$latent[1,] * dat[,'X2'] +
        m1.samp[[i]]$latent[2,] * dat[,'X3'] +
        m1.samp[[i]]$latent[3,] * dat[,'X4'] +
        m1.samp[[i]]$latent[4,] * dat[,'X5'] +
        m1.samp[[i]]$latent[5,] * dat[,'X6'] + 
        dat$reg +
        dat$dat_typ
        field_mean[,1]
      
      ##
      dens_hat[,i] <- exp(fixedeff[,i])
      pop_hat[,i] <- dens_hat[,i]*dat$bld
      
    }
    
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
  run=100
  system.time(str(sim.densb<- post_simb(mod.estb, dat.predb,Apredictionb, run)))
  write.csv(sim.densb$est_data, file=paste0(result_path,"/grid_data_ag.csv"))
  #system.time(str(sim.dens2b <- post_sim(mod.est2,datamb,Aprediction2, run)))
  

  ###--------Obtain the model validation and model fit metrics
  cpo <- mod.estb$cpo$cpo
  lcpo <- -sum(log(cpo), na.rm=T)
  dic <-  mod.estb$dic$dic
  mets <- model_metrics(dat.predb$pop, sim.densb$est_data$mean_pop_hat, sim.densb$est_data$upper_pop_hat, 
                        sim.densb$est_data$lower_pop_hat)
  
  metrics <- list(lcpo=lcpo, dic=dic, mets = mets)
  capture.output(metrics, file=paste0(result_path, "/metrics_ag.txt"))
  
  
##---Add predictions to data
  dat.fit$pdens_ag <- sim.densb$est_data$mean_dens_hat
  dat.fit$ppop_ag <- sim.densb$est_data$mean_pop_hat
  
  
  ####---Join the posterior sample to the prediction data
  data.simb <- data.frame(cbind(dat.predb[, "regionID"], sim.densb$pop_hat))
  
  names(data.simb)
  names(data.simb)[1] <- "regionID"
  
  ##--------Calculate national estimates and uncertainties (AREA2GRID)
  nat_totalb<- function(dat, run)
  {
    p_hat <- dat[,2:(run+1)]
    tots <- apply(p_hat,2, sum, na.rm=T) #Col sums
    
    tot_sd  <- sd(tots, na.rm=T)
    
    tot_mean  <- mean(tots, na.rm=T)
    
    tot_lower <- quantile(tots, probs=c(0.025))
    tot_median <- quantile(tots, probs=c(0.5))
    tot_upper <- quantile(tots, probs=c(0.975))
    
    return(estimates <- data.frame(estimates=unlist(list(total=tot_mean, lower=tot_lower, median=tot_median, upper=tot_upper))))
  }
  (nationalb <- nat_totalb(data.simb, run))
  (nationalb <- data.frame(total= nationalb[1,],
                           lower = nationalb[2,],
                           median=nationalb[3,],
                           upper=nationalb[4,]))
  
  sum(predb <- round(fitb*dat.estb$bldg2))
  write.csv(nationalb, file=paste0(result_path,"/national_estimates_ag.csv"))
  # write.csv(national, file=paste0(results_path, "/estimates/National_estimates_final.csv"))
  
  
  
  
  (national <- nat_total(data.sim, run))
  (national <- data.frame(total= national[1,],
                          lower = national[2,],
                          median=national[3,],
                          upper=national[4,]))
  
  
  
  ##---Areal level estimates (AREA2GRID)
  admin_estb <- function(datr, run)
  {
    uniR <-as.numeric(unique(datr$regionID))
    outR <- matrix(0, nrow=length(uniR), ncol=5)
    for(j in uniR)
    {
      reg <- datr[datr$regionID==j,]
      rtots <- apply(reg[,2:(1+run)], 2, sum, na.rm=T)
      #rtot_mean  <- sum(reg$pop_hat, na.rm=T)
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
    return(reg_est <- data.frame(ID = uniR,
                                 total = outR[,1],
                                 lower = outR[,2],
                                 median = outR[,3],
                                 upper = outR[,4],
                                 uncertainty = outR[,5]))
  }
  (admin.estb <- admin_estb(data.simb, run))
  sum(admin.estb$total)
  write.csv(admin.estb, file=paste0(result_path,"/area_estimates_ag.csv"))
  

##---------------------------------------
##---Repeat under various proportions of missing observations
###--------------------------------------------------------
  for(j in 1:length(cover))###----Do the following for each coverage percentage 
  {
    
    j=4
    result_path2 <- paste0(result_path,"_at_", cover[j]*100,"%_coverage")
    if (file.exists(result_path2)){
      setwd(file.path(result_path2))
    } else {
      dir.create(file.path(result_path2))
      setwd(file.path(result_path2))
    } 
    
 
##---Select the subset of the data as observed at random.
    pp8 <- cover[j]
    opts8 <- sample(nrow(coord), pp8*nrow(coord))
    dim(dat8 <-  datam[opts8,]) #---Locations with observations
    
##--------------
    coord8 <- cbind(dat8$Lon, dat8$Lat)
    dim(coord8)
    
    
##---Build mesh at grid cell level at the coverage level
    bnd8 <- inla.nonconvex.hull(coord8, -0.035, -0.04, resolution = c(100, 100))
    mesh8<- inla.mesh.2d(boundary = bnd8, 
                         offset=c(0.1, 0.15), max.edge=c(0.070920000001,0.3)) 
    #plot(mesh8)
    mesh8$n
    #plot(Area, add=T, col="light green", lwd=1.5)
    #plot(mesh8, add=T)
    #points(coord8, pch="*", col="brown", cex=1)
    #points(coord[-opts,], pch="*", col="blue", cex=1.5)
    #mesh8$n
    
###---Build projector matrix A 
    A8<-inla.spde.make.A(mesh=mesh8,loc=as.matrix(coord8));dim(A8)
    
##---Create the SPDE
    spde8 <- inla.spde2.matern(mesh8, alpha=2)
    
##----specify the observation indices for estimation 
    iset8<- inla.spde.make.index(name = "spatial.field", spde8$n.spde)
    
###----Scale the covariates for comparability
    stdize <- function(x)
    {
      stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
      return(stdz)
    }
  
    dat.est8 <- dat8
    
    dat.est8[, c("X2", "X3", "X4", "X5", "X6")] <- apply(dat.est8[, c("X2", "X3", "X4", "X5", "X6")], 2, stdize)
    
    covars.est8 <- dat.est8[,c("X2", "X3", "X4", "X5", "X6", "regionID", "dat_typ_ID")]; dim(covars.est8) ##---Population density
    
    
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
    
    
    
    ######---------INLA FORMULA   
    form.est8 <- y ~ -1 + Intercept + X2 + X3 + X4 + X5 + X6 + f(spatial.field, model=spde8)+ 
      f(dat_typ_ID, model='iid') + f(regionID, model='iid') 
    mod.est8 <-inla(form.est8, #the formula
                    data=inla.stack.data(stk.est8,spde=spde8),  #the data stack
                    family= 'gamma',   #which family the data comes from
                    control.predictor=list(A=inla.stack.A(stk.est8),compute=TRUE),  #compute gives you the marginals of the linear predictor
                    control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                    verbose = FALSE) 
    summary(mod.est8)
    
    #saveRDS(mod.est8, file=paste0(result_path2,"/grid_model.rds"))

##---Explore the posterior
    ind8 <-inla.stack.index(stk.est8, "est")$data
    fit8 <- exp(mod.est8$summary.linear.predictor[ind8,"mean"])
  

    sum(pred8 <- round(fit8*dat.est8$bld), na.rm=T)
    #(obs_pred8 <- data.frame(pop8=dat.est8$y,pred8))
    #apply(obs_pred8, 2, sum, na.rm=T)
    #plot(dat.est8$y,pred8)
    #abline(a=0, b=1)
    
    
####----View the spatial fields of the best fit model
    #looking at the spatial field and what it looks like
    gproj8 <- inla.mesh.projector(mesh8,  dims = c(300, 300))
    
    col <- viridis(100)
    
    g.mean8 <- inla.mesh.project(gproj8, mod.est8$summary.random$spatial.field$mean)
    g.sd8 <- inla.mesh.project(gproj8, mod.est8$summary.random$spatial.field$sd)
    
    
    grid.arrange(levelplot(g.mean8, scales=list(draw=F), xlab='', ylab='', cex.lab=2, main='Mean',col.regions = col),
                levelplot(g.sd8, scal=list(draw=F), xla='', yla='', main='SD' ,col.regions = col), nrow=1)
    
    dev.off()
    
    
    
####-----Make prediction at grid cell through posterior simulation
    #############################################################################################
    ###----Prediction Projection matrix 
    Aprediction8 <- inla.spde.make.A(mesh = mesh8, loc = as.matrix(coord))
    dim(Aprediction8)#--prediction projectionnmatrix
    
    dat.fit8 <- dat.fit
    dim(dat.fit8)
    names(dat.fit8)
  
##-----Extract regional effects
    dat.fit8$reg <- data_type_re(dat.fit8, mod.est8$summary.random$regionID$mean,dat.fit8$regionID)
    
##-----Extract data types effects
    dat.fit8$dat_typ <- data_type_re(dat.fit8, mod.est8$summary.random$dat_typ_ID$mean,dat.fit8$dat_typ_ID)
    
    
####-------Run posterior and grid cell prediction (GRID2GRID)
    post_sim <- function(model, dat, Aprediction, run)
    {
      fixedeff <- dens_hat <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
      inla.seed = as.integer(runif(1)*.Machine$integer.max)
      set.seed(inla.seed)
      print(inla.seed)
      m1.samp <- inla.posterior.sample(run, model, seed = inla.seed, selection=list(X2=1, X3=1, X4=1,
                                                                                    X5=1, X6=1),num.threads="1:1")
      
      sfield_nodes_mean <- model$summary.random$spatial.field['mean']
      field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
      for(i in 1:run)
      {
        fixedeff[,i] <- 
          model$summary.fixed['Intercept', 'mean'] +
          m1.samp[[i]]$latent[1,] * dat[,'X2'] +
          m1.samp[[i]]$latent[2,] * dat[,'X3'] +
          m1.samp[[i]]$latent[3,] * dat[,'X4'] +
          m1.samp[[i]]$latent[4,] * dat[,'X5'] +
          m1.samp[[i]]$latent[5,] * dat[,'X6'] + 
          dat$reg + ##--AREAL-LEVEL RANDOM EFFECT
          dat$dat_typ + #Data source random effect
          field_mean[,1]
        
        ##
        dens_hat[,i] <- exp(fixedeff[,i])
        pop_hat[,i] <- dens_hat[,i]*dat$bld
        
      }
      
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
    run=100
    system.time(str(sim.dens8 <- post_sim(mod.est8,dat.fit8,Aprediction8, run)))
    
    write.csv(sim.dens8$est_data, file=paste0(result_path2,"/grid_data_gg.csv"))
   
##---Calculate model fit and model cross-validation metrics  
    cpo <- mod.est8$cpo$cpo
    lcpo <- -sum(log(cpo), na.rm=T)
    dic <-  mod.est8$dic$dic
    mets <- model_metrics(dat.fit8$pop, sim.dens8$est_data$mean_pop_hat, sim.dens8$est_data$upper_pop_hat, 
                          sim.dens8$est_data$lower_pop_hat)
    metrics <- list(lcpo=lcpo, dic=dic, mets=mets)
    capture.output(metrics, file=paste0(result_path2, "/metrics_gg.txt"))
    
    
    
    
    
  ####---Join the posterior sample to the prediction data
    data.sim8 <- data.frame(cbind(dat.fit8[,c("regionID")], sim.dens8$pop_hat))
    
    names(data.sim8)
    
    
    ##--------Calculate and save National total with uncertainties
    nat_total <- function(dat, run)
    {
      p_hat <- dat[,2:(run+1)]
      tots <- apply(p_hat,2, sum, na.rm=T) #Col sums
      
      tot_sd  <- sd(tots, na.rm=T)
      
      tot_mean  <- mean(tots, na.rm=T)
      
      tot_lower <- quantile(tots, probs=c(0.025))
      tot_median <- quantile(tots, probs=c(0.5))
      tot_upper <- quantile(tots, probs=c(0.975))
      
      return(estimates <- data.frame(estimates=unlist(list(total=tot_mean, lower=tot_lower, median=tot_median, upper=tot_upper))))
    }
    (national8 <- nat_total(data.sim8, run))
    (national8 <- data.frame(total= national8[1,],
                             lower = national8[2,],
                             median=national8[3,],
                             upper=national8[4,]))
# write.csv(national, file=paste0(results_path, "/estimates/National_estimates_final.csv"))
    write.csv(national8, file=paste0(result_path2,"/national_estimates_gg.csv"))
    sum(predb <- round(fitb*dat.estb$bld2))
   
##---Regional estimates
    admin_est <- function(datr, run)
    {
      uniR <- unique(datr$X1)
      outR <- matrix(0, nrow=length(uniR), ncol=5)
      for(j in uniR)
      {
        reg <- datr[datr$X1==j,]
        rtots <- apply(reg[,2:(1+run)], 2, sum, na.rm=T)
        #rtot_mean  <- sum(reg$pop_hat, na.rm=T)
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
      return(reg_est <- data.frame(ClusterID = uniR,
                                   total = outR[,1],
                                   lower = outR[,2],
                                   median = outR[,3],
                                   upper = outR[,4],
                                   uncertainty = outR[,5]))
    }
    (admin.est8 <- admin_est(data.sim8, run))
    sum(admin.est8$total)
    write.csv(admin.est8, file=paste0(result_path2,"/area_estimates_gg.csv"))
    

###----------aggregate the partial observation to areal-unit level
    dat.agg81<- data.frame(datam1 <-   dat8 %>%
                            group_by(regionID, dat_typ_ID) %>%
                            summarise(bld2 = sum(bld), pop2 = sum(pop), X2 = mean(X2),
                                      X3=mean(X3), X4= mean(X4), X5=mean(X5), X6 = mean(X6)))
    
    
    dim(dat.agg81); head(dat.agg81)
    
###---
    dat.agg82 <- data.frame(datamm <-  dat.agg81  %>%
                             group_by(regionID) %>%
                             summarise(bld3 = sum(bld2), pop3= sum(pop2), X2 = mean(X2),
                                       X3=mean(X3), X4= mean(X4), X5=mean(X5), X6 = mean(X6))) 
    
    head(dat.agg82); dim(dat.agg82)
    
    locs.agg8 <- coordinates(Area)
    locs.agg8 <- data.frame(ID = 1:nrow(locs.agg8),
                            Lon=locs.agg8[,1],
                            Lat=locs.agg8[,2])
    
    ind.agg8 <- which(dat.agg82$regionID %in% locs.agg8$ID)#--select only observed areas
    dat.agg82$Lon <- locs.agg8$Lon[ind.agg8]
    dat.agg82$Lat <- locs.agg8$Lat[ind.agg8]
    ####
    
  ##---Add Lon-Lat to dat.agg1
    dat.agg81$Lon <- rep(1, nrow(dat.agg81))
    dat.agg81$Lat <- rep(1, nrow(dat.agg81))
    
    for(i in 1:nrow(dat.agg82))
    {
      dat.agg81$Lon[dat.agg81$regionID == i] = dat.agg82$Lon[dat.agg82$regionID == i]
      dat.agg81$Lat[dat.agg81$regionID == i] = dat.agg82$Lat[dat.agg82$regionID == i]
    }
  
    
  ##--Check
    par(mfrow=c(1,2))
    plot(dat.agg81$Lon, dat.agg81$Lat)
    plot(dat.agg82$Lon, dat.agg82$Lat)
    

    
    par(mfrow=c(1,1))
    dat.agg8 <- dat.agg81
    head(dat.agg8)
    #write.csv(data.area, file=paste0(result_path,"/grid_estimates.csv"))
    
    ##---
    write.csv(dat.agg8, file=paste0(result_path2,"/area_data.csv"))
    head(dat.agg8)
    
    length(unique(dat8$regionID))
    ##########
    coord8b <- cbind(dat.agg8$Lon, dat.agg8$Lat)
    
    
    #dat.agg8$pred_ga <- admin.est8$total
    
  #Some functions
    #mesh8b: fine triangulated mesh
    bnd8b <- inla.nonconvex.hull(coord8b, -0.08, -0.7, resolution = c(100, 100))
    mesh8b<- inla.mesh.2d(boundary = bnd8b, 
                          offset=c(0.1, 0.15), max.edge=c(0.07,0.3)) 
    #plot(mesh8b)
    #plot(Area, add=T, col="light green", lwd=1.5)
    #plot(mesh8b, add=T)
    #points(coord8b, pch="*", col="red", cex=1)
    #points(coord[-opts,], pch="*", col="blue", cex=1.5)
    mesh8b$n
    
    
    #--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ###---Build projector matrix A
    A8b<-inla.spde.make.A(mesh=mesh8b,loc=as.matrix(coord8b));dim(A8b)
    
    ##---Create the SPDE
    spde8b <- inla.spde2.matern(mesh8b, alpha=2)
    
    ##----specify the observation indices for estimation 
    iset8b<- inla.spde.make.index(name = "spatial.field", spde8b$n.spde)
    
    
    
  #----Scale areal-level covariates for comparability with the grid-level ones
    dat.agg8$bldg2[dat.agg8$bld2==0]=1
    dat.agg8$dens <- dat.agg8$pop2/dat.agg8$bld2
    dat.est8b <- dat.agg8
    stdize <- function(x)
    {
      stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
      return(stdz)
    }
    
    dat.est8b[, c("X2", "X3", "X4", "X5", "X6")] <- apply(dat.est8b[, c("X2", "X3", "X4", "X5", "X6")], 2, stdize)
    covars.est8b <- dat.est8b[,c("X2", "X3", "X4", "X5", "X6", "regionID", "dat_typ_ID")]; dim(covars.est8b) ##---Population density
    
    
    
    names(datam)
    #---Build the stack
    stk.est8b <- inla.stack(data=list(y=dat.est8b$dens), #the response
                            
                            A=list(A8b,1),  #the A matrix; the 1 is included to make the list(covariates)
                            
                            effects=list(c(list(Intercept=1), #the Intercept
                                           iset8b),  #the spatial index
                                         #the covariates
                                         list(covars.est8b)
                            ), 
                            #this is a quick name so you can call upon easily
                            tag='est')
    
    
    
    ######---------INLA formula   
    form.est8b <- y ~ -1 + Intercept + X2 + X3 + X4 + X5 + X6 + f(spatial.field, model=spde8b)+
      f(regionID, model='iid')+ f(dat_typ_ID, model='iid')
    mod.est8b <-inla(form.est8b, #the formula
                     data=inla.stack.data(stk.est8b,spde=spde8b),  #the data stack
                     family= 'gamma',   #which family the data comes from
                     control.predictor=list(A=inla.stack.A(stk.est8b),compute=TRUE),  #compute gives you the marginals of the linear predictor
                     control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                     verbose = FALSE) 
    summary(mod.est8b)
    
    
    #saveRDS(mod.est8b, file=paste0(result_path2,"/area_model.rds"))
    
  ##---Explore the posterior
    ind8b <-inla.stack.index(stk.est8b, "est")$data
    fit8b <- exp(mod.est8b$summary.linear.predictor[ind8b,"mean"])
    
    sum(pred8b <- round(fit8b*dat.est8b$bld2))

    #(obs_pred8b <- data.frame(popb=dat.est8b$pop2,pred8b))
    #apply(obs_pred8b, 2, sum, na.rm=T)
    #plot(dat.est8b$pop2,pred8b)
    #abline(a=0, b=1)
    
    
    
    
####-----Make prediction at grid cell level

    dat.fit8b <- dat.fit
    dim(dat.fit8b)
    
    coord88b <- cbind(dat.fit8b$Lon, dat.fit8b$Lat)
    Aprediction8b <- inla.spde.make.A(mesh = mesh8b, loc = as.matrix(coord88b))
    dim(Aprediction8b)
    
##-----Extract regional effects
    dat.fit8b$reg <- data_type_re(dat.fit8b, mod.est8b$summary.random$regionID$mean,dat.fit8b$regionID)
    
##-----Extract data types effects
    dat.fit8b$dat_typ <- data_type_re(dat.fit8b, mod.est8b$summary.random$dat_typ_ID$mean,dat.fit8b$dat_typ_ID)
    
####----Run posterior Simulation and grid cell predictions (AREA2GRID)
    post_sim <- function(model, dat, Aprediction, run)
    {
      fixedeff <- dens_hat <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
      inla.seed = as.integer(runif(1)*.Machine$integer.max)
      #inla.seed = 
      
      set.seed(inla.seed)
      print(inla.seed)
      m1.samp <- inla.posterior.sample(run, model, seed = inla.seed, selection=list(X2=1, X3=1, X4=1,
                                                                                    X5=1, X6=1),num.threads="1:1")
      
      sfield_nodes_mean <- model$summary.random$spatial.field['mean']
      field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
      for(i in 1:run)
      {
        fixedeff[,i] <- 
          model$summary.fixed['Intercept', 'mean'] +
          m1.samp[[i]]$latent[1,] * dat[,'X2'] +
          m1.samp[[i]]$latent[2,] * dat[,'X3'] +
          m1.samp[[i]]$latent[3,] * dat[,'X4'] +
          m1.samp[[i]]$latent[4,] * dat[,'X5'] +
          m1.samp[[i]]$latent[5,] * dat[,'X6'] + 
          dat$reg +
          dat$dat_typ +
          field_mean[,1]
        
        ##
        dens_hat[,i] <- exp(fixedeff[,i])
        pop_hat[,i] <- dens_hat[,i]*dat$bld
        
      }
      
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
    run=100
  system.time(str(sim.dens8b <- post_sim(mod.est8b,dat.fit8b,Aprediction8b, run)))
 
    #data_ag[,j] <- sim.dens8b$est_data$mean_pop_hat
    write.csv(sim.dens8b$est_data, file=paste0(result_path2,"/grid_data_ag.csv"))
    #write.csv(data.area, file=paste0(result_path,"/grid_estimates.csv"))
###-----Model cross-validation and fit checks
    cpo <- mod.est8b$cpo$cpo
    lcpo <- -sum(log(cpo), na.rm=T)
    dic <-  mod.est8b$dic$dic
    mets <- model_metrics(dat.fit8b$pop, sim.dens8b$est_data$mean_pop_hat, sim.dens8b$est_data$upper_pop_hat, 
                          sim.dens8b$est_data$lower_pop_hat)
    metrics <- list(lcpo=lcpo, dic=dic, mets=mets)
    capture.output(metrics, file=paste0(result_path2, "/metrics_ag.txt"))
    
  ####---Join the posterior sample to the prediction data
    data.sim8b <- data.frame(cbind(dat.fit8b[,c("regionID")], sim.dens8b$pop_hat))
    
    names(data.sim8b)
    
    
    ##--------Calculate and save National total with uncertainties
    nat_total <- function(dat, run)
    {
      p_hat <- dat[,2:(run+1)]
      tots <- apply(p_hat,2, sum, na.rm=T) #Col sums
      
      tot_sd  <- sd(tots, na.rm=T)
      
      tot_mean  <- mean(tots, na.rm=T)
      
      tot_lower <- quantile(tots, probs=c(0.025))
      tot_median <- quantile(tots, probs=c(0.5))
      tot_upper <- quantile(tots, probs=c(0.975))
      
      return(estimates <- data.frame(estimates=unlist(list(total=tot_mean, lower=tot_lower, median=tot_median, upper=tot_upper))))
    }
    (national8b <- nat_total(data.sim8b, run))
    (national8b <- data.frame(total= national8b[1,],
                              lower = national8b[2,],
                              median=national8b[3,],
                              upper=national8b[4,]))
    write.csv(national8b, file=paste0(result_path2,"/national_estimates_ag.csv"))
    #write.csv(data.area, file=paste0(result_path,"/grid_estimates.csv"))
    # write.csv(national, file=paste0(results_path, "/estimates/National_estimates_final.csv"))
    
    
    
    
    (national8 <- nat_total(data.sim8, run))
    (national8 <- data.frame(total= national8[1,],
                             lower = national8[2,],
                             median=national8[3,],
                             upper=national8[4,]))
    
    
    (national <- nat_total(data.sim, run))
    #plot(dat.fit$y,sim.dens8$est_data$mean_pop_hat)
    #cor(dat.fit$y,sim.dens8$est_data$mean_pop_hat)
    
    
    #plot(sim.dens8b$est_data$mean_pop_hat, dat.fit$y)
    #plot(sim.dens8b$est_data$mean_dens_hat,dat.fit$dens)
    #abline(a=0, b=1)
    #cor(dat.fit$y,sim.dens8b$est_data$mean_pop_hat)
    #
    ##---Regional/Areal estimates
    admin_est <- function(datr, run)
    {
      uniR <- unique(datr$X1)
      outR <- matrix(0, nrow=length(uniR), ncol=5)
      for(j in uniR)
      {
        reg <- datr[datr$X1==j,]
        rtots <- apply(reg[,2:(1+run)], 2, sum, na.rm=T)
        #rtot_mean  <- sum(reg$pop_hat, na.rm=T)
        rtot_mean  <- mean(rtots, na.rm=T)
        rtot_sd <- sd(rtots, na.rm=T)
        rtot_lower <- quantile(rtots, probs=c(0.025))
        rtot_median <- quantile(rtots, probs=c(0.5))
        rtot_upper <- quantile(rtots, probs=c(0.975))
        rtot_uncert <- (rtot_upper - rtot_lower)/rtot_mean
        
        restimates <- round(c(rtot_mean, rtot_lower, rtot_median,rtot_upper, rtot_uncert),4)
        outR[j,] <- restimates
      }
      outR <- data.frame(outR)
      return(reg_est <- data.frame(ClusterID = uniR,
                                   total = outR[,1],
                                   lower = outR[,2],
                                   median = outR[,3],
                                   upper = outR[,4],
                                   uncertainty = outR[,5]))
    }
    (admin.est8b <- admin_est(data.sim8b, run))
    sum(admin.est8b$total)
    
    write.csv(admin.est8b, file=paste0(result_path2,"/area_estimates_ag.csv"))

  }
  
}



#save.image(paste0(results_path, "/posterior_samples/CMR_MODEL_main_workspace_final.Rdata"))

