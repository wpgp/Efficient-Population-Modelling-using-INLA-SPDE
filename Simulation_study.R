####--TITLE: R SCRIPTS FOR MODELLED POPULATION ESTIMATES BASED ON R-INLA:SIMULATION STUDY----#####
####--METHODS: GEOSTATISTICAL BAYESIAN HIERARCHICAL REGRESSION MODEL---------------------#####
####--AUTHOR: DR CHIBUZOR CHRISTOPHER NNANATU ------------------------#####
####--INSTITUTION: WORLDPOP, UNIVERSITY OF SOUTHAMPTON----------------------------------##### 
####--DATE: DECEMBER 2022---------------------------------------------------------------####
####@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@####



library(INLA); library(raster); library(maptools)
library(gtools); library(sp); library(spdep)
library(fields); library(mvtnorm); library(gtools)
library(geoR)
library(actuar)
library(viridisLite)
require(grid)
require(gridExtra)
require(lattice)
require(tidyverse)

path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3/Working/CMR/Chris_N/codes/INLA_models/final_models/updated/sim_study"
file_path <- paste0(path, "/files")
out_path <- paste0(path, "/outputs")

AdmUnit <- c(8, 10, 20, 30, 40)
cover <- c(0.8, 0.6, 0.4, 0.2) 
pg <- 120# = 14400 prediction grids



####
model_metrics <- function(obs, pred, upper, lower)
{
  residual = pred - obs
  INACCURACY = mean(abs(residual), na.rm=T)#MAE
  MSE = mean(residual^2, na.rm=T)
  RMSE = sqrt(MSE)
  BIAS = mean(residual, na.rm=T)
  In_IC = mean(obs<upper & obs> lower, na.rm=T)*100
  corr = cor(obs[!is.na(obs)],pred[!is.na(obs)])
  
  output <- list(MAE  = INACCURACY ,
                 RMSE = RMSE,
                 BIAS = abs(BIAS),
                 In_IC = In_IC,
                 corr = corr)
  return(output)
}
#model_metrics(obs, pred, upper, lower)

###

for(i in AdmUnit)
{
print(i)
result_path <- paste0(out_path,"/estimates_for_", i^2,"_area_units")
if (file.exists(result_path)){
  setwd(file.path(result_path))
} else {
  dir.create(file.path(result_path))
  setwd(file.path(result_path))
}




Area <- as(raster(nrow=i, ncol=i, xmn=0, xmx=1, ymn=0, ymx=1), "SpatialPolygons")
proj4string(Area) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#
tiff(paste0(result_path,"/shape_plot.tiff"), width = 4, height = 4, units = 'in', res = 300)
plot(Area, lwd=2)
dev.off()

############
#Create raster
pgd = raster(nrow=pg, ncol=pg, xmn=0, xmx=1, ymn=0, ymx=1) #--prediction grid
pgd1 = pgd


######
coord.1 = xyFromCell(pgd1, 1:ncell(pgd1), spatial=FALSE)
coord = coord.1
#dim(coord)
#plot(coord)
#plot(Area, lwd=2, add=T)

####
spol = Area
#Generate points
c.sp=SpatialPoints(coord)
proj4string(c.sp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
sp.1 <- rep(NA, nrow(coord))
for(i in 1:length(spol)){
  sp.1[as.vector(which(!is.na(over(c.sp, spol[i]))))] <- i
}
spred <- sp.1
#table(spred)

#plot(coord)
#plot(Area, add=T, lwd=1.5)

dim(dat.all <- as.data.frame(coord))
table(dat.all$regionID <- as.factor(spred))
table(dat.all$clusterID <- as.factor(1:nrow(coord)))



dim(coord)
bnd <- inla.nonconvex.hull(coord, -0.035, -0.04, resolution = c(100, 100))
meshfit <- inla.mesh.2d(boundary = bnd, 
                        offset=c(0.1, 0.15), max.edge=c(0.06472101,0.3)) 





#--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r0 <- 0.3
nu <- 1
sigma0 <- 1
kappa0 <- sqrt(8*nu)/r0
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)

#
spde <- inla.spde2.matern(meshfit, 
                          B.tau = matrix(c(log(tau0), -1, +1), nrow=1, ncol=3),
                          B.kappa = matrix(c(log(kappa0), 0, -1), nrow=1, ncol=3),
                          theta.prior.mean = c(0,0),
                          theta.prior.prec = c(0.1, 0.1))



#
Q <- inla.spde2.precision(spde=spde, theta=c(0,0))


#---Simulate the GMRF
sam <- as.vector(inla.qsample(
  n = 1, Q = Q, seed=100))
#length(sam)



###---Build projector matrix A
A <- inla.spde.make.A(mesh=meshfit, loc=coord);dim(A)
S.pred <- as.vector(A %*% sam)
#hist(S.pred)


#####
tau2 =1
dim(Q); class(Q)
Q <- as.matrix(Q)
sig <- solve(Q)
mm <- meshfit$n
phi.a = rmvnorm(1,rep(0,mm), sig)
phi.a = as.vector(phi.a)  #for i=1,...,n

#Construct phi for the pixels
phi.pred <- numeric()
for(i in 1:length(spred)) phi.pred[i] <- phi.a[spred[i]]


##----specify the observation indices for estimation 
iset <- inla.spde.make.index(name = "spatial.field", spde$n.spde)

#length(iset)
## Nugget effect/iid term
#set.seed(20)
sigma_e <- 0.05
ee.pred <- rnorm(nrow(coord), 0, sigma_e) 

sig_reg= 0.08
reg.pred <- rnorm(nrow(coord), 0,sig_reg) 
#hist(reg.pred)

#Covariate vector 
nn <- nrow(dat.all)
covs = cbind(1,runif(nn),rnorm(nn),rpois(nn,2),rnorm(nn), runif(nn))

apply(covs, 2, mean)
dim(dat.all)

dim(ddat <- cbind(dat.all, covs))
ddat <- data.frame(ddat)


###############
#Parameters 
dim(zpred <- covs)


rr=0.01


#beta <- as.vector(round(apply(zpred, 2, mean),3))
beta <- c(1.5, 0.76, 0.2, 0.2, -0.028, 0.5)
resp <- mu <- numeric(nn) #
for (i in 1:nn)
{
  #i=5
  mu[i] <- exp(zpred[i,1]*beta[1] + zpred[i,2]*beta[2] + 
                 zpred[i,3]*beta[3] + zpred[i,4]*beta[4] + 
                 zpred[i,5]*beta[5] + zpred[i,6]*beta[6] + S.pred[i] + phi.pred[i] + ee.pred[i])
  
 resp[i] <- rgamma(1, mu[i]^2/rr,mu[i]/rr)
}

mu[is.na(mu)]

bldg <- rpois(nn,50)
##pop <- rpois(nn, resp*bldg)
#resp[resp<0.00001]=NA

pop <- rpois(nn, resp*bldg)# THERE is at least one building per unit



#####

datam <- ddat %>% 
  mutate(Lon = x, Lat=y) %>%
  select(-x, -y)
datam <- as.data.frame(datam)

class(datam)
names(datam)
dim(ddat)
###
datam$bldg <- bldg
datam$y <- pop
datam$dens <- resp
#datam$dens <- mu
hist(resp)



names(datam)
#----Declare covariates for stack
#----Declare covariates for stack
stdize <- function(x)
{
  stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  return(stdz)
}

dim(dat.fit <- datam)
dat.fit[, c("X2", "X3", "X4", "X5", "X6")] <- apply(dat.fit[, c("X2", "X3", "X4", "X5", "X6")], 2, stdize)

covars.fit <- dat.fit[,c("X1", "X2", "X3", "X4", "X5", "X6", "regionID", "clusterID")]; dim(covars.fit) ##---Population density


names(dat.fit)

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



######---------   
form.fit <- y ~ -1 + Intercept + X2 + X3 + X4 + X5 + X6 + #f(clusterID, model='iid')+ 
  f(spatial.field, model=spde) + f(regionID, model='iid')
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
fit1L <- exp(mod.fit$summary.linear.predictor[ind.fit,"0.025quant"])
fit1U <- exp(mod.fit$summary.linear.predictor[ind.fit,"0.975quant"])


sum(pred1 <- round(fit1*dat.fit$bldg))
sum(pred1L <- round(fit1L*dat.fit$bldg))
sum(pred1U <- round(fit1U*dat.fit$bldg))


dat.fit$pdens1 <-fit1 
dat.fit$pred1 <- pred1

#########################################################################################################
####----View the spatial fields of the best fit model
#looking at the spatial field and what it looks like
gproj <- inla.mesh.projector(meshfit,  dims = c(300, 300))




col <- viridis(100)

g.mean <- inla.mesh.project(gproj, mod.fit$summary.random$spatial.field$mean)
g.sd <- inla.mesh.project(gproj, mod.fit$summary.random$spatial.field$sd)



#grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', cex.lab=2, main='Mean',col.regions = col),
#             levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='SD' ,col.regions = col), nrow=1)


dat.fit2 <- dat.fit
dim(dat.fit2)
names(dat.fit2)
dat.fit2$regionID <- dat.all$regionID
########
dat.fit2$reg <- rep(1, nrow(dat.fit2))
rg <- mod.fit$summary.random$regionID$mean
dat.fit2$regionID <- as.numeric(dat.fit2$regionID)
uniq <- unique(dat.fit2$regionID)
for(i in  1:nrow(dat.fit2))
{
  
  for(j in uniq)
  {
    if(dat.fit2$regionID[i]==uniq[j]) dat.fit2$reg[i] = rg[j]
  }
  
}
dat.fit2$reg



####-------Run the uncertainty estimates 
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
      field_mean[,1]
    
    ##
    dens_hat[,i] <- exp(fixedeff[,i])
    pop_hat[,i] <- dens_hat[,i]*dat$bldg
    
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
system.time(str(sim.dens <- post_sim(mod.fit,dat.fit2,A, run)))
write.csv(sim.dens$est_data, file=paste0(result_path,"/grid_data_gg.csv"))



cpo <- mod.fit$cpo$cpo
lcpo <- -sum(log(cpo), na.rm=T)
dic <-  mod.fit$dic$dic
mets <- model_metrics(dat.fit$y, sim.dens$est_data$mean_pop_hat, sim.dens$est_data$upper_pop_hat, 
                      sim.dens$est_data$lower_pop_hat)
metrics <- list(lcpo=lcpo, dic=dic, mets=mets)
capture.output(metrics, file=paste0(result_path, "/metrics_grid.txt"))




dat.fit$pdens_gg <- sim.dens$est_data$mean_dens_hat
dat.fit$ppop_gg <- sim.dens$est_data$mean_pop_hat
mean.pred <- sim.dens$est_data$mean_pop_hat



#####
dat.fit2$pred <- mean.pred 



####---Join the posterior sample to the prediction data
data.sim <- data.frame(cbind(dat.fit2[,c("regionID")], sim.dens$pop_hat))

names(data.sim)[1] <- "regionID"


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
(national <- nat_total(data.sim, run))
(national <- data.frame(total= national[1,],
                        lower = national[2,],
                        median=national[3,],
                        upper=national[4,]))

sum(pred1 <- round(fit1*dat.fit$bldg))
write.csv(national, file=paste0(result_path,"/national_estimates_gg.csv"))
# write.csv(national, file=paste0(results_path, "/estimates/National_estimates_final.csv"))


##---Regional estimates
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
write.csv(admin.est, file=paste0(result_path,"/area_estimates_ga.csv"))




##################################################################################################
#dim(datm <- datam[opts,])
dat.agg <- data.frame(datam1 <-  datam %>%
                        group_by(regionID) %>%
                        summarise(bldg2 = sum(bldg), pop2 = sum(y), X2 = mean(X2),
                                  X3=mean(X3), X4= mean(X4), X5=mean(X5), X6 = mean(X6)))

head(dat.agg)


##---
dat.agg$Lon <- coordinates(Area)[,1]
dat.agg$Lat <- coordinates(Area)[,2]
#plot(dat.agg$Lon, dat.agg$Lat)
names(dat.agg)


##########
#par(mfrow=c(1,1))
coordb <- cbind(dat.agg$Lon, dat.agg$Lat)
#plot(coordb)

#plot(Area, add=T)



#Some functions
#meshfit: fine triangulated mesh
bndb <- inla.nonconvex.hull(coordb, -0.08, -0.7, resolution = c(100, 100))
meshb<- inla.mesh.2d(boundary = bndb, 
                     offset=c(0.1, 0.15), max.edge=c(0.07,0.3)) 

meshb$n


#--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###---Build projector matrix A
Ab<-inla.spde.make.A(mesh=meshb,loc=as.matrix(coordb));dim(Ab)

##---Create the SPDE
spdeb <- inla.spde2.matern(meshb, alpha=2)

##----specify the observation indices for estimation 
isetb<- inla.spde.make.index(name = "spatial.field", spdeb$n.spde)



#----Declare covariates for stack
dat.agg$dens <- dat.agg$pop2/dat.agg$bldg2

stdize <- function(x)
{
  stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  return(stdz)
}

dat.estb <- dat.agg
names(dat.agg)
dat.estb[, c("X2", "X3", "X4", "X5", "X6")] <- apply(dat.estb[, c("X2", "X3", "X4", "X5", "X6")], 2, stdize)

covars.estb <- dat.estb[,c("X2", "X3", "X4", "X5", "X6", "regionID")]; dim(covars.estb) ##---Population density



#---Build the stack
stk.estb <- inla.stack(data=list(y=dat.estb$dens), #the response
                       
                       A=list(Ab,1),  #the A matrix; the 1 is included to make the list(covariates)
                       
                       effects=list(c(list(Intercept=1), #the Intercept
                                      isetb),  #the spatial index
                                    #the covariates
                                    list(covars.estb)
                       ), 
                       #this is a quick name so you can call upon easily
                       tag='est')



######---------   
form.estb <- y ~ -1 + Intercept + X2 + X3 + X4 + X5 + X6 + 
  f(regionID, model='iid')+ f(spatial.field, model=spdeb)
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


sum(predb <- round(fitb*dat.estb$bldg2))


dat.agg$pred1 <- predb
dat.agg$dens1 <- fitb
dat.agg$ppop_ga <- admin.est$total

#############################################################################################
###----Prediction Projection matrix 
#dim(dat.pred <- datam[-opts,])
dim(dat.predb <- dat.fit)
#dim(coord.pred <- coord[-opts,])
dim(coord.predb <- cbind(dat.predb$Lon, dat.predb$Lat))
Apredictionb<- inla.spde.make.A(mesh = meshb, loc = coord.predb); dim(Apredictionb)
dat.fit[, c("X2", "X3", "X4", "X5", "X6")]
covars.estb <- dat.predb[,c("X1", "X2", "X3", "X4", "X5", "X6", "regionID")]; dim(covars.estb) ##---Population density



########
dat.predb$reg <- rep(1, nrow(dat.predb))
rgb <- mod.estb$summary.random$regionID$mean
dat.predb$regionID <- as.numeric(dat.predb$regionID)
for(i in  1:nrow(dat.predb))
{uniq <- unique(dat.predb$regionID)

  
  for(j in uniq)
  {
    if(dat.predb$regionID[i]==uniq[j]) dat.predb$reg[i] = rgb[j]
  }
  
}
dat.predb$reg



####-------Run the uncertainty estimates 
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

      field_mean[,1]
    
    ##
    dens_hat[,i] <- exp(fixedeff[,i])
    pop_hat[,i] <- dens_hat[,i]*dat$bldg
    
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

cpo <- mod.estb$cpo$cpo
lcpo <- -sum(log(cpo), na.rm=T)
dic <-  mod.estb$dic$dic
mets <- model_metrics(dat.predb$y, sim.densb$est_data$mean_pop_hat, sim.densb$est_data$upper_pop_hat, 
                      sim.densb$est_data$lower_pop_hat)
metrics <- list(lcpo=lcpo, dic=dic, mets = mets)
capture.output(metrics, file=paste0(result_path, "/metrics_area.txt"))


dat.fit$pdens_ag <- sim.densb$est_data$mean_dens_hat
dat.fit$ppop_ag <- sim.densb$est_data$mean_pop_hat



####---Join the posterior sample to the prediction data
data.simb <- data.frame(cbind(dat.predb[, "regionID"], sim.densb$pop_hat))

names(data.simb)
names(data.simb)[1] <- "regionID"

##--------Calculate and save National total with uncertainties
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

##---Regional estimates
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




for(j in 1:length(cover))
{
  result_path2 <- paste0(result_path,"_at_", cover[j]*100,"%_coverage")
  if (file.exists(result_path2)){
    setwd(file.path(result_path2))
  } else {
    dir.create(file.path(result_path2))
    setwd(file.path(result_path2))
  } 



pp8 <- cover[j]
opts8 <- sample(nrow(coord), pp8*nrow(coord))
dim(dat8 <-  datam[opts8,]) #---Locations with observations



coord8 <- cbind(dat8$Lon, dat8$Lat)
dim(coord8)
#plot(coord8)





#######
bnd8 <- inla.nonconvex.hull(coord8, -0.035, -0.04, resolution = c(100, 100))
mesh8<- inla.mesh.2d(boundary = bnd8, 
                     offset=c(0.1, 0.15), max.edge=c(0.070920000001,0.3)) 


#plot(mesh8)
mesh8$n



#--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###---Build projector matrix A
A8<-inla.spde.make.A(mesh=mesh8,loc=as.matrix(coord8));dim(A8)

##---Create the SPDE
spde8 <- inla.spde2.matern(mesh8, alpha=2)

##----specify the observation indices for estimation 
iset8<- inla.spde.make.index(name = "spatial.field", spde8$n.spde)



#----Declare covariates for stack
stdize <- function(x)
{
  stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  return(stdz)
}


dat.est8 <- dat8

dat.est8[, c("X2", "X3", "X4", "X5", "X6")] <- apply(dat.est8[, c("X2", "X3", "X4", "X5", "X6")], 2, stdize)

covars.est8 <- dat.est8[,c("X2", "X3", "X4", "X5", "X6", "regionID")]; dim(covars.est8) ##---Population density



names(datam)
#---Build the stack
stk.est8 <- inla.stack(data=list(y=dat.est8$dens), #the response
                      
                      A=list(A8,1),  #the A matrix; the 1 is included to make the list(covariates)
                      
                      effects=list(c(list(Intercept=1), #the Intercept
                                     iset8),  #the spatial index
                                   #the covariates
                                   list(covars.est8)
                      ), 
                      #this is a quick name so you can call upon easily
                      tag='est')



######---------   
form.est8 <- y ~ -1 + Intercept + X2 + X3 + X4 + X5 + X6 + #f(clusterID, model='iid')+ 
   f(spatial.field, model=spde8)+ f(regionID, model='iid')
mod.est8 <-inla(form.est8, #the formula
               data=inla.stack.data(stk.est8,spde=spde8),  #the data stack
               family= 'gamma',   #which family the data comes from
               control.predictor=list(A=inla.stack.A(stk.est8),compute=TRUE),  #compute gives you the marginals of the linear predictor
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
               verbose = FALSE) 
summary(mod.est8)

#saveRDS(mod.est8, file=paste0(result_path2,"/grid_model.rds"))

ind8 <-inla.stack.index(stk.est8, "est")$data
fit8 <- exp(mod.est8$summary.linear.predictor[ind8,"mean"])
fit8L <- exp(mod.est8$summary.linear.predictor[ind8,"0.025quant"])
fit8U <- exp(mod.est8$summary.linear.predictor[ind8,"0.975quant"])


sum(pred8 <- round(fit8*dat.est8$bldg))
sum(pred8L <- round(fit8L*dat.est8$bldg))
sum(pred8U <- round(fit8U*dat.est8$bldg))




#########################################################################################################
####----View the spatial fields of the best fit model
#looking at the spatial field and what it looks like
gproj8 <- inla.mesh.projector(mesh8,  dims = c(300, 300))



col <- viridis(100)

g.mean8 <- inla.mesh.project(gproj8, mod.est8$summary.random$spatial.field$mean)
g.sd8 <- inla.mesh.project(gproj8, mod.est8$summary.random$spatial.field$sd)


#grid.arrange(levelplot(g.mean8, scales=list(draw=F), xlab='', ylab='', cex.lab=2, main='Mean',col.regions = col),
 #            levelplot(g.sd8, scal=list(draw=F), xla='', yla='', main='SD' ,col.regions = col), nrow=1)
#
#dev.off()



####-----Make prediction
#############################################################################################
###----Prediction Projection matrix 
Aprediction8 <- inla.spde.make.A(mesh = mesh8, loc = as.matrix(coord))
dim(Aprediction8)

#covars.est <- dat.est[,c("X1", "X2", "X3", "X4", "X5", "X6", "clusterID")]; dim(covars.est) ##---Population density


#dat.fit8 <- dat.pred8
dat.fit8 <- dat.fit
dim(dat.fit8)
names(dat.fit8)
########
dat.fit8$reg <- rep(1, nrow(dat.fit8))
rg8 <- mod.est8$summary.random$regionID$mean
dat.fit8$regionID <- as.numeric(dat.fit8$regionID)
uniq <- unique(dat.fit8$regionID)
for(i in  1:nrow(dat.fit8))
{

  
  for(j in uniq)
  {
    if(dat.fit8$regionID[i]==uniq[j]) dat.fit8$reg[i] = rg8[j]
  }
  
}
dat.fit8$reg



####-------Run the uncertainty estimates 
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
     # rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[2]) + #---cluster level random effect
      field_mean[,1]
    
    ##
    dens_hat[,i] <- exp(fixedeff[,i])
    pop_hat[,i] <- dens_hat[,i]*dat$bldg
    
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
#hist(log(sim.dens8$est_data$mean_dens_hat))


#data_gg[,j] <- sim.dens8$est_data$mean_pop_hat

cpo <- mod.est8$cpo$cpo
lcpo <- -sum(log(cpo), na.rm=T)
dic <-  mod.est8$dic$dic
mets <- model_metrics(dat.fit8$y, sim.dens8$est_data$mean_pop_hat, sim.dens8$est_data$upper_pop_hat, 
                      sim.dens8$est_data$lower_pop_hat)
metrics <- list(lcpo=lcpo, dic=dic, mets=mets)
capture.output(metrics, file=paste0(result_path2, "/metrics_grid.txt"))





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
sum(predb <- round(fitb*dat.estb$bldg2))
# write.csv(national, file=paste0(results_path, "/estimates/National_estimates_final.csv"))




(national <- nat_total(data.sim, run))

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
write.csv(admin.est8, file=paste0(result_path2,"/area_estimates_ga.csv"))

########
##-----Agggregate





##################################################################################################

dat.agg8<- data.frame(datam1 <-  dat8 %>%
  group_by(regionID) %>%
  summarise(bldg2 = sum(bldg), pop2 = sum(y), X2 = mean(X2),
            X3=mean(X3), X4= mean(X4), X5=mean(X5), X6 = mean(X6)))


dim(dat.agg8)
write.csv(dat.agg8, file=paste0(result_path2,"/area_data.csv"))
#write.csv(data.area, file=paste0(result_path,"/grid_estimates.csv"))

##---
dat.agg8$Lon <- coordinates(Area)[,1]
dat.agg8$Lat <- coordinates(Area)[,2]



length(unique(dat8$regionID))
##########
coord8b <- cbind(dat.agg8$Lon, dat.agg8$Lat)


#dat.agg8$pred_ga <- admin.est8$total

#Some functions
#meshfit: fine triangulated mesh


bnd8b <- inla.nonconvex.hull(coord8b, -0.08, -0.7, resolution = c(100, 100))
mesh8b<- inla.mesh.2d(boundary = bnd8b, 
                        offset=c(0.1, 0.15), max.edge=c(0.07,0.3)) 

mesh8b$n


#--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###---Build projector matrix A
A8b<-inla.spde.make.A(mesh=mesh8b,loc=as.matrix(coord8b));dim(A8b)

##---Create the SPDE
spde8b <- inla.spde2.matern(mesh8b, alpha=2)

##----specify the observation indices for estimation 
iset8b<- inla.spde.make.index(name = "spatial.field", spde8b$n.spde)



#----Declare covariates for stack
dat.agg8$bldg2[dat.agg8$bldg2==0]=1
dat.agg8$dens <- dat.agg8$pop2/dat.agg8$bldg2
dat.est8b <- dat.agg8
stdize <- function(x)
{
  stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  return(stdz)
}

dat.est8b[, c("X2", "X3", "X4", "X5", "X6")] <- apply(dat.est8b[, c("X2", "X3", "X4", "X5", "X6")], 2, stdize)

covars.est8b <- dat.est8b[,c("X2", "X3", "X4", "X5", "X6", "regionID")]; dim(covars.est8b) ##---Population density



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



######---------   
form.est8b <- y ~ -1 + Intercept + X2 + X3 + X4 + X5 + X6 + 
  f(regionID, model='iid')+ f(spatial.field, model=spde8b)
mod.est8b <-inla(form.est8b, #the formula
               data=inla.stack.data(stk.est8b,spde=spde8b),  #the data stack
               family= 'gamma',   #which family the data comes from
               control.predictor=list(A=inla.stack.A(stk.est8b),compute=TRUE),  #compute gives you the marginals of the linear predictor
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
               verbose = FALSE) 
summary(mod.est8b)


#saveRDS(mod.est8b, file=paste0(result_path2,"/area_model.rds"))


ind8b <-inla.stack.index(stk.est8b, "est")$data
fit8b <- exp(mod.est8b$summary.linear.predictor[ind8b,"mean"])
fit8bL <- exp(mod.est8b$summary.linear.predictor[ind8b,"0.025quant"])
fit8bU <- exp(mod.est8b$summary.linear.predictor[ind8b,"0.975quant"])


sum(pred8b <- round(fit8b*dat.est8b$bldg2))
sum(pred8bL <- round(fit8bL*dat.est8b$bldg2))
sum(pred8bU <- round(fit8bU*dat.est8b$bldg2))




####-----Make prediction
#############################################################################################
###----Prediction Projection matrix 

Aprediction8b <- inla.spde.make.A(mesh = mesh8b, loc = as.matrix(coord))
dim(Aprediction8b)


#dat.pred8b <- dat.fit
#datamb<- datam
#dat.pred8b[, c("X2", "X3", "X4", "X5", "X6")] 


dat.fit8b <- dat.fit
dim(dat.fit8b)
names(dat.fit8b)
dat.fit8b$reg <- rep(1, nrow(dat.fit8b))
rg8b <- mod.est8b$summary.random$regionID$mean
dat.fit8b$regionID <- as.numeric(dat.fit8b$regionID)
uniq <- unique(dat.fit8b$regionID)
for(i in  1:nrow(dat.fit8b))
{
  
  for(j in uniq)
  {
    if(dat.fit8b$regionID[i]==uniq[j]) dat.fit8b$reg[i] = rg8b[j]
  }
  
}
length(dat.fit8b$reg)



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
      #rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[2]) + #---cluster level random effect
      field_mean[,1]
    
    ##
    dens_hat[,i] <- exp(fixedeff[,i])
    pop_hat[,i] <- dens_hat[,i]*dat$bldg
    
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
#system.time(str(sim.dens2b <- post_sim(mod.est2,datamb,Aprediction2, run)))

#hist(log(sim.dens8b$est_data$mean_dens_hat))


#data_ag[,j] <- sim.dens8b$est_data$mean_pop_hat
write.csv(sim.dens8b$est_data, file=paste0(result_path2,"/grid_data_ag.csv"))
#write.csv(data.area, file=paste0(result_path,"/grid_estimates.csv"))

cpo <- mod.est8b$cpo$cpo
lcpo <- -sum(log(cpo), na.rm=T)
dic <-  mod.est8b$dic$dic
mets <- model_metrics(dat.fit8b$y, sim.dens8b$est_data$mean_pop_hat, sim.dens8b$est_data$upper_pop_hat, 
                      sim.dens8b$est_data$lower_pop_hat)
metrics <- list(lcpo=lcpo, dic=dic, mets=mets)
capture.output(metrics, file=paste0(result_path2, "/metrics_area.txt"))

###############


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
(admin.est8b <- admin_est(data.sim8b, run))
sum(admin.est8b$total)

write.csv(admin.est8b, file=paste0(result_path2,"/area_estimates_aa.csv"))



}

}



#save.image(paste0(results_path, "/posterior_samples/CMR_MODEL_main_workspace_final.Rdata"))
