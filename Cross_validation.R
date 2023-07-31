

####--TITLE: R SCRIPTS FOR MODELLED POPULATION ESTIMATES FOR CAMEROON BASED ON R-INLA----#####
####--METHODS: GEOSTATISTICAL BAYESIAN HIERARCHICAL REGRESSION MODEL---------------------#####
####--AUTHOR: DR CHIBUZOR CHRISTOPHER NNANATU & DR ORTIS YANKEY-------------------------#####
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


#################################################################################
#####-------------------Cross validation
####################################################################################
#(metrics <- model_metrics(dat2$Imputed_LHHSIZE,  
#                         pred9a, pred9aU,  pred9aL))

cross_val <- function(dat, mod, mod_pred, mesh, spde,
                      A, shp, formula, k_folds)
{
  set.seed(13235)                                     
  N <- nrow(dat)
  ######
  ind_train <- factor(sample(x = rep(1:k_folds, each = floor(N / k_folds)),  # Sample IDs for training data
                             size = N))
  
  table(as.numeric(ind_train)) 
  dat$k_fold <- as.numeric(ind_train)
  coords <- cbind(dat$longitude, dat$latitude)
  
  
  k_uniq <-sort(unique(dat$k_fold))
  metrics_cv <- matrix(0, nrow=length(k_uniq), ncol=5)#5 metrics
  corr_vec <- rep(0, length(k_uniq))
  par(mfrow=c(3,2), mar=c(2,2,1,1))
  plot(dat$Imputed_LHHSIZE, mod_pred, xlab = "Observed", 
       ylab = "Predicted", col=c('dark green','orange'),
       pch=c(16,16), cex.axis=1.5)
  abline(0,1)
  legend("topleft", c("Observed", "Predicted"), col=c("dark green", "orange"), pch=c(16,16),
         bty="n", cex=1.5)
  
  cvs <- list()
  for(i in 1:length(k_uniq))
  {
    #i=1
    print(paste0("fold_", i, sep=""))
    train_ind <- which(dat$k_fold!=k_uniq[i])
    dim(train <- dat[train_ind, ])#---train set for fold i
    dim(test <- dat[-train_ind, ]) #---test set for fold i
    
    train_coords <- coords[train_ind,]
    test_coords <- coords[-train_ind,]
    
    
    ###---Create projection matrices for training and testing datasets
    Ae<-inla.spde.make.A(mesh=mesh,loc=as.matrix(train_coords));dim(Ae) #training
    
    png(paste0(results_path,"/fold_",i,"_sample points for cross-validation.png"))
    plot(shp)
    points(train_coords, pch=16,  col="blue", cex=0.6)
    points(test_coords, pch=2,  col="red", cex=0.6)
    dev.off()
    
    
    #####------------------------
    covars_train <- train[,c("x2", "x3", "x17", "x21", "x32", "x34", "x40",
                             "x42", "set_reg", "set_typ", "region", "IDsr")]; dim(covars_train)
    stk_train <- inla.stack(data=list(y=train$dens_bldg), #the response
                            
                            A=list(Ae,1),  #the A matrix; the 1 is included to make the list(covariates)
                            
                            effects=list(c(list(Intercept=1), #the Intercept
                                           iset),  #the spatial index
                                         #the covariates
                                         list(covars_train)
                            ), 
                            tag='train')
    
    ####-----
    covars_test <- test[,c("x2", "x3", "x17", "x21", "x32", "x34", "x40",
                           "x42", "set_reg", "set_typ", "region", "IDsr")]; dim(covars_test)
    
    
    
    ###---Rerun INLA for model test prediction
    model <-inla(formula, #the formula
                 data=inla.stack.data(stk_train,spde=spde),  #the data stack
                 family= 'gamma',   #which family the data comes from
                 control.predictor=list(A=inla.stack.A(stk_train),compute=TRUE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
    summary(model)
    
    
    
    sfield_nodes_mean <- mod$summary.random$spatial.field['mean']
    field_mean <- (A%*% as.data.frame(sfield_nodes_mean)[, 1])
    
    sfield_nodesL <- mod$summary.random$spatial.field['0.025quant']
    fieldL <- (A%*% as.data.frame(sfield_nodesL)[, 1])
    
    sfield_nodesU<- mod$summary.random$spatial.field['0.975quant']
    fieldU <- (A%*% as.data.frame(sfield_nodesU)[, 1])
    
    
    
    ##--------
    fixed <-  
      #m1.samp[[i]]$latent[1,] +
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
      mod$summary.random$IDsr['mean'][-train_ind,1] + #--uncorrelated spatial random effects
      
      field_mean[-train_ind,1]
    
    dens_ht <- exp(fixed)
    sum(pop_ht <- dens_ht*test$Total_Building_Count)
    
    
    ###
    ######----Lower
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
      
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$`0.025quant`[2]) + #---settlement type and region nested effects
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$`0.025quant`[5]) + #---settlement type random effect
      mod$summary.random$IDsr['0.025quant'][-train_ind,1] + #--uncorrelated spatial random effects
      
      fieldL[-train_ind,1]
    
    dens_htL <- exp(fixedL)
    sum(pop_htL <- dens_htL*test$Total_Building_Count)
    
    
    #=========Upper
    fixedU <-   
      model$summary.fixed['Intercept', '0.975quant'] +
      model$summary.fixed['x2', '0.975quant'] * test[,'x2'] +
      model$summary.fixed['x3', '0.975quant'] * test[,'x3'] +
      model$summary.fixed['x17', '0.975quant'] * test[,'x17'] +
      model$summary.fixed['x21', '0.975quant'] * test[,'x21'] +
      model$summary.fixed['x32', '0.975quant'] * test[,'x32'] +
      model$summary.fixed['x34', '0.975quant'] * test[,'x34'] + 
      model$summary.fixed['x40', '0.975quant'] * test[,'x40'] +
      model$summary.fixed['x42', '0.975quant'] * test[,'x42'] +
      
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$`0.975quant`[2]) + #---settlement type and region nested effects
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$`0.975quant`[5]) + #---settlement type random effect
      mod$summary.random$IDsr['0.975quant'][-train_ind,1] + #--uncorrelated spatial random effects
      
      fieldU[-train_ind,1]
    
    dens_htU <- exp(fixedU)
    sum(pop_htU <- dens_htU*test$Total_Building_Count)
    
    
    
    plot(test$Imputed_LHHSIZE, pop_ht, xlab = "Observed", 
         ylab = "Predicted", col=c('dark green','orange'),
         pch=c(16,16), cex.axis=1.5)
    abline(0,1)
    legend("topleft", c("Observed", "Predicted"), col=c("dark green", "orange"), pch=c(16,16),
           bty="n", cex=1.5) 
    
    
    
    corr_vec[i] <- round(cor(test$Imputed_LHHSIZE, pop_ht),4)
    (met <- model_metrics(test$Imputed_LHHSIZE,  
                          pop_ht, pop_htU,  pop_htL))
    metrics_cv[i, ] <- as.vector(unlist(met)) 
    
    cvs[[i]] <- data.frame(test=test$Imputed_LHHSIZE, pred=pop_ht, fold=rep(k_uniq[i], length(as.vector(pop_ht))))
  }
  
  metrics_cv <- data.frame(Metrics = c("MAE", "RMSE", "Abias","IMPRECISION","%Accuracy","corr"),
                           Fold_1=c(metrics_cv[1,],corr_vec[1]),
                           Fold_2 =c(metrics_cv[2,],corr_vec[2]),
                           Fold_3 =c(metrics_cv[3,],corr_vec[3]),
                           Fold_4 =c(metrics_cv[4,],corr_vec[4]),
                           Fold_5=c(metrics_cv[5,],corr_vec[5]))
  metrics_cv$mean <- c(apply(metrics_cv[,2:6], 1, mean))
  return(list(metrics=metrics_cv, data = cvs))
}

#print(nm)
#number of folds of different test data sets
#modfit: The name of the best fit model
#predfit: Predicted values based on the best fit model
#mesh: The mesh built for the entire study location
#shpf: Shapefile of the EAs
#formfit: Formula of the best fit model
#spde: spde object 
#dat2: the input dataset from where the k-folds are created
#A: Projection matrix

k_folds=5 #----
(cv <- cross_val(dat2, modfit, predfit, mesh, spde,
                 A, shpf, formfit, k_folds))




