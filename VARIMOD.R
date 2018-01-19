
Rprof("d1.prof")


######### main R script for statistical modeling tool VARIMOD

#### Caution: Do not touch this R script in normal tool use!

###############################################################################
#####
##### reset R
#####
###############################################################################
#### reset R environment
rm( list=ls() );

###############################################################################
#####
##### init VARIMOD - NOTE: This is not done in a separate function since
#####                      VARIMOD functions are not loaded yet. Instead,
#####                      they are loaded as a part of the initialization.
#####
###############################################################################

#### init -- read file created by RunVarimod.tcsh
source(".varimod_init.R");

#### open logfile
sink(file=VarimodLogFile, append=FALSE, split=TRUE);
cat("VARIMOD -- Statistical device modeling tool\n")
cat("===========================================\n");

#### echo start date
time_start <- Sys.time();
cat("Start: ",date(),"\n", sep="");
cat("===========================================\n");

#### load required R packages
cat("Initializing\n")
cat("===========================================\n");
cat("  loading required R packages\n");
library("ellipse"); 
library("Matrix");
library("BB");
library("compiler") 
library("foreach",   lib.loc="/eas/projekte/varimod-2_180564/work/R_pkgs/")
library("iterators", lib.loc="/eas/projekte/varimod-2_180564/work/R_pkgs/")
library("doParallel",lib.loc="/eas/projekte/varimod-2_180564/work/R_pkgs/")

##################################################################################
#####  This package handles running much larger chunks of computations in parallel
##################################################################################
library("parallel")                     
##################################################################################
#####detectcores() tries to determine the number of CPU cores in the machine 
#####on which R is running 
##################################################################################

no_cores <- detectCores() -1;                   


#### load R functions
cat("  loading R functions\n");
## find all .R files in VARIMOD install dir
#VarimodInstallDir <- dirname(sys.frame(1)$ofile)
Rfiles <- list.files(path=VarimodInstallDir, full.names=TRUE, pattern=".R$");

cl<-makeCluster(2)
registerDoParallel(cl)
files <-foreach(i=seq_along(Rfiles)) %do%{                     
 ## source all .R files, but skip main script
  if( grepl(pattern="VARIMOD.R$", x=Rfiles[i]) ) {
    ## nothing to be done here
  } else {
    cat("    ",sub(pattern="^.*/", replacement="", x=Rfiles[i]),"\n", sep="");
    source(Rfiles[i]);
  }
}
stopCluster(cl)
###############################################################################
#####
##### load user settings from config file
#####
###############################################################################

## note: variable VarimodRunSourceFile (defined in .varimod_init.R)
##       contains the config file name
UserSettings <- LoadEvalUserSettings(local_source_file=VarimodRunSourceFile,
                                     SETTINGS=SETTINGS,
                                     FAB_DATA=FAB_DATA,
                                     MODELS=MODELS,
                                     SIMULATIONS=SIMULATIONS,
                                     CORNERS=CORNERS)
SETTINGS <- UserSettings$SETTINGS;

###############################################################################
#####
##### read and prepare input models
#####
###############################################################################

MODELS <- ParseAndPrepareModels(Models=UserSettings$MODELS);

###############################################################################
#####
##### init SKEW_PARAMETERS - list/reactiveValues to store skew parameter values
#####
###############################################################################

SKEW_PARAMETERS <- InitSkewParameters(Models=MODELS)

###############################################################################
#####
##### import measurement results
#####
###############################################################################

FAB_DATA <- ImportMeasurementResults( FabData=UserSettings$FAB_DATA )

###############################################################################
#####
##### align SIMULATIONS and FAB_DATA
#####   --> skip simulations where measurement data is not available
#####   --> skip measurement values where simulation setup is not available
#####
###############################################################################

dummy <- AlignMeasurementsSimulations( Simulations=UserSettings$SIMULATIONS ,
                                       FabData=FAB_DATA
                                     )
## extract revised list FAB_DATA and new list SIMULATIONS
FAB_DATA <- dummy$FAB_DATA;
SIMULATIONS <- dummy$SIMULATIONS;

###############################################################################
#####
##### prepare simulations: simulations, inits, data block
#####
###############################################################################

SIM_NETLISTS <- SimulationPreparation(PCM_sims=SIMULATIONS$PCM,
                                      MAT_sims=SIMULATIONS$MATCHING,
                                      SkewPars=SKEW_PARAMETERS$golden,
                                      Models=MODELS,
                                      sim_dir=SETTINGS$simulation_dir,
                                      param_init_file=SETTINGS$param_init_file,
                                      data_block_file_pcm=SETTINGS$data_block_file_PCM,
                                      data_block_file_mat=SETTINGS$data_block_file_mat,
                                      data_block_name=SETTINGS$data_block_name,
                                      DEBUG=SETTINGS$debug
                                     );

## note: in the simulation netlists and the init file, <<skew_parameter>>.nom,
##       <<skew_parameter>>.delta_global, and <<skew_parameter>>.delta_local 
##       are available to access skew parameter values

###############################################################################
#####
##### simulate golden model
#####
###############################################################################

cat("===========================================\n");
cat("Simulating golden model\n");
cat("===========================================\n");

#### create data block
## SKEW_PARAMETERS$golden needs to be written to SPICE .datablock syntax
## adapt skew parameter names by appending ".nom"

Var_Golden <- SKEW_PARAMETERS$golden;
colnames(Var_Golden) <- paste0(names(SKEW_PARAMETERS$purpose),".nom");
SKEW_PARAMETERS$golden <- CreateDataBlock(filename=paste(SETTINGS$simulation_dir,"/",SETTINGS$data_block_file_PCM, sep=""),
                                          dataname=SETTINGS$data_block_name,
                                          analysis="points",
                                          BASE_POINT=Var_Golden,
                                          DEBUG=SETTINGS$debug
                                         );
## export skew parameters to .csv for debugging purposes
write.csv(x=SKEW_PARAMETERS$golden, 
          file=paste(SETTINGS$results_dir,"/skew_golden.csv", sep=""), row.names=FALSE);

#### simulate
## run the simulation; RunHspice returns a dataframe of simulation results
## all PCM simulation results are stored in the list PCM
PCM <- list();
PCM$golden <- RunHspice(netlist=SIM_NETLISTS$PCM,
                        performance_params=names(SIMULATIONS$PCM),
                        hspice_command <- SETTINGS$HSPICE,
	                      ECHO=SETTINGS$HSPICE_ECHO,
                        DEBUG=SETTINGS$debug
                       );
write.csv(x=PCM$golden, 
          file=paste(SETTINGS$results_dir,"/pcm_golden.csv", sep=""), row.names=FALSE);

###############################################################################
#####
##### sensitivity analysis and BPV for golden model
#####
###############################################################################

## this section examines the golden models. furthermore, it allows a short cut for VARIMOD
## a backward propagation of variance (BPV) is performed at the golden models to determine
## a set of statistical models around the golden models
## the model standard deviations, golden model parameters, and skew parameter limits are
## combined to Cpk values that provide a quick insight into the quality of golden models

cat("===========================================\n");
cat("Checking golden model\n");
cat("===========================================\n");
#### run a sensitivity analysis
## create data block & parameter deflections - use <<skew_parameter>>.delta_global
Var_Sens <- 0*SKEW_PARAMETERS$golden; 
names(Var_Sens) <- paste0(names(SKEW_PARAMETERS$purpose),".delta_global");
SKEW_PARAMETERS$sens_golden <- CreateDataBlock(filename=paste(SETTINGS$simulation_dir,"/",SETTINGS$data_block_file_PCM, sep=""),
                                               dataname=SETTINGS$data_block_name,
                                               analysis="sensitivity",
                                               BASE_POINT=Var_Sens,
                                               param_ranges=as.numeric(apply(SKEW_PARAMETERS$limits,2,diff)),   ##Updated by Nishant 
                                               DEBUG=SETTINGS$debug
);
write.csv(x=SKEW_PARAMETERS$sens_golden, 
          file=paste(SETTINGS$results_dir,"/skew_sens_golden.csv", sep=""), row.names=FALSE);

## simulate
PCM$sens_golden <- RunHspice(netlist=SIM_NETLISTS$PCM,
                             performance_params=names(SIMULATIONS$PCM),
                             hspice_command <- SETTINGS$HSPICE,
                             ECHO=SETTINGS$HSPICE_ECHO,
                             DEBUG=SETTINGS$debug
                            );
write.csv(x=PCM$sens_golden, 
          file=paste(SETTINGS$results_dir,"/pcm_sens_golden.csv", sep=""), row.names=FALSE);

## determine sensitivity matrix
SENSITIVITY <- list();
SENSITIVITY$pcm_golden <- DetermineSensitivityMatrix(X=SKEW_PARAMETERS$sens_golden,
                                                     Y=PCM$sens_golden)
## export sensitivity matrix at golden model
TMP <- data.frame(SENSITIVITY$pcm_golden);
names(TMP) <- names(SKEW_PARAMETERS$limits);
write.csv(x=TMP, 
          file=paste(SETTINGS$results_dir,"/sensitivity_pcm_golden.csv", sep=""), row.names=FAB_DATA$pcm_parameters);

## BPV
BPV_golden <- BPV(CovY=cov(FAB_DATA$PCM_values),
                  S=as.matrix(SENSITIVITY$pcm_golden),
                  yweights=FAB_DATA$pcm_weights,
                  varnames=names(SKEW_PARAMETERS$golden),
                  optim_ctrl=SETTINGS$optim_ctrl_BPV,
                  DEBUG=SETTINGS$debug
                 );
## note: BPV_golden is a list with entries,
##       BPV_golden$sd contains the skew parameter standard deviations

## check skew parameters with respect to their limits: calculate Cpk values
CPK <- list();
CPK$golden <- CheckParamLimits(means=as.numeric(SKEW_PARAMETERS$golden),
                               sdevs=BPV_golden$sd,
                               mins=as.numeric(SKEW_PARAMETERS$limits[1,]),
                               maxs=as.numeric(SKEW_PARAMETERS$limits[2,]),
                               DEBUG=SETTINGS$debug
                              );
## print results
## this section in the logfile can be used to evaluate golden model quality
cat("  Cpk values for skew parameters at golden models\n");
cl<-makePSOCKcluster(2);
registerDoParallel(cl);
quality<-foreach( isk = 1:length(CPK$golden) ) %do% {       
  cat("    ",names(SKEW_PARAMETERS$purpose)[isk],": ",round(CPK$golden[isk],2),"\n", sep="");
}
stopCluster(cl);
#### allow short cut: did user only request to check the golden models?
## if yes (SETTINGS$check_golden_model_only==TRUE), stop tool, close logfile, and quit
if( SETTINGS$check_golden_model_only==TRUE ) {
  cat("===========================================\n");
  cat("End: "); print(date());
  time_stop <- Sys.time();
  runtime <- difftime(time_stop,time_start);
  runtime_h <- floor(as.numeric(runtime, "hours"));
  runtime_m <- floor(as.numeric(runtime, "mins")) - 60*runtime_h;
  runtime_s <- round(as.numeric(runtime, "secs")) - 60*runtime_m - 3600*runtime_h;
  cat("Runtime: ",runtime_h,"h ",runtime_m,"min ",runtime_s,"s\n", sep="");
  cat("===========================================\n");
  cat("Stopping since user requested checking golden models only.\n");
  sink();
  q(save='no');
} else {
  ## continue VARIMOD
  tmp <- "dummy";
}

###############################################################################
#####
##### model centering
#####
###############################################################################

## model centering is the step that adopts skew parameter values so that PCM 
## simulations yield the mean PCM measurement results (best possible)

cat("===========================================\n");
cat("Centering models\n");
cat("--> major VARIMOD step 1/4\n");
cat("===========================================\n");

y_target <- data.frame(matrix(apply(FAB_DATA$PCM_values,2,mean), nrow=1));       
names(y_target) <- names(FAB_DATA$PCM_values);
CenteringResults <- ModelCentering(x_start=SKEW_PARAMETERS$golden,
                                   x_use=SKEW_PARAMETERS$purpose[1,],
                                   x_lower=as.numeric(SKEW_PARAMETERS$center_limits[1,]),
                                   x_upper=as.numeric(SKEW_PARAMETERS$center_limits[2,]),
                                   y_target=y_target,
                                   y_weights=FAB_DATA$pcm_weights,
                                   optim_ctrl=SETTINGS$optim_ctrl_centering,
                                   hspice_command="hspice",
                                   hspice_echo=SETTINGS$HSPICE_ECHO,
                                   netlist=SIM_NETLISTS$PCM,
                                   data_block_file=paste(SETTINGS$simulation_dir,"/",SETTINGS$data_block_file_PCM, sep=""),
                                   data_block_name=SETTINGS$data_block_name,
                                   DEBUG=SETTINGS$debug
                                  );

## note: CenteringResults is a list containing
##         $x: centered skew parameters
##         $y: corresponding PCM simulation results

## store and export centered models and corresponding PCNM values
SKEW_PARAMETERS$center <- CenteringResults$x;
PCM$center <- CenteringResults$y;
write.csv(x=SKEW_PARAMETERS$center, 
          file=paste(SETTINGS$results_dir,"/skew_center.csv", sep=""), row.names=FALSE);
write.csv(x=PCM$center, 
          file=paste(SETTINGS$results_dir,"/pcm_center.csv", sep=""), row.names=FALSE);

#### this section provides some evaluation of model centering
#### this allows a further analysis of the golden models

## how much did the skew parameters changed during centering?
## print out parameter changes during centering in .
cat("  Parameter shifts during centering\n");
param_shift_abs <- rep(0,ncol(SKEW_PARAMETERS$center));
param_shift_rel <- param_shift_abs
for( i = 1:ncol(SKEW_PARAMETERS$center) ){
  sname <- names(SKEW_PARAMETERS$center)[i];
  golden_ind <- which(sname==names(SKEW_PARAMETERS$golden))
  
  if( length(golden_ind)==0 ) { stop("INTERNAL ERROR in model centering parameter shift summary") }
  
  SK_golden <- as.numeric(SKEW_PARAMETERS$golden[1,golden_ind]);
  
  param_shift_abs[i] <- SKEW_PARAMETERS$center[1,i]-SK_golden;
  param_shift_rel[i] <- 100*param_shift_abs[i]/SK_golden
  cat("    ",sub(pattern="\\.nom", replacement="", x=names(SKEW_PARAMETERS$center)[i]),": ",
      "from  ",
      SK_golden,
      "  to  ",
      SKEW_PARAMETERS$center[1,i],
      "  by  ",
      signif(param_shift_abs[i],SETTINGS$signif_digits)," (",
      round(param_shift_rel[i],1),"%)\n", sep="");
}

#### at this point, a second short cut is possible
#### VARIMOD can be stopped due to a user request (when SETTINGS$stop_after_centering==TRUE)
#### or when a user-defined maximum parameter shift is exceeded (SETTINGS$max_param_change_centering)

## stop here due to user request?
if( SETTINGS$stop_after_centering==TRUE ) {
  cat("===========================================\n");
  cat("End: "); print(date());
  time_stop <- Sys.time();
  runtime <- difftime(time_stop,time_start);
  runtime_h <- floor(as.numeric(runtime, "hours"));
  runtime_m <- floor(as.numeric(runtime, "mins")) - 60*runtime_h;
  runtime_s <- round(as.numeric(runtime, "secs")) - 60*runtime_m - 3600*runtime_h;
  cat("Runtime: ",runtime_h,"h ",runtime_m,"min ",runtime_s,"s\n", sep="");
  cat("===========================================\n");
  cat("Stopping after model centering due to user request.\n");
  sink();
  q(save='no');
} else {
  # continue
  tmp <- "dummy";
}

## stop here since param change was too large?
if( any(SETTINGS$max_param_change_centering<abs(param_shift_rel)) ) {
  ind <- which(SETTINGS$max_param_change_centering<abs(param_shift_rel));
  cat("===========================================\n");
  cat("End: "); print(date());
  time_stop <- Sys.time();
  runtime <- difftime(time_stop,time_start);
  runtime_h <- floor(as.numeric(runtime, "hours"));
  runtime_m <- floor(as.numeric(runtime, "mins")) - 60*runtime_h;
  runtime_s <- round(as.numeric(runtime, "secs")) - 60*runtime_m - 3600*runtime_h;
  cat("Runtime: ",runtime_h,"h ",runtime_m,"min ",runtime_s,"s\n", sep="");
  cat("===========================================\n");
  cat("Stopping after model centering since parameter changes exceed limit of ",
      SETTINGS$max_param_change_centering,".:\n");
  for( i in ind ) {
    cat("         ",names(SKEW_PARAMETERS$golden)[i],"\n", sep="");
  }
  sink();
  q(save='no');
}


#### update init file for variation parameters
## replace golden model parameters by centered parameters in init file for SPICE simulations
cat("===========================================\n");
cat("Updating init file for variation parameters\n")
cat("===========================================\n");
TMP <- SKEW_PARAMETERS$golden; names(TMP) <- sub(pattern=".nom", replacement="", x=names(SKEW_PARAMETERS$golden));
for( i in 1:ncol(TMP) ) {
  ind_center <- grep(pattern=names(TMP)[i], x=names(SKEW_PARAMETERS$center));
  if( length(ind_center)>0 ) {
    TMP[1,i] <- SKEW_PARAMETERS$center[1,ind_center];
  } else {
    ## do not change entry
  }
}
TMP <- CreateSkewParamsInit(output_filename=paste(SETTINGS$simulation_dir,"/",SETTINGS$param_init_file, sep=""),
                            skew_parameters=TMP,
                            DEBUG=SETTINGS$debug
);

Rprof(NULL)
d1prof = summaryRprof("d1.prof")
head(d1prof$by.self)


###############################################################################
#####
##### sensitivity analysis and BPV for PCM -- total variation
#####
###############################################################################

cat("===========================================\n");
cat("Sensitivity analysis and BPV for PCM - total variation\n");
cat("--> major VARIMOD step 2/4\n");
cat("===========================================\n");
#### create data block & parameter deflections


###############################################################################
#####
##### Parallelization of CentreDataBlock of PCM and Local Variation starts
#####
##############################################################################
cl <- makePSOCKcluster(2); ## Cluster defining number of cores created 

ind_global_stats <- which(SKEW_PARAMETERS$purpose[2,]!=0)
TMP <- data.frame(matrix(data=0, ncol=length(ind_global_stats), nrow=1));
colnames(TMP) <- paste0(names(SKEW_PARAMETERS$purpose)[ind_global_stats],".delta_global")
param_ranges <- as.numeric(parApply(cl,SKEW_PARAMETERS$limits[,ind_global_stats],2,diff));     ##Data Values totally parallelized in 2 cores


SKEW_PARAMETERS$sens_total_center <- CreateDataBlock(filename=paste(SETTINGS$simulation_dir,"/",SETTINGS$data_block_file_PCM, sep=""),
                                                     dataname=SETTINGS$data_block_name,
                                                     analysis="sensitivity",
                                                     BASE_POINT=TMP,
                                                     param_ranges=param_ranges,
                                                     DEBUG=SETTINGS$debug
);

write.csv(x=SKEW_PARAMETERS$sens_total_center, 
          file=paste(SETTINGS$results_dir,"/skew_sens_total_center.csv", sep=""), row.names=FALSE);
if( length(FAB_DATA$matching)>0 ) {
ind_local_stats <- which(SKEW_PARAMETERS$purpose[3,]!=0)
dummy_local <- data.frame(matrix(data=0, ncol=length(ind_local_stats), nrow=1));
colnames(dummy_local) <- paste0(names(SKEW_PARAMETERS$purpose)[ind_local_stats],".delta_local")
param_ranges <- 1e-8*as.numeric(parApply(cl,SKEW_PARAMETERS$limits[,ind_local_stats],2,diff));  ##Data Values totally parallelized in 2 cores

SKEW_PARAMETERS$sens_local_center <- CreateDataBlock(filename=paste(SETTINGS$simulation_dir,"/",SETTINGS$data_block_file_mat, sep=""),
                                                     dataname=SETTINGS$data_block_name,
                                                     analysis="sensitivity",
                                                     BASE_POINT=dummy_local,
                                                     param_ranges=param_ranges,
                                                     DEBUG=SETTINGS$debug
);
write.csv(x=SKEW_PARAMETERS$sens_local_center, 
          file=paste(SETTINGS$results_dir,"/skew_sens_local_center.csv", sep=""), row.names=FALSE);
stopCluster(cl);
}else {
  ## no matching simulations -- no local variations
  cat("===========================================\n");
  cat("No matching analysis required.\n");
  cat("--> skipping major VARIMOD step 3/4\n");
  cat("===========================================\n");
  BPV_local <- list();
  MAT <- list();
  MAT$sens <- data.frame();
}
###################################################################################################
###
### Cluster for parallely executing create datablock function ends
###
###################################################################################################
#### simulate
#### simulate

MY_SIMULATION_DATA <- list(
  PCM=data.frame(netlist=SIM_NETLISTS$PCM, 
                 performance_params=paste(names(SIMULATIONS$PCM), collapse=" "),
                 hspice_command <- SETTINGS$HSPICE,
                 ECHO=SETTINGS$HSPICE_ECHO,
                 DEBUG=SETTINGS$debug),
  MATCHING=data.frame(netlist=SIM_NETLISTS$Matching,
                      performance_params=paste((FAB_DATA$matching_parameters), collapse=" "),
                      hspice_command <- SETTINGS$HSPICE,
                      ECHO=SETTINGS$HSPICE_ECHO,
                      DEBUG=SETTINGS$debug)
)

cl <- makePSOCKcluster(2)
RESULTS2 <- parLapply(cl,MY_SIMULATION_DATA, RunHspiceWrapper)
stopCluster(cl);

PCM$sens_total_center<-RESULTS2$PCM;
MAT <- list();
MAT$sens<-RESULTS2$MATCHING
#PCM$sens_total_center <- RunHspice(netlist=SIM_NETLISTS$PCM,
#                                   performance_params=names(SIMULATIONS$PCM),
#                                   hspice_command <- SETTINGS$HSPICE,
#                                   ECHO=SETTINGS$HSPICE_ECHO,
#                                   DEBUG=SETTINGS$debug
#);
write.csv(x=PCM$sens_total_center, 
          file=paste(SETTINGS$results_dir,"/pcm_sens_total_center.csv", sep=""), row.names=FALSE);

#### determine sensitivity matrix
SENSITIVITY$pcm_center <- DetermineSensitivityMatrix(X=SKEW_PARAMETERS$sens_total_center,
                                                     Y=PCM$sens_total_center)

#### export sensitivity matrix at centered model
TMP <- data.frame(SENSITIVITY$pcm_center);
names(TMP) <- sub(pattern="\\.delta_global", replacement="", x=names(SKEW_PARAMETERS$global_stats));
write.csv(x=TMP, 
          file=paste(SETTINGS$results_dir,"/sensitivity_pcm_center.csv", sep=""), row.names=FAB_DATA$pcm_parameters);

#### BPV
BPV_center <- BPV(CovY=cov(FAB_DATA$PCM_values),
                  S=as.matrix(SENSITIVITY$pcm_center),
                  yweights=FAB_DATA$pcm_weights,
                  varnames=sub(pattern="global", replacement="total", x=names(SKEW_PARAMETERS$sens_total_center)),
                  optim_ctrl=SETTINGS$optim_ctrl_BPV,
                  DEBUG=SETTINGS$debug
);

#### check skew parameters with respect to their limits
means <- c(); mins <- c(); maxs <- c();
for( i in seq_along(BPV_center$sd) ) {
  pname <- sub(pattern="\\.delta_total", replacement="", x=BPV_center$names[i]);
  ind <- grep(pattern=paste0(pname,"\\."), x=names(SKEW_PARAMETERS$center));
  if( length(ind)>0 ) {
    means <- c( means , SKEW_PARAMETERS$center[1,ind]);
  } else {
    ind <- grep(pattern=paste0(pname,"\\."), x=names(SKEW_PARAMETERS$golden));
    means <- c( means , SKEW_PARAMETERS$golden[1,ind]);
  }
  ind <- which(pname==names(SKEW_PARAMETERS$limits));
  mins <- c( mins , SKEW_PARAMETERS$limits[1,ind] );
  maxs <- c( maxs , SKEW_PARAMETERS$limits[2,ind] );
}

CPK$center <- CheckParamLimits(means=means,
                               sdevs=BPV_center$sd,
                               mins=mins,
                               maxs=maxs,
                               DEBUG=SETTINGS$debug
);

#### print results
cat("  Cpk values for skew parameters at centered models\n");
for( isk in seq_along(CPK$center) ) {
  cat("    ",
      sub(pattern="\\.delta_.*$", replacement="", x=BPV_center$names[isk]),
      ": ",round(CPK$center[isk],2),"\n");
}

#### consider total variation as global variation
## note: this part may be subject to future toiol extensions
##       the basic problem is that this approach double-counts local variations which are 
##        present in PCM results
##       however, a separation of global and local variations is not straight forward...
BPV_global <- BPV_center;
BPV_global$names <- names(SKEW_PARAMETERS$sens_total_center);

#### export global skew parameter covariance matrix to .csv
TMP <- data.frame(BPV_global$cov); names(TMP) <- BPV_global$names;
write.csv(x=TMP, 
          file=paste(SETTINGS$results_dir,"/covariance_matrix_global.csv", sep=""), row.names=BPV_global$names);

###############################################################################
#####
##### sensitivity analysis and BPV for matching -- local variation
#####
###############################################################################
#cl <- makePSOCKcluster(2); 
#### is matching data available?
#### if yes, repeat sensitivity analysis and BPV with matching netlist to extract local variations
if( length(FAB_DATA$matching)>0 ) {
  cat("===========================================\n");
  cat("Sensitivity analysis and BPV for matching - local variation\n");
  cat("--> major VARIMOD step 3/4\n");
  cat("===========================================\n");
  
  #### create data block & parameter deflections
#  ind_local_stats <- which(SKEW_PARAMETERS$purpose[3,]!=0)
#  TMP <- data.frame(matrix(data=0, ncol=length(ind_local_stats), nrow=1));
#  colnames(TMP) <- paste0(names(SKEW_PARAMETERS$purpose)[ind_local_stats],".delta_local")
#  param_ranges <- 1e-8*as.numeric(parApply(cl,SKEW_PARAMETERS$limits[,ind_local_stats],2,diff));
#  stopCluster(cl);
#  SKEW_PARAMETERS$sens_local_center <- CreateDataBlock(filename=paste(SETTINGS$simulation_dir,"/",SETTINGS$data_block_file_mat, sep=""),
#                                                       dataname=SETTINGS$data_block_name,
#                                                       analysis="sensitivity",
#                                                       BASE_POINT=TMP,
#                                                       param_ranges=param_ranges,
#                                                       DEBUG=SETTINGS$debug
#  );
# write.csv(x=SKEW_PARAMETERS$sens_local_center, 
#            file=paste(SETTINGS$results_dir,"/skew_sens_local_center.csv", sep=""), row.names=FALSE);
  
  #### simulate
 # MAT <- list();
#  MAT$sens <- RunHspice(netlist=SIM_NETLISTS$Matching,
#                        performance_params=FAB_DATA$matching_parameters,
 #                       hspice_command <- SETTINGS$HSPICE,
 #                       ECHO=SETTINGS$HSPICE_ECHO,
 #                       DEBUG=SETTINGS$debug
#  );
  write.csv(x=MAT$sens, 
            file=paste(SETTINGS$results_dir,"/mat_sens.csv", sep=""), row.names=FALSE);
  
  #### determine sensitivity matrix
  SENSITIVITY$mat <- DetermineSensitivityMatrix(X=2*sqrt(2)*SKEW_PARAMETERS$sens_local_center,
                                                Y=MAT$sens)
  
  #### BPV
  ## BPV has to be run separately for each model since local variations in different models
  ## do not interact
  BPV_local <- list();
  for( imodel in seq_along(MODELS) ) {
    ## search for model-specific entries in SKEW_PARAMETERS$sens_local_center
    ## and FAB_DATA$matching_parameters
    ind_sk <- grep(pattern=paste(names(MODELS)[imodel],"..", sep=""), x=names(SKEW_PARAMETERS$sens_local_center));
    ind_ma <- grep(pattern=paste(names(MODELS)[imodel],"..", sep=""), x=FAB_DATA$matching_parameters);
    if( any(SENSITIVITY$mat[ind_ma,ind_sk]!=0) ) {
      #### potentially, a section determining a maximum local variation 
      #### will be required in the future
      ### find minimum size device for this model in matching simulations
      #ind_sim <- grep(pattern=tolower(paste(names(MODELS)[imodel],"..", sep="")), 
      #                x=tolower(names(SIMULATIONS$MATCHING)));
      #area <- Inf
      #for( ia in ind_sim ) {
      #  le <- SIMULATIONS$MATCHING[[ia]]$params$length;
      #  le <- sub(pattern="[Nn]", replacement="e-09", x=le);
      #  le <- sub(pattern="[Uu]", replacement="e-06", x=le);
      #  le <- as.numeric(le);
      #  wi <- SIMULATIONS$MATCHING[[ia]]$params$width;
      #  wi <- sub(pattern="[Nn]", replacement="e-09", x=wi);
      #  wi <- sub(pattern="[Uu]", replacement="e-06", x=wi);
      #  wi <- as.numeric(wi);
      #  area <- min( area , le*wi );
      #}
      ### derive limits for skew parameter standard deviations
      #var_x_upper <- diag(BPV_center$cov[ind_sk,ind_sk]) * area; 
      ### run BPV
      #TMP <- BPV(CovY=cov(FAB_DATA$MATCHING_values[,ind_ma], use="complete.obs"),
      #           S=as.matrix(SENSITIVITY$mat[ind_ma,ind_sk]/2),
      #           var_x_upper=var_x_upper,
      #           varnames=names(SKEW_PARAMETERS$sens_local_center)[ind_sk],
      #           DEBUG=SETTINGS$debug
      #);
      ## run BPV
      imodel_fab_data <- which(names(FAB_DATA$matching)==names(MODELS)[imodel])
      TMP <- BPV(CovY=cov(FAB_DATA$MATCHING_values[,ind_ma], use="complete.obs"),
                 S=as.matrix(SENSITIVITY$mat[ind_ma,ind_sk]),
                 yweights=FAB_DATA$matching[[imodel_fab_data]]$weights,
                 varnames=names(SKEW_PARAMETERS$sens_local_center)[ind_sk],
                 optim_ctrl=SETTINGS$optim_ctrl_BPV,
                 DEBUG=SETTINGS$debug
      );
      ## store BPV results in list BPV_local
      ## finally, there will be model-specific list elements
      ## BPV_local$<<model_name>> with components
      ##   $cor, $cov, $sd, $var, $names
      BPV_local <- lappend( l=BPV_local, name=names(MODELS)[imodel], value=TMP );
      
      ## export local skew parameter covariance matrix to .csv
      TMP2 <- data.frame(TMP); names(TMP2) <- TMP$names;
      write.csv(x=TMP2, 
                file=paste(SETTINGS$results_dir,"/covariance_matrix_local_",names(MODELS)[imodel],".csv", sep=""), 
                row.names=TMP$names);
      ## export local sensitivity matrix to .csv
      TMP <- data.frame(SENSITIVITY$mat[ind_ma,ind_sk]/2); 
      names(TMP) <- names(SKEW_PARAMETERS$sens_local_center)[ind_sk];
      write.csv(x=TMP, 
                file=paste(SETTINGS$results_dir,"/sensitivity_matching_",names(MODELS)[imodel],".csv", sep=""), 
                row.names=FAB_DATA$matching_parameters[ind_ma]);
    } else {
      # nothing to be done here
    }
  }
  
  #### print Pelgrom coefficients
  cat("  Pelgrom coefficients for skew parameters\n");
  for( imodel in seq_along(BPV_local) ) {
    model_name <- names(BPV_local)[imodel];
    model_ind <- which(model_name==names(MODELS));
    for( ip in seq_along(BPV_local[[imodel]]$sd) ) {
      param_name <- sub(pattern=".delta_local", replacement="", x=BPV_local[[imodel]]$names[ip]);
      param_name <- sub(pattern="^.*\\.\\.", replacement="", x=param_name);
      param_ind <- which(param_name==names(MODELS[[imodel]]$skew_parameters))
      cat("      sigma( ",
          model_name,"..",param_name,
          " ) = ",
          signif(BPV_local[[imodel]]$sd[ip],SETTINGS$signif_digits),
          " * ",
          MODELS[[imodel]]$skew_parameters[[param_ind]]$dep,
          "\n", sep="");
    }
  }
  
} else {
  ## no matching simulations -- no local variations
  cat("===========================================\n");
  cat("No matching analysis required.\n");
  cat("--> skipping major VARIMOD step 3/4\n");
  cat("===========================================\n");
  BPV_local <- list();
  MAT <- list();
  MAT$sens <- data.frame();
}
###############################################################################
#####
##### derivation of corners from global variations
#####
###############################################################################

#### corner derivation from global variations
#### two approaches are implemented
####  "univariate corners": all parameters are deflected by k standard deviations
####  "multivariate corners": parameters are deflected in user-defined directions such that
####                          the corresponding point on the k-sigma ellipsoid is found

cat("===========================================\n");
cat("Deriving corners from global variations\n");
cat("===========================================\n");

## note: list CORNERS contains corner definitions
##       CORNERS$<<corner_name>> with entries
##                     $mult -- sigma multiplier - shift skew params by how many standard deviations?
##                     $<<model_name>>$skew_parameter -- direction of parameter shift:
##                                                       (-1: reduce; 1: increase; 0: keep constant)

## prepare list for corner handling
CORNER_VALUES <- list();
## corners from univariate statistics
p_mins <- c();
p_maxs <- c();
for( i in seq_along(BPV_global$sd) ) {
  pname <- sub(pattern="\\.delta_.*$", replacement="", x=BPV_global$names[i])
  ind <- grep(pattern=pname, x=names(SKEW_PARAMETERS$center));
  if( length(ind)>0 ) {
    pc <- SKEW_PARAMETERS$center[1,ind]
  } else {
    ind <- grep(pattern=pname, x=names(SKEW_PARAMETERS$golden));
    pc <- SKEW_PARAMETERS$golden[1,ind]
  }
  ind <- grep(pattern=pname, x=names(SKEW_PARAMETERS$limits));
  p_mins <- c( p_mins , as.numeric(SKEW_PARAMETERS$limits[1,ind]-pc) );
  p_maxs <- c( p_maxs , as.numeric(SKEW_PARAMETERS$limits[2,ind]-pc) );
}
CORNER_VALUES$univariate <- DeriveSkewCorners(p_cor=BPV_global$cor,
                                              p_sds=BPV_global$sd,
                                              p_min=p_mins,
                                              p_max=p_maxs,
                                              param_names=BPV_global$names,
                                              CORNERS=UserSettings$CORNERS,
                                              type="univariate",
                                              DEBUG=SETTINGS$debug
                                             );
## create data block and simulate
SKEW_PARAMETERS$univ_corners <- CreateDataBlock(filename=paste(SETTINGS$simulation_dir,"/",SETTINGS$data_block_file_PCM, sep=""),
                                                dataname=SETTINGS$data_block_name,
                                                analysis="points",
                                                BASE_POINT=CORNER_VALUES$univariate,
                                                DEBUG=SETTINGS$debug
);
write.csv(x=SKEW_PARAMETERS$univ_corners, 
          file=paste(SETTINGS$results_dir,"/skew_univ_corners.csv", sep=""), row.names=FALSE);
PCM$univ_corners <- RunHspice(netlist=SIM_NETLISTS$PCM,
                              performance_params=names(SIMULATIONS$PCM),
                              hspice_command <- SETTINGS$HSPICE,
	                            ECHO=SETTINGS$HSPICE_ECHO,
                              DEBUG=SETTINGS$debug
);
rownames(x=PCM$univ_corners) <- rownames(SKEW_PARAMETERS$univ_corners);
write.csv(x=PCM$univ_corners, 
          file=paste(SETTINGS$results_dir,"/pcm_univ_corners.csv", sep=""));

#### corners from multiariate statistics
CORNER_VALUES$multivariate <- DeriveSkewCorners(p_cor=BPV_global$cor,
                                                p_sds=BPV_global$sd,
                                                p_min=p_mins,
                                                p_max=p_maxs,
                                                param_names=BPV_global$names,
                                                CORNERS=UserSettings$CORNERS,
                                                type="multivariate",
                                                DEBUG=SETTINGS$debug
);
## create data block and simulate
SKEW_PARAMETERS$mult_corners <- CreateDataBlock(filename=paste(SETTINGS$simulation_dir,"/",SETTINGS$data_block_file_PCM, sep=""),
                                                dataname=SETTINGS$data_block_name,
                                                analysis="points",
                                                BASE_POINT=CORNER_VALUES$multivariate,
                                                DEBUG=SETTINGS$debug
);
write.csv(x=SKEW_PARAMETERS$mult_corners, 
          file=paste(SETTINGS$results_dir,"/skew_mult_corners.csv", sep=""), row.names=FALSE);
PCM$mult_corners <- RunHspice(netlist=SIM_NETLISTS$PCM,
                              performance_params=names(SIMULATIONS$PCM),
                              hspice_command <- SETTINGS$HSPICE,
                  	          ECHO=SETTINGS$HSPICE_ECHO,
                              DEBUG=SETTINGS$debug
);
rownames(x=PCM$mult_corners) <- rownames(SKEW_PARAMETERS$mult_corners);
write.csv(x=PCM$mult_corners, 
          file=paste(SETTINGS$results_dir,"/pcm_mult_corners.csv", sep=""));

###############################################################################
#####
##### conversion of analysis results into HSPICE syntax
#####
###############################################################################

cat("===========================================\n");
cat("Exporting analysis results to HSPICE models\n");
cat("===========================================\n");

model_files <- HspiceExport( output_file_prefix=SETTINGS$model_output,
                             center_parameters=SKEW_PARAMETERS$center,
                             golden_parameters=SKEW_PARAMETERS$golden,
                             BPV_GLOBAL=BPV_global,
                             CORNERS=CORNER_VALUES,
                             BPV_LOCAL=BPV_local,
                             MODELS=MODELS,
                             significant_digits=SETTINGS$signif_digits,
                             relative_variation=SETTINGS$rel_variation,
                             DEBUG=SETTINGS$debug
);

###############################################################################
#####
##### HSPICE Monte Carlo Simulation to check modeling results
#####
###############################################################################

#### if user requested MC simulation (SETTINGS$RUN_MC=TRUE), analysis results 
#### are validated by MC simulations
if( SETTINGS$RUN_MC ) {
  ## run internal or external MC simulation, i.e. VARIMOD-internal sampling 
  ## or simulator-internal sampling?
  cat("===========================================\n");
  cat("Running Monte Carlo analysis\n")
  cat("--> major VARIMOD step 4/4\n");
  cat("===========================================\n");
  
  if( SETTINGS$run_internal_MC ) {
    cat("  Running internal Monte Carlo analysis\n");
    #### internal sampling
    ## sample sets of skew parameters from global statistics
    library("MASS");
    ## multi-variate Gaussian sampling
    TMP <- mvrnorm(n=SETTINGS$mc_samples, mu=rep(0,length(BPV_global$sd)), Sigma=BPV_global$cor)
    for( i in 1:ncol(TMP) ) {
      TMP[,i] <- TMP[,i] * BPV_global$sd[i];
    }
    ## adaption of parameter names
    TMP <- data.frame(TMP); names(TMP) <- BPV_global$names;
    #
    SKEW_PARAMETERS$mc <- TMP; 
    names(SKEW_PARAMETERS$mc) <- names(SKEW_PARAMETERS$limits);
    for( i in 1:ncol(TMP) ) {
      SKEW_PARAMETERS$mc[,i] <- TMP[,i] + SKEW_PARAMETERS$center[1,i];
    }
    ## create data block
    TMP <- CreateDataBlock(filename=paste(SETTINGS$simulation_dir,"/",SETTINGS$data_block_file_PCM, sep=""),
                                          dataname=SETTINGS$data_block_name,
                                          analysis="points",
                                          BASE_POINT=TMP,
                                          DEBUG=SETTINGS$debug
    );
    ## simulate
    PCM$mc <- RunHspice(netlist=SIM_NETLISTS$PCM,
                        performance_params=names(SIMULATIONS$PCM),
                        hspice_command <- SETTINGS$HSPICE,
                        ECHO=SETTINGS$HSPICE_ECHO,
                        DEBUG=SETTINGS$debug
    );
    ## export
    write.csv(x=SKEW_PARAMETERS$mc, 
              file=paste(SETTINGS$results_dir,"/skew_mc.csv", sep=""), row.names=FALSE);
    write.csv(x=PCM$mc, 
              file=paste(SETTINGS$results_dir,"/pcm_mc.csv", sep=""), row.names=FALSE);
    
    #### Matching
    cat("  Do not analyze matching in internal MC simulation.\n");
    MAT$mc <- data.frame();
  } else {
    cat("  Running external Monte Carlo SPICE simulation\n");
    ind_model_files <- grep(pattern="corner", x=model_files)
    #### PCM
    cat("    PCM\n");
    ## prepare netlist: re-format simulation netlists to support HSPICE MC simulations
    MC_NETLISTS <- list();
    MC_NETLISTS$PCM <- SimNetlist2Mc( input_netlist=SIM_NETLISTS$PCM,
                                      model_files_to_load=model_files[-ind_model_files],
                                      MODEL_PARAMS=MODELS,
                                      mc_sample_size=SETTINGS$mc_samples,
                                      data_block_file=SETTINGS$data_block_file_PCM,
                                      param_init_file=SETTINGS$param_init_file,
                                      DEBUG=SETTINGS$debug
                                    );
    ## simulate
#    PCM$mc <- RunHspice(netlist=MC_NETLISTS$PCM,                                                          #here
#                        performance_params=names(SIMULATIONS$PCM),
#                        hspice_command <- SETTINGS$HSPICE,
#                      	ECHO=SETTINGS$HSPICE_ECHO,
#                       DEBUG=SETTINGS$debug
#    );
   
    
    #### Matching
    if( length(FAB_DATA$matching)>0 ) {
      cat("    Matching\n");
      ## prepare netlist
      MC_NETLISTS$Matching <- SimNetlist2Mc( input_netlist=SIM_NETLISTS$Matching,
                                             model_files_to_load=model_files[-ind_model_files],
                                             MODEL_PARAMS=MODELS,
                                             mc_sample_size=SETTINGS$mc_samples,
                                             data_block_file=SETTINGS$data_block_file_matching,
                                             param_init_file=SETTINGS$param_init_file,
                                             DEBUG=SETTINGS$debug
                                           );
    
      
      ## simulate
 #     MAT$mc <- RunHspice(netlist=MC_NETLISTS$Matching,                                                   #here
 #                         performance_params=FAB_DATA$matching_parameters,
 #                         hspice_command <- SETTINGS$HSPICE,
 # 	                      ECHO=SETTINGS$HSPICE_ECHO,
 #                         DEBUG=SETTINGS$debug
 #    );
      
    } else {
      ## no matching data -- no matching MC simulation
      MAT$mc <- data.frame();
    }
  }

  MY_SIMULATION_DATA <- list( 
    PCM=data.frame(netlist=MC_NETLISTS$PCM,                 
                   performance_params=paste(names(SIMULATIONS$PCM), collapse=" "),                 
                   hspice_command <- SETTINGS$HSPICE,                 
                   ECHO=SETTINGS$HSPICE_ECHO,                 
                   DEBUG=SETTINGS$debug),
    MATCHING=data.frame(netlist=MC_NETLISTS$Matching,                     
                        performance_params=paste((FAB_DATA$matching_parameters), collapse=" "),                     
                        hspice_command <- SETTINGS$HSPICE,                     
                        ECHO=SETTINGS$HSPICE_ECHO,                     
                        DEBUG=SETTINGS$debug) 
  )
  
  
  cl <- makePSOCKcluster(2)
  RESULTS3 <- parLapply(cl,MY_SIMULATION_DATA, RunHspiceWrapper)
  stopCluster(cl);
  PCM$mc<-RESULTS3$PCM;
  MAT$mc<-RESULTS3$MATCHING
  write.csv(x=PCM$mc, 
            file=paste(SETTINGS$results_dir,"/pcm_mc.csv", sep=""), row.names=FALSE);
  write.csv(x=MAT$mc, 
            file=paste(SETTINGS$results_dir,"/mat_mc.csv", sep=""), row.names=FALSE);
  
  #### MC simulation results are compared with user-defined limits for PCM parameters
  #### if PCM limits are defined in FAB_DATA$PCM_min and/or FAB_DATA$PCM_max
  if( "PCM_min" %in% names(FAB_DATA) || "PCM_max" %in% names(FAB_DATA) ) {
    cat("  Checking PCM Monte Carlo simulation results against user-defined limits.\n");
    if( "PCM_min" %in% names(FAB_DATA) ) {
      PCM_min <- FAB_DATA$PCM_min;
      for( i in 2:SETTINGS$mc_samples ) {
        PCM_min <- rbind(PCM_min, FAB_DATA$PCM_min)
      }
      ind <- which(PCM_min>PCM$mc, arr.ind=TRUE);
      if( length(ind)>0 ) {
        # report simulation runs in which PCM limits are missed
        cat("    Minimum PCM limits missed in\n");
        for( i in 1:nrow(ind) ) {
          cat("      run ",ind[i,1]," for \"",FAB_DATA$pcm_parameters[ind[i,2]],"\"\n", sep="");
        }
        # remove values -- convert to NA
        PCM$mc[ind[i,1],] <- NA;
      } else {
        cat("    all lower PCM limits met.\n");
      }
    } else {
      cat("    lower PCM limits not defined.\n")
    }
    
    if( "PCM_max" %in% names(FAB_DATA) ) {
      PCM_max <- FAB_DATA$PCM_max;
      for( i in 2:SETTINGS$mc_samples ) {
        PCM_max <- rbind(PCM_max, FAB_DATA$PCM_max)
      }
      ind <- which(PCM_max<PCM$mc, arr.ind=TRUE);
      if( length(ind)>0 ) {
        # report
        cat("    Maximum PCM limits missed in\n");
        for( i in 1:nrow(ind) ) {
          cat("      run ",ind[i,1]," for \"",FAB_DATA$pcm_parameters[ind[i,2]],"\"\n", sep="");
        }
        # remove values
        PCM$mc[ind[i,1],] <- NA;
      } else {
        cat("    all upper PCM limits met.\n");
      }
    } else {
      cat("    upper PCM limits not defined.\n")
    }
  } else {
    # do not check PCM limits
    tmp <- "dummy";
  }
} else {
  cat("===========================================\n");
  cat("Skipping validation by MC simulations due to user request by \"SETTINGS$RUN_MC\".\n");
  cat("--> skipping major VARIMOD step 4/4\n");
  cat("===========================================\n");
}

###############################################################################
#####
##### Visualization of results
#####
###############################################################################

if( SETTINGS$RUN_MC ) {
  #### visualize VARIMOD results only if MC simulations were run
  cat("===========================================\n");
  cat("Visualizing results \n");
  cat("===========================================\n");
  ## PCM
  cat("  PCM:");
  for( i in seq_along(SETTINGS$plot_format) ) {          #Modified by Nishant from 1:length() to seq_along()
    cat(" ",SETTINGS$plot_format[i], sep="");
    TMP <- ScatterplotPCM(MEAS=FAB_DATA$PCM_values,
                          SIM=PCM$mc,
                          SIM_SETUP=SIMULATIONS$PCM,
                          CORNERS=PCM$univ_corners,
                          corners2show=SETTINGS$corners_to_plot,
                          CENTER=PCM$center,
                          GOLDEN=PCM$golden,
                          pcm_vs_pcm=SETTINGS$plot_pcm_vs_pcm,
                          file_format=SETTINGS$plot_format[i],
                          show.meas.points=SETTINGS$plot_meas_points,
                          show.meas.ellipses=SETTINGS$plot_meas_ellipse,
                          show.sim.points=SETTINGS$plot_sim_points,
                          show.sim.ellipses=SETTINGS$plot_sim_ellipse,
                          plot.cdf=SETTINGS$plot_cdf,
                          show.golden=SETTINGS$plot_golden_models,
                          output_dir=paste(SETTINGS$results_dir,"/pics", sep=""),
                          DEBUG=FALSE);
  }
  cat("\n");
  
  ## global skew parameter variations
  cat("  skew parameters:");
  for( i in seq_along(SETTINGS$plot_format) ) {          #Modified by Nishant from 1:length() to seq_along() 
    cat(" ",SETTINGS$plot_format[i], sep="");
    TMP <- VisualizeSkewsGlobal(CovX=BPV_global$cov, 
                                CovX_param_names=BPV_global$names, 
                                CENTER=SKEW_PARAMETERS$center, 
                                GOLDEN=SKEW_PARAMETERS$golden, 
                                CORNERS=SKEW_PARAMETERS$univ_corners, 
                                corners2show=SETTINGS$corners_to_plot, 
                                file_format=SETTINGS$plot_format[i], 
                                show.golden=SETTINGS$plot_golden_models, 
                                output_dir=paste(SETTINGS$results_dir,"/pics", sep=""),
                                DEBUG=FALSE)
  }
  cat("\n");
  
  if( SETTINGS$visualize_matching && any(grepl("MATCHING_values", names(FAB_DATA))) ) {
    cat("  matching:");
    ## matching
    if( any(grepl(pattern="MATCHING_values", x=names(FAB_DATA))) ) {
      ## visualization of matching is adapted to X-FAB matching measurements
      ## matching visualization is, therefore, not applicable to PCM example
      ## hence, matching visualization can be disabled by user (SETTINGS$visualize_matching=FALSE)
      for( i in seq_along(SETTINGS$plot_format) ) {       #Modified by Nishant from 1:length() to seq_along()
        cat(" ",SETTINGS$plot_format[i], sep="");
        tmp <- VisualizeMatching(MEAS=FAB_DATA$MATCHING_values, 
                                 SIM=MAT$mc, 
                                 SIM_SETTINGS=SIMULATIONS, 
                                 gate_overdrives=SETTINGS$gate_overdrives,
                                 file_format=SETTINGS$plot_format[i], 
                                 output_dir=paste(SETTINGS$results_dir,"/pics/", sep=""), 
                                 DEBUG=FALSE);
      }
      cat("\n");
    }
  } else {
    ## do not visualize matching results
    tmp <- "dummy";
  }
  
  ## create LaTeX files for documentation of VARIMOD results
  ## only if first plot format is NOT "svg" - since svg is not suitable to LaTeX
  if( SETTINGS$plot_format[1]!="svg" ) {
    tmp <- LatexExport( figdir=paste(SETTINGS$results_dir,"/pics", sep=""),
                        fig_format=SETTINGS$plot_format[1],
                        texdir=paste(SETTINGS$results_dir,"/tex", sep=""),
                        details_two_column=SETTINGS$latex_2_columns,
                        DEBUG=SETTINGS$debug
                      );
  }
} else {
  # do not visualize and export to LaTeX
  tmp <- "dummy";
}  

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               ###############################################################################
#####
##### that's it
#####
###############################################################################

cat("===========================================\n");
cat("Stop: ",date(),"\n", sep=""); 
time_stop <- Sys.time();
runtime <- difftime(time_stop,time_start);
runtime_h <- floor(as.numeric(runtime, "hours"));
runtime_m <- floor(as.numeric(runtime, "mins")) - 60*runtime_h;
runtime_s <- round(as.numeric(runtime, "secs")) - 60*runtime_m - 3600*runtime_h;
cat("Runtime: ",runtime_h,"h ",runtime_m,"min ",runtime_s,"s\n", sep="");
cat("===========================================\n");
cat("Done.\n\n");
sink();
