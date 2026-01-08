library(rtop)
library(deSolve)
library(rootSolve)
library(FME)
library(rsq)
library(foreach)
library(doParallel)
library(zoo)
library(minpack.lm) 

source("/media/DATADRIVE1/Model/Millennial/calibration_mm2/d02_script/functions/MM_v3.1.R")

ParaName    = c('param_pi'  ,'param_pa'  ,'kaff_pl' ,'alpha_pl','eact_pl' ,'rate_pa',
                'rate_break','rate_leach','kaff_des','param_p1','param_p2','kaff_lb',
                'alpha_lb'  ,'eact_lb'   ,'rate_bd' ,'rate_ma' ,'cue_ref' ,'cue_t'  ,
                'tae_ref'   ,'matpot'    ,'lambda'  ,'porosity','kamin'   ,'param_pb')

ParaNameNew = c('param_pi'  ,'param_pa'  ,'kaff_pl' ,'alpha_pl','eact_pl' ,'rate_pa',
                'rate_break','rate_leach','kaff_des','param_p1','param_p2','kaff_lb',
                'alpha_lb'  ,'eact_lb'   ,'rate_bd' ,'rate_ma' ,'cue_ref' ,'cue_t'  ,
                'tae_ref'   ,'matpot'    ,'lambda'  ,'porosity','kamin'   ,'param_pb', 'param_qmax','param_pH')

nUpper  = c(1,            1,      15000,    3.75e12,   70000,    0.85,
            0.1,     0.00225,      0.036,        0.5,   0.324,    1000,
            3.9e12,   75000,     0.0054,        0.1,     0.9,   0.018,
            22.5,        22.5,    3.15e-4,        0.9,       0.3,       1)


nLower  = c(0,            0,       5000,    1.25e12,   50000,     0.0,
            0,       0.00075,      0.012,       0.05,   0.108,     100,
            1.0e12,   50000,     0.0018,      0.001,     0.3,   0.006,
            7.5,         7.5,    1.05e-4,        0.3,     0.1,       0)

par0    = c(0.66,      0.4,      10000,     2.6e+12,   63339,    0.012,
            0.026,   0.0015,       0.02,      0.078,   0.216,    710.8,
            1.2e+12,  60428,     0.0044,     0.0052,    0.53,    0.012,
            15,          15,    0.00021,       0.62,     0.2,     0.52)


logor   = c(FALSE, FALSE,  TRUE,   TRUE,  TRUE,  FALSE,
            FALSE, FALSE,  FALSE,  FALSE, FALSE, TRUE,
            TRUE,  TRUE,   FALSE,  FALSE, FALSE, FALSE,
            FALSE, FALSE,  FALSE,  FALSE, FALSE, FALSE)

params_df <- data.frame(
  Parameter = ParaName,
  Lower = nLower,
  Initial = par0,
  Upper = nUpper,
  Logor = logor
)

important_params_set6 <- c('param_pa', 'rate_pa', 'rate_break', 'rate_leach',
                           'rate_ma',  'matpot',  'eact_pl',    'eact_lb',
                           'cue_ref',  'tae_ref', 'kaff_des',   'param_p1',
                           'param_pi')

selected_params <- params_df[params_df$Parameter %in% important_params_set6, ]


derivs.rangeland.MV3 <- function(step.num,state,params_sceua,forc_st,forc_sw,forc_npp) {
  
  MV3_params <- list(
    param_pi   = 0.66,
    param_pa   = 0.40,
    kaff_pl    = 10000,
    alpha_pl   = 2.6e+12,
    eact_pl    = 63339,
    rate_pa    = 0.012,
    rate_break = 0.026,
    rate_leach = 0.0015,
    kaff_des   = 0.02,
    param_p1   = 0.078,
    param_p2   = 0.216,
    kaff_lb    = 710.8,
    alpha_lb   = 1.2e+12,
    eact_lb    = 60428,
    rate_bd    = 0.0044,
    rate_ma    = 0.0052,
    cue_ref    = 0.53,
    cue_t      = 0.012,
    tae_ref    = 15,
    matpot     = 15,
    lambda     = 0.00021,
    porosity   = 0.62,
    kamin      = 0.2,
    param_pb   = 0.52,
    param_qmax = 6000,
    param_pH   = 6.2
  )
  
  parameters <- modifyList(MV3_params, params_sceua)
  
  with(as.list(c(state,parameters)), {
    
    #POM <- max(state[1], 0)
    #LMWC <- max(state[2], 0)
    #AGG <- max(state[3], 0)
    #MIC <- max(state[4], 0)
    #MAOM <- max(state[5], 0)
    # Soil type properties  
    #Equation 10
    kaff_lm = exp(-parameters$param_p1 * parameters$param_pH - parameters$param_p2) * parameters$kaff_des
    
    # Hydrological properties
    
    #Equation 4
    scalar_wd = (forc_sw(step.num) / parameters$porosity)^0.5
    
    #Equation 15
    scalar_wb = exp(parameters$lambda * -parameters$matpot) * (parameters$kamin + (1 - parameters$kamin) * ((parameters$porosity - forc_sw(step.num)) / parameters$porosity)^0.5) * scalar_wd
    
    # Decomposition
    
    gas_const <- 8.31446
    
    #Equation 3
    vmax_pl = parameters$alpha_pl * exp(-parameters$eact_pl / (gas_const * (forc_st(step.num) + 273.15)))
    
    #Equation 2
    # POM -> LMWC
    if(POM>0 && MIC>0){
      f_PO_LM = vmax_pl * scalar_wd * POM * MIC / (parameters$kaff_pl + MIC)
    }else{
      f_PO_LM=0
    }
    
    #Equation 5
    # POM -> AGG
    if(POM>0){
      f_PO_AG = parameters$rate_pa * scalar_wd * POM
    }else{
      f_PO_AG=0
    }
    
    #Equation 6
    # AGG -> MAOM + POM
    if(AGG>0){
      f_AG_break = parameters$rate_break * scalar_wd * AGG
    }else{
      f_AG_break=0
    }
    
    #Equation 8
    # LMWC -> out of system leaching
    if(LMWC>0){
      f_LM_leach = parameters$rate_leach * scalar_wd * LMWC
    }else{
      f_LM_leach=0
    }
    
    #Equation 9
    # LMWC -> MAOM
    if(LMWC>0 && MAOM>0){
      f_LM_MA = scalar_wd * kaff_lm * LMWC * (1 - MAOM / parameters$param_qmax)
    }else{
      f_LM_MA=0
    }
    
    #Equation 12
    # MAOM -> LMWC
    if(MAOM>0){
      f_MA_LM = parameters$kaff_des * MAOM / parameters$param_qmax
    }else{
      f_MA_LM=0
    }
    
    #Equation 4
    vmax_lb = parameters$alpha_lb * exp(-parameters$eact_lb / (gas_const * (forc_st(step.num) + 273.15)))
    
    #Equation 13
    # LMWC -> MIC
    if(LMWC>0 && MIC>0){
      f_LM_MB = vmax_lb * scalar_wb * MIC * LMWC / (parameters$kaff_lb + LMWC)
    }else{
      f_LM_MB=0
    }
    
    #Equation 16
    # MIC -> MAOM + LMWC
    if(MIC>0){
      f_MB_turn = parameters$rate_bd * MIC^2.0
    }else{
      f_MB_turn=0
    }
    
    #Equation 18
    # MAOM -> AGG
    if(MAOM>0){  
      f_MA_AG = parameters$rate_ma * scalar_wd * MAOM
    }else{
      f_MA_AG=0
    }
    
    #Equation 22
    # microbial growth flux, but is not used in mass balance
    
    #Equation 21
    # MIC -> atmosphere
    if(MIC>0 && LMWC>0){
      
      #f_MB_atm = f_LM_MB * (1 - (parameters$cue_ref - parameters$cue_t * (forc_st(step.num) - parameters$tae_ref) ) )
      #f_MB_gw  = f_LM_MB * (parameters$cue_ref - parameters$cue_t * (forc_st(step.num) - parameters$tae_ref))
      CUE_tmp = parameters$cue_ref - parameters$cue_t * (forc_st(step.num) - parameters$tae_ref)
      
      f_MB_atm = f_LM_MB * (1 - (CUE_tmp + CUE_corr) )
      f_MB_gw  = f_LM_MB * (CUE_tmp + CUE_corr)
      
      
      }else{
      f_MB_atm = 0
      f_MB_gw  = 0
    }
    
    # Update state variables
    
    #Equation 1
    dPOM = forc_npp(step.num) * parameters$param_pi + f_AG_break * parameters$param_pa - f_PO_AG - f_PO_LM
    
    #Equation 7
    dLMWC = forc_npp(step.num) * (1. - parameters$param_pi) - f_LM_leach + f_PO_LM - f_LM_MA - f_LM_MB + f_MB_turn * (1. - parameters$param_pb) + f_MA_LM
    
    #Equation 17
    dAGG = f_MA_AG + f_PO_AG - f_AG_break
    
    #Equation 20
    dMIC = f_LM_MB - f_MB_turn - f_MB_atm
    
    #Equation 19
    dMAOM = f_LM_MA - f_MA_LM + f_MB_turn * parameters$param_pb - f_MA_AG + f_AG_break * (1. - parameters$param_pa)
    
    #if(POM <= 1e-6 && dPOM < 0) {dPOM <- 0; POM <- 0.001}
    #if(LMWC <= 1e-6 && dLMWC < 0) {dLMWC <- 0; LMWC <- 0.001}
    #if(AGG <= 1e-6 && dAGG < 0) {dAGG <- 0; AGG <- 0.001}
    #if(MIC <= 1e-6 && dMIC < 0) {dMIC <- 0; MIC <- 0.001}
    #if(MAOM <= 1e-6 && dMAOM < 0) {dMAOM <- 0; MAOM <- 0.001}    # if(MAOM <= 1e-6 && dMAOM < 0) {dMAOM <- 0; MAOM <- 0.001}
    
    #outputs_delta <- c(dPOM, dLMWC, dAGG, dMIC, dMAOM)
    
    #if (any(abs(unlist(outputs_delta)) > 1000)){
    #  stop("something is wrong!")
    #}    # }
    
    # new_values <- state + outputs_delta
    # 
    # if (any(unlist(new_values) < 0)){
    #   #stop("something is wrong!")
    #   if(state[1] <=1e-6) state[1] <- 0
    #   if(state[2] <=1e-6) state[2] <- 0
    #   if(state[3] <=1e-6) state[3] <- 0
    #   if(state[4] <=1e-6) state[4] <- 0
    #   if(state[5] <=1e-6) state[5] <- 0
    # }
    
    f_mic_resp = f_MB_atm + f_MB_gw
     
    if (!is.na(f_mic_resp) && f_mic_resp > 0){
      rCUE = f_MB_gw/(f_MB_atm+f_MB_gw)
    } else {
      rCUE = 0
    }
    
    return(list(c(dPOM, dLMWC, dAGG, dMIC, dMAOM),'CUE'=rCUE))
    
  })
}

derivs.rangeland.MV3.adj <- function(step.num,state,params_sceua,forc_st,forc_sw,forc_npp) {
  
  MV3_params <- list(
    param_pi   = 0.66,
    param_pa   = 0.40,
    kaff_pl    = 10000,
    alpha_pl   = 2.6e+12,
    eact_pl    = 63339,
    rate_pa    = 0.012,
    rate_break = 0.026,
    rate_leach = 0.0015,
    kaff_des   = 0.02,
    param_p1   = 0.078,
    param_p2   = 0.216,
    kaff_lb    = 710.8,
    alpha_lb   = 1.2e+12,
    eact_lb    = 60428,
    rate_bd    = 0.0044,
    rate_ma    = 0.0052,
    cue_ref    = 0.53,
    cue_t      = 0.012,
    tae_ref    = 15,
    matpot     = 15,
    lambda     = 0.00021,
    porosity   = 0.62,
    kamin      = 0.2,
    param_pb   = 0.52,
    param_qmax = 6000,
    param_pH   = 6.2,
    param_corr = 0.2)
  
  parameters <- modifyList(MV3_params, params_sceua)
  
  with(as.list(c(state,parameters)), {
    
    #POM <- max(state[1], 0)
    #LMWC <- max(state[2], 0)
    #AGG <- max(state[3], 0)
    #MIC <- max(state[4], 0)
    #MAOM <- max(state[5], 0)
    # Soil type properties  
    #Equation 10
    kaff_lm = exp(-parameters$param_p1 * parameters$param_pH - parameters$param_p2) * parameters$kaff_des
    
    # Hydrological properties
    
    #Equation 4
    scalar_wd = (forc_sw(step.num) / parameters$porosity)^0.5
    
    #Equation 15
    scalar_wb = exp(parameters$lambda * -parameters$matpot) * (parameters$kamin + (1 - parameters$kamin) * ((parameters$porosity - forc_sw(step.num)) / parameters$porosity)^0.5) * scalar_wd
    
    # Decomposition
    
    gas_const <- 8.31446
    
    #Equation 3
    vmax_pl = parameters$alpha_pl * exp(-parameters$eact_pl / (gas_const * (forc_st(step.num) + 273.15)))
    
    #Equation 2
    # POM -> LMWC
    if(POM>0 && MIC>0){
      f_PO_LM = vmax_pl * scalar_wd * POM * MIC / (parameters$kaff_pl + MIC)
    }else{
      f_PO_LM=0
    }
    
    #Equation 5
    # POM -> AGG
    if(POM>0){
      f_PO_AG = parameters$rate_pa * scalar_wd * POM
    }else{
      f_PO_AG=0
    }
    
    #Equation 6
    # AGG -> MAOM + POM
    if(AGG>0){
      f_AG_break = parameters$rate_break * scalar_wd * AGG
    }else{
      f_AG_break=0
    }
    
    #Equation 8
    # LMWC -> out of system leaching
    if(LMWC>0){
      f_LM_leach = parameters$rate_leach * scalar_wd * LMWC
    }else{
      f_LM_leach=0
    }
    
    #Equation 9
    # LMWC -> MAOM
    if(LMWC>0 && MAOM>0){
      f_LM_MA = scalar_wd * kaff_lm * LMWC * (1 - MAOM / parameters$param_qmax)
    }else{
      f_LM_MA=0
    }
    
    #Equation 12
    # MAOM -> LMWC
    if(MAOM>0){
      f_MA_LM = parameters$kaff_des * MAOM / parameters$param_qmax
    }else{
      f_MA_LM=0
    }
    
    #Equation 4
    vmax_lb = parameters$alpha_lb * exp(-parameters$eact_lb / (gas_const * (forc_st(step.num) + 273.15)))
    
    #Equation 13
    # LMWC -> MIC
    if(LMWC>0 && MIC>0){
      f_LM_MB = vmax_lb * scalar_wb * MIC * LMWC / (parameters$kaff_lb + LMWC)
    }else{
      f_LM_MB=0
    }
    
    #Equation 16
    # MIC -> MAOM + LMWC
    if(MIC>0){
      f_MB_turn = parameters$rate_bd * MIC^2.0
    }else{
      f_MB_turn=0
    }
    
    #Equation 18
    # MAOM -> AGG
    if(MAOM>0){  
      f_MA_AG = parameters$rate_ma * scalar_wd * MAOM
    }else{
      f_MA_AG=0
    }
    
    #Equation 22
    # microbial growth flux, but is not used in mass balance
    
    #Equation 21
    # MIC -> atmosphere
    if(MIC>0 && LMWC>0){
      
      #f_MB_atm = f_LM_MB * (1 - (parameters$cue_ref - parameters$cue_t * (forc_st(step.num) - parameters$tae_ref) ) )
      #f_MB_gw  = f_LM_MB * (parameters$cue_ref - parameters$cue_t * (forc_st(step.num) - parameters$tae_ref))
      CUE_tmp = parameters$cue_ref - parameters$cue_t * (forc_st(step.num) - parameters$tae_ref)
      
      f_MB_atm = f_LM_MB * (1 - (CUE_tmp + parameters$param_corr) )
      f_MB_gw  = f_LM_MB * (CUE_tmp + parameters$param_corr)
      
      
    }else{
      f_MB_atm = 0
      f_MB_gw  = 0
    }
    
    # Update state variables
    
    #Equation 1
    dPOM = forc_npp(step.num) * parameters$param_pi + f_AG_break * parameters$param_pa - f_PO_AG - f_PO_LM
    
    #Equation 7
    dLMWC = forc_npp(step.num) * (1. - parameters$param_pi) - f_LM_leach + f_PO_LM - f_LM_MA - f_LM_MB + f_MB_turn * (1. - parameters$param_pb) + f_MA_LM
    
    #Equation 17
    dAGG = f_MA_AG + f_PO_AG - f_AG_break
    
    #Equation 20
    dMIC = f_LM_MB - f_MB_turn - f_MB_atm
    
    #Equation 19
    dMAOM = f_LM_MA - f_MA_LM + f_MB_turn * parameters$param_pb - f_MA_AG + f_AG_break * (1. - parameters$param_pa)
    
    #if(POM <= 1e-6 && dPOM < 0) {dPOM <- 0; POM <- 0.001}
    #if(LMWC <= 1e-6 && dLMWC < 0) {dLMWC <- 0; LMWC <- 0.001}
    #if(AGG <= 1e-6 && dAGG < 0) {dAGG <- 0; AGG <- 0.001}
    #if(MIC <= 1e-6 && dMIC < 0) {dMIC <- 0; MIC <- 0.001}
    #if(MAOM <= 1e-6 && dMAOM < 0) {dMAOM <- 0; MAOM <- 0.001}    # if(MAOM <= 1e-6 && dMAOM < 0) {dMAOM <- 0; MAOM <- 0.001}
    
    #outputs_delta <- c(dPOM, dLMWC, dAGG, dMIC, dMAOM)
    
    #if (any(abs(unlist(outputs_delta)) > 100)){
    #  stop("something is wrong!")
    #}    # }
    
    #new_values <- state + outputs_delta
    # 
    # if (any(unlist(new_values) < 0)){
    #   #stop("something is wrong!")
    #   if(state[1] <=1e-6) state[1] <- 0
    #   if(state[2] <=1e-6) state[2] <- 0
    #   if(state[3] <=1e-6) state[3] <- 0
    #   if(state[4] <=1e-6) state[4] <- 0
    #   if(state[5] <=1e-6) state[5] <- 0
    # }
    #if (any(unlist(new_values) < 0)) {
    #   stop("something is wrong!")
    #}
    
    f_mic_resp = f_MB_atm + f_MB_gw
    
    if (!is.na(f_mic_resp) && f_mic_resp > 0){
      rCUE = f_MB_gw/(f_MB_atm+f_MB_gw)
    } else {
      rCUE = 0
    }
    
    return(list(c(dPOM, dLMWC, dAGG, dMIC, dMAOM),'CUE'=rCUE))
    
  })
}

derivs_Aus_wrapper <- function(step.num,state,params_sceua,forc_st,forc_sw,forc_npp) {
  
  #adjusted_state <- state
  #adjusted_state[adjusted_state <= 1e-6] <- 0.001
  
  #output <- derivs.rangeland.MV3.adj(step.num,adjusted_state,params_sceua,forc_st,forc_sw,forc_npp)
  output <- derivs.rangeland.MV3.adj(step.num,state,params_sceua,forc_st,forc_sw,forc_npp)
  
  return(list(output[[1]],output[[2]]))
  
}

derivs_Aus_wrapper_adj <- function(step.num,state,params_sceua,forc_st,forc_sw,forc_npp) {
  #output <- derivs.rangeland.MV3(step.num,state,params_sceua,forc_st,forc_sw,forc_npp)
  #new_state <- state + output[[1]] * 1
  #new_state[new_state <= 1e-6] <- 0.001
  #return(list(output[[1]],output[[2]]))
  
    # Check if any state values are too small and adjust them
    adjusted_state <- state
    adjusted_state[adjusted_state <= 1e-6] <- 0.001
    
    # Use the adjusted state for calculations
    output <- derivs.rangeland.MV3(step.num, adjusted_state, params_sceua, forc_st, forc_sw, forc_npp)
    
    # Get derivatives
    derivs <- output[[1]]
    
    # Prevent derivatives from making small values go negative
    #for(i in 1:length(state)) {
    #  if(state[i] < 0.01 && derivs[i] < 0) {  # Wider threshold
    #    derivs[i] <- derivs[i] * 0.1  # Reduce negative derivative rather than zero it
    #  }
    #  if(state[i] <= 1e-6 && derivs[i] < 0) {
    #    derivs[i] <- 0
    #  }
    #}    # }
    
    return(list(derivs, output[[2]]))
  
}

SStime     <- 36500
run.steps  <- seq(1,36500)

compile_muresk_csat <- read.csv("/media/DATADRIVE1/Project/muresk/model/d02_data/for_model/compile_muresk_csat.txt")
compile_muresk_csat <- read.csv("/media/DATADRIVE1/Project/muresk/model/d02_data/for_model/compile_muresk_csat_all.txt")
compile_muresk_csat <- compile_muresk_csat[,c(1:2,4:19)]
compile_muresk_csat$pocmac030 <- with(compile_muresk_csat,POCmac_30/100*BD_30*30*10000)
compile_muresk_csat$pocmic030 <- with(compile_muresk_csat,POCmic_30/100*BD_30*30*10000)
compile_muresk_csat$maoc030   <- with(compile_muresk_csat,MAOC_30/100*BD_30*30*10000)
compile_muresk_csat$toc030    <- with(compile_muresk_csat,toc_30/100*BD_30*30*10000)
compile_muresk_csat$qmax030   <- with(compile_muresk_csat,Csat/100*BD_30*30*10000)
compile_muresk_csat <- compile_muresk_csat[,c(2:20)]
compile_muresk_csat <- compile_muresk_csat[,c(1,15,16,17,18,13,19,8)]
colnames(compile_muresk_csat) <- c("site","POM","AGG","MAOM","SOM","pH","qmax","BD")


senst_ode_100 <- function(params,sitename,dat,pH_sens) {
  
  obs_row  <- compile_muresk_csat[which(compile_muresk_csat$site == sitename),]
  
  initial_values <- c(POM = obs_row$POM, LMWC = 6, AGG = obs_row$AGG, MIC = 6, MAOM=obs_row$MAOM)  
  
  #100years
  forc_st_fun  <- approxfun(1:SStime, rep(dat$forc_st, 100))
  forc_sw_fun  <- approxfun(1:SStime, rep(dat$forc_sw, 100))
  forc_npp_fun <- approxfun(1:SStime, rep(dat$forc_npp,100))
  
  param_qmax  <- obs_row$qmax
  param_pH    <- pH_sens
  
  params_sceua <- c(as.list(setNames(params, selected_params$Parameter)),
                    param_qmax = param_qmax,
                    param_pH = param_pH)
  
  
  ODE.MMRMM <- ode(y = initial_values, times=run.steps,  func = derivs_Aus_wrapper, parms = params_sceua, forc_st=forc_st_fun, forc_sw=forc_sw_fun, forc_npp=forc_npp_fun, method="rk4")
  #ODE.MMRMM <- forcing_delta_MM(initial_values,params_sceua,forcing_dat)
  
  ODE.MMRMM <- as.data.frame(ODE.MMRMM)
  # get the mean of equlibrium state, suggesting 120 rows
  ODE.MMRMM$SOM <- with(ODE.MMRMM,POM+LMWC+AGG+MIC+MAOM)
  
  #ODE.MMRMM_tail <- tail(ODE.MMRMM, 365)
  
  #row.names(ODE.MMRMM_tail) <- dat$doy
  
  #ODE.MMRMM_tail$pH <- pH_sens
  sim_update <- rollapply(ODE.MMRMM, width = 365, by = 365, FUN = mean)
  sim_update <- as.data.frame(sim_update)
  
  names(sim_update)[7] <- c("CUE")
  
  POM10000       <- mean(tail(sim_update$POM, 1))
  LMWC10000      <- mean(tail(sim_update$LMWC, 1))
  AGG10000       <- mean(tail(sim_update$AGG, 1))
  MIC10000       <- mean(tail(sim_update$MIC, 1))
  MAOM10000      <- mean(tail(sim_update$MAOM,1))
  SOM10000       <- mean(tail(sim_update$SOM,1))
  CUE10000       <- mean(tail(sim_update$CUE,1))
  
  #POM_dev     <- abs(POM10000 - obs_row$POM)
  #AGG_dev     <- abs(AGG10000 - obs_row$AGG)
  #MAOM_dev    <- abs(MAOM10000 - obs_row$MAOM)
  #SOM_dev     <- abs(SOM10000  - obs_row$SOM)
  
  return(c(SOM10000, MAOM10000, POM10000,AGG10000,MIC10000, LMWC10000, CUE10000, pH_sens))
}

safe_ode_100_fn <- function(params,sitename,dat,pH_sens) {
  J_result <- tryCatch({
    # Replace 'actual_objective_function' with your actual function
    senst_ode_100(params,sitename,dat,pH_sens)
  }, error = function(e) {
    # Handle the error, e.g., by returning a large penalty value or NA
    return(Inf)  # or NA, depending on how sceua() handles it
  })
  
  return(J_result)
}

senst_ode_100 <- function(params,sitename,dat,qmax_sens) {
  
  obs_row  <- compile_muresk_csat[which(compile_muresk_csat$site == sitename),]
  
  initial_values <- c(POM = obs_row$POM, LMWC = 6, AGG = obs_row$AGG, MIC = 6, MAOM=obs_row$MAOM)  
  
  #100years
  forc_st_fun  <- approxfun(1:SStime, rep(dat$forc_st, 100))
  forc_sw_fun  <- approxfun(1:SStime, rep(dat$forc_sw, 100))
  forc_npp_fun <- approxfun(1:SStime, rep(dat$forc_npp,100))
  
  param_qmax  <- qmax_sens/100*obs_row$BD*30*10000
  param_pH    <- obs_row$pH
  
  params_sceua <- c(as.list(setNames(params, selected_params$Parameter)),
                    param_qmax = param_qmax,
                    param_pH = param_pH)
  
  
  ODE.MMRMM <- ode(y = initial_values, times=run.steps,  func = derivs_Aus_wrapper, parms = params_sceua, forc_st=forc_st_fun, forc_sw=forc_sw_fun, forc_npp=forc_npp_fun, method="rk4")
  #ODE.MMRMM <- forcing_delta_MM(initial_values,params_sceua,forcing_dat)
  
  ODE.MMRMM <- as.data.frame(ODE.MMRMM)
  # get the mean of equlibrium state, suggesting 120 rows
  ODE.MMRMM$SOM <- with(ODE.MMRMM,POM+LMWC+AGG+MIC+MAOM)
  
  #ODE.MMRMM_tail <- tail(ODE.MMRMM, 365)
  
  #row.names(ODE.MMRMM_tail) <- dat$doy
  
  #ODE.MMRMM_tail$pH <- pH_sens
  sim_update <- rollapply(ODE.MMRMM, width = 365, by = 365, FUN = mean)
  sim_update <- as.data.frame(sim_update)
  
  names(sim_update)[7] <- c("CUE")
  
  POM10000       <- mean(tail(sim_update$POM, 1))
  LMWC10000      <- mean(tail(sim_update$LMWC, 1))
  AGG10000       <- mean(tail(sim_update$AGG, 1))
  MIC10000       <- mean(tail(sim_update$MIC, 1))
  MAOM10000      <- mean(tail(sim_update$MAOM,1))
  SOM10000       <- mean(tail(sim_update$SOM,1))
  CUE10000       <- mean(tail(sim_update$CUE,1))
  
  #POM_dev     <- abs(POM10000 - obs_row$POM)
  #AGG_dev     <- abs(AGG10000 - obs_row$AGG)
  #MAOM_dev    <- abs(MAOM10000 - obs_row$MAOM)
  #SOM_dev     <- abs(SOM10000  - obs_row$SOM)
  
  return(c(SOM10000, MAOM10000, POM10000,AGG10000,MIC10000, LMWC10000, CUE10000, param_qmax))
}

safe_ode_100_fn <- function(params,sitename,dat,qmax_sens) {
  J_result <- tryCatch({
    # Replace 'actual_objective_function' with your actual function
    senst_ode_100(params,sitename,dat,qmax_sens)
  }, error = function(e) {
    # Handle the error, e.g., by returning a large penalty value or NA
    return(Inf)  # or NA, depending on how sceua() handles it
  })
  
  return(J_result)
}


# Event function to detect when states go negative
eventfun <- function(t, y, parms, ...) {
  return(y - 1e-6)  # Trigger when any state approaches 1e-6
}

# Root function to reset states
rootfun <- function(t, y, parms, ...) {
  y[y < 1e-6] <- 0.001
  return(y)
}


senst_ode_cue_100 <- function(params, sitename, dat, CUE_corr) {
  
  obs_row  <- compile_muresk_csat[which(compile_muresk_csat$site == sitename),]
  
  initial_values <- c(POM = obs_row$POM, LMWC = 6, AGG = obs_row$AGG, MIC = 6, MAOM=obs_row$MAOM)  
  
  #100years
  forc_st_fun  <- approxfun(1:SStime, rep(dat$forc_st, 100))
  forc_sw_fun  <- approxfun(1:SStime, rep(dat$forc_sw, 100))
  forc_npp_fun <- approxfun(1:SStime, rep(dat$forc_npp,100))
  
  param_qmax  <- obs_row$qmax
  param_pH    <- obs_row$pH
  param_corr  <- CUE_corr
  
  params_sceua <- c(as.list(setNames(params, selected_params$Parameter)),
                    param_qmax = param_qmax,
                    param_pH   = param_pH,
                    param_corr = CUE_corr)
  
  
  ODE.MMRMM <- ode(y = initial_values, times=run.steps,  func = derivs_Aus_wrapper, parms = params_sceua, forc_st=forc_st_fun, forc_sw=forc_sw_fun, forc_npp=forc_npp_fun, method="rk4")
  #ODE.MMRMM <- forcing_delta_MM(initial_values,params_sceua,forcing_dat)
  
  ODE.MMRMM <- as.data.frame(ODE.MMRMM)
  # get the mean of equlibrium state, suggesting 120 rows
  ODE.MMRMM$SOM <- with(ODE.MMRMM,POM+LMWC+AGG+MIC+MAOM)
  
  #ODE.MMRMM_tail <- tail(ODE.MMRMM, 365)
  
  #row.names(ODE.MMRMM_tail) <- dat$doy
  
  #ODE.MMRMM_tail$pH <- pH_sens
  sim_update <- rollapply(ODE.MMRMM, width = 365, by = 365, FUN = mean)
  sim_update <- as.data.frame(sim_update)
  
  names(sim_update)[7] <- c("CUE")
  
  POM10000       <- mean(tail(sim_update$POM, 1))
  LMWC10000      <- mean(tail(sim_update$LMWC, 1))
  AGG10000       <- mean(tail(sim_update$AGG, 1))
  MIC10000       <- mean(tail(sim_update$MIC, 1))
  MAOM10000      <- mean(tail(sim_update$MAOM,1))
  SOM10000       <- mean(tail(sim_update$SOM,1))
  CUE10000       <- mean(tail(sim_update$CUE,1))
  
  #POM_dev     <- abs(POM10000 - obs_row$POM)
  #AGG_dev     <- abs(AGG10000 - obs_row$AGG)
  #MAOM_dev    <- abs(MAOM10000 - obs_row$MAOM)
  #SOM_dev     <- abs(SOM10000  - obs_row$SOM)
  
  return(c(SOM10000, MAOM10000, POM10000,AGG10000,MIC10000, LMWC10000, CUE10000))
}

safe_ode_100_cue_fn <- function(params, sitename, dat, CUE_corr) {
  J_result <- tryCatch({
    # Replace 'actual_objective_function' with your actual function
    senst_ode_cue_100(params,sitename,dat,CUE_corr)
  }, error = function(e) {
    # Handle the error, e.g., by returning a large penalty value or NA
    return(Inf)  # or NA, depending on how sceua() handles it
  })
  
  return(J_result)
}

senmix_ode_100 <- function(params,sitename,dat,qmax_sens,pH_sens) {
  
  obs_row  <- compile_muresk_csat[which(compile_muresk_csat$site == sitename),]
  
  initial_values <- c(POM = obs_row$POM, LMWC = 6, AGG = obs_row$AGG, MIC = 6, MAOM=obs_row$MAOM)  
  
  #100years
  forc_st_fun  <- approxfun(1:SStime, rep(dat$forc_st, 100))
  forc_sw_fun  <- approxfun(1:SStime, rep(dat$forc_sw, 100))
  forc_npp_fun <- approxfun(1:SStime, rep(dat$forc_npp,100))
  
  param_qmax  <- qmax_sens/100*obs_row$BD*30*10000
  param_pH    <- pH_sens
  
  params_sceua <- c(as.list(setNames(params, selected_params$Parameter)),
                    param_qmax = param_qmax,
                    param_pH = param_pH)
  
  
  ODE.MMRMM <- ode(y = initial_values, times=run.steps,  func = derivs_Aus_wrapper, parms = params_sceua, forc_st=forc_st_fun, forc_sw=forc_sw_fun, forc_npp=forc_npp_fun, method="rk4")
  #ODE.MMRMM <- forcing_delta_MM(initial_values,params_sceua,forcing_dat)
  
  ODE.MMRMM <- as.data.frame(ODE.MMRMM)
  # get the mean of equlibrium state, suggesting 120 rows
  ODE.MMRMM$SOM <- with(ODE.MMRMM,POM+LMWC+AGG+MIC+MAOM)
  
  #ODE.MMRMM_tail <- tail(ODE.MMRMM, 365)
  
  #row.names(ODE.MMRMM_tail) <- dat$doy
  
  #ODE.MMRMM_tail$pH <- pH_sens
  sim_update <- rollapply(ODE.MMRMM, width = 365, by = 365, FUN = mean)
  sim_update <- as.data.frame(sim_update)
  
  names(sim_update)[7] <- c("CUE")
  
  POM10000       <- mean(tail(sim_update$POM, 1))
  LMWC10000      <- mean(tail(sim_update$LMWC, 1))
  AGG10000       <- mean(tail(sim_update$AGG, 1))
  MIC10000       <- mean(tail(sim_update$MIC, 1))
  MAOM10000      <- mean(tail(sim_update$MAOM,1))
  SOM10000       <- mean(tail(sim_update$SOM,1))
  CUE10000       <- mean(tail(sim_update$CUE,1))
  
  #POM_dev     <- abs(POM10000 - obs_row$POM)
  #AGG_dev     <- abs(AGG10000 - obs_row$AGG)
  #MAOM_dev    <- abs(MAOM10000 - obs_row$MAOM)
  #SOM_dev     <- abs(SOM10000  - obs_row$SOM)
  
  return(c(SOM10000, MAOM10000, POM10000,AGG10000,MIC10000, LMWC10000, CUE10000, param_qmax, param_pH))
}

safe_ode_mix_fn <- function(params,sitename,dat,qmax_sens,pH_sens) {
  J_result <- tryCatch({
    # Replace 'actual_objective_function' with your actual function
    senmix_ode_100(params,sitename,dat,qmax_sens,pH_sens)
  }, error = function(e) {
    # Handle the error, e.g., by returning a large penalty value or NA
    return(Inf)  # or NA, depending on how sceua() handles it
  })
  
  return(J_result)
}


set.seed(123)

flist <- list.files("/media/DATADRIVE1/Data/muresk/data/d02_data/muresk_daily_driving/",full.names = TRUE,recursive = TRUE)
flist <- flist[!flist %in% c("/media/DATADRIVE1/Data/muresk/data/d02_data/muresk_daily_driving//5644.csv")]
flist <- flist[!flist %in% c("/media/DATADRIVE1/Data/muresk/data/d02_data/muresk_daily_driving//5636.csv")]
#flist <- flist[!flist %in% c("/media/DATADRIVE1/Data/muresk/data/d02_data/muresk_daily_driving//5610.csv")]
pars_set <- readRDS("/media/DATADRIVE1/Data/muresk/data/d04_output/sitebysite_benchmark_muresk_ode_100y_0515_init.rds")
pars_set <- pars_set[,1:13]
###########################################################
#ode version pH sensitivity
cl <- parallel::makeCluster(15)
doParallel::registerDoParallel(cl)


sim_col <- foreach(
  i            = seq_along(flist),
  .packages    = c("doParallel","deSolve","zoo")
  # if safe_ode_100_fn or pars_set live in your global env:
  # for reproducible *and distinct* streams, uncomment:
  # .options.RNG = 123,
  ) %dopar% {
  
  ## (1) No set.seed() here, or if you need it,
  ##     use doRNG/.options.RNG so each i gets its own stream.
  set.seed(123)
  
  driving_input <- read.csv(flist[i])
  dates <- as.Date(with(driving_input, ISOdate(year,month,day)), tz="UTC")
  doy   <- as.POSIXlt(dates)$yday + 1
  keep  <- format(dates, "%m-%d") != "02-29"
  agg_NPP   <- tapply(driving_input$Npp[keep],         doy[keep], mean, na.rm=TRUE)
  agg_moist <- tapply(driving_input$SoilMoisture[keep],doy[keep], mean, na.rm=TRUE)
  agg_temp  <- tapply(driving_input$SoilTemp[keep],    doy[keep], mean, na.rm=TRUE)
  
  agg <- data.frame(
    doy      = as.integer(names(agg_NPP)),
    forc_npp = as.numeric(agg_NPP),
    forc_sw  = as.numeric(agg_moist),
    forc_st  = as.numeric(agg_temp)
  )
  agg$avg_date <- as.Date(agg$doy, origin = as.Date("2000-12-31"))
  agg <- agg[1:365, ]
  
  agg$forc_npp <- with(agg,forc_npp*1.5)
  sitename      <- sub("\\.csv$", "", basename(flist[i]))
  selected_row  <- pars_set[sitename, ]
  
  ## (2) Build up a *local* data.frame with <-, not <<-
  sim_output_init <- NULL
  for (ph in seq(4, 6.5, 0.5)) {
    out_pH <- safe_ode_100_fn(selected_row, sitename, agg, ph)
    sim_output_init <- rbind(sim_output_init, out_pH)
  }
  
  ## scale & clean
  sim_output_init <- na.omit(sim_output_init)
  sens <- scale(sim_output_init[,1:6])
  sens <- as.data.frame(sens)
  sens$pH  <- seq(4, 6.5, 0.5)
  sens$CUE <- sim_output_init[,7]
  colnames(sens)[1:6] <- c('SOM1000','MAOM1000','POM1000',
                           'AGG1000','MIC1000','LMWC1000')
  ## add a column so you know which site this came from
  sens$site <- sitename
  sens$NPP_ratio  <- 1.5
  sens
}

stopCluster(cl)

new_sim_outs <- do.call(rbind, sim_col)

new_sim_outs_npp15 <- new_sim_outs

rm(new_sim_outs_npp05)

colnames(new_sim_outs) <- c('SOM1000', 'MAOM1000', 'POM1000', 'AGG1000', 'MIC1000', 'LMWC1000',"pH")

new_sim_outs <- as.data.frame(new_sim_outs)

saveRDS(new_sim_outs_npp15,"/media/DATADRIVE1/Data/muresk/data/d04_output/sitebysite_benchmark_muresk_ode_100y_0521_npp15_sens.rds")


new_sim_outs_ph_4 <- new_sim_outs_npp05[which(new_sim_outs_npp05$pH == 4.0),]

new_sim_outs_ph_4 <- new_sim_outs_ph_4[order(new_sim_outs_ph_4$CUE),]

new_sim_outs_sort <- new_sim_outs_npp15[order(new_sim_outs_npp15$CUE),]

plot(new_sim_outs_sort$pH,new_sim_outs_sort$SOM1000)
#########################################################
#clay

outline_df <-readRDS("/media/DATADRIVE1/Data/muresk/data/d02_data/for_model/muresk_outline.RDS")

get_csat_frombound <- function(boundary_df, csat_df) {
  # a vector to store deficit
  csat_df["ycsat"] <- NA
  
  # x step in the boundary_df
  xbound_step <- 0.5
  
  # get boundary at xi
  xi <- csat_df[1, "x"]
  if (xi %in% boundary_df$x) {
    ybound_i <- boundary_df[boundary_df$x == xi, "y"]
  }
  else {
    xi_left <- floor(xi / xbound_step) * 0.5
    xi_right <- floor(xi / xbound_step) * 0.5 + 0.5
    
    yi_left <- boundary_df[boundary_df$x == xi_left, "y"]
    yi_right <- boundary_df[boundary_df$x == xi_right, "y"]
    
    xi_slope <- (yi_right - yi_left) / xbound_step
    ybound_i <- yi_left + (xi - xi_left) * xi_slope
  }
  
  # calculate deficit
  csat_df[1, "ycsat"] <- ybound_i
  
  
  return(csat_df)
  
}

cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)

sim_col <- foreach(
  i            = seq_along(flist),
  .packages    = c("doParallel","deSolve","zoo")
  # if safe_ode_100_fn or pars_set live in your global env:
  # for reproducible *and distinct* streams, uncomment:
  # .options.RNG = 123,
) %dopar% {
  
  ## (1) No set.seed() here, or if you need it,
  ##     use doRNG/.options.RNG so each i gets its own stream.
  set.seed(123)
  
  driving_input <- read.csv(flist[i])
  dates <- as.Date(with(driving_input, ISOdate(year,month,day)), tz="UTC")
  doy   <- as.POSIXlt(dates)$yday + 1
  keep  <- format(dates, "%m-%d") != "02-29"
  agg_NPP   <- tapply(driving_input$Npp[keep],         doy[keep], mean, na.rm=TRUE)
  agg_moist <- tapply(driving_input$SoilMoisture[keep],doy[keep], mean, na.rm=TRUE)
  agg_temp  <- tapply(driving_input$SoilTemp[keep],    doy[keep], mean, na.rm=TRUE)
  
  agg <- data.frame(
    doy      = as.integer(names(agg_NPP)),
    forc_npp = as.numeric(agg_NPP),
    forc_sw  = as.numeric(agg_moist),
    forc_st  = as.numeric(agg_temp)
  )
  agg$avg_date <- as.Date(agg$doy, origin = as.Date("2000-12-31"))
  agg <- agg[1:365, ]
  
  agg$forc_npp <- with(agg,forc_npp*1.5)
  sitename      <- sub("\\.csv$", "", basename(flist[i]))
  selected_row  <- pars_set[sitename, ]
  
  ## (2) Build up a *local* data.frame with <-, not <<-
  sim_output_init <- NULL
  
  for (clay_silt in seq(16, 36, 4)) {
    
        csat_df <- as.data.frame(list("x"  = clay_silt,"y"  = 1.8))
        
        csat_df_target <- get_csat_frombound(outline_df,csat_df)
        
        out_qmax <- safe_ode_100_fn(selected_row, sitename, agg, csat_df_target$ycsat)
        
        out_qmax <- c(out_qmax, Npp_ratio = 1.5)
        
        out_qmax <- c(out_qmax, clay_silt = clay_silt)
        
        sim_output_init <- rbind(sim_output_init, out_qmax)
        
        csat_df_target = NULL
         
  
    
  }
  
  ## scale & clean
  sim_output_init <- na.omit(sim_output_init)
  sens <- scale(sim_output_init[,1:6])
  sens <- as.data.frame(sens)
  #sens$clay_silt  <- seq(16, 36, 4)
  #sens$pH         <- sim_output_init[,8]
  sens$CUE        <- sim_output_init[,7]
  sens$clay_silt  <- sim_output_init[,10]
  sens$Npp_ratio  <- sim_output_init[,9]
  colnames(sens)[1:6] <- c('SOM1000','MAOM1000','POM1000',
                           'AGG1000','MIC1000','LMWC1000')
  ## add a column so you know which site this came from
  sens$site <- sitename
  sens
}

stopCluster(cl)

new_sim_outs <- do.call(rbind, sim_col)

#new_sim_outs_sort <- new_sim_outs[order(new_sim_outs$clay_silt),]

plot(new_sim_outs$clay_silt,new_sim_outs$SOM1000)
#colnames(new_sim_outs) <- c('SOM1000', 'MAOM1000', 'POM1000', 'AGG1000', 'MIC1000', 'LMWC1000',"pH")

new_sim_outs <- as.data.frame(new_sim_outs)

saveRDS(new_sim_outs,"/media/DATADRIVE1/Data/muresk/data/d04_output/sitebysite_benchmark_muresk_ode_100y_0521_npp15_claysiltsens.rds")

claysilt_sim_outs <- readRDS("/media/DATADRIVE1/Data/muresk/data/d04_output/sitebysite_benchmark_muresk_ode_100y_0521_npp10_claysiltsens.rds")

claysilt_sim_outs$Npp_ratio  <- 1.0

saveRDS(new_sim_outs,"/media/DATADRIVE1/Data/muresk/data/d04_output/sitebysite_benchmark_muresk_ode_100y_0523_mixedsens.rds")

ggplot(new_sim_outs_comp, aes(x = clay_silt, y = SOM1000, fill=Npp_ratio)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "lightgray") +
  # stat_summary(fun = mean,
  #              geom = "point",
  #              shape = 23, size = 3, fill = "white") +
  labs(
    title = "Distribution of Value by Class (Group A)",
    x     = "Class",
    y     = "Value"
  ) +
  theme_classic()


##############################################

library(dplyr)
library(ggplot2)

# 1) compute all the stats you want
stats_df <- new_sim_outs_comp %>%
  group_by(clay_silt, Npp_ratio) %>%
  summarise(
    q5     = quantile(SOM1000, 0.05, na.rm = TRUE),
    q25    = quantile(SOM1000, 0.25, na.rm = TRUE),
    q50    = quantile(SOM1000, 0.50, na.rm = TRUE),
    q75    = quantile(SOM1000, 0.75, na.rm = TRUE),
    q95    = quantile(SOM1000, 0.95, na.rm = TRUE),
    mean_v = mean(SOM1000,        na.rm = TRUE),
    .groups = "drop"
  )

# 2) draw “boxplots” from those numbers, then add the mean
ggplot(stats_df, aes(x = clay_silt, fill = Npp_ratio)) +
  geom_boxplot(
    aes(ymin   = q5,
        lower  = q25,
        middle = q50,
        upper  = q75,
        ymax   = q95),
    stat  = "identity",
    width = 0.7,
    alpha = 0.6
  ) +
  scale_fill_discrete(name = "C input factor")+
  labs(
    title = "",
    x     = "Clay–Silt (%)",
    y     = "Standardized SOC"
  ) +
  theme_classic()

#########################################################
#pH+NPP+clay
cl <- parallel::makeCluster(15)
doParallel::registerDoParallel(cl)



sim_col <- foreach(
  i            = seq_along(flist),
  .packages    = c("doParallel","deSolve","zoo")
  # if safe_ode_100_fn or pars_set live in your global env:
  # for reproducible *and distinct* streams, uncomment:
  # .options.RNG = 123,
) %dopar% {
  
  ## (1) No set.seed() here, or if you need it,
  ##     use doRNG/.options.RNG so each i gets its own stream.
  set.seed(123)
  
  driving_input <- read.csv(flist[i])
  dates <- as.Date(with(driving_input, ISOdate(year,month,day)), tz="UTC")
  doy   <- as.POSIXlt(dates)$yday + 1
  keep  <- format(dates, "%m-%d") != "02-29"
  agg_NPP   <- tapply(driving_input$Npp[keep],         doy[keep], mean, na.rm=TRUE)
  agg_moist <- tapply(driving_input$SoilMoisture[keep],doy[keep], mean, na.rm=TRUE)
  agg_temp  <- tapply(driving_input$SoilTemp[keep],    doy[keep], mean, na.rm=TRUE)
  
  agg <- data.frame(
    doy      = as.integer(names(agg_NPP)),
    forc_npp = as.numeric(agg_NPP),
    forc_sw  = as.numeric(agg_moist),
    forc_st  = as.numeric(agg_temp)
  )
  agg$avg_date <- as.Date(agg$doy, origin = as.Date("2000-12-31"))
  agg <- agg[1:365, ]
  
  #agg$forc_npp <- with(agg,forc_npp*1.5)
  sitename      <- sub("\\.csv$", "", basename(flist[i]))
  selected_row  <- pars_set[sitename, ]
  
  ## (2) Build up a *local* data.frame with <-, not <<-
  sim_output_init <- NULL
  
  for (clay_silt in seq(16, 36, 4)) {
    
    for (ph_sens in seq(4, 6.5, 0.5)) {
      
      for (npp_ratio in seq(0.5,1.5,0.5)) {
        
          set.seed(123)
    
          csat_df <- as.data.frame(list("x"  = clay_silt,"y"  = 1.8))
          
          csat_df_target <- get_csat_frombound(outline_df,csat_df)
          
          temp_npp     <- agg
          
          temp_npp$forc_npp <- with(temp_npp,forc_npp*npp_ratio)
          
          out_mix <- safe_ode_mix_fn(selected_row, sitename, temp_npp, csat_df_target$ycsat,ph_sens)
          
          out_mix <- c(out_mix, Npp_ratio = npp_ratio)
          
          out_mix <- c(out_mix, clay_silt = clay_silt)
          
          sim_output_init <- rbind(sim_output_init, out_mix)
          
          csat_df_target = NULL
          
    
        }
      
    }
    
  }
  
  ## scale & clean
  sim_output_init <- na.omit(sim_output_init)
  sens <- scale(sim_output_init[,1:6])
  sens <- as.data.frame(sens)
  sens$CUE        <- sim_output_init[,7]
  sens$pH         <- sim_output_init[,9]
  sens$Npp_ratio  <- sim_output_init[,10]
  sens$clay_silt  <- sim_output_init[,11]
  
  colnames(sens)[1:6] <- c('SOM1000','MAOM1000','POM1000',
                           'AGG1000','MIC1000','LMWC1000')
  ## add a column so you know which site this came from
  sens$site <- sitename
  sens

}

stopCluster(cl)

new_sim_outs <- do.call(rbind, sim_col)

saveRDS(new_sim_outs,"/media/DATADRIVE1/Data/muresk/data/d04_output/sitebysite_benchmark_muresk_ode_100y_0530_mixed_sens.rds")

################################################################################

################################################################################
#CUE+NPP
cl <- parallel::makeCluster(15)
doParallel::registerDoParallel(cl)

fourier_temp <- function(t, a0, a1, b1, a2, b2, a3, b3) {
  omega <- 2 * pi / 365  # yearly frequency in days
  a0 + 
    a1 * cos(omega * t) + b1 * sin(omega * t) +
    a2 * cos(2 * omega * t) + b2 * sin(2 * omega * t) +
    a3 * cos(3 * omega * t) + b3 * sin(3 * omega * t)
}


# Create the transformed function

sim_col <- foreach(
  i            = seq_along(flist),
  .packages    = c("doParallel","deSolve","zoo","minpack.lm")
  # if safe_ode_100_fn or pars_set live in your global env:
  # for reproducible *and distinct* streams, uncomment:
  # .options.RNG = 123,
) %dopar% {
  
  ## (1) No set.seed() here, or if you need it,
  ##     use doRNG/.options.RNG so each i gets its own stream.
  set.seed(123)
  
  driving_input <- read.csv(flist[i])
  dates <- as.Date(with(driving_input, ISOdate(year,month,day)), tz="UTC")
  doy   <- as.POSIXlt(dates)$yday + 1
  keep  <- format(dates, "%m-%d") != "02-29"
  agg_NPP   <- tapply(driving_input$Npp[keep],         doy[keep], mean, na.rm=TRUE)
  agg_moist <- tapply(driving_input$SoilMoisture[keep],doy[keep], mean, na.rm=TRUE)
  agg_temp  <- tapply(driving_input$SoilTemp[keep],    doy[keep], mean, na.rm=TRUE)
  
  agg <- data.frame(
    doy      = as.integer(names(agg_NPP)),
    forc_npp = as.numeric(agg_NPP),
    forc_sw  = as.numeric(agg_moist),
    forc_st  = as.numeric(agg_temp)
  )
  agg$avg_date <- as.Date(agg$doy, origin = as.Date("2000-12-31"))
  agg <- agg[1:365, ]
  

  #agg$forc_npp <- with(agg,forc_npp*1.5)
  sitename      <- sub("\\.csv$", "", basename(flist[i]))
  selected_row  <- pars_set[sitename, ]
  #temperatures = agg$forc_st
  #days = agg$doy
  #fit <- nlsLM(temperatures~ fourier_temp(days, a0, a1, b1, a2, b2, a3, b3),
  #             start = list(a0 = mean(temperatures), a1 = 1, b1 = 1, 
  #                          a2 = 0.5, b2 = 0.5, a3 = 0.1, b3 = 0.1))
  #
  ## Extract coefficients
  #coeffss <- coef(fit)
  #
  #y_function <- function(day,coeffss,selected_row) {
  #  
  #  coeffs_fit = coeffss
  #  
  #  temp <- fourier_temp(day, coeffs_fit[1], coeffs_fit[2], coeffs_fit[3],
  #                       coeffs_fit[4], coeffs_fit[5], coeffs_fit[6], coeffs_fit[7])
  #  y <- selected_row$cue_ref - 0.012 * (temp - selected_row$tae_ref)
  #  
  #  return(y)
  #}
  
  #y_values <- sapply(1:365, function(day) y_function(day, coeffss,selected_row))
  #average_y <- mean(y_values)
  average_temp <- mean(agg$forc_st)
  
  ## (2) Build up a *local* data.frame with <-, not <<-
  sim_output_init <- NULL
  
  for (CUE in seq(0.2,0.7,0.1)){
  
    for (npp_ratio in seq(1,1.5,0.5)) {
        
        set.seed(123)
        
        temp_npp     <- agg
        
        temp_npp$forc_npp <- with(temp_npp,forc_npp*npp_ratio)
        
        #temp_npp$forc_st <- with(temp_npp,forc_st+(average_y-CUE)/0.012)
        CUE_corr <- CUE - selected_row$cue_ref - 0.012*selected_row$tae_ref + 0.012*average_temp
        
        out_mix <- senst_ode_cue_100(selected_row, sitename, temp_npp, CUE_corr)
        
        out_mix <- c(out_mix, Npp_ratio = npp_ratio)
        
        out_mix <- c(out_mix, CUE_corr = CUE_corr)
        
        sim_output_init <- rbind(sim_output_init, out_mix)
        
    
   }
    
  }
  
  ## scale & clean
  sim_output_init <- na.omit(sim_output_init)
  sens <- scale(sim_output_init[,1:6])
  sens <- as.data.frame(sens)
  sens$CUE        <- sim_output_init[,7]
  sens$Npp_ratio  <- sim_output_init[,8]
  sens$CUE_corr   <- sim_output_init[,9]
  colnames(sens)[1:6] <- c('SOM1000','MAOM1000','POM1000',
                           'AGG1000','MIC1000','LMWC1000')
  ## add a column so you know which site this came from
  sens$site <- sitename
  sens
  
}

stopCluster(cl)

new_sim_outs <- do.call(rbind, sim_col)

saveRDS(new_sim_outs,"/media/DATADRIVE1/Data/muresk/data/d04_output/sitebysite_benchmark_muresk_ode_100y_0603_cue_sens.rds")






################################################################################
A = c(16,20,24,28,32,36)

B = c(4.0,4.5,5.0,5.5,6.0,6.5)

C = c(0.5,1.0,1.5)

sim_output_out <- NULL

sitelists <- unique(new_sim_outs$site)

for (sitename in sitelists){

site_sim_outs <- new_sim_outs[which(new_sim_outs$site == sitename),]

for (i in seq(1:6)) {
  
  for (j in seq(1:6)) {
    
    for (k in seq(1:3)) {
      
        index = 18*(i-1)+3*(j-1)+k
        
        sim_output_index <- site_sim_outs[index,]
          
        sim_output_index <- cbind(sim_output_index, Npp_ratio = C[k])
          
        sim_output_index <- cbind(sim_output_index, clay_silt = A[i])
          
        sim_output_out   <- rbind(sim_output_out,sim_output_index)
        
      }
        
    }
    
  }
  
}

sim_output_out_npp10 <- sim_output_out[which(sim_output_out$Npp_ratio == 1.0),]
sim_output_out_npp05 <- sim_output_out[which(sim_output_out$Npp_ratio == 0.5),]
sim_output_out_npp15 <- sim_output_out[which(sim_output_out$Npp_ratio == 1.5),]

plot(sim_output_out_npp05$CUE,sim_output_out_npp05$SOM1000)
plot(sim_output_out_npp10$CUE,sim_output_out_npp10$SOM1000)
plot(sim_output_out_npp15$CUE,sim_output_out_npp15$SOM1000)


sim_output_out_npp10_ph65 <- sim_output_out[which(sim_output_out$Npp_ratio == 1.0 & sim_output_out$pH == 6.5),]
sim_output_out_npp10_ph60 <- sim_output_out[which(sim_output_out$Npp_ratio == 1.0 & sim_output_out$pH == 6.0),]
sim_output_out_npp10_ph55 <- sim_output_out[which(sim_output_out$Npp_ratio == 1.0 & sim_output_out$pH == 5.5),]
sim_output_out_npp10_ph50 <- sim_output_out[which(sim_output_out$Npp_ratio == 1.0 & sim_output_out$pH == 5.0),]
sim_output_out_npp10_ph45 <- sim_output_out[which(sim_output_out$Npp_ratio == 1.0 & sim_output_out$pH == 4.5),]
sim_output_out_npp10_ph40 <- sim_output_out[which(sim_output_out$Npp_ratio == 1.0 & sim_output_out$pH == 4.0),]


sim_output_out_npp10_clay16 <- sim_output_out[which(sim_output_out$Npp_ratio == 1.0 & sim_output_out$clay_silt == 16),]
sim_output_out_npp10_clay20 <- sim_output_out[which(sim_output_out$Npp_ratio == 1.0 & sim_output_out$clay_silt == 20),]
sim_output_out_npp10_clay24 <- sim_output_out[which(sim_output_out$Npp_ratio == 1.0 & sim_output_out$clay_silt == 24),]
sim_output_out_npp10_clay28 <- sim_output_out[which(sim_output_out$Npp_ratio == 1.0 & sim_output_out$clay_silt == 28),]
sim_output_out_npp10_clay32 <- sim_output_out[which(sim_output_out$Npp_ratio == 1.0 & sim_output_out$clay_silt == 32),]
sim_output_out_npp10_clay36 <- sim_output_out[which(sim_output_out$Npp_ratio == 1.0 & sim_output_out$clay_silt == 36),]

sim_output_out_npp10_clay16_pH65 <- sim_output_out[which(sim_output_out$Npp_ratio == 1.0 & sim_output_out$clay_silt == 16 & sim_output_out$pH == 6.5),]

plot(sim_output_out_npp10_clay16_pH65$CUE,sim_output_out_npp10_clay16_pH65$SOM1000)

#summary(lm(sim_output_out_npp10_clay16_pH65$CUE ~ sim_output_out_npp10_clay16_pH65$SOM1000))

plot(sim_output_out_npp10_ph65$clay_silt,sim_output_out_npp10_ph65$MAOM1000)
plot(sim_output_out_npp10_ph60$clay_silt,sim_output_out_npp10_ph60$MAOM1000)
plot(sim_output_out_npp10_ph55$clay_silt,sim_output_out_npp10_ph55$MAOM1000)
plot(sim_output_out_npp10_ph50$clay_silt,sim_output_out_npp10_ph50$MAOM1000)
plot(sim_output_out_npp10_ph45$clay_silt,sim_output_out_npp10_ph45$MAOM1000)
plot(sim_output_out_npp10_ph40$clay_silt,sim_output_out_npp10_ph40$MAOM1000)

plot(sim_output_out_npp10_clay16$pH,sim_output_out_npp10_clay16$SOM1000)
plot(sim_output_out_npp10_clay20$pH,sim_output_out_npp10_clay20$SOM1000)
plot(sim_output_out_npp10_clay24$pH,sim_output_out_npp10_clay24$SOM1000)
plot(sim_output_out_npp10_clay28$pH,sim_output_out_npp10_clay28$SOM1000)
plot(sim_output_out_npp10_clay32$pH,sim_output_out_npp10_clay32$SOM1000)
plot(sim_output_out_npp10_clay36$pH,sim_output_out_npp10_clay36$SOM1000)


saveRDS(sim_output_out,"/media/DATADRIVE1/Data/muresk/data/d04_output/sitebysite_benchmark_muresk_ode_100y_0522_mixed_sens.rds")



#########################################################
#full ode version pH

sim_col <- foreach(i=seq(flist),.packages = c("rtop","doParallel","deSolve","zoo")) %dopar% {
  
  set.seed(123)
  
  driving_input <- read.csv(flist[i])
  
  dates <- as.Date(with(driving_input, ISOdate(year, month, day)), tz = "UTC")
  # compute doy and drop Feb 29
  doy <- as.POSIXlt(dates)$yday + 1
  keep <- !(format(dates, "%m-%d")=="02-29")
  # aggregate each variable
  agg_NPP  <- tapply(driving_input$Npp[keep],  doy[keep], mean, na.rm=TRUE)
  agg_moist<- tapply(driving_input$SoilMoisture[keep], doy[keep], mean, na.rm=TRUE)
  agg_temp <- tapply(driving_input$SoilTemp[keep],doy[keep], mean, na.rm=TRUE)
  # combine into a data.frame
  # "forc_st", "forc_sw", "forc_npp"
  agg <- data.frame(
    doy             = as.integer(names(agg_NPP)),
    forc_npp        = as.numeric(agg_NPP),
    forc_sw         = as.numeric(agg_moist),
    forc_st         = as.numeric(agg_temp)
  )
  agg$avg_date <- as.Date(agg$doy, origin = as.Date("2000-12-31"))
  agg <- agg[1:365,]
  
  sitename <- gsub("*.csv","",basename(flist[i]))
  
  selected_row <- pars_set[sitename, ]
  
  sim_output_init <- NULL
  
  for (ph in seq(4,6.5,0.5)){
    
    sim_output_pH <-  safe_ode_100_fn(selected_row,sitename,agg,ph)
    
    sim_output_fullpH <-  safe_ode_fn(selected_row,sitename,agg,ph)
    #rownames(sim_output_pH) <- sitename
    
    sim_output_init <<- rbind(sim_output_init, sim_output_pH)
  }
  
  sim_output_init_sens <- scale(sim_output_init[ , 1:6], center = TRUE, scale = TRUE)
  sim_output_init_sens <- as.data.frame(sim_output_init_sens)
  sim_output_init_sens$pH <- c(seq(4,6.5,0.5))
  rownames(sim_output_init_sens) <- NULL
  colnames(sim_output_init_sens) <- c('SOM1000', 'MAOM1000', 'POM1000', 'AGG1000', 'MIC1000', 'LMWC1000',"pH")
  
  return(sim_output_init_sens)
  
}

new_sim_outs <- do.call(rbind, sim_col)

colnames(new_sim_outs) <- c('SOM1000', 'MAOM1000', 'POM1000', 'AGG1000', 'MIC1000', 'LMWC1000',"pH")

new_sim_outs <- as.data.frame(new_sim_outs)















