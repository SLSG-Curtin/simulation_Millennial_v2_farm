library(rtop)
library(deSolve)
library(FME)
library(rsq)
library(foreach)
library(doParallel)
library(tools)
library(dplyr)
library(lubridate)
library(zoo)

set.seed(123)

#
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
            22.5,        22.5,    3.15e-4,        0.9,       0.3,    1)


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

#excluded_params <- c('param_pi','rate_leach','kaff_pl', 'param_p2', 'kaff_lb','matpot', 'kamin', 'lambda')

important_params_set6 <- c('param_pa', 'rate_pa', 'rate_break', 'rate_leach',
                           'rate_ma',  'matpot',  'eact_pl',    'eact_lb',
                           'cue_ref',  'tae_ref', 'kaff_des',   'param_p1',
                           'param_pi')

#selected_params <- params_df[!params_df$Parameter %in% excluded_params, ]
selected_params <- params_df[params_df$Parameter %in% important_params_set6, ]


flist <- list.files("/media/DATADRIVE1/Project/muresk/driving/forcing_inputs_new/",full.names = TRUE,recursive = TRUE)
flists <- sub(".csv","",basename(flist))

litter_soil_CNratio <- read.csv("/home/286159b/Desktop/muresk/CNratio/litter_soil_CNratio.csv")
litter_soil_CNratio <- litter_soil_CNratio[,c(2,15)]
colnames(litter_soil_CNratio) <- c("site","CNratio")

compile_muresk_csat <- read.csv("/media/DATADRIVE1/Project/muresk/model/d02_data/for_model/compile_muresk_csat.txt")

#cfractions in gC/m2
compile_muresk_csat$pocmac030 <- with(compile_muresk_csat,POCmac_30/100*BD_30*30*10000)
compile_muresk_csat$pocmic030 <- with(compile_muresk_csat,POCmic_30/100*BD_30*30*10000)
compile_muresk_csat$maoc030   <- with(compile_muresk_csat,MAOC_30/100*BD_30*30*10000)
compile_muresk_csat$toc030    <- with(compile_muresk_csat,toc_30/100*BD_30*30*10000)
compile_muresk_csat$qmax030   <- with(compile_muresk_csat,Csat/100*BD_30*30*10000)

# compile_muresk_csat <- compile_muresk_csat[,c(2:20)]
compile_muresk_csat <- compile_muresk_csat[,c(2,16,17,18,19,14,20,9)]
colnames(compile_muresk_csat) <- c("site","POM","AGG","MAOM","SOM","pH","qmax","BD")
compile_muresk_csat <- merge(compile_muresk_csat,litter_soil_CNratio,by="site")

# parameter sets
pars_set <- readRDS("/home/286159b/Desktop/muresk/sitebysite_muresk_threefrac_ode_100y_1113_init.rds")
pars_set <- pars_set[apply(pars_set, 1, function(x) all(is.finite(x))), ]
pars_set <- pars_set[,1:13]
#5594,5658

derivs.rangeland.MV6 <- function(step.num,state,params_sceua,forc_st,forc_sw,forc_npp) {
  
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
    param_CN   = 1.0
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
      
      if (CUE_tmp > 0.8) {
        stop("CUE exceeds 0.8 - invalid parameter combination")
      }
      
      Adj_CNratio = parameters$param_CN^(-0.6)
      
      #Adj_CNratio = pmax(0,parameters$param_CN^(-0.6))
      
      CUE_realized = CUE_tmp*Adj_CNratio
      #
      if (CUE_realized > 0.8) {
        stop("CUE exceeds 0.8 - invalid parameter combination")
      }
      
      #CUE_realized = pmin(1.0, pmax(0.0, CUE_realized))
      
      f_MB_atm = f_LM_MB * (1 - CUE_realized)
      f_MB_gw  = f_LM_MB * CUE_realized
      
    } else {
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
    
    f_mic_resp = f_MB_atm + f_MB_gw
    
    if (!is.na(f_mic_resp) && f_mic_resp > 0){
      rCUE = f_MB_gw/(f_MB_atm+f_MB_gw)
    } else {
      rCUE = 0
    }
    
    return(list(c(dPOM, dLMWC, dAGG, dMIC, dMAOM),'CUE'=rCUE))
    
  })
}

derivs.rangeland.MV6.CNcorr <- function(step.num,state,params_sceua,forc_st,forc_sw,forc_npp) {
  
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
    param_CN   = 1,
    param_CUEcorr = 0.0,
    param_CNcorr = 0.0,
    use_CN_adj = FALSE
  )
  
  parameters <- modifyList(MV3_params, params_sceua)
  
  with(as.list(c(state,parameters)), {
    
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
    
    if(MIC>0 && LMWC>0){
      
      #f_MB_atm = f_LM_MB * (1 - (parameters$cue_ref - parameters$cue_t * (forc_st(step.num) - parameters$tae_ref) ) )
      #f_MB_gw  = f_LM_MB * (parameters$cue_ref - parameters$cue_t * (forc_st(step.num) - parameters$tae_ref))
      CUE_tmp = parameters$cue_ref - parameters$cue_t * (forc_st(step.num) - parameters$tae_ref)
      
      CUE_adjusted = CUE_tmp + parameters$param_CUEcorr
      
      #f (CUE_adjusted > 0.8) {
      # stop("CUE exceeds 0.8 - invalid parameter combination")
      #
      
      #use_CN_adj_flag <- if (!is.null(parameters$use_CN_adj)) isTRUE(parameters$use_CN_adj) else TRUE
      
      use_CN_adj_flag <- parameters$use_CN_adj
      
      if (use_CN_adj_flag) {
        CNratio_effective = parameters$param_CN + parameters$param_CNcorr
        Adj_CNratio = CNratio_effective^(-0.6)
      } else {
        Adj_CNratio = 1
      }
      
      CUE_realized = CUE_adjusted*Adj_CNratio
      #
      #if (CUE_realized > 0.8) {
      #  stop("CUE exceeds 0.8 - invalid parameter combination")
      #}
      
      f_MB_atm = f_LM_MB * (1 - CUE_adjusted*Adj_CNratio)
      f_MB_gw  = f_LM_MB * CUE_adjusted*Adj_CNratio
      
    } else {
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
    
    f_mic_resp = f_MB_atm + f_MB_gw
    
    if (!is.na(f_mic_resp) && f_mic_resp > 0){
      rCUE = f_MB_gw/(f_MB_atm+f_MB_gw)
    } else {
      rCUE = 0
    }
    
    return(list(c(dPOM, dLMWC, dAGG, dMIC, dMAOM),'CUE'=rCUE))
    
  })
}

derivs.rangeland.MV6.All <- function(step.num,state,params_sceua,forc_st,forc_sw,forc_npp) {
  
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
    param_CN   = 1,
    param_CUEcorr = 0.0,
    param_CNcorr = 0.0,
    use_CN_adj = FALSE
  )
  
  parameters <- modifyList(MV3_params, params_sceua)
  
  with(as.list(c(state,parameters)), {
    
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
    
    if(MIC>0 && LMWC>0){
      
      #f_MB_atm = f_LM_MB * (1 - (parameters$cue_ref - parameters$cue_t * (forc_st(step.num) - parameters$tae_ref) ) )
      #f_MB_gw  = f_LM_MB * (parameters$cue_ref - parameters$cue_t * (forc_st(step.num) - parameters$tae_ref))
      CUE_tmp = parameters$cue_ref - parameters$cue_t * (forc_st(step.num) - parameters$tae_ref)
      
      CUE_adjusted = CUE_tmp + parameters$param_CUEcorr
      
      #f (CUE_adjusted > 0.8) {
      # stop("CUE exceeds 0.8 - invalid parameter combination")
      #
      
      #use_CN_adj_flag <- if (!is.null(parameters$use_CN_adj)) isTRUE(parameters$use_CN_adj) else TRUE
      
      use_CN_adj_flag <- parameters$use_CN_adj
      
      if (use_CN_adj_flag) {
        CNratio_effective = parameters$param_CN + parameters$param_CNcorr
        Adj_CNratio = CNratio_effective^(-0.6)
      } else {
        Adj_CNratio = 1
      }
      
      CUE_realized = CUE_adjusted*Adj_CNratio
      #
      #if (CUE_realized > 0.8) {
      #  stop("CUE exceeds 0.8 - invalid parameter combination")
      #}
      
      f_MB_atm = f_LM_MB * (1 - CUE_adjusted*Adj_CNratio)
      f_MB_gw  = f_LM_MB * CUE_adjusted*Adj_CNratio
      
    } else {
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
  
  #output <- derivs.rangeland.MV6.CNcorr(step.num,state,params_sceua,forc_st,forc_sw,forc_npp)
  output <- derivs.rangeland.MV6.All(step.num,state,params_sceua,forc_st,forc_sw,forc_npp)
  return(list(output[[1]],output[[2]]))
}


senmix_ode_100 <- function(params,sitename,dat,CUE_corr,CNratio_corr,use_CN_adj) {
  
  obs_row  <- compile_muresk_csat[which(compile_muresk_csat$site == sitename),]
  
  initial_values <- c(POM = obs_row$POM, LMWC = 6, AGG = obs_row$AGG, MIC = 6, MAOM=obs_row$MAOM)  
  
  #100years
  forc_st_fun  <- approxfun(1:SStime, rep(dat$forc_st, 100))
  forc_sw_fun  <- approxfun(1:SStime, rep(dat$forc_sw, 100))
  forc_npp_fun <- approxfun(1:SStime, rep(dat$forc_npp,100))
  
  param_qmax  <- obs_row$qmax
  param_pH    <- obs_row$pH
  #param_qmax    <- qmax_sens/100*obs_row$BD*30*10000
  #param_pH      <- pH_sens
  param_CN    <- obs_row$CNratio
  
  params_sceua <- c(as.list(setNames(params, selected_params$Parameter)),
                    param_qmax    = param_qmax,
                    param_pH      = param_pH,
                    param_CN      = param_CN,
                    param_CUEcorr = CUE_corr,
                    param_CNcorr  = CNratio_corr,
                    use_CN_adj    = use_CN_adj)
  
  
  ODE.MMRMM <- ode(y = initial_values, times=run.steps,  func = derivs_Aus_wrapper, parms = params_sceua, forc_st=forc_st_fun, forc_sw=forc_sw_fun, forc_npp=forc_npp_fun, method="rk4")
  #ODE.MMRMM <- forcing_delta_MM(initial_values,params_sceua,forcing_dat)
  
  ODE.MMRMM <- as.data.frame(ODE.MMRMM)
  # get the mean of equlibrium state, suggesting 120 rows
  ODE.MMRMM$SOM <- with(ODE.MMRMM,POM+LMWC+AGG+MIC+MAOM)
  
  tail_idx <- tail(seq_len(nrow(ODE.MMRMM)), 365)
  ode.df <- as.data.frame(ODE.MMRMM)
  names(ode.df)[7] <- c("CUE")
  sim_update <- colMeans(ode.df[tail_idx, c("POM", "LMWC", "AGG", "MIC", "MAOM", "SOM", "CUE")])
  
  return(c(sim_update))
  
}

safe_ode_mix_fn <- function(params,sitename,dat,CUE_corr,CNratio_corr,use_CN_adj) {
  J_result <- tryCatch({
    # Replace 'actual_objective_function' with your actual function
    senmix_ode_100(params,sitename,dat,CUE_corr,CNratio_corr,use_CN_adj)
  }, error = function(e) {
    # Handle the error, e.g., by returning a large penalty value or NA
    return(Inf)  # or NA, depending on how sceua() handles it
  })
  
  return(J_result)
}

senmix_ode_all <- function(params,sitename,dat,qmax_sens,pH_sens,CUE_corr,CNratio_corr,use_CN_adj) {
  
  obs_row  <- compile_muresk_csat[which(compile_muresk_csat$site == sitename),]
  
  initial_values <- c(POM = obs_row$POM, LMWC = 6, AGG = obs_row$AGG, MIC = 6, MAOM=obs_row$MAOM)  
  
  #100years
  forc_st_fun  <- approxfun(1:SStime, rep(dat$forc_st, 100))
  forc_sw_fun  <- approxfun(1:SStime, rep(dat$forc_sw, 100))
  forc_npp_fun <- approxfun(1:SStime, rep(dat$forc_npp,100))
  
  param_qmax    <- qmax_sens/100*obs_row$BD*30*10000
  param_pH      <- pH_sens
  param_CUEcorr <- CUE_corr      # CUE correction to shift base level
  param_CNcorr  <- CNratio_corr  # CNratio correction for substrate quality
  param_CN      <- obs_row$CNratio
  
  params_sceua <- c(as.list(setNames(params, selected_params$Parameter)),
                    param_qmax = param_qmax,
                    param_pH   = param_pH,
                    param_CUEcorr = param_CUEcorr,  # CUE correction
                    param_CNcorr = param_CNcorr,     # CNratio correction
                    param_CN   = param_CN)
  
  
  ODE.MMRMM <- ode(y = initial_values, times=run.steps,  func = derivs_Aus_wrapper, parms = params_sceua, forc_st=forc_st_fun, forc_sw=forc_sw_fun, forc_npp=forc_npp_fun, method="rk4")
  #ODE.MMRMM <- forcing_delta_MM(initial_values,params_sceua,forcing_dat)
  
  ODE.MMRMM <- as.data.frame(ODE.MMRMM)
  # get the mean of equlibrium state, suggesting 120 rows
  ODE.MMRMM$SOM <- with(ODE.MMRMM,POM+LMWC+AGG+MIC+MAOM)
  
  #ODE.MMRMM_tail <- tail(ODE.MMRMM, 365)
  
  #row.names(ODE.MMRMM_tail) <- dat$doy
  
  #ODE.MMRMM_tail$pH <- pH_sens
  #sim_update <- rollapply(ODE.MMRMM, width = 365, by = 365, FUN = mean)
  #sim_update <- as.data.frame(sim_update)
  
  tail_idx <- tail(seq_len(nrow(ODE.MMRMM)), 365)
  ode.df <- as.data.frame(ODE.MMRMM)
  names(ode.df)[7] <- c("CUE")
  sim_update <- colMeans(ode.df[tail_idx, c("POM", "LMWC", "AGG", "MIC", "MAOM", "SOM", "CUE")])
  
  #POM10000       <- mean(tail(sim_update$POM, 1))
  #LMWC10000      <- mean(tail(sim_update$LMWC, 1))
  #AGG10000       <- mean(tail(sim_update$AGG, 1))
  #MIC10000       <- mean(tail(sim_update$MIC, 1))
  #MAOM10000      <- mean(tail(sim_update$MAOM,1))
  #SOM10000       <- mean(tail(sim_update$SOM,1))
  #CUE10000       <- mean(tail(sim_update$CUE,1))
  
  #POM_dev     <- abs(POM10000 - obs_row$POM)
  #AGG_dev     <- abs(AGG10000 - obs_row$AGG)
  #MAOM_dev    <- abs(MAOM10000 - obs_row$MAOM)
  #SOM_dev     <- abs(SOM10000  - obs_row$SOM)
  
  #return(c(SOM10000, MAOM10000, POM10000,AGG10000,MIC10000, LMWC10000, CUE10000, Qmax = param_qmax, pH = param_pH))
  return(c(sim_update, Qmax = param_qmax, pH = param_pH))
}

safe_ode_all_fn <- function(params,sitename,dat,qmax_sens,pH_sens,CUE_corr,CNratio_corr,use_CN_adj) {
  J_result <- tryCatch({
    # Replace 'actual_objective_function' with your actual function
    senmix_ode_all(params,sitename,dat,qmax_sens,pH_sens,CUE_corr,CNratio_corr,use_CN_adj)
  }, error = function(e) {
    # Handle the error, e.g., by returning a large penalty value or NA
    return(Inf)  # or NA, depending on how sceua() handles it
  })
  
  return(J_result)
}


SStime     <- 36500
run.steps  <- seq(1,36500)


set.seed(123)

#new quality forcing
#flist <- list.files("/media/DATADRIVE1/Project/muresk/driving/forcing_inputs_new/",full.names = TRUE,recursive = TRUE)
#flists <- sub(".csv","",basename(flist))

fnames <- sub(".csv","",basename(flist))

#CUE only sensitivity
flist <- list.files("/media/DATADRIVE1/Project/muresk/driving/forcing_inputs_new/",full.names = TRUE,recursive = TRUE)
flist <- flist[!file_path_sans_ext(basename(flist)) %in% c("5594", "5658")]
flists <- sub(".csv","",basename(flist))

fnames <- sub(".csv","",basename(flist))

outline_df <-readRDS("/media/DATADRIVE1/Project/muresk/model/d02_data/for_model/muresk_outline.RDS")

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
  .packages    = c("doParallel","deSolve","zoo","minpack.lm")
  
) %dopar% {
  
  
  set.seed(123)
  
  driving_input <- read.csv(flist[i])
  dates <- as.Date(with(driving_input, ISOdate(Year, Month, Day)), tz = "UTC")
  # compute doy and drop Feb 29
  doy <- as.POSIXlt(dates)$yday + 1
  keep <- !(format(dates, "%m-%d")=="02-29")
  
  names(driving_input)[5:7] <- c("Npp","SoilMoisture","SoilTemp")
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
  
  sitename      <- sub("\\.csv$", "", basename(flist[i]))
  
  selected_row  <- pars_set[sitename, ]
  
  average_temp  <- mean(agg$forc_st)
  
  sim_output_org       <- NULL
  sim_output_init      <- NULL
  sim_output_init_noCN <- NULL
  
for (npp_ratio in seq(0.5,1.5,0.5)) {
  
  for (CUE in seq(0.2,0.8,0.1)){
    
    for (clay_silt in seq(16, 36, 4)){
      
      for (pH_sens in seq(4, 6.5, 0.5)){
        
        set.seed(123)
        
        temp_npp     <- agg
        
        temp_npp$forc_npp <- with(temp_npp,forc_npp*npp_ratio)
        
        csat_df <- as.data.frame(list("x"  = clay_silt,"y"  = 1.8))
        
        csat_df_target <- get_csat_frombound(outline_df,csat_df)
        
        average_temp <- mean(agg$forc_st)
        min_temp     <- min(agg$forc_st)  # Coldest day - highest CUE_tmp
            
        
        CUE_tmp_avg = selected_row$cue_ref - 0.012 * (average_temp - selected_row$tae_ref)
        #CUE_tmp_avg = 7.765605e-01 - 0.012*(average_temp - 9.985720e0)
        
        # Maximum daily CUE_tmp (at minimum temperature, since cue_t > 0)
        CUE_tmp_max = selected_row$cue_ref - 0.012 * (min_temp - selected_row$tae_ref)
        #CUE_tmp_max = 7.765605e-01 - 0.012* (min_temp - 9.985720e0)
        
        param_CNratio <- compile_muresk_csat[which(compile_muresk_csat$site == sitename),"CNratio"]
        
        #Adj_CN_base <- pmin(1, param_CNratio^(-0.6))
        Adj_CN_base <- param_CNratio^(-0.6)
        delta_tmp   <- CUE_tmp_max - CUE_tmp_avg
        # Always set CUE correction to hit target
        #CUE_corr_org = CUE - CUE_tmp_avg
        
        # Check daily constraint at worst case (coldest day)
        denom_tmp <- delta_tmp + CUE  # = CUE_tmp_max + CUE_corr
        
        CUE_realised <- denom_tmp * Adj_CN_base
        
        if (denom_tmp > 0.8) {
          if (CUE_realised <= 0.8) {
            Adj_CN_target <- CUE/CUE_tmp_avg
            CNratio_effective <- Adj_CN_target^(-1 / 0.6)
            CNratio_corr <- CNratio_effective - param_CNratio
            CUE_corr     <- 0
            CUE_capped <- TRUE
          } else {
            Adj_CN_limit <- 0.8 / CUE_tmp_avg
            CNratio_effective <- Adj_CN_limit^(-1 / 0.6)
            CNratio_corr <- CNratio_effective - param_CNratio
            CUE_corr     <- 0
            CUE_capped <- TRUE
          }
        } else if (denom_tmp > 0) {
          Adj_CN_target <- CUE/CUE_tmp_avg
          CNratio_effective <- Adj_CN_target^(-1 / 0.6)
          CNratio_corr <- CNratio_effective - param_CNratio
          CUE_corr     <- 0
          CUE_capped <- TRUE
        } else {
          CNratio_corr <- 0
          CUE_corr     <- 0
          CUE_capped <- FALSE
        }
        
        #out_mix_org <- safe_ode_mix_fn(selected_row, sitename, temp_npp, CUE_corr = CUE_corr_org, CNratio_corr, use_CN_adj = FALSE)
        #  
        #out_mix_org <- c(out_mix_org, Npp_ratio = npp_ratio)
        #out_mix_org <- c(out_mix_org, CUE_corr = CUE_corr_org)      # CUE correction applied
        #out_mix_org <- c(out_mix_org, CNratio_corr = CNratio_corr)  # CNratio correction applied
        #out_mix_org <- c(out_mix_org, CUE_target = CUE)             # Target CUE from loop
        #out_mix_org <- c(out_mix_org, CUE_capped = CUE_capped)      # Was CUE capped at 0.8?
        #  
        #sim_output_org <- rbind(sim_output_org, out_mix_org)
        #  
        #out_mix <- safe_ode_mix_fn(selected_row, sitename, temp_npp, CUE_corr, CNratio_corr, use_CN_adj = TRUE)
        out_mix <- safe_ode_all_fn(selected_row, sitename, temp_npp, qmax_sens = csat_df_target$ycsat, pH_sens = pH_sens, CUE_corr, CNratio_corr, use_CN_adj = TRUE)
        
        #  
        out_mix <- c(out_mix, Npp_ratio = npp_ratio)
        out_mix <- c(out_mix, CUE_corr = CUE_corr)           # CUE correction applied
        out_mix <- c(out_mix, CNratio_corr = CNratio_corr)   # CNratio correction applied
        out_mix <- c(out_mix, ratio_constraints = (CNratio_corr + param_CNratio)^(-0.6))      # Was CNratio real?
        #out_mix <- c(out_mix, CUE_target = CUE)             # Target CUE from loop
        out_mix <- c(out_mix, CUE_capped = CUE_capped)       # Was CUE capped at 0.8?
        #out_mix <- c(out_mix, pH = pH_sens)                  # pH
        out_mix <- c(out_mix, ClaySilt = clay_silt)          # clay_silt
        #  
        sim_output_init <- rbind(sim_output_init, out_mix)
      
        # Run sensitivity without CN adjustments (CNratio_corr forced to 0, use_CN_adj = FALSE)
        
        #CUE_corr = CUE/Adj_CN_base - CUE_tmp_avg
        
        #out_mix_noCN <- safe_ode_mix_fn(selected_row, sitename, temp_npp, CUE_corr = CUE_corr, CNratio_corr = 0, use_CN_adj = TRUE)
        
        #  out_mix_noCN <- c(out_mix_noCN, Npp_ratio = npp_ratio)
        #  out_mix_noCN <- c(out_mix_noCN, CUE_corr = CUE_corr)
        #  out_mix_noCN <- c(out_mix_noCN, CNratio_corr = 0)
        #  out_mix_noCN <- c(out_mix_noCN, CUE_target = CUE)
        #  out_mix_noCN <- c(out_mix_noCN, CUE_capped = FALSE)
        
        #sim_output_init_noCN <- rbind(sim_output_init_noCN, out_mix_noCN)
      
      
        }
    
      }
  
    }
  
  }
  #sim_output_org <- na.omit(sim_output_org)
  #sens_org <- data.frame(scale(sim_output_org[,1:6]))
  #sens_org$CUE            <- sim_output_org[,7]   # Realized CUE from model
  #sens_org$Npp_ratio      <- sim_output_org[,8]
  #sens_org$CUE_corr       <- sim_output_org[,9]   # CUE correction applied
  #sens_org$CNratio_corr   <- sim_output_org[,10]  # CNratio correction applied
  #sens_org$CUE_target     <- sim_output_org[,11]  # Target CUE from loop
  #sens_org$CUE_capped     <- sim_output_org[,12]  # Was CUE capped at 0.8?
  #colnames(sens_org)[1:6] <- c('SOM1000','MAOM1000','POM1000',
  #                         'AGG1000','MIC1000','LMWC1000')
  #sens_org$site <- sitename
  #sens_org$analysis <- "Org"
  
  
  ## scale & clean
  sim_output_init <- na.omit(sim_output_init)
  #sens <- data.frame(scale(sim_output_init[,1:6]))
  #sens$CUE            <- sim_output_init[,7]   # Realized CUE from model
  #sens$Qmax           <- sim_output_init[,8]
  #sens$pH             <- sim_output_init[,9]
  #sens$Npp_ratio      <- sim_output_init[,10]
  #sens$CUE_corr       <- sim_output_init[,11]   # CUE correction applied
  #sens$CNratio_corr   <- sim_output_init[,12]  # CNratio correction applied
  #sens$ratio_constraints <- sim_output_init[,13]
  ##sens$CUE_target     <- sim_output_init[,11]  # Target CUE from loop
  #sens$CUE_capped     <- sim_output_init[,14]  # Was CUE capped at 0.8?
  #sens$ClaySilt       <- sim_output_init[,15]  # Was CUE capped at 0.8?
  #colnames(sens)[1:6] <- c('SOM1000','MAOM1000','POM1000',
  #                         'AGG1000','MIC1000','LMWC1000')
  #sens$site <- sitename
  #sens$analysis <- "All_adjusted"
  sim_output_init          <- data.frame(sim_output_init)
  sim_output_init$site     <- sitename
  sim_output_init$analysis <- "All_adjusted"
  
  #sim_output_init_noCN <- na.omit(sim_output_init_noCN)
  #sens_noCN <- data.frame(scale(sim_output_init_noCN[,1:6]))
  #sens_noCN$CUE            <- sim_output_init_noCN[,7]
  #sens_noCN$Npp_ratio      <- sim_output_init_noCN[,8]
  #sens_noCN$CUE_corr       <- sim_output_init_noCN[,9]
  #sens_noCN$CNratio_corr   <- sim_output_init_noCN[,10]
  #sens_noCN$CUE_target     <- sim_output_init_noCN[,11]
  #sens_noCN$CUE_capped     <- sim_output_init_noCN[,12]
  #colnames(sens_noCN)[1:6] <- c('SOM1000','MAOM1000','POM1000',
  #                              'AGG1000','MIC1000','LMWC1000')
  #sens_noCN$site <- sitename
  #sens_noCN$analysis <- "CN_fixed"
  
  #sens_test <- rbind(sens_org, sens, sens_noCN)
  #rbind(sens, sens_noCN)
  #return(sens)
  return(sim_output_init)
  
}

stopCluster(cl)

new_sim_outs <- do.call(rbind, sim_col)

# Convenience data frames for downstream analysis
#new_sim_outs_CN_fixed <- subset(new_sim_outs, analysis == "CN_fixed")
#new_sim_outs_all <- subset(new_sim_outs, analysis == "All_adjusted")

saveRDS(new_sim_outs,"/media/DATADRIVE1/Project/muresk/model/d04_output/All_interactions_1118.rds")