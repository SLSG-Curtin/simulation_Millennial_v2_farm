delta_MM_v2 <- function(fraction_values,params_sceua,forc_st,forc_sw,forc_npp) {
  
  default_params <- list(
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
    param_qmax = 20000,
    param_pH   = 6
  )
  
  params_sceua <- modifyList(default_params, params_sceua)
  
  #with(as.list(c(fraction_values,params_sceua)), {
  with(c(fraction_values,params_sceua), {
    # Soil type properties  
    #Equation 11
    #Mayes 2012, SSAJ
    kaff_lm = exp(-param_p1 * param_pH - param_p2) * kaff_des
    
    #Equation 12
    #Georgiou in review
    #param_qmax = param_bulkd * param_c1 * param_claysilt 
    param_qmax  = param_qmax
    
    # Hydrological properties
    
    #Equation 5
    scalar_wd = (forc_sw / porosity)^0.5
    
    #Equation 4
    
    if(porosity > forc_sw){
      scalar_wb = exp(lambda * -matpot) * (kamin + (1 - kamin) * ((porosity - forc_sw) / porosity)^0.5) * scalar_wd
    }else{
      scalar_wb = exp(lambda * -matpot) * (kamin + (1 - kamin) * 0.5) * scalar_wd
    }
    
    
    # Decomposition
    
    gas_const <- 8.31446
    
    #Equation 3
    vmax_pl = alpha_pl * exp(-eact_pl / (gas_const * (forc_st + 273.15)))
    
    #Equation 2
    # POM -> LMWC
    if(POM>0 && MIC>0){
      f_PO_LM = vmax_pl * scalar_wd * POM * MIC / (kaff_pl + MIC)
    }else{
      f_PO_LM=0
    }
    
    #Equation 6
    # POM -> AGG
    if(POM>0){
      f_PO_AG = rate_pa * scalar_wd * POM
    }else{
      f_PO_AG=0
    }
    
    #Equation 7
    # AGG -> MAOM
    if(AGG>0){
      f_AG_break = rate_break * scalar_wd * AGG
    }else{
      f_AG_break=0
    }
    
    #Equation 9
    # LMWC -> out of system leaching
    if(LMWC>0){
      f_LM_leach = rate_leach * scalar_wd * LMWC
    }else{
      f_LM_leach=0
    }
    
    #Equation 10
    # LMWC -> MAOM
    if(LMWC>0 && MAOM>0){
      f_LM_MA = scalar_wd * kaff_lm * LMWC * (1 - MAOM / param_qmax)
    }else{
      f_LM_MA=0
    }
    
    #Equation 13
    # MAOM -> LMWC
    if(MAOM>0){
      f_MA_LM = kaff_des * MAOM / param_qmax
    }else{
      f_MA_LM=0
    }
    
    #Equation 15
    vmax_lb = alpha_lb * exp(-eact_lb / (gas_const * (forc_st + 273.15)))
    
    #Equation 14
    # LMWC -> MIC
    if(LMWC>0 && MIC>0){
      f_LM_MB = vmax_lb * scalar_wb * MIC * LMWC / (kaff_lb + LMWC)
    }else{
      f_LM_MB=0
    }
    
    #Equation 16
    # MIC -> MAOM/LMWC
    if(MIC>0){
      f_MB_turn = rate_bd * MIC^2.0
    }else{
      f_MB_turn=0
    }
    
    #Equation 18
    # MAOM -> AGG
    if(MAOM>0){  
      f_MA_AG = rate_ma * scalar_wd * MAOM
    }else{
      f_MA_AG=0
    }
    
    #Equation 22
    # microbial growth flux, but is not used in mass balance
    
    #Equation 21
    # MIC -> atmosphere
    if(MIC>0 && LMWC>0){ 
      f_MB_atm = f_LM_MB * (1 - (cue_ref - cue_t * (forc_st - tae_ref) ) )
    }else{
      f_MB_atm=0
    }
    
    # Update state variables
    
    #Equation 1
    dPOM = forc_npp * param_pi + f_AG_break * param_pa - f_PO_AG - f_PO_LM
    
    #Equation 8
    dLMWC = forc_npp * (1. - param_pi) - f_LM_leach + f_PO_LM - f_LM_MA - f_LM_MB + f_MB_turn * (1. - param_pb) + f_MA_LM
    
    #Equation 17
    dAGG = f_MA_AG + f_PO_AG - f_AG_break
    
    #Equation 19
    dMAOM = f_LM_MA - f_MA_LM + f_MB_turn * param_pb - f_MA_AG + f_AG_break * (1. - param_pa)
    
    #Equation 20
    dMIC = f_LM_MB - f_MB_turn - f_MB_atm
    
    return(list(c(dPOM, dLMWC, dAGG, dMIC, dMAOM)))
  })
}
