#' Get a texture class based on clay content
#'
#' Three texture classes are available: coarse, medium, and fine textures.
#' Select a texture class based on clay content.
#'
#' @param x The clay content of soil (0-1)
#' @export
get_texture <- function(x) {
  texture <- ifelse(x < 0.2, "coarse", ifelse(x >= 0.2 & x < 0.5, "medium", "fine"))
  return(texture)
}

#' Calculate lower limit water (mm)
#'
#' Calculate lower limit water (mm) in the top 30 cm or user-defined layer.
#'
#' @param x The clay content of soil (0-1)
#' @param d A depth (cm), 30 cm by default
#' @export
calc_lower_limit <- function(x, d = 30) {
  (0.0041*x*100 + 0.0445)*d
}

#' Calculate the maximum TSMD under a crop
#'
#' Calculate the maximum topsoil mositure deficity (TSMD) under a crop.
#'
#' @param x The clay content of soil (0-1)
#' @param d A depth (cm), 30 cm by default
#' @export
calc_TSMD_max <- function(x, d = 30) {
  (20 + 130*x - 100*x*x)*d/23
}

#' Calculate a saturated water content
#'
#' Calculate a saturated water content.
#'
#' @param d A depth (cm)
#' @param BD The bulk density of soil (g per cubic cm)
#' @param lim The lower limit water of soil (mm)
#' @param TSMD Topsoil mositure deficity
#' @export
calc_sat_water <- function(d = 30, BD, lim, TSMD) {
  sat_water <- ifelse((1 - BD/2.65)*d < TSMD + lim, lim + TSMD, (1 - BD/2.65)*d)
  return(sat_water)
}

#' Calculate the maximum infiliteration (mm/d) of soil
#'
#' Calculate the maximum infiliteration (mm/d) for soil.
#'
#' @param x A texture class based on clay content
#' @export
get_infil_max <- function(x) {
  fine_infil_max <- 24 # by default (mm/d)
  medium_infil_max <- 120
  coarse_infil_max <- 200
  infil_max <- ifelse(x == "fine", fine_infil_max,
                      ifelse(x == "medium", medium_infil_max,
                             coarse_infil_max)) # if "coarse"
  return(infil_max)
}

#' A fallow category based on previous 7 day's rainfall (mm)
#'
#' Determine a fallow category based on previous 7 day's rainfall (mm).
#'
#' @param x Previous 7 day's rainfall (mm)
#' @export
get_fallow_cat <- function(x) {
  fallow_cat <- ifelse(x > 25, 1,
                       ifelse(x <= 25 & x > 10, 2,
                              ifelse(x <= 10 & x > 1, 3,
                                     4)))
  return(fallow_cat)
}

#' Get bare soil evaporation
#'
#' Get bare soil evaporation based on a fallow category ("in_general").
#'
#' @param x A fallow category that should be 1, 2, 3 or 4
#' @export
get_bare_soil_evap <- function(x) {
  # fraction of pan_evap achieved
  if (x == 1) bs_evap <- 0.4 # >25 mm rain 7 days
  if (x == 2) bs_evap <- 0.3 # from 10 to 25 mm rain 7 days
  if (x == 3) bs_evap <- 0.2 # 1 to 10 mm rain 7 days
  if (x == 4) bs_evap <- 0.1 # <1 mm rain 7 days
  return(bs_evap)
}

#' Read all crop parameters by crop
#'
#' Read all crop parameters - not calibrated yet. The data include 56 different crops.
#'
#' @param x A specific crop name
#' @export
get_crop_params <- function(x) {
  params <- rothc::crop_params_all
  crop <- params[params$crop == x, ] # from "crop_params.csv"
  return(crop)
}

#' Scale a selected crop parameters
#'
#' Scale a selected crop parameters. Do not use yet!
#'
#' @param x A specific crop name
#' @export
check_crop_scaleoutput <- function(x) {
  params <- rothc::crop_scale_all
  scale <- params[params$crop == x, "scaleoutput"] # from "crop_scaleoutput.csv"
  crop_scale <- as.character(scale) # 'y' or 'n'
  return(crop_scale)
}

#' Set automated irrigation
#'
#' @param crop A list of crop parameters
#' @param site A site
#' @export
auto_irrigation <- function(crop, site) {

  IR_min_PAWC_frac <- crop$IR_min_PAWC_frac # (6)
  IR_max_PAWC_frac <- crop$IR_max_PAWC_frac # (7)
  IR_irr_PAWC_frac <- IR_max_PAWC_frac - IR_min_PAWC_frac
  IR_irr_amount <- IR_irr_PAWC_frac*(site$awc*crop$plant_PAWC_frac) # based on plant available water capacity to 1 m

  texture <- get_texture(x = site$cly) # this is clay fraction
  infil_max <- get_infil_max(x = texture) # maximum infiltration for selected soil (mm/d)

  # where maximum infiltration is low and the difference between the minimum and maximum PAWC fraction is high,
  # it may be necessary to irrigate on more than 5 subsequent days!
  irr_d <- floor(IR_irr_amount/infil_max + 1)

  no <- NULL
  for (i in seq(irr_d)) {
    s <- paste("D", i, " irrigation", sep = "")
    no <- c(no, s)
  }

  no <- as.list(c(no, "no irrigation"))
  event <- as.list(c(seq(1:irr_d), 0)) # event #
  amount <- as.list(rep(0, irr_d + 1)) # irrigation (mm)
  autoirri <- do.call(rbind, Map(data.frame, "no" = no, "event" = event, "amount" = amount))

  if (IR_irr_amount < infil_max) {
    autoirri$amount[1] <- IR_irr_amount
  } else { autoirri$amount[1] <- infil_max }

  irr_sum <- autoirri$amount[1]
  for (i in 2:5) {
    if (irr_sum < IR_irr_amount) {
      if (IR_irr_amount - irr_sum < infil_max) {
        autoirri$amount[i] <- IR_irr_amount - irr_sum
      } else {autoirri$amount[i] <- infil_max}
      irr_sum <- irr_sum + autoirri$amount[i]
    } else { break }
  }

  return(autoirri)
}

#' Set an estimate of the decomposability of the incoming plant material
#'
#' Incoming plant C is split between DPM and RPM based on a DPM/RPM ratio.
#' Three values of the ratio are permitted in the model version 26.3:
#' 1.44 for most agricultural crops and improved grassland,
#' 0.67 for unimproved grassland and scrub, including Savanna, and
#' 0.25 for a deciduous or tropical woodland.
#'
#' @param x A specific land use
#' @export
set_ratio_DPM_to_RPM <- function(x) {
    if (x == 0) DPM_RPM_ratio <- 0.67 # 0 = other (quick and dirty)
    if (x == 1) DPM_RPM_ratio <- 0.25 # 1 = native woodland or forest
    if (x == 2) DPM_RPM_ratio <- 0.67 # 2 = unimproved pastures
    if (x == 3) DPM_RPM_ratio <- 1.44 # 3 = improved pastures
    if (x == 4) DPM_RPM_ratio <- 0.25 # 4 = production woodland or forest
    if (x == 5) DPM_RPM_ratio <- 1.44 # 5 = cropping, except for grazing modified pastures
    return(DPM_RPM_ratio)
  }

#' The month by day
#'
#' Find the month by day of normal year
#'
#' @param x A day of the year
#' @export
which_month <- function(x) {
  mon <- NULL
  if (x <= 31) mon <- 1
  if (x > 31 & x <= 59) mon <- 2 # 60 for leap year
  if (x > 59 & x <= 90) mon <- 3
  if (x > 90 & x <= 120) mon <- 4
  if (x > 120 & x <= 151) mon <- 5
  if (x > 151 & x <= 181) mon <- 6
  if (x > 181 & x <= 212) mon <- 7
  if (x > 212 & x <= 243) mon <- 8
  if (x > 243 & x <= 273) mon <- 9
  if (x > 273 & x <= 304) mon <- 10
  if (x > 304 & x <= 334) mon <- 11
  if (x > 334) mon <- 12
  mon
}

#' Calculate the sum for previous k rows
#'
#' Calculate the sum for previous k rows.
#'
#' @param x A list of values
#' @param k Rows
#' @export
sum2 <- function(x, k) {
  v <- NULL
  n <- length(x)
  for (i in 1:n) {
    if (i <= k) {
      s <- sum(x[1:i])
    } else { # if i > k
      s <- sum(x[(i - k):(i - 1)])
    }
    v <- c(v, s)
  }
  return(v)
}

#' Create a directory recursively
#'
#' This function is re-written from "dir.create" to avoid any permission issue
#' on the network when recursively creating the directory.
#'
#' @param x A path
#' @export
dir.create.rev <- function(x) {

  subdir <- sub(".*/Soil_MMM-VISCAR-SE06315", "", x)
  home <- sub(subdir, "", x)

  d <- unlist(strsplit(subdir, "/"))
  for (i in seq(d)) {
    newsubdir <- paste0(home, d[i], "/")
    home <- newsubdir
    dir.create(newsubdir)
  }
}


# Get a temperature index for the selected perennial crop
get_pe_TI <- function(x, crop) {
  pe_TI_b_sub <- crop$pe_TI_b_sub     # changes inflexion point of TI curve (27)
  pe_TI_b_supra <- crop$pe_TI_b_supra # (28)
  
  if (any(x < 0.5)) {
    pe_TI <- 1 - ((2*x)^pe_TI_b_sub/2)
  } else {
    pe_TI <- 0.5*(2*(1 - x))^pe_TI_b_supra
  }
  pe_TI
}

