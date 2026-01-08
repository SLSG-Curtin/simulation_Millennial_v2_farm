# predict soil temperature at depth

# step0: read in the weather


xy_new_cp <- read.csv("/media/DATADRIVE1/Project/muresk/croptype/geo_croptype.csv")

weatherID <- xy_new_cp[2, c("SILO")]

weather <- read.csv(paste0("/media/DATADRIVE1/Project/muresk/model/d02_data/SILO/", weatherID, ".csv"))

# step1: Find the appropriate parameters in Table 3 for the soil  temperature to be measured,
# for air minimum, air maximum, and log(rain + 1).

# step2: Calculate alpha_i, beta_i, and rho_i using the formulae in Table 2 for each of the models specified in Table 3.

# step3: Calculate Nd , Md , and Rd using the Eqns 5, 6, and 7

# step4: Calculate the coefficients using Table 4 and Eqn 9.

# step5: Calculate the soil temperature using Eqn 8

# non-linear model a
# model_A = (a/depth)*(((b*depth)/(1+c*depth))^Lag_days)

# non-linear model b
# model_B = a*(((b*depth)/(1+c*depth))^Lag_days)

# non-linear model c
# model_C = (a/depth)*((c+b*depth)^Lag_days)

# non-linear model d
# model_D = a*((c+b*depth)^Lag_days)


# Eqn 5
calculate_weightings_mintemp_with_yearly_substitution <- function(varible_df, row, a, b, c, depth) {
  # varible_df should have columns: date, temp
  # date should be in Date format

  varible_df$date <- as.Date(varible_df$date)
  total_days <- nrow(varible_df)
  N_days <- numeric(total_days)

  # Helper function to find substitute dates for missing historical data
  find_substitute_min <- function(target_date, varible_df) {
    # Extract month and day from target date
    target_month_day <- format(target_date, "%m-%d")

    # Find all dates in the dataset with the same month-day
    same_month_day <- varible_df[format(varible_df$date, "%m-%d") == target_month_day, ]

    if (nrow(same_month_day) > 0) {
      # Return the mean of all available years for this date
      return(mean(same_month_day$min, na.rm = TRUE))
    } else {
      return(NA)
    }
  }


  sum_upper <- 0
  sum_lower <- 0
  days <- varible_df[row, "julianday"]
  current_date <- varible_df$date[days]


  for (j in 0:28) {
    # model B
    alpha_j <- a * (((b * depth) / (1 + c * depth))^j)

    # Calculate target date (i days before current date)
    target_date <- current_date - j

    # First, try to find exact match in the dataset
    target_row <- which(varible_df$date == target_date)

    if (length(target_row) > 0) {
      # Direct match found
      AirTemp_days_j <- varible_df$min[target_row[1]]
    } else {
      # No direct match - use substitute from corresponding dates in other years
      AirTemp_days_j <- find_substitute_min(target_date, varible_df)
    }

    # Only include if we have valid temperature data
    if (!is.na(AirTemp_days_j)) {
      sum_upper <- sum_upper + alpha_j * AirTemp_days_j
      sum_lower <- sum_lower + alpha_j
    }
  }

  # Calculate weighted average
  if (sum_lower > 0) {
    N_days <- sum_upper / sum_lower
  } else {
    N_days <- NA
  }

  return(N_days)
}

# Eqn 6
calculate_weightings_maxtemp_with_yearly_substitution <- function(varible_df, row, a, b, c, depth) {
  # varible_df should have columns: date, temp
  # date should be in Date format

  varible_df$date <- as.Date(varible_df$date)
  total_days <- nrow(varible_df)
  M_days <- numeric(total_days)

  # Helper function to find substitute dates for missing historical data
  find_substitute_max <- function(target_date, varible_df) {
    # Extract month and day from target date
    target_month_day <- format(target_date, "%m-%d")

    # Find all dates in the dataset with the same month-day
    same_month_day <- varible_df[format(varible_df$date, "%m-%d") == target_month_day, ]

    if (nrow(same_month_day) > 0) {
      # Return the mean of all available years for this date
      return(mean(same_month_day$max, na.rm = TRUE))
    } else {
      return(NA)
    }
  }

  sum_upper <- 0
  sum_lower <- 0
  days <- varible_df[row, "julianday"]
  current_date <- varible_df$date[days]

  for (j in 0:28) {
    # model C
    beta_j <- (a / depth) * ((c + b * depth)^j)

    # Calculate target date (i days before current date)
    target_date <- current_date - j

    # First, try to find exact match in the dataset
    target_row <- which(varible_df$date == target_date)

    if (length(target_row) > 0) {
      # Direct match found
      AirTemp_days_j <- varible_df$max[target_row[1]]
    } else {
      # No direct match - use substitute from corresponding dates in other years
      AirTemp_days_j <- find_substitute_max(target_date, varible_df)
    }

    # Only include if we have valid temperature data
    if (!is.na(AirTemp_days_j)) {
      sum_upper <- sum_upper + beta_j * AirTemp_days_j
      sum_lower <- sum_lower + beta_j
    }
  }

  # Calculate weighted average
  if (sum_lower > 0) {
    M_days <- sum_upper / sum_lower
  } else {
    M_days <- NA
  }

  return(M_days)
}

# Eqn 7
calculate_weightings_rain_with_yearly_substitution <- function(varible_df, row, a, b, c, depth) {
  # varible_df should have columns: date, rain
  # date should be in Date format

  varible_df$date <- as.Date(varible_df$date)
  total_days <- nrow(varible_df)
  R_days <- numeric(total_days)

  # Helper function to find substitute dates for missing historical data
  find_substitute_rain <- function(target_date, varible_df) {
    # Extract month and day from target date
    target_month_day <- format(target_date, "%m-%d")

    # Find all dates in the dataset with the same month-day
    same_month_day <- varible_df[format(varible_df$date, "%m-%d") == target_month_day, ]

    if (nrow(same_month_day) > 0) {
      # Return the mean of all available years for this date
      return(mean(same_month_day$rainfall, na.rm = TRUE))
    } else {
      return(NA)
    }
  }

  # for (days in 1:total_days) {

  sum_upper <- 0
  sum_lower <- 0
  days <- varible_df[row, "julianday"]
  current_date <- varible_df$date[days]

  for (j in 0:28) {
    # model D
    rho_j <- a * ((c + b * depth)^j)

    # Calculate target date (i days before current date)
    target_date <- current_date - j

    # First, try to find exact match in the dataset
    target_row <- which(varible_df$date == target_date)

    if (length(target_row) > 0) {
      # Direct match found
      Rain_days_j <- varible_df$rainfall[target_row[1]]
    } else {
      # No direct match - use substitute from corresponding dates in other years
      Rain_days_j <- find_substitute_rain(target_date, varible_df)
    }

    # Only include if we have valid temperature data
    if (!is.na(Rain_days_j)) {
      sum_upper <- sum_upper + rho_j * log(Rain_days_j + 1)
      sum_lower <- sum_lower + rho_j
    }
  }

  # Calculate weighted average
  R_days <- sum_upper / sum_lower

  return(R_days)
}

# improved version

# Define coefficient tables
# Table 1: Coefficients for N_d, M_d, R_d calculations
weighting_coefficients <- data.frame(
  Variable = c("Air_minimum", "Air_maximum", "Rain"),
  Model = c("Average", "Average", "Average"),
  a = c(0.087, 1.648, -0.244),
  b = c(0.691, 0.0056, 0.0002),
  c = c(0.959, 0.432, 0.879)
)

# Table 2: Coefficients for Equation 8
equation8_coefficients <- data.frame(
  Parameter = c(
    "Intercept", "Latitude", "Interaction", "Rain",
    "Air_maximum", "Air_minimum", "Phase", "Seasonality"
  ),
  A = c(3.6586, 0.1451, -0.3313, 6.2933, 0.8751, 0.0927, 78.8887, 4.0215),
  B = c(-0.0637, -0.0007, 0.0004, 0.0182, 0.0036, -0.0046, -0.4984, -0.0121),
  C = c(4.1486, 0.0501, 0.0305, -1.1385, -0.1608, 0.1215, 2.0673, -0.1989)
)

# Initialize results vector
T_days_depth <- numeric(total_days)

# get_weighting_coeffs <- function(variable_name) {
#  coeffs <- weighting_coefficients[weighting_coefficients$Variable == variable_name, ]
#  return(list(a = coeffs$a, b = coeffs$b, c = coeffs$c))
# }
#
# get_equation8_coeff <- function(parameter_name, depth) {
#  coeffs <- equation8_coefficients[equation8_coefficients$Parameter == parameter_name, ]
#  return(coeffs$A + depth * coeffs$B + coeffs$C * log(depth))
# }

# Main calculation loop


# variable_df <- weather[,c(1:7)]
# varible_df$date <- with(varible_df,ISOdate(year, month, day))


xy_new_cp <- read.csv("/media/DATADRIVE1/Project/muresk/croptype/geo_croptype.csv")

num_sites <- nrow(xy_new_cp)

library(doParallel)

cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)


foreach::foreach(i = seq(num_sites)) %dopar% {
  site_lat <- xy_new_cp[i, c("y")]
  weatherID <- xy_new_cp[i, c("SILO")]
  site_id <- xy_new_cp[i, c("site")]
  weather <- read.csv(paste0("/media/DATADRIVE1/Project/muresk/model/d02_data/SILO/", weatherID, ".csv"))
  varible_df <- weather[, c(1:7)]
  varible_df$date <- with(varible_df, ISOdate(year, month, day))


  for (depth in c(5, 10, 20, 50, 100)) {
    T_days_depth <- numeric(9132)

    for (row in 1:nrow(varible_df)) {
      # Get coefficients for N_d (Air minimum temperature)
      air_min_coeffs <- weighting_coefficients[weighting_coefficients$Variable == "Air_minimum", ]
      a <- air_min_coeffs$a
      b <- air_min_coeffs$b
      c <- air_min_coeffs$c

      N_d <- calculate_weightings_mintemp_with_yearly_substitution(varible_df, row, a, b, c, depth)

      # Get coefficients for M_d (Air maximum temperature)
      air_max_coeffs <- weighting_coefficients[weighting_coefficients$Variable == "Air_maximum", ]
      a <- air_max_coeffs$a
      b <- air_max_coeffs$b
      c <- air_max_coeffs$c

      M_d <- calculate_weightings_maxtemp_with_yearly_substitution(varible_df, row, a, b, c, depth)

      # Get coefficients for R_d (Rain)
      rain_coeffs <- weighting_coefficients[weighting_coefficients$Variable == "Rain", ]
      a <- rain_coeffs$a
      b <- rain_coeffs$b
      c <- rain_coeffs$c

      R_d <- calculate_weightings_rain_with_yearly_substitution(varible_df, row, a, b, c, depth)

      # Calculate coefficients for Equation 8 using Table 2
      # Extract coefficients for each parameter
      intercept_coeffs <- equation8_coefficients[equation8_coefficients$Parameter == "Intercept", ]
      latitude_coeffs <- equation8_coefficients[equation8_coefficients$Parameter == "Latitude", ]
      interaction_coeffs <- equation8_coefficients[equation8_coefficients$Parameter == "Interaction", ]
      rain_coeffs_eq8 <- equation8_coefficients[equation8_coefficients$Parameter == "Rain", ]
      air_max_coeffs_eq8 <- equation8_coefficients[equation8_coefficients$Parameter == "Air_maximum", ]
      air_min_coeffs_eq8 <- equation8_coefficients[equation8_coefficients$Parameter == "Air_minimum", ]
      phase_coeffs <- equation8_coefficients[equation8_coefficients$Parameter == "Phase", ]
      seasonality_coeffs <- equation8_coefficients[equation8_coefficients$Parameter == "Seasonality", ]

      # Equation 9: Calculate coefficients for Equation 8
      Intercept <- intercept_coeffs$A + depth * intercept_coeffs$B + intercept_coeffs$C * log(depth)
      Latitude <- latitude_coeffs$A + depth * latitude_coeffs$B + latitude_coeffs$C * log(depth)
      Interaction <- interaction_coeffs$A + depth * interaction_coeffs$B + interaction_coeffs$C * log(depth)
      Rain <- rain_coeffs_eq8$A + depth * rain_coeffs_eq8$B + rain_coeffs_eq8$C * log(depth)
      Air_max <- air_max_coeffs_eq8$A + depth * air_max_coeffs_eq8$B + air_max_coeffs_eq8$C * log(depth)
      Air_min <- air_min_coeffs_eq8$A + depth * air_min_coeffs_eq8$B + air_min_coeffs_eq8$C * log(depth)
      Phase <- phase_coeffs$A + depth * phase_coeffs$B + phase_coeffs$C * log(depth)
      Seasonality <- seasonality_coeffs$A + depth * seasonality_coeffs$B + seasonality_coeffs$C * log(depth)

      # Equation 8: Calculate soil temperature

      days <- varible_df[row, "julianday"]

      T_days_depth[row] <- Intercept + Latitude * site_lat + sin(days * Seasonality + Phase) +
        Air_min * N_d + Air_max * M_d + Rain * R_d + Interaction * M_d * R_d
    }

    newcolumn <- paste0("soiltemp_", depth)

    varible_df[newcolumn] <- T_days_depth
  }

  varible_df <- varible_df[, c(1, 2, 3, 4, 9, 10, 11, 12, 13)]
  write.csv(varible_df, paste0("/media/DATADRIVE1/Project/muresk/soiltemp/Inputs/", site_id, ".csv"))
}



calc_avg_soil_temp <- function(temps) {
  depths <- c(5, 10, 20, 50)

  # Extrapolate T(0)
  slope_0_5 <- (temps[2] - temps[1]) / (depths[2] - depths[1])
  T0 <- temps[1] - slope_0_5 * depths[1]

  # Interpolate T(30)
  slope_20_50 <- (temps[4] - temps[3]) / (depths[4] - depths[3])
  T30 <- temps[3] + slope_20_50 * (30 - depths[3])

  # Trapezoidal integral
  integral <- (5 * (T0 + temps[1]) / 2) +
    (5 * (temps[1] + temps[2]) / 2) +
    (10 * (temps[2] + temps[3]) / 2) +
    (10 * (temps[3] + T30) / 2)
  return(integral / 30)
}

f <- list.files(path = "/media/DATADRIVE1/Project/muresk/soiltemp/Inputs/", pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

foreach::foreach(i = seq(f), .combine = "cbind", .packages = "raster") %dopar% {
  b <- read.csv(f[i])
  b_base <- basename(f[i])
  temp_data <- b[, c(6:9)]
  temp_data <- data.frame(temp_data)
  temp_data$soiltemp_030 <- apply(temp_data, 1, calc_avg_soil_temp)
  temp_outs <- cbind(b[, c(2:5)], temp_data)
  write.csv(temp_outs, paste0("/media/DATADRIVE1/Project/muresk/soiltemp/Outputs/", b_base))
}


stopCluster(cl)


# Soil Moisture
################################################################################

calc_daily_sw_update <- function(x, crop, site) {
  # from the site soil data
  depth <- 30 # depth of soil layer (cm)
  bulk_density <- site$BD
  clay <- site$clay # clay fraction - can be used to calculate plant available water content (PAWC)
  soil_PAWC <- site$PAWC # PAWC to 1 m depth (mm), from the map

  texture <- get_texture(x = clay)
  lower_limit <- calc_lower_limit(x = clay, d = depth) # lower limit water (mm)
  TSMD_init <- TSMD_max <- calc_TSMD_max(x = clay, d = depth) # set to be equal to maximum top soil mositure deficity
  # assuming that the soil drys out during summer
  sat_water <- calc_sat_water(d = depth, BD = bulk_density, lim = lower_limit, TSMD = TSMD_max) # in mm
  pot_drain <- sat_water - lower_limit - TSMD_max # in mm
  infil_max <- get_infil_max(x = texture) # maximum infiltration for selected soil (mm/d)

  # from the crop parameters
  crop_type <- crop$season # season type - summer, winter or perennial (2)
  veg_type <- crop$lifeform # vegetation type - annual or perennial	(3)
  Eo_frac <- crop$Eo_frac # maximum Et as fraction of pan evaporation (8)
  plant_PAWC_frac <- crop$plant_PAWC_frac # fraction (0-1) (14)
  PAWC <- soil_PAWC * plant_PAWC_frac # use clay content to calculate it?
  pot_deep_PAW <- (PAWC - TSMD_max)
  if (pot_deep_PAW < 0) {
    pot_deep_PAW <- 0
  }

  swc_init <- 10 # in mm new for soil-water calculation

  sow_start <- crop$sow_start
  if (is.na(sow_start)) sow_start <- 0 # (15)
  harvest <- crop$harvest
  if (is.na(harvest)) harvest <- 0 # (16) should check

  SW <- x # start from the daily weather
  year <- SW$year
  julianday <- SW$julianday
  panevap <- SW$panevap
  rainfall <- SW$rainfall

  vars <- c(
    "winter_crop_r", "summer_crop_r", "crop_r", "irrevent_r", "irrigation_r",
    "TSMD_init_r", "TSMD_mod_r", "fallow_cat_r", "ET_actual_r", "net_water_r",
    "infiltration_r", "TSMD_acc_r", "excess_water_r", "drainage_r",
    "runoff_r", "drainage_frac_r", "TSMD_final_r", "PAW_shallow_r",
    "fallow_deep_r", "swc_init_r", "PAWC_frac_r", "swc_final_r", "deep_drain_r"
  )
  SW[vars] <- 0 # initialize

  if (crop_type != "winter") {
    SW$winter_crop_r <- 0
  } else {
    SW$winter_crop_r <- ifelse(julianday < sow_start | julianday > harvest, 0, 1)
  }
  if (crop_type != "summer") {
    SW$summer_crop_r <- 0
  } else {
    # SW$summer_crop_r <- ifelse(julianday < sow_start, 0, 1)
    SW$summer_crop_r <- ifelse(julianday < sow_start & julianday + 500 > harvest, 0, 1)
  }
  SW$crop_r <- 1 # if perennial
  if (crop_type == "winter") {
    SW$crop_r <- SW$winter_crop_r
  }
  if (crop_type == "summer") {
    SW$crop_r <- SW$summer_crop_r
  }
  SW$crop_r <- with(SW, ifelse(year == 1 & julianday < 8 & veg_type == "annual", 0, crop_r))
  SW$crop_r <- with(SW, ifelse(year == 1 & julianday == 1, 0, crop_r))

  SW[c("winter_crop_r", "summer_crop_r", "crop_r")] <- sapply(SW[c("winter_crop_r", "summer_crop_r", "crop_r")], as.integer)

  # categories based on previous 7 day's rainfall
  # used to dictate what fraction of pan_evap is lost during fallow
  # 1: >25 mm, 2: >10 mm and <=25 mm, 3: > 1 mm and <=10 mm, 4: <=1 mm
  SW$rain_7days <- sum2(rainfall, 7)
  SW$fallow_cat_r <- as.integer(with(SW, get_fallow_cat(x = rain_7days)))
  SW$rain_7days <- NULL # remove

  # vectorize
  winter_crop_r <- SW$winter_crop_r
  summer_crop_r <- SW$summer_crop_r
  crop_r <- SW$crop_r # 0 or 1
  irrevent_r <- SW$irrevent_r # no irrigation considered yet
  irrigation_r <- SW$irrigation_r # no irrigation considered yet
  TSMD_init_r <- SW$TSMD_init_r
  TSMD_mod_r <- SW$TSMD_mod_r
  fallow_cat_r <- SW$fallow_cat_r
  ET_actual_r <- SW$ET_actual_r
  net_water_r <- SW$net_water_r
  infiltration_r <- SW$infiltration_r
  TSMD_acc_r <- SW$TSMD_acc_r
  excess_water_r <- SW$excess_water_r
  drainage_r <- SW$drainage_r
  runoff_r <- SW$runoff_r
  drainage_frac_r <- SW$drainage_frac_r
  TSMD_final_r <- SW$TSMD_final_r
  PAW_shallow_r <- SW$PAW_shallow_r
  fallow_deep_r <- SW$fallow_deep_r
  swc_init_r <- SW$swc_init_r
  PAWC_frac_r <- SW$PAWC_frac_r
  swc_final_r <- SW$swc_final_r
  deep_drain_r <- SW$deep_drain_r

  # initial default values for day 1:
  ET_init <- 2
  net_water_init <- rainfall[1] + irrigation_r[1] - ET_init
  infiltration_init <- pmin(net_water_init, infil_max)
  swc_init_start <- pmax(0, pmin(PAWC, infiltration_init + swc_init))

  for (i in 1:nrow(SW)) { # have to run in sequence

    # initial top soil moisture deficit (TSMD) for the day (mm)
    TSMD_init_r[i] <- ifelse(i == 1, TSMD_init, TSMD_final_r[i - 1])

    # a scaled value (0-1) with 0 equivalent to maximum TSMD and 1 equivalent to TSMD = 0
    # used to run nitrification and denitrification equations, modify ET, and calculate decomposition
    TSMD_mod_r[i] <- 1 - TSMD_init_r[i] / TSMD_max

    # adds infiltration to previous day's final soil water content, limits to PAWC (mm)
    swc_init_r[i] <- ifelse(i == 1, swc_init_start, swc_final_r[i - 1])

    # a fraction (0-1) of SWC/PAWC - used to modify ET when a crop is present
    PAWC_frac_r[i] <- swc_init_r[i] / PAWC

    # calculates evapotranspiration (mm) from pan evaporation
    ET_actual_r[i] <- ifelse(i == 1, ET_init, Eo_frac * panevap[i] * pmin(TSMD_mod_r[i], PAWC_frac_r[i]))

    # rainfall minus actual ET (mm)
    net_water_r[i] <- rainfall[i] + irrigation_r[i] - ET_actual_r[i]

    # infilteration rate
    infiltration_r[i] <- pmin(net_water_r[i], infil_max)

    # interim accumulated TSMD
    # can go negative, in which case there is excess water displayed in the next column
    # that can go to drainage or runoff if the top 30 cm is full
    TSMD_acc_r[i] <- pmin(TSMD_init_r[i] - infiltration_r[i], TSMD_max)

    # water in excess of the amount that can exist in the top 30 cm as defined by TSMD_max
    excess_water_r[i] <- pmax(0, -TSMD_acc_r[i])

    # limits excess water that can enter the deep soil to the maximum amount that can be held by the top 30 cm (pot_drainage)
    drainage_r[i] <- pmin(excess_water_r[i], pot_drain)

    # runoff (mm)
    runoff_r[i] <- pmax(net_water_r[i] - infiltration_r[i], 0) + pmax(excess_water_r[i] - pot_drain, 0)

    # fraction of existing water in the top 30 cm (equivalent to maximum TSMD when drainage occurs) that drains
    # used for leaching nitrate
    drainage_frac_r[i] <- pmin(pmax(drainage_r[i] / TSMD_max, 0), 1)

    # final TSMD, ranging between 0 and maximum top soil moisture deficit
    TSMD_final_r[i] <- pmax(0, pmin(TSMD_max, TSMD_acc_r[i]))

    # amount (mm) of water existing in the top 30 cm
    PAW_shallow_r[i] <- TSMD_max - TSMD_final_r[i]

    # during fallow, picks up any drainage that occurs, and evaporates it off if the top 30 cm of soil is dry
    # limits to potential fallow deepwater (calculated as PAWC - TSMD_max).
    if (i == 1) {
      fallow_deep_r[i] <- 0
    } else {
      if (crop_r[i] != 0) {
        temp <- 0
      } else {
        add_drain <- drainage_r[i] + fallow_deep_r[i - 1]
        if (TSMD_final_r[i] < TSMD_max) {
          temp <- add_drain
        } else if (add_drain > 0) {
          evap_bare <- get_bare_soil_evap(x = fallow_cat_r[i]) * panevap[i]
          temp <- add_drain - evap_bare
        } else {
          temp <- 0
        }
      }
      fallow_deep_r[i] <- pmin(pmax(temp, 0), pot_deep_PAW)
    }

    # final soil water content for day by subtracting ET from initial soil water content (infiltration + previous day's final SWC)
    # limits to 0 if soil dries out
    swc_final_r[i] <- pmax(0, pmin(PAWC, swc_init_r[i] + infiltration_r[i]))

    # drainage from the rooting zone (mm)
    deep_drain_r[i] <- ifelse(i == 1, 0, pmax(swc_init_r[i] + infiltration_r[i] - PAWC, 0))
  }

  # update
  SW$TSMD_init_r <- TSMD_init_r
  SW$TSMD_mod_r <- TSMD_mod_r
  SW$swc_init_r <- swc_init_r
  SW$PAWC_frac_r <- PAWC_frac_r
  SW$ET_actual_r <- ET_actual_r
  SW$net_water_r <- net_water_r
  SW$infiltration_r <- infiltration_r
  SW$TSMD_acc_r <- TSMD_acc_r
  SW$excess_water_r <- excess_water_r
  SW$drainage_r <- drainage_r
  SW$runoff_r <- runoff_r
  SW$drainage_frac_r <- drainage_frac_r
  SW$TSMD_final_r <- TSMD_final_r
  SW$PAW_shallow_r <- PAW_shallow_r
  SW$fallow_deep_r <- fallow_deep_r
  SW$swc_final_r <- swc_final_r
  SW$deep_drain_r <- deep_drain_r
  return(SW)
}


calc_daily_sw_grok <- function(x, crop, site) {
  # from the site soil data
  depth <- 30 # depth of soil layer (cm)
  bulk_density <- site$BD
  clay <- site$clay # clay fraction - can be used to calculate plant available water content (PAWC)
  soil_PAWC <- site$PAWC # PAWC to 1 m depth (mm), from the map

  texture <- get_texture(x = clay)
  lower_limit <- calc_lower_limit(x = clay, d = depth) # lower limit water (mm)
  TSMD_init <- TSMD_max <- calc_TSMD_max(x = clay, d = depth) # set to be equal to maximum top soil mositure deficity
  # assuming that the soil drys out during summer
  sat_water <- calc_sat_water(d = depth, BD = bulk_density, lim = lower_limit, TSMD = TSMD_max) # in mm
  pot_drain <- sat_water - lower_limit - TSMD_max # in mm
  infil_max <- get_infil_max(x = texture) # maximum infiltration for selected soil (mm/d)

  # from the crop parameters
  crop_type <- crop$season # season type - summer, winter or perennial (2)
  veg_type <- crop$lifeform # vegetation type - annual or perennial	(3)
  Eo_frac <- crop$Eo_frac # maximum Et as fraction of pan evaporation (8)
  plant_PAWC_frac <- crop$plant_PAWC_frac # fraction (0-1) (14)
  PAWC <- soil_PAWC * plant_PAWC_frac # use clay content to calculate it?
  pot_deep_PAW <- (PAWC - TSMD_max)
  if (pot_deep_PAW < 0) {
    pot_deep_PAW <- 0
  }

  swc_init <- 10 # in mm new for soil-water calculation

  sow_start <- crop$sow_start
  if (is.na(sow_start)) sow_start <- 0 # (15)
  harvest <- crop$harvest
  if (is.na(harvest)) harvest <- 0 # (16) should check

  SW <- x # start from the daily weather
  year <- SW$year
  julianday <- SW$julianday
  panevap <- SW$panevap
  rainfall <- SW$rainfall

  vars <- c(
    "winter_crop_r", "summer_crop_r", "crop_r", "irrevent_r", "irrigation_r",
    "TSMD_init_r", "TSMD_mod_r", "fallow_cat_r", "ET_actual_r", "net_water_r",
    "infiltration_r", "TSMD_acc_r", "excess_water_r", "drainage_r",
    "runoff_r", "drainage_frac_r", "TSMD_final_r", "PAW_shallow_r",
    "fallow_deep_r", "swc_init_r", "PAWC_frac_r", "swc_final_r", "deep_drain_r"
  )
  SW[vars] <- 0 # initialize

  if (crop_type != "winter") {
    SW$winter_crop_r <- 0
  } else {
    SW$winter_crop_r <- ifelse(julianday < sow_start | julianday > harvest, 0, 1)
  }
  if (crop_type != "summer") {
    SW$summer_crop_r <- 0
  } else {
    # SW$summer_crop_r <- ifelse(julianday < sow_start, 0, 1)
    SW$summer_crop_r <- ifelse(julianday < sow_start & julianday + 500 > harvest, 0, 1)
  }
  SW$crop_r <- 1 # if perennial
  if (crop_type == "winter") {
    SW$crop_r <- SW$winter_crop_r
  }
  if (crop_type == "summer") {
    SW$crop_r <- SW$summer_crop_r
  }
  SW$crop_r <- with(SW, ifelse(year == 1 & julianday < 8 & veg_type == "annual", 0, crop_r))
  SW$crop_r <- with(SW, ifelse(year == 1 & julianday == 1, 0, crop_r))

  SW[c("winter_crop_r", "summer_crop_r", "crop_r")] <- sapply(SW[c("winter_crop_r", "summer_crop_r", "crop_r")], as.integer)

  # categories based on previous 7 day's rainfall
  # used to dictate what fraction of pan_evap is lost during fallow
  # 1: >25 mm, 2: >10 mm and <=25 mm, 3: > 1 mm and <=10 mm, 4: <=1 mm
  SW$rain_7days <- sum2(rainfall, 7)
  SW$fallow_cat_r <- as.integer(with(SW, get_fallow_cat(x = rain_7days)))
  SW$rain_7days <- NULL # remove

  # vectorize
  winter_crop_r <- SW$winter_crop_r
  summer_crop_r <- SW$summer_crop_r
  crop_r <- SW$crop_r # 0 or 1
  irrevent_r <- SW$irrevent_r # no irrigation considered yet
  irrigation_r <- SW$irrigation_r # no irrigation considered yet
  TSMD_init_r <- SW$TSMD_init_r
  TSMD_mod_r <- SW$TSMD_mod_r
  fallow_cat_r <- SW$fallow_cat_r
  ET_actual_r <- SW$ET_actual_r
  net_water_r <- SW$net_water_r
  infiltration_r <- SW$infiltration_r
  TSMD_acc_r <- SW$TSMD_acc_r
  excess_water_r <- SW$excess_water_r
  drainage_r <- SW$drainage_r
  runoff_r <- SW$runoff_r
  drainage_frac_r <- SW$drainage_frac_r
  TSMD_final_r <- SW$TSMD_final_r
  PAW_shallow_r <- SW$PAW_shallow_r
  fallow_deep_r <- SW$fallow_deep_r
  swc_init_r <- SW$swc_init_r
  PAWC_frac_r <- SW$PAWC_frac_r
  swc_final_r <- SW$swc_final_r
  deep_drain_r <- SW$deep_drain_r

  # initial default values for day 1:
  ET_init <- 2
  net_water_init <- rainfall[1] + irrigation_r[1] - ET_init
  infiltration_init <- pmin(net_water_init, infil_max)
  swc_init_start <- pmax(0, pmin(PAWC, infiltration_init + swc_init))

  for (i in 1:nrow(SW)) { # have to run in sequence

    # initial top soil moisture deficit (TSMD) for the day (mm)
    TSMD_init_r[i] <- ifelse(i == 1, TSMD_init, TSMD_final_r[i - 1])

    # a scaled value (0-1) with 0 equivalent to maximum TSMD and 1 equivalent to TSMD = 0
    # used to run nitrification and denitrification equations, modify ET, and calculate decomposition
    TSMD_mod_r[i] <- 1 - TSMD_init_r[i] / TSMD_max

    # adds infiltration to previous day's final soil water content, limits to PAWC (mm)
    swc_init_r[i] <- ifelse(i == 1, swc_init_start, swc_final_r[i - 1])

    # a fraction (0-1) of SWC/PAWC - used to modify ET when a crop is present
    PAWC_frac_r[i] <- swc_init_r[i] / PAWC

    # calculates evapotranspiration (mm) from pan evaporation
    ET_actual_r[i] <- ifelse(i == 1, ET_init, Eo_frac * panevap[i] * pmin(TSMD_mod_r[i], PAWC_frac_r[i]))

    # rainfall minus actual ET (mm)
    net_water_r[i] <- rainfall[i] + irrigation_r[i] - ET_actual_r[i]

    # infilteration rate
    infiltration_r[i] <- pmin(net_water_r[i], infil_max)

    # interim accumulated TSMD
    # can go negative, in which case there is excess water displayed in the next column
    # that can go to drainage or runoff if the top 30 cm is full
    TSMD_acc_r[i] <- pmin(TSMD_init_r[i] - infiltration_r[i], TSMD_max)

    # water in excess of the amount that can exist in the top 30 cm as defined by TSMD_max
    excess_water_r[i] <- pmax(0, -TSMD_acc_r[i])

    # limits excess water that can enter the deep soil to the maximum amount that can be held by the top 30 cm (pot_drainage)
    drainage_r[i] <- pmin(excess_water_r[i], pot_drain)

    # runoff (mm)
    runoff_r[i] <- pmax(net_water_r[i] - infiltration_r[i], 0) + pmax(excess_water_r[i] - pot_drain, 0)

    # fraction of existing water in the top 30 cm (equivalent to maximum TSMD when drainage occurs) that drains
    # used for leaching nitrate
    drainage_frac_r[i] <- pmin(pmax(drainage_r[i] / TSMD_max, 0), 1)

    # final TSMD, ranging between 0 and maximum top soil moisture deficit
    TSMD_final_r[i] <- pmax(0, pmin(TSMD_max, TSMD_acc_r[i]))

    # amount (mm) of water existing in the top 30 cm
    PAW_shallow_r[i] <- TSMD_max - TSMD_final_r[i]

    # during fallow, picks up any drainage that occurs, and evaporates it off if the top 30 cm of soil is dry
    # limits to potential fallow deepwater (calculated as PAWC - TSMD_max).
    if (i == 1) {
      fallow_deep_r[i] <- 0
    } else {
      if (crop_r[i] != 0) {
        temp <- 0
      } else {
        add_drain <- drainage_r[i] + fallow_deep_r[i - 1]
        if (TSMD_final_r[i] < TSMD_max) {
          temp <- add_drain
        } else if (add_drain > 0) {
          evap_bare <- get_bare_soil_evap(x = fallow_cat_r[i]) * panevap[i]
          temp <- add_drain - evap_bare
        } else {
          temp <- 0
        }
      }
      fallow_deep_r[i] <- pmin(pmax(temp, 0), pot_deep_PAW)
    }

    # final soil water content for day by subtracting ET from initial soil water content (infiltration + previous day's final SWC)
    # limits to 0 if soil dries out
    swc_final_r[i] <- pmax(0, pmin(PAWC, swc_init_r[i] + infiltration_r[i]))

    # drainage from the rooting zone (mm)
    deep_drain_r[i] <- ifelse(i == 1, 0, pmax(swc_init_r[i] + infiltration_r[i] - PAWC, 0))
  }

  # update
  SW$TSMD_init_r <- TSMD_init_r
  SW$TSMD_mod_r <- TSMD_mod_r
  SW$swc_init_r <- swc_init_r
  SW$PAWC_frac_r <- PAWC_frac_r
  SW$ET_actual_r <- ET_actual_r
  SW$net_water_r <- net_water_r
  SW$infiltration_r <- infiltration_r
  SW$TSMD_acc_r <- TSMD_acc_r
  SW$excess_water_r <- excess_water_r
  SW$drainage_r <- drainage_r
  SW$runoff_r <- runoff_r
  SW$drainage_frac_r <- drainage_frac_r
  SW$TSMD_final_r <- TSMD_final_r
  SW$PAW_shallow_r <- PAW_shallow_r
  SW$fallow_deep_r <- fallow_deep_r
  SW$swc_final_r <- swc_final_r
  SW$deep_drain_r <- deep_drain_r
  return(SW)
}

source("/media/DATADRIVE1/Project/muresk/model/d03_script/functions_rev_update.R")
source("/media/DATADRIVE1/Project/muresk/model/d03_script/calc_soil_water.R")

f <- list.files(path = "/media/DATADRIVE1/Project/muresk/soiltemp/Inputs/", pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

xy_new_cp <- read.csv("/media/DATADRIVE1/Project/muresk/croptype/geo_croptype.csv")
xy_new_cp$PAWC <- with(xy_new_cp, as.integer(AWC / 100 * 1000))

num_sites <- nrow(xy_new_cp)

for (i in 1:num_sites) {
  # Print progress (optional, but helpful for long loops)
  cat("Processing site:", i, "of", num_sites, "\n")

  ################################################################################
  # 1. Get crop parameters for the current row

  # Using [i,] to get data from the current row in the loop
  crop <- get_crop_params(xy_new_cp[i, c("cropname")])

  ################################################################################
  # 2. Read in weather data for the current row

  weatherID <- xy_new_cp[i, c("SILO")]
  weather_path <- paste0("/media/DATADRIVE1/Project/muresk/model/d02_data/SILO/", weatherID, ".csv")

  # It's good practice to check if the file exists before trying to read it
  if (!file.exists(weather_path)) {
    warning(paste("Weather file not found for site", i, ":", weather_path, "- Skipping this site."))
    next # Skips to the next iteration of the loop
  }

  weather <- read.csv(weather_path)

  ################################################################################
  # 3. Perform soil-water calculations and C inputs for the current row

  site <- xy_new_cp[i, c("BD", "clay", "PAWC")]
  weather$temp <- with(weather, (max + min) / 2)
  weather <- weather[, c(1:4, 9, 7, 8)]
  weather$date <- with(weather, as.Date(paste(year, month, day, sep = "-")))

  # SW <- calc_daily_sw(x = weather, crop = crop, site = site)
  SW_update <- calc_daily_sw_grok(x = weather, crop = crop, site = site)

  # SW_update$soilmoisture <- with(SW_update,swc_final_r/300)

  SW_update$soilmoisture <- round(with(SW_update, swc_final_r / 300), 2)

  SW_update <- SW_update[, c(1, 2, 3, 4, 8, 32)]

  siteID <- xy_new_cp[i, c("site")]

  write.csv(SW_update, paste0("/media/DATADRIVE1/Project/muresk/driving/SoilMoisture_New/", siteID, ".csv"))
}

################################################################################

f.st <- list.files(path = "/media/DATADRIVE1/Project/muresk/driving/SoilTemp_Outputs/", pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
f.sm <- list.files(path = "/media/DATADRIVE1/Project/muresk/driving/SoilMoisture_New/", pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
f.cp <- list.files(path = "/media/DATADRIVE1/Project/muresk/model/d04_output/cinputs_update/", pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)

cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)

library(doParallel)

foreach::foreach(i = seq(f.st), .combine = "cbind", .packages = "raster") %dopar% {
  st <- read.csv(f.st[i])
  sm <- read.csv(f.sm[i])
  cp <- readRDS(f.cp[i])

  f_base <- basename(f.st[i])
  temp_outs <- cbind(st[, c(2, 3, 4, 10)], sm[, c(7)], cp[, c(8)])
  temp_outs <- temp_outs[, c(1, 2, 3, 6, 5, 4)]
  names(temp_outs) <- c("year", "month", "day", "Cinputs", "soilmoisture", "soiltemp_030")
  temp_outs[, c(4)] <- with(temp_outs, sprintf("%.2f", Cinputs * 0.42 / 10))
  temp_outs$soilmoisture <- with(temp_outs, sprintf("%.2f", soilmoisture))
  temp_outs$soiltemp_030 <- with(temp_outs, sprintf("%.2f", soiltemp_030))
  colnames(temp_outs) <- c("Year", "Month", "Day", "forc_npp", "forc_sw", "forc_st")
  temp_outs <- temp_outs[367:9132, ]
  write.csv(temp_outs, paste0("/media/DATADRIVE1/Project/muresk/driving/forcing_inputs_new/", f_base))
}

stopCluster(cl)
