# explainable soil water balance model
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

# annual crop
calc_annuals_new <- function(crop, by) { # for annuals ("annual_calc")

    an_TE <- crop$TE # transpiration efficiency of system (9)
    transpiration_frac_s <- crop$transpiration_frac_s # (10)
    transpiration_frac_d <- crop$transpiration_frac_d # (11)

    # sow_rain <- 0                       # rain to sow (mm/7 days)
    # sow_PAWC_frac <- crop$sow_PAWC_frac # fraction of plant available water capacity to sow (13)

    # for annuals - growth function
    an_max_DM <- crop$max_DM # maximum permissible total cumulative shoot DM in any year (kg/ha) (17)
    an_sow_day <- crop$sow # (30)
    an_harvest <- crop$harvest.1 # (31)
    an_growth_a <- crop$grwth_a # (33)
    an_growth_b <- crop$grwth_b # (34)

    an_max_days <- crop$harvest.1 - crop$sow
    if (an_max_days < 0) { # for a crop harvested in the second year of sowing
        an_max_days <- (365 - crop$sow) + crop$harvest.1
    }

    shoot_root_ratio <- ifelse(!is.null(crop$shoot_root_ratio_sampled), crop$shoot_root_ratio_sampled, crop$shoot_root_ratio) # shoot-to-root ratio (35)
    min_shoot_DM <- ifelse(!is.null(crop$min_shoot_DM_sampled), crop$min_shoot_DM_sampled, crop$min_shoot_DM) # minimum permissible shoot DM (kg/ha) (36)

    grazing_frac <- ifelse(!is.null(crop$grazing_frac_sampled), crop$grazing_frac_sampled, crop$grazing_frac) # fraction of todays growth grazed (37)
    if (is.na(grazing_frac)) {
        grazing_frac <- 0 # assume no grazing if missing
    }

    # here, we name shed as root exudation
    grazed_root_shed <- ifelse(!is.null(crop$grazed_root_shed_sampled), crop$grazed_root_shed_sampled, crop$grazed_root_shed) # fractional root response to shoot grazing (38)

    if (grazed_root_shed == 0) {
        # root exudation is 0.5, as fraction of DM allocated to roots lost
        # as root exudation is 0.5, ranging between 0.4 and 0.7
        grazed_root_shed <- 0.5 # assume no grazing if missing
        # if the live root biomass turnover rate is 0.65 per year, the daily rate is 0.0018
        # grazed_root_shed <- 0.0032 # assume no grazing if missing
    }

    offtake_as_dung <- crop$offtake_as_dung # dung returns (39)

    # non-legume %N function
    # an_nonleg_pcnt_N_a	<- 3.35
    # an_nonleg_pcnt_N_b	<- 0.65
    # an_nonleg_pcnt_N_c	<- 0.13
    # an_nonleg_pcnt_N_frac <- 0.5
    # an_nonleg_harvest_N_pcnt <- 1
    # an_nonleg_residue_N_pcnt <- 0.3
    # an_nonleg_max_grain_N_pcnt	<- 3
    # an_nonleg_root_N_frac	<- 0.4
    # an_nonleg_root_harvest_N_pcnt	<- 0.5
    # an_nonleg_shoot_C_frac	<- 0.42
    # an_nonleg_root_C_frac	<- 0.42
    # an_nonleg_HI	<- 0.30

    # by site weather and soil-water conditions
    AN <- by[c("date", "julianday", "month", "year", "crop_r", "ET_actual_r", "fallow_deep_r")]

    vars <- c(
        "an_season_day", "an_WLDM_increment", "an_sum_DM", "max_DM", "daily_DM",
        "daily_growth", "root_growth", "standing_shoot_DM_start", "living_root_DM_start",
        "grazing_offtake", "roots_shed", "shoots_post_grazing", "roots_post_grazing",
        "dung_returns"
    )
    AN[vars] <- 0 # initialize

    # set scalars
    input_multiplier <- 1 # input scenarios?
    fertility_scalar <- 1 # site or management specific scalar to increase or decrease potential pasture or crop growth

    # the potential water limited crop DM is calculated from available water and ET
    AN$an_season_day <- with(AN, stats::ave(AN$crop_r, cumsum(AN$crop_r == 0), FUN = cumsum))

    AN$an_WLDM_increment <- with(AN, ifelse(crop_r == 0, 0,
        ((ET_actual_r * transpiration_frac_s) + (fallow_deep_r * transpiration_frac_d)) * an_TE * fertility_scalar
    ))

    # cumulative sum of crop DM for any season type (kg DM/ha)
    AN$an_sum_DM <- stats::ave(AN$an_WLDM_increment, cumsum(AN$an_WLDM_increment == 0), FUN = cumsum)

    AN$max_DM <- with(AN, ifelse(crop_r == 0, 0, stats::ave(an_sum_DM, cumsum(an_sum_DM == 0), FUN = max))) # total_DM

    # check if maximum cumulative DM is within the maximum permissible total cumulative DM
    AN$max_DM <- with(AN, ifelse(max_DM > an_max_DM, an_max_DM, max_DM))

    AN$daily_DM <- with(AN, ifelse(crop_r == 0, 0,
        max_DM / (1 + exp(-(an_season_day - an_growth_a * an_max_days) / ((an_growth_b * an_sow_day) * an_max_days)))
    ))

    AN$daily_growth[1] <- 0 # assume zero initial growth

    for (i in 2:nrow(AN)) {
        AN$daily_growth[i] <- ifelse(AN$crop_r[i] == 0, 0, AN$daily_DM[i] - AN$daily_DM[i - 1])
    }

    AN$root_growth <- with(AN, daily_growth * shoot_root_ratio)

    AN$grazing_offtake <- with(AN, ifelse(daily_DM < min_shoot_DM, 0, daily_growth * grazing_frac))

    # abover-ground residues means what reminder is left on the soil surface
    # we treat the shoot and root decay rate during growth stage as negligible

    AN$roots_shed <- with(AN, ifelse(daily_DM < min_shoot_DM, 0, root_growth * grazing_frac * grazed_root_shed))

    AN$dung_returns <- with(AN, grazing_offtake * offtake_as_dung)

    # also assume no initial shoot and root
    AN$standing_shoot_DM_start[1] <- 0
    AN$living_root_DM_start[1] <- 0
    AN$shoots_post_grazing[1] <- 0
    AN$roots_post_grazing[1] <- 0

    for (i in 2:nrow(AN)) {
        AN$standing_shoot_DM_start[i] <- ifelse(AN$crop_r[i] == 0, 0, # harvested
            AN$daily_growth[i] + AN$shoots_post_grazing[i - 1]
        )

        AN$living_root_DM_start[i] <- ifelse(AN$crop_r[i] == 0, 0,
            AN$root_growth[i] + AN$roots_post_grazing[i - 1]
        )

        AN$shoots_post_grazing[i] <- AN$standing_shoot_DM_start[i] - AN$grazing_offtake[i]

        AN$roots_post_grazing[i] <- AN$living_root_DM_start[i] - AN$roots_shed[i]
    }

    AN
}

# perennial grass
calc_perennials_v3 <- function(crop, by) {
    # for perennials ("perennial_calc")

    veg_type <- crop$lifeform # vegetation type - annual or perennial (3)
    T_limit <- ifelse(veg_type == "perennial", "yes", "no") # temperature limitation to plant production

    pe_nonleg_DM_init <- 50 # initial perennial non-legume dry matter (kg DM/ha)
    pe_leg_DM_init <- 50 # initail perennial legume dry matter (kg DM/ha)
    # pe_leg_DM_frac <- 0.5 # fraction of the system that is legume
    pe_TE <- crop$TE # transpiration efficiency of the grass species (kg DM/mm) (9)
    pe_TEc <- 1

    critical_SWI <- crop$critical_SWI # fraction of PAWC at which sensecence occurs (18)
    tree_cover <- crop$tree_cover # tree cover fraction (19)
    pe_max_growth <- crop$max_growth # kg (20)
    die_off <- crop$die_off # the amount of grass biomass remaining after die off due to senescence
    # that occurs when critical soil water or transpiration is reached (21)
    pe_critical_transp <- crop$critical_T # the transpiration rate at which senescence occurs (mm) (22)
    # the proportion of the land surface covered by vegetation (grass and trees)
    # is used to calculate grass transpiration after taking tree cover into account
    veg_cover <- crop$veg_cover # fraction (23)

    # temperature index constants (C3 crop species)
    pe_TI_opt <- crop$pe_TI_opt # optimum temperature for growth (24)
    pe_TI_min <- crop$pe_TI_min # minimum temperature for growth (25)
    pe_TI_max <- crop$pe_TI_max # maximum temperature for growth (26)
    pe_TI_b_sub <- crop$pe_TI_b_sub # changes inflexion point of TI curve (27)
    pe_TI_b_supra <- crop$pe_TI_b_supra # (28)

    # fixed parameter value
    shoot_root_ratio <- ifelse(!is.null(crop$shoot_root_ratio_sampled), crop$shoot_root_ratio_sampled, crop$shoot_root_ratio) # shoot-to-root ratio (35)
    min_shoot_DM <- ifelse(!is.null(crop$min_shoot_DM_sampled), crop$min_shoot_DM_sampled, crop$min_shoot_DM) # mimum permissible shoot DM (kg/ha) (36)
    grazing_frac <- ifelse(!is.null(crop$grazing_frac_sampled), crop$grazing_frac_sampled, crop$grazing_frac) # fraction of today's growth grazed (37)
    grazed_root_shed <- ifelse(!is.null(crop$grazed_root_shed_sampled), crop$grazed_root_shed_sampled, crop$grazed_root_shed) # fractional root response to shoot grazing (38)
    offtake_as_dung <- crop$offtake_as_dung # dung returns (39), can be setted as one-third, see "Carbon sequestration under subtropical perennial pastures II: Carbon dynamics"
    root_die_off <- crop$root_die_off # (40)


    # fire
    pe_critical_fireload <- 20 # kg DM
    pe_fire_residual <- 75 # kg/ha
    pe_burn <- "no" # yes or no
    pe_burn_DOY <- 200 # day

    # by site weather and soil-water variables
    PE <- by[c("date", "julianday", "month", "year", "temp", "PAWC_frac_r", "ET_actual_r")]

    vars <- c(
        "pe_WL_transp", "pe_TI_x", "pe_TI", "pe_WL_growth", "root_growth",
        "standing_shoot_DM_start", "living_root_DM_start", "grazing_offtake",
        "roots_shed", "shoots_post_grazing", "roots_post_grazing", "pe_dieoff",
        "root_die_off", "shoot_DM_end", "root_DM_end", "shoot_and_root_residues",
        "dung_returns"
    )

    PE[vars] <- 0 # initialize

    # set scalars
    input_multiplier <- 1 # input scenarios?
    fertility_scalar <- 1 # site or management specific scalar to increase or decreas potential pasture or crop growth

    temp <- PE$temp
    PAWC_frac_r <- PE$PAWC_frac_r
    ET_actual_r <- PE$ET_actual_r

    for (i in 1:nrow(PE)) { # first calculate temparature and water limited growth

        # water limited transpriation (mm)
        PE$pe_WL_transp[i] <- ET_actual_r[i] * veg_cover * (1 - tree_cover)

        # Nix (1981) temperature index (TI) function
        PE$pe_TI_x[i] <- abs(ifelse(temp[i] < pe_TI_min, (pe_TI_opt - temp[i]) / (pe_TI_opt - pe_TI_min),
            (temp[i] - pe_TI_opt) / (pe_TI_max - pe_TI_opt)
        ))

        PE$pe_TI[i] <- ifelse(temp[i] < pe_TI_min | temp[i] > pe_TI_max, 0,
            ifelse(T_limit == "no", 1, get_pe_TI(x = PE$pe_TI_x[i], crop = crop))
        )

        # water limited (WL) growth
        PE$pe_WL_growth[i] <- ifelse(pe_max_growth < PE$pe_WL_transp[i] * pe_TE * pe_TEc * PE$pe_TI[i] * fertility_scalar,
            pe_max_growth, PE$pe_WL_transp[i] * pe_TE * pe_TEc * PE$pe_TI[i] * fertility_scalar
        )

        PE$root_growth[i] <- PE$pe_WL_growth[i] * shoot_root_ratio
    }

    # initialise (for i = 1)
    # in our cases, the shoot DM is less than 1.5 t/ha
    # min_shoot_DM <- 1500
    PE$standing_shoot_DM_start[1] <- min_shoot_DM
    PE$living_root_DM_start[1] <- min_shoot_DM * shoot_root_ratio
    PE$grazing_offtake[1] <- 0
    PE$roots_shed[1] <- 0
    PE$shoots_post_grazing[1] <- min_shoot_DM
    PE$roots_post_grazing[1] <- min_shoot_DM * shoot_root_ratio
    PE$pe_dieoff[1] <- 0
    PE$root_die_off[1] <- 0
    PE$shoot_DM_end[1] <- min_shoot_DM
    PE$root_DM_end[1] <- min_shoot_DM * shoot_root_ratio
    PE$shoot_and_root_residues[1] <- 0
    PE$dung_returns[1] <- 0
    PE$shoot_turnover[1] <- 0
    PE$shoot_stress_loss[1] <- 0
    PE$root_stress_loss[1] <- 0
    PE$root_dieoff[1] <- 0

    # Use sampled shoot_turnover_rate if available, otherwise use default
    shoot_turnover_rate <- ifelse(!is.null(crop$shoot_turnover_rate_sampled),
        crop$shoot_turnover_rate_sampled, 0.001
    )

    for (i in 2:nrow(PE)) {
        PE$standing_shoot_DM_start[i] <- PE$pe_WL_growth[i] + PE$shoot_DM_end[i - 1]

        PE$living_root_DM_start[i] <- PE$root_growth[i] + PE$root_DM_end[i - 1]

        PE$grazing_offtake[i] <- ifelse(PE$shoot_DM_end[i - 1] < min_shoot_DM, 0,
            PE$pe_WL_growth[i] * grazing_frac
        )

        # post-grazing biomass
        PE$shoots_post_grazing[i] <- PE$standing_shoot_DM_start[i] - PE$grazing_offtake[i]

        is_growing <- PE$shoots_post_grazing[i] > 0

        # root exudates only occur during active growth
        PE$roots_shed[i] <- ifelse(is_growing || PE$shoot_DM_end[i - 1] >= min_shoot_DM,
            PE$root_growth[i] * grazed_root_shed, 0
        )

        PE$roots_post_grazing[i] <- PE$living_root_DM_start[i] - PE$roots_shed[i]

        # daily shoot turnover
        PE$shoot_turnover[i] <- ifelse(is_growing, PE$shoots_post_grazing[i] * shoot_turnover_rate, 0)

        # senescence condition
        stress_condition <- PE$pe_WL_transp[i] < pe_critical_transp | PE$PAWC_frac_r[i] < critical_SWI

        # stress related shoot and root losses
        PE$shoot_stress_loss[i] <- ifelse(stress_condition,
            pmax(PE$shoots_post_grazing[i] - PE$shoot_turnover[i], 0) * (1 - die_off), 0
        )

        PE$root_stress_loss[i] <- ifelse(stress_condition,
            PE$roots_post_grazing[i] * (1 - die_off) * root_die_off, 0
        )

        PE$pe_dieoff[i] <- PE$shoot_turnover[i] + PE$shoot_stress_loss[i]

        PE$root_dieoff[i] <- PE$roots_shed[i] + PE$root_stress_loss[i]

        # check fire
        PE$shoot_DM_end[i] <- ifelse(pe_burn == "yes" & PE$julianday[i] == pe_burn_DOY & PE$standing_shoot_DM_start[i] > pe_critical_fireload,
            pe_fire_residual,
            PE$shoots_post_grazing[i] - PE$pe_dieoff[i]
        )

        PE$root_DM_end[i] <- PE$roots_post_grazing[i] - PE$root_dieoff[i]

        PE$shoot_and_root_residues[i] <- (PE$root_dieoff[i] + PE$pe_dieoff[i]) * input_multiplier

        PE$dung_returns[i] <- (offtake_as_dung * PE$grazing_offtake[i]) * input_multiplier
    }

    PE
}

# daily crop residues
get_AN_daily_v3 <- function(x, crop, shoot_pool_init = 0, root_pool_init = 0,
                            HI = NULL, shoot_turnover_rate = NULL,
                            shoot_decay_rate = NULL, root_decay_rate = NULL) {
    AN_d <- x
    AN_d$s <- 0
    AN_d$plant_residue_returns <- 0
    AN_d$shoot_residue_pool <- 0
    AN_d$root_residue_pool <- 0
    AN_d$root_residue_pool_old <- 0
    AN_d$root_residue_pool_new <- 0

    # Use provided HI or get from crop parameters
    if (is.null(HI)) {
        HI <- crop$harvest_index
        if (is.na(HI) || length(HI) == 0) {
            HI <- 0.1
        }
    }

    # Use provided rates or defaults
    if (is.null(shoot_turnover_rate)) shoot_turnover_rate <- 0.001
    if (is.null(shoot_decay_rate)) shoot_decay_rate <- 0.0043
    if (is.null(root_decay_rate)) root_decay_rate <- 0.005

    n_days <- nrow(AN_d)
    if (n_days == 0) {
        return(AN_d)
    }

    AN_d$shoot_residue_pool[1] <- shoot_pool_init
    AN_d$root_residue_pool_old[1] <- root_pool_init

    for (i in seq_len(n_days)) {
        if (i > 1) {
            AN_d$shoot_residue_pool[i] <- AN_d$shoot_residue_pool[i - 1]
            AN_d$root_residue_pool_old[i] <- AN_d$root_residue_pool_old[i - 1]
            AN_d$root_residue_pool_new[i] <- AN_d$root_residue_pool_new[i - 1]
        }

        is_growing <- AN_d$shoots_post_grazing[i] > 0
        will_be_growing <- i < n_days && AN_d$shoots_post_grazing[i + 1] > 0
        is_harvest_day <- is_growing && !will_be_growing

        shoot_turnover <- if (is_growing) AN_d$shoots_post_grazing[i] * shoot_turnover_rate else 0
        root_shed <- if (is_growing) AN_d$roots_shed[i] else 0

        # 1. decay existing pools before adding today's inputs
        shoot_decomp <- AN_d$shoot_residue_pool[i] * shoot_decay_rate
        root_decomp_old <- AN_d$root_residue_pool_old[i] * root_decay_rate
        root_decomp_new <- if (!is_growing) AN_d$root_residue_pool_new[i] * root_decay_rate else 0

        AN_d$shoot_residue_pool[i] <- AN_d$shoot_residue_pool[i] - shoot_decomp
        AN_d$root_residue_pool_old[i] <- AN_d$root_residue_pool_old[i] - root_decomp_old
        AN_d$root_residue_pool_new[i] <- AN_d$root_residue_pool_new[i] - root_decomp_new

        # 2. add today's fresh inputs after decay
        AN_d$s[i] <- shoot_turnover
        AN_d$shoot_residue_pool[i] <- AN_d$shoot_residue_pool[i] + shoot_turnover
        AN_d$root_residue_pool_new[i] <- AN_d$root_residue_pool_new[i] + root_shed

        unharvested_shoots <- 0
        harvest_root_add <- 0
        if (is_harvest_day) {
            unharvested_shoots <- max(AN_d$shoots_post_grazing[i] - shoot_turnover, 0) * (1 - HI)
            harvest_root_add <- AN_d$roots_post_grazing[i]
            AN_d$shoot_residue_pool[i] <- AN_d$shoot_residue_pool[i] + unharvested_shoots
            AN_d$root_residue_pool_old[i] <- AN_d$root_residue_pool_old[i] + AN_d$root_residue_pool_new[i] + harvest_root_add
            AN_d$root_residue_pool_new[i] <- 0
        } else if (!is_growing) {
            AN_d$root_residue_pool_old[i] <- AN_d$root_residue_pool_old[i] + AN_d$root_residue_pool_new[i]
            AN_d$root_residue_pool_new[i] <- 0
        }

        AN_d$root_residue_pool[i] <- AN_d$root_residue_pool_old[i] + AN_d$root_residue_pool_new[i]

        decomp_flux <- shoot_decomp + root_decomp_old + root_decomp_new
        daily_inputs <- shoot_turnover + root_shed

        if (is_growing) {
            AN_d$plant_residue_returns[i] <- decomp_flux + daily_inputs
        } else {
            AN_d$plant_residue_returns[i] <- decomp_flux
        }
    }

    attr(AN_d, "shoot_residue_pool_final") <- AN_d$shoot_residue_pool[n_days]
    attr(AN_d, "root_residue_pool_final") <- AN_d$root_residue_pool[n_days]

    AN_d
}

################################################################################
# Function to validate root exudation constraint
# Checks if sum(roots_shed) / sum(daily_growth) is between 7% and 15%
validate_root_exudation_constraint <- function(result_data, min_ratio = 0.07, max_ratio = 0.15) {
    # Only check during growth periods (crop_r == 1 for annuals, active growth for perennials)
    if ("crop_r" %in% names(result_data)) {
        # Annual crops
        growth_data <- result_data[result_data$crop_r == 1, ]
    } else {
        # For perennials, check when shoots_post_grazing > 0
        growth_data <- result_data[result_data$shoots_post_grazing > 0 | result_data$pe_WL_growth > 0, ]
    }

    if (nrow(growth_data) == 0) {
        return(TRUE) # No growth data, skip constraint
    }

    # Calculate sums
    if ("daily_growth" %in% names(growth_data)) {
        # Annual crops
        total_shoot_growth <- sum(growth_data$daily_growth, na.rm = TRUE)
        total_root_exudation <- sum(growth_data$roots_shed, na.rm = TRUE)
    } else if ("pe_WL_growth" %in% names(growth_data)) {
        # Perennial crops
        total_shoot_growth <- sum(growth_data$pe_WL_growth, na.rm = TRUE)
        total_root_exudation <- sum(growth_data$roots_shed, na.rm = TRUE)
    } else {
        return(TRUE) # Cannot validate, skip
    }

    if (total_shoot_growth <= 0) {
        return(TRUE) # Avoid division by zero
    }

    # Calculate actual ratio: root exudation as fraction of shoot growth
    actual_ratio <- total_root_exudation / total_shoot_growth

    # Check if within bounds
    is_valid <- actual_ratio >= min_ratio && actual_ratio <= max_ratio

    return(is_valid)
}

################################################################################
# Function to sample uncertain parameters for Monte Carlo simulation
# with rejection sampling to enforce root exudation constraint
sample_uncertain_params_constrained <- function(crop_base, SW_update,
                                                min_ratio = 0.07, max_ratio = 0.15,
                                                max_attempts = 50) {
    attempt <- 0

    while (attempt < max_attempts) {
        attempt <- attempt + 1

        # Sample parameters
        mc_params <- list(
            HI = runif(1, min = 0.2, max = 0.6), # Harvest index
            shoot_turnover_rate = runif(1, min = 0.001, max = 0.005), # Shoot turnover rate
            shoot_decay_rate = runif(1, min = 0.001, max = 0.005), # Shoot decay rate
            root_decay_rate = runif(1, min = 0.001, max = 0.005), # Root decay rate
            grazed_root_shed = runif(1, min = 0.4, max = 0.7), # Root exudation (grazed_root_shed)
            grazing_frac = runif(1, min = 0.0, max = 0.5), # Grazing fraction
            min_shoot_DM = runif(1, min = 1200, max = 2000), # Minimum shoot DM
            shoot_root_ratio = runif(1, min = 0.1, max = 0.5) # Shoot-to-root ratio
        )

        # Create test crop with sampled parameters
        crop_test <- crop_base
        crop_test$grazed_root_shed_sampled <- mc_params$grazed_root_shed
        crop_test$grazing_frac_sampled <- mc_params$grazing_frac
        crop_test$min_shoot_DM_sampled <- mc_params$min_shoot_DM
        crop_test$shoot_turnover_rate_sampled <- mc_params$shoot_turnover_rate
        crop_test$shoot_root_ratio_sampled <- mc_params$shoot_root_ratio

        # Run a quick test simulation
        if (crop_test$lifeform == "annual") {
            test_result <- calc_annuals_new(crop = crop_test, by = SW_update)
        } else {
            test_result <- calc_perennials_v3(crop = crop_test, by = SW_update)
        }

        # Validate constraint
        if (validate_root_exudation_constraint(test_result, min_ratio, max_ratio)) {
            # Add constraint info to parameters
            mc_params$constraint_attempts <- attempt
            mc_params$constraint_satisfied <- TRUE
            return(mc_params)
        }
    }

    # If max attempts reached, return parameters anyway with flag
    # This prevents infinite loops but allows tracking of failures
    mc_params$constraint_attempts <- max_attempts
    mc_params$constraint_satisfied <- FALSE
    warning(paste("Could not satisfy root exudation constraint after", max_attempts, "attempts"))

    return(mc_params)
}

################################################################################
# Simple unconstrained sampling (for backward compatibility or when constraint disabled)
sample_uncertain_params <- function() {
    list(
        HI = runif(1, min = 0.2, max = 0.6), # Harvest index
        shoot_turnover_rate = runif(1, min = 0.001, max = 0.005), # Shoot turnover rate
        shoot_decay_rate = runif(1, min = 0.001, max = 0.005), # Shoot decay rate
        root_decay_rate = runif(1, min = 0.001, max = 0.005), # Root decay rate
        grazed_root_shed = runif(1, min = 0.4, max = 0.7), # Root exudation (grazed_root_shed)
        grazing_frac = runif(1, min = 0.0, max = 0.5), # Grazing fraction
        min_shoot_DM = runif(1, min = 1200, max = 2000), # Minimum shoot DM
        shoot_root_ratio = runif(1, min = 0.1, max = 0.5), # Shoot-to-root ratio
        constraint_satisfied = NA # Not checked
    )
}

################################################################################
# Output variables to retain when saving results
annual_output_vars <- c(
    "plant_residue_returns", "shoot_residue_pool", "root_residue_pool",
    "daily_DM", "daily_growth", "root_growth", "roots_shed"
)

perennial_output_vars <- c(
    "shoot_and_root_residues", "pe_WL_growth", "root_growth",
    "root_dieoff", "pe_dieoff", "shoot_DM_end", "root_DM_end"
)

output_variables_to_keep <- unique(c(annual_output_vars, perennial_output_vars))

source("/media/DATADRIVE1/Project/muresk/model/d03_script/functions_rev_update.R")
source("/media/DATADRIVE1/Project/muresk/model/d03_script/calc_soil_water.R")


xy_new_cp <- read.csv("/media/DATADRIVE1/Project/muresk/croptype/geo_croptype.csv")
xy_new_cp$PAWC <- with(xy_new_cp, as.integer(AWC / 100 * 1000))

# Get the total number of rows to loop through
num_sites <- nrow(xy_new_cp)

# Set number of Monte Carlo iterations
n_mc_iterations <- 100 # Adjust this number based on computational resources

################################################################################
# CONSTRAINT CONFIGURATION
################################################################################
# Enable/disable root exudation constraint
ENABLE_ROOT_EXUDATION_CONSTRAINT <- TRUE # Set to FALSE to disable constraint

# Root exudation ratio bounds (as fraction of shoot growth)
ROOT_EXUDATION_MIN_RATIO <- 0.07 # Minimum: root exudation = 7% of shoot growth
ROOT_EXUDATION_MAX_RATIO <- 0.15 # Maximum: root exudation = 15% of shoot growth

# Maximum attempts for rejection sampling
MAX_CONSTRAINT_ATTEMPTS <- 50 # Increase if many rejections occur

cat("\n=== Constraint Configuration ===\n")
cat("Root exudation constraint:", ifelse(ENABLE_ROOT_EXUDATION_CONSTRAINT, "ENABLED", "DISABLED"), "\n")
if (ENABLE_ROOT_EXUDATION_CONSTRAINT) {
    cat("  Acceptable range: ", ROOT_EXUDATION_MIN_RATIO * 100, "% to ", ROOT_EXUDATION_MAX_RATIO * 100, "%\n", sep = "")
    cat("  Max sampling attempts:", MAX_CONSTRAINT_ATTEMPTS, "\n")
}
cat("================================\n\n")

################################################################################
# Sequential loop to process each site
cat("Starting sequential Monte Carlo simulation across", num_sites, "sites\n")
cat("Each site will run", n_mc_iterations, "Monte Carlo iterations\n")
cat("Total simulations:", num_sites * n_mc_iterations, "\n\n")

all_results_full <- list()
output_dir <- "/media/DATADRIVE1/Project/muresk/model/d04_output"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

for (i in seq_len(num_sites)) {
    cat("Processing site", i, "of", num_sites, ":", xy_new_cp$site[i], "\n")

    site_results <- tryCatch(
        {
            ################################################################################
            # 1. Get crop parameters for the current row

            crop_base <- get_crop_params(xy_new_cp[i, c("cropname")])

            ################################################################################
            # 2. Read in weather data for the current row

            weatherID <- xy_new_cp[i, c("SILO")]
            weather_path <- paste0("/media/DATADRIVE1/Project/muresk/model/d02_data/SILO/", weatherID, ".csv")

            # Check if the file exists before trying to read it
            if (!file.exists(weather_path)) {
                warning(paste("Weather file not found for site", i, ":", weather_path, "- Skipping this site."))
                return(list())
            }

            weather <- read.csv(weather_path)

            ################################################################################
            # 3. Perform soil-water calculations (only once per site)

            site <- xy_new_cp[i, c("BD", "clay", "PAWC")]
            weather$temp <- with(weather, (max + min) / 2)
            weather <- weather[, c(1:4, 9, 7, 8)]
            weather$date <- with(weather, as.Date(paste(year, month, day, sep = "-")))

            SW_update <- calc_daily_sw_update(x = weather, crop = crop_base, site = site)

            ################################################################################
            # 4. Monte Carlo iterations with parameter uncertainty

            site_results_inner <- vector("list", n_mc_iterations)

            for (mc in seq_len(n_mc_iterations)) {
                # Sample uncertain parameters (with or without constraint)
                if (ENABLE_ROOT_EXUDATION_CONSTRAINT) {
                    mc_params <- sample_uncertain_params_constrained(
                        crop_base = crop_base,
                        SW_update = SW_update,
                        min_ratio = ROOT_EXUDATION_MIN_RATIO,
                        max_ratio = ROOT_EXUDATION_MAX_RATIO,
                        max_attempts = MAX_CONSTRAINT_ATTEMPTS
                    )
                } else {
                    mc_params <- sample_uncertain_params()
                }

                # Create crop parameters with sampled values
                crop <- crop_base
                crop$grazed_root_shed_sampled <- mc_params$grazed_root_shed
                crop$grazing_frac_sampled <- mc_params$grazing_frac
                crop$min_shoot_DM_sampled <- mc_params$min_shoot_DM
                crop$shoot_turnover_rate_sampled <- mc_params$shoot_turnover_rate
                crop$shoot_root_ratio_sampled <- mc_params$shoot_root_ratio

                # Run model with sampled parameters
                if (crop$lifeform == "annual") {
                    # Calculate daily growth
                    AN_update <- calc_annuals_new(crop = crop, by = SW_update)

                    # Change the monthly and annual function to daily function with MC parameters
                    AN_con <- get_AN_daily_v3(
                        x = AN_update,
                        crop = crop,
                        shoot_pool_init = 1500,
                        root_pool_init = 900,
                        HI = mc_params$HI,
                        shoot_turnover_rate = mc_params$shoot_turnover_rate,
                        shoot_decay_rate = mc_params$shoot_decay_rate,
                        root_decay_rate = mc_params$root_decay_rate
                    )

                    LN_d <- AN_con
                } else {
                    # For perennials
                    PE <- calc_perennials_v3(crop = crop, by = SW_update)

                    LN_d <- PE
                }

                # Add site and Monte Carlo iteration information
                LN_d$lifeform <- crop_base$lifeform
                LN_d$site <- xy_new_cp[i, c("site")]
                LN_d$landuse <- xy_new_cp[i, c("land_use")]
                LN_d$cropname <- xy_new_cp[i, c("cropname")]
                LN_d$mc_iteration <- mc
                LN_d$HI_sampled <- mc_params$HI
                LN_d$shoot_turnover_rate_sampled <- mc_params$shoot_turnover_rate
                LN_d$shoot_decay_rate_sampled <- mc_params$shoot_decay_rate
                LN_d$root_decay_rate_sampled <- mc_params$root_decay_rate
                LN_d$grazed_root_shed_sampled <- mc_params$grazed_root_shed
                LN_d$grazing_frac_sampled <- mc_params$grazing_frac
                LN_d$min_shoot_DM_sampled <- mc_params$min_shoot_DM
                LN_d$shoot_root_ratio_sampled <- mc_params$shoot_root_ratio

                # Add constraint validation info
                if (!is.null(mc_params$constraint_satisfied)) {
                    LN_d$constraint_satisfied <- mc_params$constraint_satisfied
                    LN_d$constraint_attempts <- mc_params$constraint_attempts
                } else {
                    LN_d$constraint_satisfied <- NA
                    LN_d$constraint_attempts <- NA
                }

                ################################################################################
                # 5. Store the result in the site results list

                site_results_inner[[mc]] <- LN_d
            } # End of Monte Carlo loop

            site_results_inner
        },
        error = function(e) {
            warning(paste("Error processing site", i, ":", e$message))
            list()
        }
    )

    if (length(site_results) == 0) {
        next
    }

    site_results <- Filter(function(x) !is.null(x) && nrow(x) > 0, site_results)
    if (length(site_results) == 0) {
        next
    }

    site_results_table_full <- do.call(rbind, site_results)

    all_results_full[[length(all_results_full) + 1]] <- site_results_table_full

    available_vars <- intersect(names(site_results_table_full), output_variables_to_keep)
    results_list_table <- site_results_table_full[, available_vars, drop = FALSE]

    site_name_raw <- xy_new_cp$site[i]
    if (is.null(site_name_raw) || is.na(site_name_raw) || identical(site_name_raw, "")) {
        site_name_raw <- paste0("site_", i)
    }
    site_file_name <- paste0(gsub("[^A-Za-z0-9_]", "_", site_name_raw), ".rds")

    saveRDS(results_list_table, file.path(output_dir, site_file_name))
    cat("  Saved site results to", site_file_name, "\n")
}

cat("\nSequential processing complete!\n")

################################################################################
# Calculate summary statistics for each site
cat("\nCalculating summary statistics...\n")

library(dplyr)

results_list_table <- readRDS("/media/DATADRIVE1/Project/muresk/model/d04_output/cinputs/5542.rds")
# Process annual crop summary if exists

annual_results <- results_list_table

annual_vars <- annual_output_vars[annual_output_vars %in% names(annual_results)]

if (length(annual_vars) > 0) {
        summary_annual <- annual_results %>%
            group_by(site, landuse, cropname, date, julianday, month, year) %>%
            summarise(
                across(all_of(annual_vars),
                    list(
                        mean = ~ mean(., na.rm = TRUE),
                        sd = ~ sd(., na.rm = TRUE),
                        q025 = ~ quantile(., 0.025, na.rm = TRUE),
                        q25 = ~ quantile(., 0.25, na.rm = TRUE),
                        median = ~ median(., na.rm = TRUE),
                        q75 = ~ quantile(., 0.75, na.rm = TRUE),
                        q975 = ~ quantile(., 0.975, na.rm = TRUE)
                    ),
                    .names = "{.col}_{.fn}"
                ),
                .groups = "drop"
            )

    saveRDS(summary_annual, "/media/DATADRIVE1/Project/muresk/model/d04_output/muresk_crop_Cinputs_MC_summary_annual.rds")
    cat("Annual crop summary statistics saved.\n")
}



# Process perennial crop summary if exists

# Filter to existing columns
perennial_vars <- perennial_output_vars[perennial_output_vars %in% names(perennial_results)]

if (length(perennial_vars) > 0) {
        summary_perennial <- perennial_results %>%
            group_by(site, landuse, cropname, date, julianday, month, year) %>%
            summarise(
                across(all_of(perennial_vars),
                    list(
                        mean = ~ mean(., na.rm = TRUE),
                        sd = ~ sd(., na.rm = TRUE),
                        q025 = ~ quantile(., 0.025, na.rm = TRUE),
                        q25 = ~ quantile(., 0.25, na.rm = TRUE),
                        median = ~ median(., na.rm = TRUE),
                        q75 = ~ quantile(., 0.75, na.rm = TRUE),
                        q975 = ~ quantile(., 0.975, na.rm = TRUE)
                    ),
                    .names = "{.col}_{.fn}"
                ),
                .groups = "drop"
            )

    saveRDS(summary_perennial, "/media/DATADRIVE1/Project/muresk/model/d04_output/muresk_crop_Cinputs_MC_summary_perennial.rds")
    cat("Perennial crop summary statistics saved.\n")
}


annual_summary <- readRDS("/media/DATADRIVE1/Project/muresk/model/d04_output/cinputs_update/5542.rds")

perennial_summary <- readRDS("/media/DATADRIVE1/Project/muresk/model/d04_output/cinputs_update/5652.rds")


library(ggplot2)
library(dplyr)
library(tidyr)

plot_annual_uncertainty <- function(annual_summary, site_name, variable = "plant_residue_returns") {
    site_data <- annual_summary %>%
        filter(site == site_name)

    if (nrow(site_data) == 0) {
        stop(paste("No data found for site:", site_name))
    }

    mean_col <- paste0(variable, "_mean")
    q025_col <- paste0(variable, "_q025")
    q975_col <- paste0(variable, "_q975")
    q25_col <- paste0(variable, "_q25")
    q75_col <- paste0(variable, "_q75")

    p <- ggplot(site_data, aes(x = date)) +
        geom_ribbon(aes(ymin = .data[[q025_col]], ymax = .data[[q975_col]]),
            fill = "lightblue", alpha = 0.3
        ) +
        geom_ribbon(aes(ymin = .data[[q25_col]], ymax = .data[[q75_col]]),
            fill = "blue", alpha = 0.3
        ) +
        geom_line(aes(y = .data[[mean_col]]), color = "darkblue", size = 1) +
        theme_minimal() +
        labs(
            title = paste("Annual Crop:", variable, "-", "5542"),
            #title = paste("Perennial grass:", variable, "-", "5622"),
            subtitle = "Mean (blue), IQR (dark blue), 95% CI (light blue)",
            x = "Date",
            y = variable
        )

    return(p)
}


plot_annual_uncertainty(annual_summary[367:9132,],"5542")


# Also create a combined summary for common variables
cat("Creating combined summary statistics...\n")
common_vars <- c("dung_returns", "root_growth", "grazing_offtake")
common_vars <- common_vars[common_vars %in% names(combined_results_full)]

if (nrow(combined_results_full) > 0 && length(common_vars) > 0) {
    summary_combined <- combined_results_full %>%
        group_by(site, landuse, cropname, date, julianday, month, year) %>%
        summarise(
            across(all_of(common_vars),
                list(
                    mean = ~ mean(., na.rm = TRUE),
                    sd = ~ sd(., na.rm = TRUE),
                    q025 = ~ quantile(., 0.025, na.rm = TRUE),
                    q25 = ~ quantile(., 0.25, na.rm = TRUE),
                    median = ~ median(., na.rm = TRUE),
                    q75 = ~ quantile(., 0.75, na.rm = TRUE),
                    q975 = ~ quantile(., 0.975, na.rm = TRUE)
                ),
                .names = "{.col}_{.fn}"
            ),
            .groups = "drop"
        )

    saveRDS(summary_combined, "/media/DATADRIVE1/Project/muresk/model/d04_output/muresk_crop_Cinputs_MC_summary.rds")
    cat("Combined summary statistics saved.\n")
} else {
    cat("No combined summary statistics generated.\n")
}

num_sites_processed <- if (nrow(combined_results_full) > 0 && "site" %in% names(combined_results_full)) {
    length(unique(combined_results_full$site))
} else {
    0
}

num_iterations <- if (nrow(combined_results_full) > 0 && "mc_iteration" %in% names(combined_results_full)) {
    length(unique(combined_results_full$mc_iteration))
} else {
    0
}

cat("\nNumber of sites processed:", num_sites_processed, "\n")
cat("Number of MC iterations:", num_iterations, "\n")
cat("Monte Carlo analysis complete!\n")
