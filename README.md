# Simulation of Millennial Model v2 at Muresk Farm Scale

> A comprehensive guide to the simulation workflow, model structure, and parameters for the Millennial Model v2 at Muresk Farm.

---

## 📚 Table of Contents

- [Key Adaptations](#-key-adaptations)
- [Simulation Workflow](#-simulation-workflow)
- [Input Data & Forcing](#-input-data--forcing)
- [Model Calibration](#-model-calibration)
- [Sensitivity Analysis](#-sensitivity-analysis)
- [SEM Analysis](#-sem-analysis)
- [Data Summary](#-data-summary)
- [References](#-references)
- [Appendix](#-appendix)

---

## 💡 Key Adaptations

This simulation workflow builds upon the calibration framework established in rangelands calibration. Key adaptations for the Muresk Farm scale include:

| Aspect | Australian Rangelands | Muresk Farm |
| :--- | :--- | :--- |
| **Measured Fractions** | Two fractions (POM.AGG, MAOM) | Three fractions (POM, AGG, MAOM) |
| **Forcing Data Source** | ERA5 reanalysis | SILO weather station data |
| **Soil Temperature** | ERA5 soil temperature layers | Calculated from weather station data |
| **Soil Moisture** | ERA5 soil moisture layers | Soil water balance model |
| **C Inputs** | NPP from satellite products | Plant growth model with Monte Carlo |
| **Calibration Strategy** | Global, Bioregional, Site-specific | Direct site-specific calibration |
| **Sensitive Parameters** | 13 PAWN-identified parameters | Same 13 parameters reused |

---

## 🌍 Simulation Workflow

The workflow for Muresk Farm involves site-specific calibration using measured fractions and tailored forcing inputs.

```text
══════════════════════════════════════════════════════════════════════════════════
                    MURESK FARM SIMULATION WORKFLOW
══════════════════════════════════════════════════════════════════════════════════

┌─────────────────────────────────────────────────────────────────────────────────┐
│                           INPUT DATA PREPARATION                                │
│                                                                                 │
│   ┌───────────────────┐   ┌──────────────────────┐   ┌──────────────────────┐   │
│   │ Measured Fractions│   │   Soil Properties    │   │ Historical Crop Types│   │
│   │ (POCmac, POCmic,  │   │  (Clay+Silt, pH, BD, │   │  • Annuals           │   │
│   │  MAOC)            │   │   C:N ratios)        │   │  • Perennials        │   │
│   └─────────┬─────────┘   └──────────┬───────────┘   │  • Crop Params (TE,  │   │
│             │                        │               │    MaxDM, Sow/Harv)  │   │
│             │                        │               └──────────┬───────────┘   │
└─────────────┼────────────────────────┼──────────────────────────┼───────────────┘
              │                        │                          │
              ▼                        ▼                          ▼
┌────────────────────────┐  ┌────────────────────────┐  ┌────────────────────────┐
│   Daily Soil Forcing   │  │   Daily C Input Forcing│  │  Monte Carlo C Inputs  │
│                        │  │                        │  │  (Uncertainty)         │
│                        │  │                        │  │                        │
│ ┌────────────────────┐ │  │ ┌────────────────────┐ │  │ ┌────────────────────┐ │
│ │ • Soil Temp        │ │  │ │ • Annual Growth    │ │◄─┼─┤ • HI (0.2-0.6)     │ │
│ │ • Soil Moist       │ │  │ │ • Perennial Growth │ │  │ │ • Turnover Rates   │ │
│ └────────────────────┘ │  │ └────────────────────┘ │  │ │ • Grazing Frac     │ │
└─────────────┬──────────┘  └──────────┬─────────────┘  │ └────────────────────┘ │
              │                        │                └────────────────────────┘
              ▼                        ▼
      ┌────────────────────────────────────────────────────────┐
      │               COMBINED FORCING INPUTS                  │
      │   (Daily Soil Temp + Soil Moisture + C Inputs)         │
      └────────────────────────────┬───────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                        SITE-SPECIFIC CALIBRATION                                │
│  • Uses 13 PAWN-identified important parameters                                 │
│  • Three measured C fractions: POM, AGG, MAOM                                   │
│  • ODE solver with SCE-UA optimization                                          │
└─────────────────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                         SENSITIVITY ANALYSIS                                    │
│  • Analyzes effects of NPP, CUE, Clay+Silt, pH                                  │
│  • Three NPP regimes: 0.5×, 1.0×, 1.5×                                          │
│  • Parameter Ranges: Clay+Silt (16-36%), pH (4.0-6.4), CUE (0.2-0.8)            │
└─────────────────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                            SEM ANALYSIS                                         │
│  • Latent Factors:                                                              │
│     - Efficiency (CUE, ratio_constraints)                                       │
│     - Sorption (MAOM_Qmax, ClaySilt, pH)                                        │
│  • Multi-group analysis by NPP regime                                           │
└─────────────────────────────────────────────────────────────────────────────────┘

══════════════════════════════════════════════════════════════════════════════════
```

---

## 📊 Input Data & Forcing

### 1. Soil Data Sources
Unlike the regional rangelands model (two fractions), the Muresk farm model relies on **three measured C fractions** for initialization and parameterization:

| Variable | Description | Unit | Source |
| :--- | :--- | :--- | :--- |
| **POCmac** | Macro particulate organic carbon (>250 µm) | % | Measured (Range: 0.017-0.748, Mean: 0.090) |
| **POCmic** | Micro particulate organic carbon (53-250 µm) | % | Measured (Range: 0.037-0.951, Mean: 0.185) |
| **MAOC** | Mineral-associated organic carbon (<53 µm) | % | Measured (Range: 0.178-1.480, Mean: 0.583) |
| **Clay + Silt** | Texture fraction  | % | Measured (Range: 16.6-35.6, Mean: 23.0) |
| **pH** | Soil pH (CaCl₂) | --- | Measured (Range: 3.95-6.43, Mean: 5.16) |
| **Litter C:N** | Litter carbon to nitrogen ratio | --- | Measured |
| **Total Soil C:N** | Soil organic matter C:N ratio | --- | Measured |
| **BD** | Bulk density | mg cm -3 | Measured |
| **Csat (Qmax)** | MAOC sorption capacity | g C m⁻² | Viscarra Rossel et al. (2024) |

**Unit Conversion**: All fractions are converted to g C m⁻² using:
```r
# Example conversion for 0-30cm depth
compile_muresk_csat$pocmac030 <- with(compile_muresk_csat, POCmac_30/100 * BD_30 * 30 * 10000)
compile_muresk_csat$maoc030   <- with(compile_muresk_csat, MAOC_30/100 * BD_30 * 30 * 10000)
compile_muresk_csat$qmax030   <- with(compile_muresk_csat, Csat/100 * BD_30 * 30 * 10000)
```

### 2. Forcing Inputs Generation
**Script**: `d01_src/soiltempmoisture_daily.R`

The model uses **SILO weather station data** (minimum air temperature, maximum air temperature, rainfall) to calculate daily soil drivers at a depth of 0-30cm.

#### 2.1 Soil Temperature Calculation

The daily soil temperature ($T_{soil}$) is calculated using the following core function (Horton and Corkrey, 2011):

$$T_{soil}(d, z) = \beta_0 + \beta_{lat} \times Lat + \sin(d \times \beta_{seas} + \beta_{phase}) + \beta_{min} \times N_d + \beta_{max} \times M_d + \beta_{rain} \times R_d + \beta_{int} \times M_d \times R_d$$

#### 2.2 Soil Moisture Calculation

Calculated using a **Soil Water Balance Model** with function `calc_daily_sw_update` (Script: `d01_src/soiltempmoisture_daily.R`).

```text
══════════════════════════════════════════════════════════════════════════════════
                    DAILY SOIL WATER BALANCE MODEL (0-30 cm)
══════════════════════════════════════════════════════════════════════════════════

┌─────────────────────────────────────────────────────────────────────────────────┐
│                           1. INITIALIZATION (Day 0)                             │
├─────────────────────────────────────────────────────────────────────────────────┤
│                                                                                 │
│   ┌─────────────────────┐           ┌─────────────────────────────────────┐     │
│   │     SITE INPUTS     │           │         HYDRAULIC LIMITS            │     │
│   │ • Clay fraction     │──────────▶│ • TSMD_max   (Max Deficit Capacity) │     │
│   │ • Depths (0-30cm)   │           │ • Infil_max  (Max Infiltration Rate)│     │
│   │ • BD, PAWC (Map)    │           │ • Pot_drain  (Drainage Capacity)    │     │
│   └─────────────────────┘           │ • PAWC       (Plant Available Water)│     │
│                                     └─────────────────────────────────────┘     │
│                                                                                 │
│   INITIAL STATE: • Soil is dry (TSMD = TSMD_max)                                │
│                  • Initial water content ~ 10mm                                 │
└─────────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                           2. DAILY LOOP (Day i)                                 │
├─────────────────────────────────────────────────────────────────────────────────┤
│                                                                                 │
│   Yesterday's State (i-1)                             Weather Forcing (i)       │
│   (TSMD, SWC)                                         (Rain, Evap)              │
│         │                                                    │                  │
│         ▼                                                    ▼                  │
│   ┌───────────────────────┐                    ┌────────────────────────────┐   │
│   │   STATE UPDATE        │                    │      WATER BALANCE         │   │
│   │ TSMD_init ← Prev TSMD │                    │ ET_act ← Evap × StressFac  │   │
│   │ StressFac ← 1-TSMD/Max│                    │ Net ← Rain + Irrig - ET    │   │
│   └───────────┬───────────┘                    │ infil ← min(Net, MaxInfil) │   │
│               │                                └─────────────┬──────────────┘   │
│               └──────────────────────┬───────────────────────┘                  │
│                                      │                                          │
│                                      ▼                                          │
│                       ┌─────────────────────────────┐                           │
│                       │   UPDATE DEFICIT (TSMD)     │                           │
│                       │ TSMD ← TSMD_init - infil    │                           │
│                       └──────────────┬──────────────┘                           │
│                                      │                                          │
│                                      ▼                                          │
│                      ┌───────────────────────────────┐                          │
│                      │      IS SOIL SATURATED?       │                          │
│                      │         (TSMD < 0)            │                          │
│                      └───────┬───────────────┬───────┘                          │
│                         YES  │               │  NO                              │
│         (Surplus Water)      ▼               ▼      (Deficit Remains)           │
│       ┌────────────────────────┐           ┌────────────────────────┐           │
│       │ Drainage ← surplus     │           │ Drainage ← 0           │           │
│       │ Runoff   ← overflow    │           │ TSMD_final ← TSMD      │           │
│       │ TSMD_final ← 0         │           └───────────┬────────────┘           │
│       └──────────────┬─────────┘                       │                        │
│                      │                                 │                        │
│                      └───────────────┬─────────────────┘                        │
│                                      │                                          │
│                                      ▼                                          │
│                          ┌───────────────────────┐                              │
│                          │  Next Day Iteration   │                              │
│                          └───────────────────────┘                              │
└─────────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                           3. OUTPUT VARIABLES                                   │
├─────────────────────────────────────────────────────────────────────────────────┤
│                                                                                 │
│   PRIMARY OUTPUT (Model Input)              SECONDARY OUTPUTS (Plant Growth)    │
│   ┌───────────────────────────────┐         ┌─────────────────────────────────┐ │
│   │ SWC_final_r (mm)              │         │ • ET_actual_r (mm)              │ │
│   │         │                     │         │   (Actual Evapotranspiration)   │ │
│   │         ▼                     │         │                                 │ │
│   │ forc_sw ← SWC_final_r / 300   │         │ • PAW_shallow_r (mm)            │ │
│   │ (Volumetric Water m³/m³)      │         │   (Plant Available Water)       │ │
│   └───────────────────────────────┘         │                                 │ │
│                                             │ • Fallow_deep_r (mm)            │ │
│                                             │   (Deep soil storage)           │ │
│                                             └─────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────────┘

══════════════════════════════════════════════════════════════════════════════════
```

**Key Abbreviations**:
| Symbol | Description | Unit |
| :--- | :--- | :--- |
| TSMD | Top Soil Moisture Deficit | mm |
| PAWC | Plant Available Water Capacity | mm |
| SWC | Soil Water Content | mm |
| ET | Evapotranspiration | mm |
| Eo_frac | Crop coefficient (ET/Pan ratio) | - |

#### 2.3 Daily C input calculation

**Script**: `d01_src/cinputs_daily.R`

Carbon inputs from plants are simulated using a plant growth model based on transpiration efficiency (Lee et al., 2021; Unkovich et al., 2018) that accounts for uncertainties through Monte Carlo sampling.

##### 2.3.1 Plant Growth Models

For annual crops (e.g., wheat, barley), growth follows a sigmoidal function:

```r
# Water-limited dry matter increment
an_WLDM_increment <- ((ET_actual * transpiration_frac_s) + 
                      (fallow_deep * transpiration_frac_d)) * 
                      TE * fertility_scalar

# Cumulative DM with sigmoidal growth curve
daily_DM <- max_DM / (1 + exp(-(season_day - a * max_days) / ((b * sow_day) * max_days)))

# Daily increment
daily_growth <- daily_DM[i] - daily_DM[i-1]

# Root growth based on shoot-to-root ratio
root_growth <- daily_growth * shoot_root_ratio

# Grazing offtake and root exudation
grazing_offtake <- ifelse(daily_DM < min_shoot_DM, 0, daily_growth * grazing_frac)
roots_shed <- ifelse(daily_DM < min_shoot_DM, 0, root_growth * grazing_frac * grazed_root_shed)
```

For perennial vegetation, growth is controlled by water and temperature:

```r
# Water-limited transpiration
pe_WL_transp <- ET_actual * veg_cover * (1 - tree_cover)

# Nix (1981) Temperature Index
pe_TI_x <- abs(ifelse(temp < TI_min, (TI_opt - temp) / (TI_opt - TI_min),
                      (temp - TI_opt) / (TI_max - TI_opt)))
pe_TI <- ifelse(temp < TI_min | temp > TI_max, 0, get_pe_TI(pe_TI_x))

# Water-limited growth
pe_WL_growth <- min(pe_max_growth, pe_WL_transp * TE * TEc * pe_TI * fertility_scalar)

# Senescence under water stress
stress_condition <- pe_WL_transp < critical_transp | PAWC_frac < critical_SWI
shoot_stress_loss <- ifelse(stress_condition, shoots_post_grazing * (1 - die_off), 0)
```

Tracks shoot and root residue pools with decay:

```r
# Pool decay rates
shoot_decomp <- shoot_residue_pool * shoot_decay_rate
root_decomp <- root_residue_pool * root_decay_rate

# Daily turnover during growth
shoot_turnover <- ifelse(is_growing, shoots_post_grazing * shoot_turnover_rate, 0)

# Harvest day residue additions
if (is_harvest_day) {
    unharvested_shoots <- max(shoots_post_grazing - shoot_turnover, 0) * (1 - HI)
    shoot_residue_pool <- shoot_residue_pool + unharvested_shoots
    root_residue_pool <- root_residue_pool + roots_post_grazing
}

# Plant residue returns to soil (C inputs)
plant_residue_returns <- decomp_flux + daily_inputs
```

##### 2.3.2 Monte Carlo Simulation

To account for parameter uncertainty, **100 Monte Carlo iterations** are performed per site.

**Uncertain Parameters**:

| Parameter | Description | Range | Units |
| :--- | :--- | :--- | :--- |
| `HI` | Harvest Index | 0.2 - 0.6 | --- |
| `shoot_turnover_rate` | Daily shoot turnover rate | 0.001 - 0.005 | d⁻¹ |
| `shoot_decay_rate` | Shoot residue decay rate | 0.001 - 0.005 | d⁻¹ |
| `root_decay_rate` | Root residue decay rate | 0.001 - 0.005 | d⁻¹ |
| `grazed_root_shed` | Root exudation fraction | 0.4 - 0.7 | --- |
| `grazing_frac` | Fraction of growth grazed | 0.0 - 0.5 | --- |
| `min_shoot_DM` | Minimum shoot dry matter | 1200 - 2000 | kg ha⁻¹ |
| `shoot_root_ratio` | Shoot-to-root allocation | 0.1 - 0.5 | --- |

**Root Exudation Constraint**:

A key physiological constraint (Sanderman et al., 2009) ensures realistic root exudation rates:

```r
# Constraint: Root exudation = 7-15% of total shoot growth
validate_root_exudation_constraint <- function(result_data, min_ratio = 0.07, max_ratio = 0.15) {
    total_shoot_growth <- sum(growth_data$daily_growth, na.rm = TRUE)
    total_root_exudation <- sum(growth_data$roots_shed, na.rm = TRUE)
    actual_ratio <- total_root_exudation / total_shoot_growth
    return(actual_ratio >= min_ratio && actual_ratio <= max_ratio)
}
```

##### 2.3.3 Combining Forcing Inputs

The final forcing inputs file combines:
1. Daily Soil Temperature (from Horton model)
2. Daily Soil Moisture (from water balance model)
3. Daily Plant C Inputs (Monte Carlo mean or specific iteration)

---

## 🎯 Model Calibration

**Script**: `d01_src/muresk_calib_site.R`

We perform **direct site-specific calibration** using the 13 important parameters identified from PAWN sensitivity analysis.

### 1. Calibration Workflow

#### 1.1 Model Function

The Millennial v2 model is implemented as an ODE system. For the Muresk farm scale, the model incorporates a Carbon Use Efficiency (CUE) adjustment based on the soil C:N ratio to better represent nitrogen limitations on microbial activity.

**Key Model Function (`derivs.rangeland.MV3.adj`)**:

```r
derivs.rangeland.MV3.adj <- function(step.num, state, params_sceua, forc_st, forc_sw, forc_npp) {
    # ... environmental modifiers (moisture, temperature) ...
    
    # Temperature-dependent potential CUE
    CUE_tmp = cue_ref - cue_t * (forc_st(step.num) - tae_ref)
    
    # C:N ratio adjustment (Improved Mechanism)
    Adj_CUE = pmax(0, parameters$param_CN^(-0.6))
    CUE_realized = CUE_tmp * Adj_CUE
    
    # Thermodynamic constraint
    if (CUE_realized > 0.8) {
        stop("CUE exceeds 0.8 - invalid parameter combination")
    }
    
    # ... derivative calculations ...
    return(list(c(dPOM, dLMWC, dAGG, dMIC, dMAOM), 'CUE' = CUE_realized))
}
```

The realized CUE is calculated as:
$$CUE_{realized} = [CUE_{ref} - CUE_t \times (T_{soil} - T_{ref})] \times CN_{ratio}^{-0.6}$$
where $CN_{ratio}$ is the measured soil total C:N ratio.

#### 1.2 Initialization

Standard initialization uses measured soil carbon fractions for the three target pools, while microbial and labile pools are set to baseline values:

```r
# Site-specific initialization from measurements
#   POM  = measured POCmac (macro-particulate organic carbon, >250 µm)
#   AGG  = measured POCmic (micro-particulate organic carbon, 53-250 µm)
#   MAOM = measured MAOC   (mineral-associated organic carbon, <53 µm)
#   LMWC = 6 g C m⁻² (fixed initial value)
#   MIC  = 6 g C m⁻² (fixed initial value)
```

#### 1.3 Optimization

- **Algorithm**: SCE-UA (Shuffled Complex Evolution)
- **simulation**: 100-year simulation (36,500 days) with repeated annual forcing cycles to reach equilibrium state.
- **Cost Function**: Minimizes the sum of absolute deviations for individual pools and total SOC:
  $$J = |POM_{sim} - POM_{obs}| + |AGG_{sim} - AGG_{obs}| + |MAOM_{sim} - MAOM_{obs}| + |SOM_{sim} - SOM_{obs}|$$
- **Constraints**:
  - $MAOM / SOM \ge 0.4$ (Ensures realistic representation of mineral-associated carbon)
  - $LMWC / SOM \le 0.03$
  - $MIC / SOM \le 0.03$
  - $CUE \le 0.8$

### 4. Output

- Calibrated parameters per site
- Daily simulated C pools
- Realized CUE time series

---

## 📉 Sensitivity Analysis

**Script**: `d01_src/muresk_sensitivity_all.R`

Following calibration, we conduct sensitivity analysis to understand the influence of key drivers on MAOM dynamics.

### 1. Analyzed Drivers

| Driver | Description | Variation |
| :--- | :--- | :--- |
| **NPP** | Net Primary Productivity | 0.5×, 1.0×, 1.5× baseline |
| **CUE** | Carbon Use Efficiency | Model-derived from calibration |
| **Clay+Silt** | Texture fraction | Measured Range: 16.6 - 35.6 % (Mean: 23.0 %) |
| **pH** | Soil acidity | Measured Range: 3.95 - 6.43 (Mean: 5.16) |
| **CUE** | Carbon Use Efficiency | Model-derived (Constraint ≤ 0.8, Ref ~0.53) |


### 2. Output Variables

For each site: simulated soil C fractions at equilibrium state

---

## 🧩 SEM Analysis

**Script**: `d01_src/maom_influence_sem_analysis.R`

Structural Equation Modeling (SEM) is used to analyze the complex causal relationships between validated drivers and MAOM (Mineral-Associated Organic Matter).

### 1. Data Preparation

```r
# Per-site z-score normalization
zscore_cols <- c("MAOM", "MAOM_Qmax")
normalized_df <- maom_base %>%
    group_by(site) %>%
    mutate(across(all_of(zscore_cols), safe_z)) %>%
    ungroup()
```

### 2. Candidate SEM Models

Multiple model structures are tested to identify the best representation:

#### Option A3 (Preferred Model):
```r
model_A3 <- '
    Sorption =~ MAOM_Qmax + ClaySilt + pH
    CUE ~ ratio_constraints
    MAOM ~ Sorption + CUE
    Sorption ~~ CUE
'
```

#### Option B (Path Model):
```r
model_B <- '
    MAOM ~ pH + ClaySilt + CUE + ratio_constraints + MAOM_Qmax
'
```

### 3. Latent Factors

| Latent Factor | Indicators | Interpretation |
| :--- | :--- | :--- |
| **Efficiency** | CUE, ratio_constraints | Microbial use efficiency |
| **Sorption** | MAOM_Qmax, ClaySilt, pH | Mineral-Organo interaction |
| **Desorption** | 1 * pH | pH-dependent release |

### 4. Model Fitting

```r
# Multi-group SEM by NPP regime
fit <- lavaan::sem(model, 
                   data = sem_df, 
                   group = 'Npp_regime',  # npp_0.5, npp_1.0, npp_1.5
                   estimator = 'MLR',      # Robust ML estimator
                   missing = 'fiml')       # Full information ML for missing data
```

### 5. Fit Statistics Evaluated

| Statistic | Good Fit Criterion |
| :--- | :--- |
| CFI | > 0.95 |
| TLI | > 0.95 |
| RMSEA | < 0.06 |
| SRMR | < 0.08 |
| AIC/BIC | Lower is better |

### 6. Output

- Standardized path coefficients
- Direct and indirect effects on MAOM
- Comparison table across candidate models
- Visualization of significant pathways by NPP regime

---

## 📂 Data Summary

### Input Data Structure

| File/Directory | Description |
| :--- | :--- |
| `d02_data/compile_muresk_csat_all.txt` | Compiled soil C fractions and properties |
| `d02_data/forcing_inputs/` | Processed forcing input files |
| `d02_data/geo_croptype.csv` | Site locations and historical crop types |

### Output Data Structure

| File | Description |
| :--- | :--- |
| `d04_report/sitebysite_muresk_threefrac_ode_100y_*.rds` | Calibrated parameters per site |

---

## 🔎 References

- Abramoff, Rose Z., et al. "Improved global-scale predictions of soil carbon stocks with Millennial Version 2." Soil Biology and Biochemistry 164 (2022): 108466.
- Horton, Brian, and Ross Corkrey. "A weighted coefficient model for estimation of Australian daily soil temperature at depths of 5 cm to 100 cm based on air temperature and rainfall." Soil Research 49.4 (2011): 305-314.
- Lee, J., Viscarra Rossel, R.A., Zhang, M., Luo, Z., & Wang, Y.P. (2021). Assessing the response of soil carbon in Australia to changing inputs and climate using a consistent modelling framework.* Biogeosciences.
- Pianosi, F., & Wagener, T. (2018). "Distribution-based sensitivity analysis from a generic input-output sample." Environmental Modelling & Software.
- Sanderman, Jonathan, Ryan Farquharson, and Jeffrey Baldock. "Soil carbon sequestration potential: a review for Australian agriculture." (2009): viii+-80.
- Unkovich, M., Baldock, J., & Farquharson, R. (2018). Field measurements of bare soil evaporation and crop transpiration, and transpiration efficiency, for rainfed grain crops in Australia–A review.* Agricultural Water Management, 205, 72-80.
- Viscarra Rossel, R. A., et al. (2024). "How much organic carbon could the soil store? The carbon sequestration potential of Australian soil." Global Change Biology.
- Jeffrey, S.J., et al. (2001). "Using spatial interpolation to construct a comprehensive archive of Australian climate data." Environmental Modelling & Software, 16(4), 309-330. (SILO)

---

## 📂 Appendix

### Millennial v2 Model Equations

See "Calibration_Millennial_v2_Australian_Rangelands.md" Appendix for the complete model structure diagram and equations.
