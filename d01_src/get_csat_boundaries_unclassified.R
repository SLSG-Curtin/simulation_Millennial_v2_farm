# c sat on unclassified data
# using: carbon_finefraction_covarites.csv
# 
# Outline:
#   - Load data
#   - get boostrap
#   - Fit boundary lines
#   - Get csat 
# Notes:
#   - fit.boundary speed, restrict size of Xy and X_fit

# clear env
rm(list = ls())

# libs
library(ggplot2)
library(doParallel)
library(snfa)
# setup
pdir <- file.path("/media/DATADRIVE1/Data/muresk/data/")
sdir <- file.path(pdir, "d03_script")
ddir <- file.path(pdir, "d02_data/for_model")
odir <- file.path(pdir, "d03_output")

# source fucntions
source(file.path(sdir,"csat_functions.R"))

print("Preparing data ...")
#

df <- read.csv(file.path(ddir, "compile_muresk_csat_all.txt"))
df_muresk <- df[,c(2,9,11,12)]
df_muresk$clay_silt <- with(df_muresk,clay_30+silt_30)

names(df_muresk)[c(1,2,5)] <- c("site","MAOC_perc","fine_frac")
# remove samples without soil type
#df <- df[!(df$ASC_Soil_type == ""), ]
# need to reset rownames
# rownames(df) <- seq(1, nrow(df))
# 81,172
# sort df based on fine_frac, low to high
# df <- df[order(df$fine_frac), ]

# get boostrap
n_boots <- 30  # 30, 50
bstrp_rnames <- bootrownames(df_muresk, b = n_boots)

#
# X <- as.matrix(df$fine_frac)
# rownames(X) <- rownames(df)
# y <- df$MAOC_perc

# at what X values to get y
# X_fit <- seq(floor(min(X)), ceiling(max(X)), by = 0.5)
# X_fit <- seq(floor(min(df$fine_frac)), ceiling(max(df$fine_frac)), by = 0.5)
X_fit <- seq(7.5, 77, by = 0.5)
X_fit <- as.matrix(X_fit)

# get boundaries for bootstraps
cores <- n_boots + 2
clstr <- makeCluster(10, outfile = "")
registerDoParallel(clstr)

print("Fitting boundaries ...")

# get y_fit matrix
y_fit_mat <- foreach(
    i = seq_len(n_boots), .combine = cbind, .packages = c("snfa")
) %dopar% {
    # get bootstrap X and y
    bstrp_indices <- bstrp_rnames[[i]]$bootin
    df_boot <- df_muresk[bstrp_indices, ]
    X_boot <- df_boot$fine_frac
    y_boot <- df_boot$MAOC_perc

    # compress Xy
    Xy_boot_compr <- compress_Xy(X_boot, y_boot)
    X_boot_compr <- Xy_boot_compr$X
    y_boot_compr <- Xy_boot_compr$y

    print(
        paste0(
            "boost: ", i, " X dim: ", dim(X_boot_compr), " y length: ",
            length(y_boot_compr)
        )
    )

    # get boundary df
    y_fit <- get_csat_boundary_yfit(X_boot_compr, y_boot_compr, X_fit)

    return(y_fit)
}

stopCluster(clstr)
registerDoSEQ()

# put X_fit and y_fit together
X_fitdf <- data.frame(X_fit = X_fit)
y_fitdf <- as.data.frame(y_fit_mat)
colnames(y_fitdf) <- paste0("y_fit_boot", seq_len(n_boots))

# boundary
boundary_df <- cbind(X_fitdf, y_fitdf)
print("boundary dataframe generated!")

y_fit_mean <- apply(y_fitdf, FUN = mean, MARGIN = 1)
y_fit_min <- apply(y_fitdf, FUN = min, MARGIN = 1)
y_fit_max <- apply(y_fitdf, FUN = max, MARGIN = 1)
y_fit_sd <- apply(y_fitdf, FUN = sd, MARGIN = 1)
boundary_pdf <- data.frame(
    X_fit = X_fit,
    y_fit_min = y_fit_min,
    y_fit_mean = y_fit_mean,
    y_fit_max = y_fit_max,
    y_fit_sd = y_fit_sd
)

print("Making plots ...")

df_bd <- df[,c(2,10)]

names(df_bd)[1] <- "site"

df_muresk <- merge(df_muresk,df_bd,by="site")

#qmax.gC.kgsoil
df_muresk$Qmax <- with(df_muresk,fine_frac*0.86)
#qmax.gC.m2
df_muresk$Qmax_perc <- with(df_muresk,Qmax*BD_30*1000*0.3/(BD_30*30*100))

# make some plots
p <- ggplot() +
    geom_point(
        data = data.frame(x = df_muresk$fine_frac, y = df_muresk$MAOC_perc),
        mapping = aes(x, y),
        colour = "grey"
    ) +
    geom_point(
    data = data.frame(x = df_muresk$fine_frac, y = df_muresk$Qmax_perc),
    mapping = aes(x, y),
    colour = "blue"
    ) +
    geom_line(
        data = boundary_pdf,
        mapping = aes(
            x = X_fit,
            # y = y_fit_min
            y = y_fit_mean - 2 * y_fit_sd
        ),
        colour = "red", linetype = "dotted"
    ) +
    geom_line(
        data = boundary_pdf,
        mapping = aes(x = X_fit, y = y_fit_mean),
        colour = "red"
    ) +
    geom_line(
        data = boundary_pdf,
        # mapping = aes(x = X_fit, y = y_fit_max),
        mapping = aes(
            x = X_fit,
            # y = y_fit_max
            y = y_fit_mean + 2 * y_fit_sd
        ),
        colour = "red", linetype = "dotted"
    ) +
    labs(
        x = "Clay-silt (%)",
        y = "MAOC, %"
    ) +
    theme_bw() +
    theme(
        panel.grid = element_blank()
    )

p

#get satuation maoc
csat_df <- df_muresk[,c(5,2)]
names(csat_df) <- c("x","y")

boundary_pdf$y_fit_outline <- with(boundary_pdf,y_fit_mean+2 * y_fit_sd)

outline_df <- boundary_pdf[,c(1,6)]

names(outline_df) <- c("x","y")

saveRDS(outline_df,"/media/DATADRIVE1/Data/muresk/data/d02_data/for_model/muresk_outline.RDS")

csat_df <- as.data.frame(list("x"  = 16,"y"  = 1.8))

get_csat_bound <- function(boundary_df, csat_df) {
  # a vector to store deficit
  csat_df["ycsat"] <- NA
  
  # x step in the boundary_df
  xbound_step <- 0.5
  
  rnames <- rownames(csat_df)
  
  for (rname in rnames) {
    # get boundary at xi
    xi <- csat_df[rname, "x"]
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
    csat_df[rname, "ycsat"] <- ybound_i
    
  }
  
  return(csat_df)
  
}

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

csat_df_target <- get_csat_frombound(outline_df,csat_df)


csat_df_new <- get_csat_bound(outline_df,csat_df)

names(csat_df_new) <- c("clay_silt","MAOC_30","Csat")

df$clay_silt <- with(df,clay_30+silt_30)

df_csat <- merge(df[,c(2:5)],csat_df_new,by="clay_silt")

write.csv(df_csat,file.path(odir,"compiled_csat.csv"))

ggsave(file.path("/media/DATADRIVE1/Data/muresk/data//d05_fig/Csat_rplot.png"),p)

