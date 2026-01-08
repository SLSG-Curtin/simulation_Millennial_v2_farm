# helper functions for calculating boundires and c saturation deficits
#

library(snfa)

# get quartiles, %
q10 <- function(x) {
    quantile(x, 0.1, na.rm = TRUE)
}
q25 <- function(x) {
    quantile(x, 1 / 4, na.rm = TRUE)
}
q50 <- function(x) {
    quantile(x, 0.5, na.rm = TRUE)
}
q75 <- function(x) {
    quantile(x, 3 / 4, na.rm = TRUE)
}
q90 <- function(x) {
    quantile(x, 0.9, na.rm = TRUE)
}

#' apply boostratpping on a dataframe and return indices
#' df, dataframe
#' b, number of boottraps
#' 
bootrownames <- function(df, b) {
    # get the rownames of the dataframe, will be used in bootstrapping
    df_rnames <- rownames(df)

    bstrp_rnames <- list()

    for (i in seq(b)) {
        set.seed(i * 0019071970)
        bstrp_rnames_in <- sample(df_rnames, length(df_rnames), replace = T)
        bstrp_rnames_out <- df_rnames[!(df_rnames %in% bstrp_rnames_in)]

        bstrp_rnames[[paste0("boot", i)]] <- list(
            bootin = bstrp_rnames_in,
            bootout = bstrp_rnames_out
        )
    }

    return(bstrp_rnames)
}

#' fit a csat boundary based on X and y 
get_csat_boundary_yfit <- function(X, y, X_fit){
    # derive a 'reflected' data set to augment the data for fitting
    reflected_data <- snfa::reflect.data(X, y)
    X_eval <- reflected_data$X
    y_eval <- reflected_data$y

    # fitting of boundarty line
    frontier <- snfa::fit.boundary(
        X_eval, y_eval, X.bounded = X, y.bounded = y, X.constrained = X_fit,
        X.fit = X_fit, method = "mc", scale.constraints = T)

    # csat_boundary
    y_fit <- frontier$y.fit

    return(y_fit)
}

#' remove the points that are far from the boundary
#'
compress_Xy <- function(X, y, xstep = 0.1) {
    df <- data.frame(x = X, y = y)
    xbreaks <- seq(6.5, 77.0, by = xstep)

    X2 <- c()
    y2 <- c()

    for (i in seq_len(length(xbreaks))) {
        x <- xbreaks[i]
        xnext <- xbreaks[i + 1]

        subdf <- df[df$x >= x & df$x < xnext, ]

        if (nrow(subdf) != 0) {
            # keep only the top 3
            boundary_point <- subdf[subdf$y == max(subdf$y), ]
            # select first in case of multiple same values
            X2 <- c(X2, boundary_point$x[1])
            y2 <- c(y2, boundary_point$y[1])
        }
    }

    X2 <- as.matrix(X2)

    return(list(X = X2, y = y2))
}

#' calculate deficit
#' boundary_df, first column is x, second is y
#' csat_df, carbon saturation dataframe, first column (x, fine_fraction), second
#'  column (y, observed MAOC content)
get_csat_deficit <- function(boundary_df, csat_df) {
    # a vector to store deficit
    csat_df["ydefi"] <- NA

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
        csat_df[rname, "ydefi"] <- ybound_i - csat_df[rname, "y"]

    }

    return(csat_df)

}

#' fill dataframe na using modelling
#' treat complete comlumns as predictors, columns with missing values as response
#' df, a dataframe with na values
#' cols_complete, columns without missing values
#' cols_missing, columns with missing values

fill_missing <- function(
    df, cols_missing = NULL, cols_complete = NULL,  ncores = 10) {

    # set up for model tunning
    fit_control <- caret::trainControl(
        method = "cv",
        number = 10,
        savePredictions = F,
        verboseIter = F,
        allowParallel = T
    )
    hyper_grid <- expand.grid(
        committees = seq(1, 20, 2),
        neighbors = seq(2, 9, 2)
    )

    col_count <- 1
    for (coli in cols_missing) {
        print(
            paste0(
                "Imputing column (", col_count, "/", length(cols_missing),
                "): ", coli
            )
        )

        # get training data
        df_train <- df[
            complete.cases(df), c(coli, cols_complete)
        ]

        # train model on complete data
        clstr <- parallel::makeCluster(ncores)
        doParallel::registerDoParallel(clstr)
        model <- caret::train(
            x = df_train[cols_complete],
            y = df_train[, coli],
            "cubist",
            tuneGrid = hyper_grid,
            trControl = fit_control
        )
        parallel::stopCluster(clstr)
        foreach::registerDoSEQ()

        # save predictions
        x_pred <- df[is.na(df[coli]), cols_complete]
        # - used for assign predictions
        na_rownames <- rownames(x_pred)

        print("Before ...")
        print(df[na_rownames, c(coli, cols_complete[seq(min(5, length(cols_complete)))])])

        # predict missing values and fill the dataframe
        y_pred <- predict(model, x_pred)
        df[na_rownames, coli] <- y_pred

        print("")
        print(
            model$results[
                model$results$committees == model$bestTune$committees &
                model$results$neighbors == model$bestTune$neighbors,
            ]
        )
        print("")
        print("After ...")
        print(df[na_rownames, c(coli, cols_complete[seq(min(5, length(cols_complete)))])])

        print(paste0(rep("-", 150), collapse = ""))

        col_count <- col_count + 1
    }

    return(df)
}

# get complete column names from dataframe
complete_cols <- function(df) {
    cols_complete <- c()
    for (colname in colnames(df)) {
        if (sum(is.na(df[colname])) == 0) {
            cols_complete <- c(cols_complete, colname)
        }
    }

    return(cols_complete)
}