create_datasets <-
  function(y, country, start, end, train_percentage) {
    # Subset the data for the given country and time window
    y_subset <- window(y, start = start, end = end)[, country]
    
    # Calculate the number of observations for the training set
    n_train <- floor(length(y_subset) * train_percentage)
    
    # Calculate the number of observations for the test set (round up if necessary)
    n_test <- ceiling(length(y_subset) * (1 - train_percentage))
    
    # Create the training and test sets
    train <- head(y_subset, n_train)
    test <- tail(y_subset, n_test)
    
    # Return the training and test sets
    return(list(train = train, test = test))
  }

explore_data <- function(train) {
  # Apply the centered moving average
  cma <- cmav(train, outplot = TRUE)
  
  # Plot the seasonal subseries
  seasplot(train)
  seasplot(train, outplot = 2)
  seasplot(train, outplot = 3)
  seasplot(head(train, 4 * 4))
  seasplot(tail(train, 4 * 4))
  
  # Decompose the time series
  dc <- decomp(train, outplot = TRUE)
  
  # Since all the plots are displayed in the console, there's no need to return anything
}

auto_model <- function(train) {
  # Fit the model
  fit <- ets(train)
  
  # Extract the model name and whether it's damped from fit$components
  model_name <- paste(fit$components[1:3], collapse = "")
  is_damped <- as.logical(fit$components[4])
  
  model <- list(model = model_name, damped = is_damped)
  
  # Return the model name and whether it's damped
  return(model)
}

aicc_model <- function(train) {
  # Define the models to fit
  models <- list(
    list(model = "ANN", damped = FALSE),
    list(model = "ANA", damped = FALSE),
    list(model = "AAN", damped = FALSE),
    list(model = "AAN", damped = TRUE),
    list(model = "AAA", damped = FALSE),
    list(model = "AAA", damped = TRUE)
  )
  
  # Initialize a vector to store the AICc of each model
  aicc <- numeric(length(models))
  
  # Loop over the models
  for (i in seq_along(models)) {
    # Fit the model
    fit <-
      ets(train,
          model = models[[i]]$model,
          damped = models[[i]]$damped)
    
    # Store the AICc of the model
    aicc[i] <- fit$aicc
  }
  
  names(aicc) <- construct_model_names(models)
  
  print(aicc)
  
  # Find the model with the minimum AICc
  min_aicc_index <- which.min(aicc)
  
  model <-
    list(model = models[[min_aicc_index]]$model, damped = models[[min_aicc_index]]$damped)
  
  # Return the model name and whether it's damped
  return(model)
}




mape_model <- function(train, ins_percentage, selected_model) {
  # Calculate the number of observations for the in-sample set
  n_ins <- round(length(train) * ins_percentage)
  
  # Create the in-sample and validation sets
  ins <- head(train, n_ins)
  val <- tail(train,-n_ins)
  
  # Define the models to fit
  models <- list(
    list(model = "ANN", damped = FALSE),
    list(model = "ANA", damped = FALSE),
    list(model = "AAN", damped = FALSE),
    list(model = "AAN", damped = TRUE),
    list(model = "AAA", damped = FALSE),
    list(model = "AAA", damped = TRUE),
    selected_model
  )
  
  # Initialize a vector to store the MAPE of each model
  mape <- numeric(length(models))
  
  # Loop over the models
  for (i in seq_along(models)) {
    # Fit the model
    fit <-
      ets(ins, model = models[[i]]$model, damped = models[[i]]$damped)
    
    # Forecast the validation set
    frc <- forecast(fit, h = length(val))
    
    # Calculate the MAPE
    mape[i] <- mean(abs((frc$mean - val) / val)) * 100
  }
  names(mape) <- construct_model_names(models)
  
  print(mape)
  
  # Find the model with the minimum MAPE
  min_mape_index <- which.min(mape)
  
  # Return the model name and whether it's damped
  return(list(model = models[[min_mape_index]]$model, damped = models[[min_mape_index]]$damped))
}

ROCV_Model <- function(train, ins_percentage, selected_model, h) {
  
  # Check if selected_model is a list and has the necessary components
  if (is.null(selected_model) || !is.list(selected_model)) {
    stop("selected_model must be a list with components 'model' and 'damped'")
  }
  
  # Calculate the number of observations for the in-sample set
  n_ins <- round(length(train) * ins_percentage)
  
  # Create the in-sample and validation sets
  ins <- head(train, n_ins)
  val <- tail(train, -n_ins)
  
  # Define the models to fit
  models <- list(
    list(model = "ANN", damped = FALSE),
    list(model = "ANA", damped = FALSE),
    list(model = "AAN", damped = FALSE),
    list(model = "AAN", damped = TRUE),
    list(model = "AAA", damped = FALSE),
    list(model = "AAA", damped = TRUE)
  )
  
  
  for (i in seq_along(selected_model)) {
    new_model <- selected_model[[i]]
    
    # Check if the model already exists in models
    if (!any(sapply(models, function(x) x$model == new_model$model && x$damped == new_model$damped))) {
      models <- c(models, list(new_model))
    }
  }
  
  # Initialize a matrix to store the MAPE of each model at each forecast origin
  err <- matrix(NA, nrow = length(val) - h + 1, ncol = length(models)+ 1)
  
  # Loop over each forecast origin
  for (o in 1:(length(val) - h + 1)) {
    # Update the in-sample and validation sets
    ins_o <- head(train, length(ins) - 1 + o)
    val_o <- tail(train, length(val)- o + 1)
    
    # Loop over the models
    for (m in seq_along(models)) {
      # Fit the model
      fit <- ets(ins_o, model = models[[m]]$model, damped = models[[m]]$damped)
      
      # Forecast the validation set
      frc <- forecast(fit, h = h)$mean
      
      # Calculate the MAPE
      err[o, m] <- mean(abs((frc - head(val_o, h)) / head(val_o, h))) * 100
    }
    
    # Calculate the MAPE for the naive model
    frc_naive <- tail(ins_o, h)[1:h]
    err[o, length(models)+1] <- mean(abs((frc_naive - head(val_o, h)[1:h]) / head(val_o, h)[1:h]))
  }
  
  # Calculate the average MAPE for each model
  err_mean <- colMeans(err, na.rm = TRUE)
  names(err_mean) <- c(construct_model_names(models),"Naive")
  
  print(err_mean)
  boxplot(err_mean)
  
  if (length(models) <= which.min(err_mean)) {
    # Return an empty list
    return(list())
  } else {
    # Find the model with the minimum average MAPE
    min_err_mean_index <- which.min(err_mean)
    
    # Return the model name and whether it's damped
    return(list(model = models[[min_err_mean_index]]$model, damped = models[[min_err_mean_index]]$damped))
  }
}




ROCV_Model_Out <- function(train, train_percentage, selected_model, h) {
  n_ins <- round(length(train) * train_percentage)
  ins <- head(train, n_ins)
  val <- tail(train, -n_ins)
  
  omaxTest <- length(val) - h + 1
  errTest <- array(NA, c(omaxTest, length(selected_model) + 3))
  frcsTest <- array(NA, c(h, length(selected_model) + 3))
  
  for (o in 1:omaxTest) {
    ins_o <- head(train, length(ins) - 1 + o)
    val_o <- head(val, h + o - 1)
    
    # Fit and forecast the selected models
    for (m in 1:length(selected_model)) {
      fitTemp <- ets(ins_o, model = selected_model[[m]]$model, damped = selected_model[[m]]$damped)
      frcsTest[, m] <- forecast(fitTemp, h = h)$mean
      errTest[o, m] <- round(MAPE(frcsTest[, m], val_o[1:h]), digits = 3) * 100
      
    }
    
    # Forecast using the seasonal naive
    frcsTest[, length(selected_model) + 1] <- tail(ins_o, h)[1:h]
    errTest[o, length(selected_model) + 1] <- round(MAPE(frcsTest[, length(selected_model) + 1], val_o[1:h]), digits = 3) * 100
    
    # Combinations: Mean
    frcsTest[, length(selected_model) + 2] <- apply(frcsTest[, 1:(length(selected_model) + 1)], 1, mean)
    errTest[o, length(selected_model) + 2] <- round(MAPE(frcsTest[, length(selected_model) + 2], val_o[1:h]), digits = 3) * 100
    
    # Combinations: Median
    frcsTest[, length(selected_model) + 3] <- apply(frcsTest[, 1:(length(selected_model) + 1)], 1, median)
    errTest[o, length(selected_model) + 3] <- round(MAPE(frcsTest[, length(selected_model) + 3], val_o[1:h]), digits = 3) * 100
  }
  
  colnames(errTest) <- c(construct_model_names(selected_model), "Naive", "Comb.Mean", "Comb.Median")
  errTestMean <- round(colMeans(errTest), digits = 3)
  print(errTestMean)
  boxplot(errTest)
}

find_weighted_combination <- function(train, selected_model, h) {
  fit <- list()
  frc <- matrix(NA, nrow = h, ncol = length(selected_model), dimnames = list(NULL, selected_model))
  
  for (i in seq_along(selected_model)) {
    model <- selected_model[[i]]
    fit[[i]] <- ets(train, model = model$model, damped = model$damped)
    frc[, i] <- forecast(fit[[i]], h = h)$mean
  }
  
  AIC <- unlist(lapply(fit, function(x) x$aic))
  dAIC <- AIC - min(AIC)
  dAIC <- exp(-0.5 * dAIC)
  waic <- dAIC / sum(dAIC)
  
  frcComb <- frc %*% waic
  frcComb <- cbind(frcComb, rowMeans(frc))
  frcComb <- cbind(frcComb, apply(frc, 1, median))
  frcComb <- cbind(frcComb, frc[, which.min(AIC)])
  colnames(frcComb) <- c("Comb.AIC", "Comb.Mean", "Comb.Median", "Selection")
  
  return(frcComb)
}

calculate_MAPE <- function(forecasts, test_data) {
  err <- forecasts - matrix(rep(test_data, ncol(forecasts)), ncol = ncol(forecasts))
  MAPE <- colMeans(abs(err / matrix(rep(test_data, ncol(forecasts)), ncol = ncol(forecasts))) * 100)
  return(round(MAPE, 2))
}

# Define a custom function to append a model to selected_model if it doesn't exist
append_unique_model <- function(new_model, selected_model) {
  
  if (length(selected_model) == 0) {
    # If selected_model is empty, simply assign the new_model
    selected_model <- list(new_model)
  } else if (length(new_model) > 0) {
    # Check if the combination already exists in selected_model
    combination_exists <- FALSE
    for (i in seq_along(selected_model)) {
      existing_model <- selected_model[[i]]
      if (identical(existing_model$model, new_model$model) &&
          identical(existing_model$damped, new_model$damped)) {
        combination_exists <- TRUE
        break
      }
    }
    
    # If the combination doesn't exist, append it to selected_model
    if (!combination_exists) {
      selected_model <- c(selected_model, list(new_model))
    }
  }
  
  return(selected_model)
}


# Function to construct model names
construct_model_names <- function(models) {
  model_names <- character(length(models))
  
  for (i in seq_along(models)) {
    model <- models[[i]]
    if (model$damped) {
      model_name <-
        paste0(substr(model$model, 1, 2),
               "d",
               substr(model$model, 3, nchar(model$model)))
    } else {
      model_name <- model$model
    }
    model_names[i] <- model_name
  }
  
  return(model_names)
}
