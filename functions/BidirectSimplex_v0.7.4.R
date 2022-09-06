####
#### Function: Bidirectional simplex projection
#### 2018.12.23 Ushio
#### R 3.4.4
####

require(rEDM)

# Define function
bidirect_simplex <- function(time_series, E=18,
                             complete_case_only = FALSE, # if TRUE, averaged predictions will be calculated only for complete pairs
                             lib = c(1, NROW(time_series)), pred = lib,
                             tau = 1,
                             num_neighbors = "e+1",
                             stats_only = TRUE,
                             silent = TRUE,
                             exclusion_radius = NULL,
                             epsilon = NULL){
  # Perform forward and backward simplex projection
  simplex_forward <- simplex(time_series, E = E, tp = 1, lib = c(1, NROW(time_series)), pred = lib, tau = tau, num_neighbors = num_neighbors, stats_only = F, silent = silent)
  simplex_backward <- simplex(time_series, E = E, tp = -E, lib = c(1, NROW(time_series)), pred = lib, tau = tau, num_neighbors = num_neighbors, stats_only = F, silent = silent)
  
  # Match time indices of the model output
  time_shared <- as.numeric(na.omit(intersect(simplex_forward$model_output[[1]]$time,
                                              simplex_backward$model_output[[1]]$time)))
  # Extract model output
  time_forward <- simplex_forward$model_output[[1]][match(time_shared, simplex_forward$model_output[[1]]$time),]
  time_backward <- simplex_backward$model_output[[1]][match(time_shared, simplex_backward$model_output[[1]]$time),]
  
  # Averaging forward and backward predictions
  if(complete_case_only){
    averaged_prediction <- (time_forward$pred + time_backward$pred)/2
  }else{
    averaged_prediction <- apply(cbind(time_forward$pred, time_backward$pred), 1, function(x) mean(x, na.rm = T))
  }
  simplex_bidirect <- data.frame(time = time_forward$time,
                                 obs = time_forward$obs,
                                 pred = averaged_prediction)
  
  # Collecting parameters
  params_forward <- simplex_forward[1:16]
  params_backward <- simplex_backward[1:16]
  params_all <- cbind(data.frame(method = c("forward_simplex","backward_simplex")),
                     rbind(params_forward, params_backward))
  
  # Calculate prediction accuracy
  stats_forward <- simplex_forward[,c("num_pred", "rho", "mae", "rmse")]
  stats_backward <- simplex_backward[,c("num_pred", "rho", "mae", "rmse")]
  stats_bidirection <- compute_stats(simplex_bidirect$obs, simplex_bidirect$pred)
  stats_all <- cbind(data.frame(method = c("bidirect_simplex","forward_simplex","backward_simplex"), E = c(E, E, E)),
                     rbind(stats_bidirection, stats_forward, stats_backward))      

  # Return output
  if(stats_only){
    simplex_bidirect_all <- data.frame(stats_all)
  }else{
    simplex_bidirect_all <- list(params = params_all, model_output = simplex_bidirect, stats = stats_all)
  }
  return(simplex_bidirect_all)
}

compute_stats <- function(obs, pred){
  # computes performance metrics for how well predictions match observations
  # obs = vector of observations
  # pred = vector of prediction
  
  num_pred <- sum(is.finite(obs) & is.finite(pred))
  rho <- cor(obs, pred, use = "pairwise.complete.obs")
  mae <- mean(abs(obs-pred), na.rm = TRUE)
  rmse <- sqrt(mean((obs-pred)^2, na.rm = TRUE))
  return(data.frame(num_pred = num_pred, rho = rho, mae = mae, rmse = rmse))
}

# Calculate best E using bidirectional simplex projection
bestE_bidirect <- function(time_series, E_range, lib = c(1, NROW(time_series)), complete_case_only = F, criteria = "rmse", show_fig = F, save_stats = F){
  simp_res <- bidirect_simplex(time_series, E = 1, lib = lib, stats_only = T, complete_case_only = complete_case_only)[1,]
  for(i in E_range[!E_range == 1]) simp_res[i,] <- bidirect_simplex(time_series, E = i, lib = lib, stats_only = T, complete_case_only = complete_case_only)[1,]
  best_E <- simp_res[which.min(simp_res[,criteria]), "E"]
  
  if(show_fig) plot(simp_res$E, simp_res[,criteria], type = "b", ylab = criteria, xlab = "E") 
  if(save_stats){
    all_res <- list(E = best_E, stats = simp_res)
  }else{
    all_res <- best_E
  }
  return(all_res)
}