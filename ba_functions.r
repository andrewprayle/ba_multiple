get_anova_table <- function(x, g, verbose = F){
    ## x is the values
    ## g is the groups
    require(dplyr)
    if (length(x) != length(g)){
        return(" length(x) != length(g) ")
    }
    #########################################################
    # make the anova_table
    anova_table <- matrix(data = NA, ncol = 4, nrow = 3)
    anova_table <- as.data.frame(anova_table)
    colnames(anova_table) <- c("df", "SS", "MS", "F")
    rownames(anova_table) <- c("Treatment", 
                               "Residuals",
                               "Total")
    ########################################################
    # import the table and make the group_table
    
    fdata <- as.data.frame(cbind(x, g))
    if (verbose) print(table(fdata$g))
    group_table <- get_group_table(x, g, verbose = verbose)
    group_table$A_i <- group_table$group_mean - mean(fdata$x)
    
    ########################################################
    # calculate the df
    anova_table$df[1] <- length(unique(fdata$g)) - 1
    anova_table$df[3] <- nrow(fdata) - 1
    anova_table$df[2] <- anova_table$df[3] - anova_table$df[1]
    
    if (verbose) print(anova_table)
    
    ########################################################
    anova_table$SS[1] <- sum((group_table$A_i)^2 * group_table$m)
    fdata$yij_yi <- NA
    # if (verbose) print(fdata)
    for (i in 1:nrow(fdata)){
        # get the squared residual
        fdata$yij_yi[i] <- (fdata$x[i] - group_table$group_mean[
                            which(group_table$g == fdata$g[i])])
    }
    anova_table$SS[2] <- sum(fdata$yij_yi^2)
    if (verbose) print(fdata)
    
    anova_table$SS[3] <- sum((fdata$x - mean(fdata$x))^2)
    
    anova_table$MS <- anova_table$SS / anova_table$df
    
    anova_table$F[1] <- anova_table$MS[1] /anova_table$MS[2]
    
    if (verbose) print(anova_table)
    return(anova_table)
}

get_group_table <- function(x, g, verbose = F){
  require(dplyr)
  ## x is the values
  ## g is the groups)
  ## returns a tibble of the group table
  
  ########################################################
  # import the table and make the group_table
  
  fdata <- as.data.frame(cbind(x, g))
  group_table <- fdata %>% group_by(g) %>% 
    summarize(group_mean = mean(x),
              m = length(x),
              group_sd = sd(x))
  group_table$A_i <- group_table$group_mean - mean(fdata$x)
  group_table$m_2 <- group_table$m^2
  if(verbose) print(group_table)
  return(group_table)
}



get_ba_result <- function(a, b, g, verbose = F, plotit = T){
  # master function, loads the data, calls the functions
  fdata <- as.data.frame(cbind(a, b, g))
  fdata$diff <- fdata$a - fdata$b
  
  # do the ba plot
  if (plotit)  plot(I(a - b) ~ I((a + b)/2), data = fdata,
       col = "white",
       xlab = "a - b",
       ylab = "average of a and b")
  if (plotit) text(x = (fdata$a + fdata$b)/2, 
       y = fdata$a - fdata$b, labels = fdata$g, adj = NULL,
       pos = NULL, offset = 0.5, vfont = NULL,
       cex = 0.5, col = NULL, font = NULL)  
  if (plotit) abline(h = 0)
  
  ## get the tables
  anova_table <- get_anova_table(fdata$diff, fdata$g, verbose = verbose)
  group_data <- get_group_table(fdata$diff, fdata$g, verbose = verbose)
  
  var_within_subject <- anova_table$MS[2]
  var_between_subject <- anova_table$MS[1] - anova_table$MS[2]
  divisor <- ((sum(group_data$m))^2 - sum(group_data$m_2))/
    ( (nrow(group_data) - 1) * sum(group_data$m) )
  
  est_component_of_variance <- var_between_subject / divisor
  est_component_of_variance
  
  tot_variance_single_diff <- est_component_of_variance + var_within_subject
  tot_variance_single_diff
  sd_single_diff <- sqrt(tot_variance_single_diff)
  sd_single_diff
  estimated_bias <- mean(data1$diff)
  estimated_bias
  
  upper_95_limit <- estimated_bias + 1.96 * sd_single_diff
  lower_95_limit <- estimated_bias - 1.96 * sd_single_diff
  
  abline(h = upper_95_limit)
  abline(h = lower_95_limit)
  
  results_table <- as.data.frame(matrix(data = NA, ncol = 2, nrow = 9))
  colnames(results_table) <- c("parameter", "value")
  results_table[1, 1] <- "Estimated variance of multiple between methods differences for the same subject"
  results_table[1, 2] <- var_within_subject
  results_table[2, 1] <- "Estimated variance for differences between the average difference across subjects"
  results_table[2, 2] <- var_between_subject
  results_table[3, 1] <- "Divisor"
  results_table[3, 2] <- divisor
  results_table[4, 1] <- "Estimated component of variance"
  results_table[4, 2] <- est_component_of_variance
  results_table[5, 1] <- "Estimated total variance for single differences on different subjects"
  results_table[5, 2] <- tot_variance_single_diff
  results_table[6, 1] <- "Estgimated sd for single difference on different subjects"
  results_table[6, 2] <- sd_single_diff
  results_table[7, 1] <- "Estimated bias"
  results_table[7, 2] <- estimated_bias
  results_table[8, 1] <- "upper 95% limit of agreement"
  results_table[8, 2] <- upper_95_limit
  results_table[9, 1] <- "lower 95% limit of agreement"
  results_table[9, 2] <- lower_95_limit
  
  return(results_table)
  
}

plot_sd <- function(a, b, g, verbose = F, plotit = T){
  fdata <- as.data.frame(cbind(a, b, g))
  
  fdata$diff <- fdata$a - fdata$b
  
  group_data <- get_group_table(fdata$diff, fdata$g)
  group_data$size <- 1 / max(group_data$m) * group_data$m
  
  # do the ba plot
  if (plotit)  plot(group_sd ~ group_mean, data = group_data,
                    col = "black",
                    xlab = "subject_mean_of_means",
                    ylab = "subject_sd_difference",
                    cex = group_data$size * 2)
  if (plotit) abline(h = 0)
}

plot_subject_averages <- function(a, b, g, verbose = F, plotit = T){
  fdata <- as.data.frame(cbind(a, b, g))
  fdata$x <- fdata$a - fdata$b
  
  group_data <- fdata %>% group_by(g) %>% 
    summarize(group_mean_diff  = mean(x),
              m                = length(x),
              group_sd         = sd(x),
              group_mean_means = mean((a + b)/2)
              )
  if (verbose) print(group_data)
  
  group_data$size <- 1 / max(group_data$m) * group_data$m
  
  if (plotit)  plot(group_mean_diff ~ group_mean_means, data = group_data,
                    col = "black",
                    xlab = "subject_mean_of_means",
                    ylab = "subject_mean_difference",
                    cex = group_data$size * 2)
  if (plotit) abline(h = 0)
}

standard_ba <- function(a, b, verbose = F, plotit = T){
  fdata <- as.data.frame(cbind(a, b))
  fdata$diff <- fdata$a - fdata$b
  fdata$mean <- (fdata$a + fdata$b)/2
  
  plot(diff ~ mean, data = fdata,
       xlab = "mean",
       ylab = "diff")
  
  bias <- mean(fdata$diff)
  upper <- bias + 1.96 * sd(fdata$diff)
  lower <- bias - 1.96 * sd(fdata$diff)
  results_table <- as.data.frame(matrix(data = NA, ncol = 2, nrow = 3))
  colnames(results_table) <- c("Parameter", "Result")
  results_table[1, 1] <- "Estimated bias"
  results_table[2, 1] <- "upper 95% limit of agreement"
  results_table[3, 1] <- "lower 95% limit of agreement "
  results_table[1, 2] <- bias
  results_table[2, 2] <- upper
  results_table[3, 2] <- lower
  
  return(results_table)
}
