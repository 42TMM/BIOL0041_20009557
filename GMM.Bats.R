####################################################################################################################################################################################################################################


      #Running Generalised Mixed Effects Models on Bat Activity# 


####################################################################################################################################################################################################################################


library(dplyr)
library(lme4)
library(glmmTMB)
library(ggplot2)
library(DHARMa)
library(tidyr)
library(stringr)

bat_data <- read.csv("activity_df_60s.csv") %>%
  select(-duration)
head(bat_data)

###Removing 2 and 4####

bat_data <- bat_data %>%
  filter(!Cluster_32PCs %in% c("2", "4"))


#add conservancy column
bat_data <- bat_data %>%
  mutate(conservancy = str_extract(site, "^[A-Za-z]+"))

#Calculating lunar illumination for each site
library(lunar)
bat_data$date <- as.Date(bat_data$date, format="%Y-%m-%d")
bat_data$lunar_illumination <- lunar.illumination(bat_data$date)

# Step 2: Extract call count
bat_call_counts <- bat_data %>%
  group_by(site, date) %>%
  summarise(Call_Count = n(), .groups = "drop")

# Step 3: Merge with covariates 
covariates <- read.csv("AM_covs.csv")
covs_numeric<-covariates[ , c(5, 14, 19, 34, 35, 36)]
covs_numeric$cattle_30min_event_rate[is.na(covs_numeric$cattle_30min_event_rate)] <- 0
covs_numeric$shoat_30min_event_rate[is.na(covs_numeric$shoat_30min_event_rate)] <- 0
covs_numeric <- as.data.frame(scale(covs_numeric))
col_1 <- covariates[, 1, drop = FALSE]  
covs_numeric <- cbind(col_1, covs_numeric)
head(covs_numeric)


bat_data <- left_join(bat_data, covs_numeric, by = c("site" = "site")) %>%
  drop_na()

#Here we lose all data from site '2019'# 

bat_call_counts <- left_join(bat_call_counts, covs_numeric, by = c("site" = "site")) %>%
  drop_na()

head(bat_call_counts)

#add conservancy and lunar illumination
bat_call_counts <- bat_call_counts %>%
  mutate(conservancy = str_extract(site, "^[A-Za-z]+"))

#Calculating lunar illumination for each site
library(lunar)
bat_call_counts$date <- as.Date(bat_call_counts$date, format="%Y-%m-%d")
bat_call_counts$lunar_illumination <- lunar.illumination(bat_call_counts$date)



# Step 4: Fit Poisson GLMM #REMOVED CONSERVANCY AS LIKELY MASKING PATTERNS
poisson_model <- glmer(Call_Count ~ cattle_30min_event_rate + shoat_30min_event_rate + Mean.savi +
                         propopen500m + waterdist_short + humdist_short + lunar_illumination + 
                         (1|site) + (1|date),
                       data = bat_call_counts, family = poisson)
summary(poisson_model)

# Step 5: Fit Negative Binomial GLMM #REMOVED CONSERVANCY AS LIKELY MASKING PATTERNS
nb_model <- glmmTMB(Call_Count ~ cattle_30min_event_rate + shoat_30min_event_rate + Mean.savi +
                      propopen500m + waterdist_short + humdist_short + lunar_illumination + 
                      (1|site) + (1|date), 
                    data = bat_call_counts, family = nbinom2)
summary(nb_model)



AIC(poisson_model, nb_model)

# Step 6: Model diagnostics
simulation_poisson <- simulateResiduals(poisson_model)
simulation_nb <- simulateResiduals(nb_model)

plot(simulation_poisson)  # Check Poisson model residuals

plot(simulation_nb)  # Check NB model residuals

#Having added in lunar illumination and conservancy as a fixed factor the 
#binomial model is now very nice...


#NB model selected because: 
#All_bats = good 
#Cluster_0 is good
#Cluster 1 is good
#Cluster 2 is good
#Cluster 3 is good 
#Cluster 4 is good 
#Cluster 5 is good 
#Cluster 6 is good 
#but in poisson many are less good 

#################################################################################################################

#Now Run it For Clusters / Guilds swap 'guild' for 'Cluster_32PCs' 

bat_data_species <- filter(bat_data, Cluster_32PCs == "6")


# Step 2: Extract call count
bat_call_counts <- bat_data_species %>%
  group_by(site, date) %>%
  summarise(Call_Count = n(), .groups = "drop")

# Step 3: Merge with covariates 
covariates <- read.csv("covariates_2018_livestockrates_min2crops_30min_min3_threshold09_md09.csv")
covs_numeric<-covariates[ , c(7, 15, 19, 37, 43, 46)]
covs_numeric$cattle_30min_event_rate[is.na(covs_numeric$cattle_30min_event_rate)] <- 0
covs_numeric$shoat_30min_event_rate[is.na(covs_numeric$shoat_30min_event_rate)] <- 0
covs_numeric <- as.data.frame(scale(covs_numeric))
col_2 <- covariates[, 2, drop = FALSE]  
covs_numeric <- cbind(col_2, covs_numeric)
head(covs_numeric)


bat_call_counts <- left_join(bat_call_counts, covs_numeric, by = c("site" = "CT_site")) %>%
  drop_na()

#add conservancy and lunar illumination
bat_call_counts <- bat_call_counts %>%
  mutate(conservancy = str_extract(site, "^[A-Za-z]+"))

#Calculating lunar illumination for each site
library(lunar)
bat_call_counts$date <- as.Date(bat_call_counts$date, format="%Y-%m-%d")
bat_call_counts$lunar_illumination <- lunar.illumination(bat_call_counts$date)


# Step 4: Fit Poisson GLMM
poisson_model <- glmer(Call_Count ~ (1 | site) + cattle_30min_event_rate + shoat_30min_event_rate + Mean.savi + 
                         propopen500m + waterdist_short + humdist_short + lunar_illumination +
                         (1|date), 
                       data = bat_call_counts, family = poisson)
summary(poisson_model)

# Step 5: Fit Negative Binomial GLMM
nb_model <- glmmTMB(Call_Count ~ cattle_30min_event_rate + shoat_30min_event_rate + Mean.savi +
                      propopen500m + waterdist_short + humdist_short + lunar_illumination + factor(conservancy) +
                      (1|site) + (1|date), 
                    data = bat_call_counts, family = nbinom2)
summary(nb_model)


 AIC(poisson_model, nb_model)

# Step 6: Model diagnostics 
simulation_poisson <- simulateResiduals(poisson_model)
simulation_nb <- simulateResiduals(nb_model)

plot(simulation_poisson)  # Check Poisson model residuals

plot(simulation_nb)

#EXTRACT THE SUMMARY TABLE LATER#


#################################################################################################################

#For cluster 0 activity increases with humdist_short although these error messages:
  #Warning message:
    #In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
     #Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')

  #boundary (singular) fit: see help('isSingular')

#For cluster 1 no significant results except site, no error messages

#For cluster 2 activity increases with cattle grazing and distance to water although this error message: 
  #boundary (singular) fit: see help('isSingular')

#For cluster 3 activity increases with waterdist_short, no error messages

#For cluster 4 activity negatively correlated with shoat grazing and positively with humdist_short although this error message:
  #boundary (singular) fit: see help('isSingular')

#For cluster 5 activity negatively correlated with shoat grazing and positively with humdist_short although this error message:
  #boundary (singular) fit: see help('isSingular')

    #dropping columns from rank-deficient conditional model: cattle_30min_event_rate, shoat_30min_event_rate, propopen500m, waterdist_short, humdist_short
    #Warning message:
      #In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
      #Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')

#################################################################################################################

#Plotting effects

library(grid)
library(sjPlot)
library(gridExtra)
library(ggplot2)

summary(nb_model)

pred_structdist <- plot_model(nb_model, type = "pred", terms = "humdist_short", title = "", bias_correction = TRUE) +
  labs(y = "Bat activity (calls per night)", x = "Distance to human structure") +
  theme_classic() + 
  theme(text = element_text(size = 14))

pred_Meansavi <- plot_model(nb_model, type = "pred", terms = "Mean.savi", title = "", bias_correction = TRUE) +
  labs(y = "", x = "SAVI") +
  theme_classic() + 
  theme(text = element_text(size = 14))

pred_waterdist <- plot_model(nb_model, type = "pred", terms = "waterdist_short", title = "", bias_correction = TRUE) +
  labs(y = "", x = "Distance to water") +
  theme_classic() + 
  theme(text = element_text(size = 14))

pred_propopen <- plot_model(nb_model, type = "pred", terms = "propopen500m", title = "", bias_correction = TRUE) +
  labs(y = "", x = "Proportion of open habitat") +
  theme_classic() + 
  theme(text = element_text(size = 14))

pred_cattle_grazing <- plot_model(nb_model, type = "pred", terms = "cattle_30min_event_rate", title = "", bias_correction = TRUE) +
  labs(y = "Bat activity (calls per night)", x = "Cattle grazing") +
  theme_classic() + 
  theme(text = element_text(size = 14))

pred_shoat_grazing <- plot_model(nb_model, type = "pred", terms = "shoat_30min_event_rate", title = "", bias_correction = TRUE) +
  labs(y = "", x = "Shoat grazing") +
  theme_classic() + 
  theme(text = element_text(size = 14))

# Get the y-axis range from all plots
all_plots <- list(pred_structdist, pred_Meansavi, pred_waterdist, 
                  pred_propopen, pred_cattle_grazing, pred_shoat_grazing)

# Extract y ranges from all plots to find the overall min and max
y_ranges <- lapply(all_plots, function(p) {
  if (!is.null(ggplot_build(p)$layout$panel_params[[1]]$y.range)) {
    return(ggplot_build(p)$layout$panel_params[[1]]$y.range)
  } else {
    # Alternative approach for newer ggplot2 versions
    return(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  }
})

y_min <- min(sapply(y_ranges, min))
y_max <- max(sapply(y_ranges, max))

# Add a bit of padding
y_min <- y_min * 0.95
y_max <- y_max * 1.05

# Now apply consistent y-axis limits to all plots
pred_structdist <- pred_structdist + ylim(y_min, y_max)
pred_Meansavi <- pred_Meansavi + ylim(y_min, y_max)
pred_waterdist <- pred_waterdist + ylim(y_min, y_max)
pred_propopen <- pred_propopen + ylim(y_min, y_max)
pred_cattle_grazing <- pred_cattle_grazing + ylim(y_min, y_max)
pred_shoat_grazing <- pred_shoat_grazing + ylim(y_min, y_max)

# Arrange plots in a grid
grid.arrange(pred_structdist, pred_waterdist, pred_propopen, 
             pred_cattle_grazing, pred_shoat_grazing, pred_Meansavi,
             ncol = 3, top = textGrob("All Bats", gp = gpar(fontsize = 18)))



              #Printing nice table
              
              # Get the model summary
              summ <- summary(nb_model)
              
              # Extract the fixed effect coefficients
              coef_table <- as.data.frame(summ$coefficients$cond)
              
              # Move row names into a new column called 'Variable'
              coef_table$Variable <- rownames(coef_table)
              
              # Reorder so 'Variable' is the first column
              coef_table <- coef_table[, c("Variable", "Estimate", "Std. Error", "z value", "Pr(>|z|)")]
              
              coef_table$Estimate <- round(coef_table$Estimate, 3)
              coef_table$`Std. Error` <- round(coef_table$`Std. Error`, 3)
              coef_table$`z value` <- round(coef_table$`z value`, 3)
              coef_table$`Pr(>|z|)` <- round(coef_table$`Pr(>|z|)`, 5)
              
              # Add significance asterisks based on p-value thresholds
              coef_table$Significant <- ""
              coef_table$Significant[coef_table$`Pr(>|z|)` < 0.05] <- "*"
              coef_table$Significant[coef_table$`Pr(>|z|)` < 0.01] <- "**"
              coef_table$Significant[coef_table$`Pr(>|z|)` < 0.001] <- "***"
              
              library(flextable)
              # Create the flextable
              flextable(coef_table)
              

#################################################################################################################
#################################################################################################################

#This does exaclty what I want, I just need to understand and simplify the code a bit more
#To make it repeatable
              
              ###NOTE TO change between guild and clusters 'guild' for 'Cluster_32PCs' or vice versa, also...
              ###change 0/1/3/5/6 to Edge/Clutter/Open

              # List of predictors
              predictors <- c("cattle_30min_event_rate", "shoat_30min_event_rate", 
                              "propopen500m", "waterdist_short", "humdist_short", "Mean.savi")
              
              # Create empty lists to store models and predictions
              models_list <- list()
              predictions_list <- list()
              
              # Loop through guilds Edge, Clutter, Open
              for(guild_id in c("Edge", "Clutter", "Open")) {
                cat(paste("Processing guild", guild_id, "\n"))
                
                # Step 1: Filter data for the current guild
                bat_data_species <- filter(bat_data, guild == as.character(guild_id))
                
                # Step 2: Extract call count
                bat_call_counts <- bat_data_species %>%
                  group_by(site, date) %>%
                  summarise(Call_Count = n(), .groups = "drop")
                
                # Step 3: Merge with covariates
                covariates <- read.csv("covariates_2018_livestockrates_min2crops_30min_min3_threshold09_md09.csv")
                covs_numeric <- covariates[, c(7, 15, 19, 37, 43, 46)]  # Selected columns including Mean.savi
                
                # Handle NAs in livestock rates
                covs_numeric$cattle_30min_event_rate[is.na(covs_numeric$cattle_30min_event_rate)] <- 0
                covs_numeric$shoat_30min_event_rate[is.na(covs_numeric$shoat_30min_event_rate)] <- 0
                
                # Scale numeric covariates
                covs_numeric <- as.data.frame(scale(covs_numeric))
                
                # Add site column back
                covs_numeric <- cbind(covariates[, 2, drop = FALSE], covs_numeric)
                
                # Join with bat data and remove NAs
                bat_call_counts <- left_join(bat_call_counts, covs_numeric, by = c("site" = "CT_site")) %>%
                  drop_na()
                
                # Add conservancy and lunar illumination
                bat_call_counts <- bat_call_counts %>%
                  mutate(conservancy = str_extract(site, "^[A-Za-z]+"))
                
                bat_call_counts$date <- as.Date(bat_call_counts$date, format="%Y-%m-%d")
                bat_call_counts$lunar_illumination <- lunar.illumination(bat_call_counts$date)
                
                # Skip if insufficient data
                if(nrow(bat_call_counts) < 10) {
                  cat(paste("  Insufficient data for guild", guild_id, "- skipping\n"))
                  next
                }
                
                # Fit model with error handling
                tryCatch({
                  # Fit Negative Binomial GLMM
                  nb_model <- glmmTMB(Call_Count ~ cattle_30min_event_rate + shoat_30min_event_rate +
                                        propopen500m + waterdist_short + humdist_short + Mean.savi + 
                                        lunar_illumination + (1|site) + (1|date), 
                                      data = bat_call_counts, family = nbinom2)
                  
                  # Store the model
                  models_list[[as.character(guild_id)]] <- nb_model
                  
                  # Get predictions for each variable
                  for(pred_var in predictors) {
                    pred_data <- plot_model(nb_model, type = "pred", terms = pred_var, bias_correction = TRUE)
                    pred_df <- as.data.frame(pred_data$data)
                    pred_df$Cluster_32PCs <- as.character(guild_id)
                    
                    # Store predictions
                    if(is.null(predictions_list[[pred_var]])) {
                      predictions_list[[pred_var]] <- pred_df
                    } else {
                      predictions_list[[pred_var]] <- rbind(predictions_list[[pred_var]], pred_df)
                    }
                  }
                  
                  cat(paste("  Successfully processed guild", guild_id, "\n"))
                  
                }, error = function(e) {
                  cat(paste("  Error in guild", guild_id, ":", e$message, "\n"))
                })
              }
              
              # Nice labels for the variables
              var_labels <- c(
                "cattle_30min_event_rate" = "Cattle grazing",
                "shoat_30min_event_rate" = "Shoat grazing",
                "propopen500m" = "Proportion of open habitat",
                "waterdist_short" = "Distance to water",
                "humdist_short" = "Distance to human structure",
                "Mean.savi" = "SAVI"
              )
              
              # Define the order of plots
              plot_order <- c("humdist_short", "waterdist_short", "propopen500m", 
                              "cattle_30min_event_rate", "shoat_30min_event_rate", "Mean.savi")
              
              # Create a color palette for the guilds
              guild_colors <- c(
                "Edge" = "forestgreen",
                "Clutter" = "purple",
                "Open" = "orange"
                #"5" = "dodgerblue",
                #"6" = "brown"
              )
              
              # Set y-axis limits for consistency across plots
              y_min <- 0
              y_max <- 100
              
              # Create individual plots
              plot_list <- list()
              
              for(pred_var in plot_order) {
                if(!is.null(predictions_list[[pred_var]])) {
                  # Rename columns for easier access
                  colnames(predictions_list[[pred_var]])[colnames(predictions_list[[pred_var]]) == "x"] <- pred_var
                  colnames(predictions_list[[pred_var]])[colnames(predictions_list[[pred_var]]) == "predicted"] <- "predicted_calls"
                  
                  # Create plot
                  p <- ggplot(predictions_list[[pred_var]], 
                              aes_string(x = pred_var, y = "predicted_calls", color = "Cluster_32PCs", 
                                         fill = "Cluster_32PCs")) +
                    geom_ribbon(aes_string(ymin = "conf.low", ymax = "conf.high"), alpha = 0.2, color = NA) +
                    geom_line(linewidth = 1) +
                    scale_color_manual(values = guild_colors, name = "Guild") +
                    scale_fill_manual(values = guild_colors, name = "Guild") +
                    ylim(y_min, y_max) +
                    labs(y = if(pred_var %in% c("humdist_short", "cattle_30min_event_rate")) "Bat activity (calls per night)" else "",
                         x = var_labels[pred_var]) +
                    theme_classic() +
                    theme(text = element_text(size = 14),
                          legend.position = "none")  # Remove legend from all plots
                  
                  plot_list[[pred_var]] <- p
                }
              }
              
              # Create a function to extract legend
              get_legend <- function(a_gplot){
                tmp <- ggplot_gtable(ggplot_build(a_gplot))
                leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
                legend <- tmp$grobs[[leg]]
                return(legend)
              }
              
              # Create a simple plot just for legend extraction
              legend_plot <- ggplot(predictions_list[[plot_order[1]]], 
                                    aes_string(x = plot_order[1], y = "predicted_calls", color = "Cluster_32PCs")) +
                geom_line() +
                scale_color_manual(values = guild_colors, name = "Cluster") +
                theme_classic() +
                theme(legend.position = "right")
              
              # Extract legend
              legend <- get_legend(legend_plot)
              
              # Arrange plots in a grid
              combined_plots <- arrangeGrob(
                grobs = plot_list[plot_order],
                ncol = 3,
                nrow = 2
              )
              
              # Combine plots with legend on the right
              final_figure <- grid.arrange(
                combined_plots, 
                legend, 
                ncol = 2, 
                widths = c(10, 1),
                top = textGrob("Guilds", gp = gpar(fontsize = 18))
              )
              
              
              # Save the figure
              ggsave("all_predictor_effects_by_guild.png", final_figure, width = 13, height = 8, dpi = 300)
              
              # Print summary of models
              cat("\nSuccessfully fit models for guilds:", paste(names(models_list), collapse = ", "), "\n")
              
              # Print model summaries
              for(guild_id in names(models_list)) {
                cat(paste("\n\n----- MODEL SUMMARY FOR Guild", guild_id, "-----\n"))
                print(summary(models_list[[guild_id]]))
              }
              
              
              

              
              
              
############################################################################################################################################
