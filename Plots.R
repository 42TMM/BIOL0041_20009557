###########################################################################################################################################


#CREATING SAVI MAP WITH PIE CHART#


###########################################################################################################################################

library(terra)
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(scatterpie)  
library(flextable)
library(webshot2)
library(pdftools)
library(magick)
library(dplyr)
library(flextable)
library(corrplot)
library(ggplot2)
library(sf)
library(akima)
library(gstat)  
library(sp)


# Load map data
map_data <- st_read("study_area_shapefile_OMCcombined.shp")

# Load Landsat bands (red and NIR)
red <- rast("LC08_L2SP_169061_20191022_20200825_02_T1_SR_B4.TIF")
nir <- rast("LC08_L2SP_169061_20191022_20200825_02_T1_SR_B5.TIF")

# Reflectance scaling
red <- red / 10000
nir <- nir / 10000

# Reproject map data to match the CRS of the raster
map_data_proj <- st_transform(map_data, crs(red))

# Check if the map and raster overlap
# Print extents to debug
print(ext(red))
print(st_bbox(map_data_proj))

# Masking and cropping based on the shapefile
# Add error checking
tryCatch({
  red_crop <- crop(red, vect(map_data_proj))
  red_mask <- mask(red_crop, vect(map_data_proj))
  nir_crop <- crop(nir, vect(map_data_proj))
  nir_mask <- mask(nir_crop, vect(map_data_proj))
  
  # Print info about masked rasters
  print(paste("Red masked cells:", global(red_mask, "notNA", na.rm=TRUE)))
  print(paste("NIR masked cells:", global(nir_mask, "notNA", na.rm=TRUE)))
  
}, error = function(e) {
  print(paste("Error in masking/cropping:", e$message))
})

# Calculate SAVI
L <- 0.5
savi <- ((nir_mask - red_mask) / (nir_mask + red_mask + L)) * (1 + L)

# Check if SAVI has data
print(paste("SAVI has data:", !is.null(savi) && !all(is.na(values(savi)))))

# Convert SAVI raster to a data frame
savi_df <- as.data.frame(savi, xy = TRUE, na.rm = TRUE)
colnames(savi_df) <- c("longitude", "latitude", "SAVI")

# Check if savi_df has rows
print(paste("SAVI data frame rows:", nrow(savi_df)))

# Load AM site coordinates
covariates <- read.csv("covariates_2018_livestockrates_min2crops_30min_min3_threshold09_md09.csv")

site_coords <- covariates[, c(2, 5, 6)]  # Assumes columns: CT_site, longitude, latitude

AM_sites <- read.csv("site_names.csv")

AM_site_coords <- site_coords %>% filter(CT_site %in% AM_sites$site)

# Transform AM sites to the CRS of the map data
AM_site_coords_sf <- st_as_sf(AM_site_coords, coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = st_crs(map_data_proj))

# Keep as SF object for plotting with geom_sf instead
# Alternatively, extract coordinates properly if using geom_point
AM_site_coords_points <- st_coordinates(AM_site_coords_sf) %>%
  as.data.frame() %>%
  rename(longitude = X, latitude = Y)

# Determine plot limits from the data
bbox <- st_bbox(map_data_proj)
xlim <- c(bbox["xmin"], bbox["xmax"])

ylim <- c(bbox["ymin"], bbox["ymax"])

# Plotting
savi_map <- ggplot() +
  geom_raster(data = savi_df, aes(x = longitude, y = latitude, fill = SAVI)) +
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  geom_sf(data = map_data_proj, fill = NA, color = "black", size = 0.5) +
  # Use geom_sf for the points to maintain correct projection
  geom_sf(data = AM_site_coords_sf, color = "black", size = 1, shape = 21, fill = "white") +
  # Let the coord_sf determine limits based on the data if your manual limits aren't working
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(title = "SAVI Index with AudioMoth Stations", fill = "SAVI")

savi_map


####################################################
#Creating my covariates matrix
####################################################


# Loading in AM sites
AM_sites <- read.csv("site_names.csv")

# Loading in old covariates (with updated grazing data)
covariates <- read.csv("covariates_2018_livestockrates_min2crops_30min_min3_threshold09_md09.csv")

# Loading in covariates csv
new_covs <- read.csv("site_cov_incomplete.csv")

# Load map data
map_data <- st_read("study_area_shapefile_OMCcombined.shp")

# Load Landsat bands (red and NIR)
red <- rast("LC08_L2SP_169061_20191022_20200825_02_T1_SR_B4.TIF")
nir <- rast("LC08_L2SP_169061_20191022_20200825_02_T1_SR_B5.TIF")

# Reflectance scaling
red <- red / 10000
nir <- nir / 10000

# Reproject map data to match the CRS of the raster
map_data_proj <- st_transform(map_data, crs(red))

# Adding new cattle_30min_event_rate and shoat_30min_event_rate calculated using CT data
cols_to_merge <- c("CT_site", "cattle_30min_event_rate", "shoat_30min_event_rate")
covs_subset <- covariates[, cols_to_merge]
names(covs_subset)[names(covs_subset) == "CT_site"] <- "site"
new_covs <- merge(new_covs, covs_subset, by="site", all.x=TRUE)

# Filter so it only has the AM_sites 
AM_covs <- new_covs %>% filter(site %in% AM_sites$site)

# Convert to sf object for spatial operations
AM_covs_sf <- st_as_sf(AM_covs, coords = c("longitude", "latitude"), crs = 4326)
# Transform to projected CRS if necessary (assuming map_data_proj exists)
AM_covs_sf <- st_transform(AM_covs_sf, crs = st_crs(map_data_proj))

# Create 500m buffers around each site
AM_buffers <- st_buffer(AM_covs_sf, dist = 500)

# Calculate SAVI using the formula adapted for Landsat 8
savi <- (nir - red) / (nir + red + 0.428) * (1.428)

# Convert buffers to terra vect format for extraction
buffer_vect <- vect(AM_buffers)

# Extract mean, min, and max SAVI for each buffer
savi_means <- terra::extract(savi, buffer_vect, fun = mean)
savi_mins <- terra::extract(savi, buffer_vect, fun = min)
savi_maxs <- terra::extract(savi, buffer_vect, fun = max)

# Add SAVI values to your dataframe
AM_covs$Mean.savi <- savi_means[, 2]
AM_covs$min_SAVI <- savi_mins[, 2]
AM_covs$max_SAVI <- savi_maxs[, 2]

# View the results
print(head(AM_covs))

# Save the updated covariates
write.csv(AM_covs, "AM_covs.csv", row.names = FALSE)





####################################################
#Adding Call counts at each AM to the map, coloured by guild
####################################################

bat_data <- read.csv("activity_df_60s.csv") %>%
  select(-duration)
head(bat_data)

###Removing 2 and 4####

bat_data <- bat_data %>%
  filter(!Cluster_32PCs %in% c("2", "4"))


#add conservancy column
bat_data <- bat_data %>%
  mutate(conservancy = str_extract(site, "^[A-Za-z]+"))

#Extract call count
bat_call_counts <- bat_data %>%
  group_by(site, date) %>%
  summarise(Call_Count = n(), .groups = "drop")

#add covariates
covariates <- read.csv("AM_covs.csv")
covs_numeric<-covariates[ , c(5, 14, 19, 34, 35, 36)]
covs_numeric$cattle_30min_event_rate[is.na(covs_numeric$cattle_30min_event_rate)] <- 0
covs_numeric$shoat_30min_event_rate[is.na(covs_numeric$shoat_30min_event_rate)] <- 0
covs_numeric <- as.data.frame(scale(covs_numeric))
col_1 <- covariates[, 1, drop = FALSE]
col_2 <- covariates[, 2, drop = FALSE]
col_3 <- covariates[, 3, drop = FALSE]
covs_numeric <- cbind(col_1, col_2, col_3, covs_numeric)
head(covs_numeric)
bat_data <- left_join(bat_data, covs_numeric, by = c("site" = "site")) %>%
  drop_na()


head(bat_data)



# Aggregate counts by site and guild
site_guild_counts <- bat_data %>%
  group_by(site, guild) %>%
  summarize(count = n(), .groups = "drop") %>%
  left_join(bat_data %>% 
              select(site, latitude, longitude) %>% 
              distinct(site, .keep_all = TRUE),
            by = "site")


# Pivot data for scatterpie
pie_data <- site_guild_counts %>%
  pivot_wider(id_cols = c(site, latitude, longitude),
              names_from = guild, values_from = count, values_fill = 0) %>%
  left_join(site_guild_counts %>%
              group_by(site) %>%
              summarize(total_count = sum(count), .groups = "drop"),
            by = "site")


#######
# First handle the SAVI map
map_data_proj <- st_transform(map_data, crs(red))

# Process the SAVI data as before
red_crop <- crop(red, vect(map_data_proj))
red_mask <- mask(red_crop, vect(map_data_proj))
nir_crop <- crop(nir, vect(map_data_proj))
nir_mask <- mask(nir_crop, vect(map_data_proj))

# Calculate SAVI
L <- 0.5
savi <- ((nir_mask - red_mask) / (nir_mask + red_mask + L)) * (1 + L)
savi_df <- as.data.frame(savi, xy = TRUE, na.rm = TRUE)
colnames(savi_df) <- c("longitude", "latitude", "SAVI")

# Get exact bounding box from the projected map data
bbox <- st_bbox(map_data_proj)
xlim <- c(bbox["xmin"], bbox["xmax"])
ylim <- c(bbox["ymin"], bbox["ymax"])

# Transform AM site coordinates to match the projection
AM_site_coords_sf <- st_as_sf(AM_site_coords, coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = st_crs(map_data_proj))

# SAVI Map - Using coord_sf with exact limits
savi_map <- ggplot() +
  geom_raster(data = savi_df, aes(x = longitude, y = latitude, fill = SAVI)) +
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  geom_sf(data = map_data_proj, fill = NA, color = "black", size = 0.5) +
  geom_sf(data = AM_site_coords_sf, color = "black", size = 1, shape = 21, fill = "white") +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  labs(fill = "SAVI", x = "Longitude", y = "Latitude") +
  theme(legend.position = "none")

# Now for the pie chart map
# First make sure the bat site data is properly transformed to the same projection
pie_data_sf <- st_as_sf(pie_data, coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = st_crs(map_data_proj))

# Extract transformed coordinates for scatterpie
pie_coords <- st_coordinates(pie_data_sf)
pie_data$plot_x <- pie_coords[, 1]
pie_data$plot_y <- pie_coords[, 2]

# Pie Chart Map - Using the SAME coord_sf with EXACT same limits
pie_map <- ggplot() +
  geom_sf(data = map_data_proj, fill = NA, color = "black", size = 0.5) +
  geom_scatterpie(data = pie_data, 
                  aes(x = plot_x, y = plot_y, r = sqrt(total_count) * 50), 
                  cols = unique(bat_data$guild), alpha = 0.8) +
  scale_fill_brewer(palette = "Set1", name = "Guild") +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = "none")

# Combine the plots with cowplot
combined_plots <- ggdraw() +
  draw_plot(savi_map) +  # Draw SAVI map first (bottom layer)
  draw_plot(pie_map) 

# Display the combined plot
combined_plots




####################################################################################################################################################################################################################################


#SITE SUMMARY TABLE# 


####################################################################################################################################################################################################################################


summary_table <- read.csv("site_summary_table.csv")
sum(summary_table$detections)


ft <- flextable(summary_table)


# Save as HTML first
save_as_html(ft, path = "my_flextable.html")

# Convert to PDF
webshot("my_flextable.html", "my_flextable.pdf")


# Path to the PDF
pdf_path <- "my_flextable.pdf"

# Convert PDF pages to images
imgs <- image_read_pdf(pdf_path, density = 300)  # 300 DPI for better quality

# Save each page as a separate PNG file
for (i in seq_along(imgs)) {
  image_write(imgs[i], path = paste0("page_", i, ".png"))
}


####################################################################################################################################################################################################################################


#SITE SUMMARY TABLE# 


####################################################################################################################################################################################################################################
#Other plots - HOW DO I INCORPORATE THESE INTO THE MAPPING FILE?? 

####AM SITES###
site_coords <- covariates[, c(1, 2, 3)]
AM_sites <- read.csv("site_names.csv")

AM_site_coords <- site_coords %>%
  filter(site %in% AM_sites$site)

sites_no_covs <- AM_sites %>%
  anti_join(AM_site_coords, by = c("site" = "site"))

# View the result
sites_no_covs

####################################################################################################################################################################################################################################


#Correlation Matrix#


####################################################################################################################################################################################################################################

# Load necessary package

# Subset data with selected covariates
covariate_data <- covs_numeric[, c("propopen500m", "waterdist_short", 
                                   "humdist_short", "Mean.savi", 
                                   "cattle_30min_event_rate", 
                                   "shoat_30min_event_rate")]

# Compute correlation matrix
cor_matrix <- cor(covariate_data, use = "pairwise.complete.obs")

# Print correlation matrix
print(cor_matrix)

# Optional: Visualize the correlation matrix
corrplot(cor_matrix, method = "color", type = "upper", tl.col = "black", tl.srt = 45)


####################################################################################################################################################################################################################################


#Detection Counts# 
#For All Analysis swap Cluster_32PCs for guild to see difference 
#Between Sonotype and Guild analyses


####################################################################################################################################################################################################################################

data <- bat_data

#Plotting Detection counts etc
head(data)

# Assuming your data is in a dataframe called `data`
# Convert `date` and `time` to a proper datetime format if needed
data$time_padded <- sprintf("%06d", as.numeric(data$time))
data$date_time <- as.POSIXct(paste(data$date, data$time_padded), format="%Y-%m-%d %H%M%S")

data["Presence"] = 1

####################################################################################################################################################################################################################################
#Detections by site
####################################################################################################################################################################################################################################

all_sites <- data %>%
  mutate(conservancy = str_extract(site, "^[A-Za-z]+")) %>%
  filter(!is.na(conservancy)) %>%
  distinct(site, conservancy)

detection_counts <- data %>%
  filter(!is.na(site)) %>%
  group_by(site) %>%
  summarise(detection_count = sum(Presence == 1, na.rm = TRUE), .groups = "drop")

detection_site <- left_join(all_sites, detection_counts, by = "site") %>%
  mutate(detection_count = replace_na(detection_count, 0))

ggplot(detection_site, aes(x = site, y = detection_count, fill = conservancy)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Detection Count Per Site",
       x = "Site",
       y = "Detection Count",
       fill = "Conservancy") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

####################################################################################################################################################################################################################################
#Calls at each conservancy#
####################################################################################################################################################################################################################################

#Table of Detections

#Create the base table with call counts per conservancy (only where presence == 1)
calls_per_conservancy <- data %>%
  filter(Presence == 1) %>%  # Filter for presence == 1
  mutate(conservancy = str_extract(site, "^[A-Za-z]+")) %>%
  filter(!is.na(conservancy)) %>%
  group_by(conservancy) %>%
  summarise(call_count = n(), .groups = "drop")

#Calculate the call count per Cluster_32PCs and Species (guild), and the percentages
calls_with_percentages <- data %>%
  filter(Presence == 1) %>%  # Filter for presence == 1
  mutate(conservancy = str_extract(site, "^[A-Za-z]+")) %>%
  filter(!is.na(conservancy)) %>%
  group_by(conservancy, Cluster_32PCs, guild) %>%
  summarise(call_count_per_cluster_guild = n(), .groups = "drop") %>%
  left_join(calls_per_conservancy, by = "conservancy") %>%
  mutate(percentage = (call_count_per_cluster_guild / call_count) * 100)

#Aggregate the percentages by guild (Clutter, Edge, Open)
guild_percentages <- calls_with_percentages %>%
  group_by(conservancy, guild) %>%
  summarise(guild_percentage = sum(percentage), .groups = "drop") %>%
  spread(key = guild, value = guild_percentage, fill = 0) %>%
  select(conservancy, Clutter, Edge, Open)

head(calls_with_percentages)

# Aggregate the cluster counts (0, 1, 2, ..., 6) and convert to percentages
cluster_counts <- calls_with_percentages %>%
  group_by(conservancy, Cluster_32PCs) %>%
  summarise(cluster_count = sum(call_count_per_cluster_guild) / first(call_count) *100, .groups = "drop") %>%
  spread(key = Cluster_32PCs, value = cluster_count, fill = 0) %>%
  select(conservancy, '0', '1', '3', '5', '6')

#Combine guild percentages and cluster counts into a single table
combined_table <- left_join(calls_per_conservancy, guild_percentages, by = "conservancy") %>%
  left_join(cluster_counts, by = "conservancy") %>%
  mutate(across(Clutter:`6`, ~ round(.x, 2)))  # Optionally round percentages and counts

#Convert to flextable and format the columns
calls_per_conservancy_flex <- flextable(combined_table)

# Renaming the columns for clarity
calls_per_conservancy_flex <- set_header_labels(calls_per_conservancy_flex, 
                                                conservancy = "Conservancy", 
                                                call_count = "Call Count", 
                                                Clutter = "Clutter %", 
                                                Edge = "Edge %", 
                                                Open = "Open %", 
                                                `0` = "Cluster 0 %", 
                                                `1` = "Cluster 1 %", 
                                                `2` = "Cluster 2 %", 
                                                `3` = "Cluster 3 %", 
                                                `4` = "Cluster 4 %", 
                                                `5` = "Cluster 5 %", 
                                                `6` = "Cluster 6 %")

#Set the column widths (adjust these as needed)
calls_per_conservancy_flex <- set_table_properties(calls_per_conservancy_flex, 
                                                   width = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

#Display the updated table
calls_per_conservancy_flex


####################################################################################################################################################################################################################################
#Diversity at each site
####################################################################################################################################################################################################################################
species_site <- data %>%
  filter(!is.na(Cluster_32PCs) & Cluster_32PCs != "") %>%  # <- exclude empty strings
  distinct(site, Cluster_32PCs) %>%
  mutate(conservancy = str_extract(site, "^[A-Za-z]+")) %>%  
  filter(!is.na(conservancy)) %>%  
  group_by(site, Cluster_32PCs, conservancy) %>%
  summarise(species_count = 1, .groups = "drop") 

ggplot(species_site, aes(x = site, y = species_count, fill = as.factor(Cluster_32PCs))) +
  geom_bar(stat = "identity") +
  facet_wrap(~ conservancy, scales = "free_x") +  
  theme_minimal() +
  labs(title = "Number of Sonotypes Per Site",
       x = "Site",
       y = "Sonotype Count",
       fill = "Sonotypes") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 10))

####################################################################################################################################################################################################################################
#Raw counts of each 'species'
####################################################################################################################################################################################################################################
df <- read.csv('merged_df.csv')

# Filter data based on time and site
df <- df %>%
  filter((time >= 170000) | (time <= 74500)) %>%
  filter(site != "2019") %>%
  filter(!Cluster_32PCs %in% c("2", "4"))

# Calculate species counts
species_counts <- df %>% 
  filter(Cluster_32PCs %in% c("0", "1", "3", "5", "6")) %>%
  group_by(Cluster_32PCs) %>%
  summarize(count = n())

#species_counts <- df %>% 
  #filter(guild %in% c("Clutter", "Edge", "Open")) %>%
  #group_by(guild) %>%
  #summarize(count = n())



ggplot(species_counts, aes(x = "", y = count, fill = as.factor(Cluster_32PCs))) +
  geom_bar(stat = "identity", width = 1, color = "white") +  
  coord_polar(theta = "y") +  
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), color = "white") +  
  theme_void() +  
  labs(title = "Sonotype Counts",
       fill = "Sontoypes")

####################################################################################################################################################################################################################################
#Time series of detections, coloured by 'species'
####################################################################################################################################################################################################################################
data$time <- sprintf("%06d", as.numeric(data$time))  # Ensure 6-digit times
data$date_time <- as.POSIXct(paste(data$date, data$time), format="%Y-%m-%d %H%M%S")

# Create all possible 15-min time intervals in a day
all_times <- format(seq(
  from = as.POSIXct("00:00:00", format="%H:%M:%S"),
  to   = as.POSIXct("23:45:00", format="%H:%M:%S"),
  by   = "15 min"
), "%H%M%S")

# Get unique species
all_species <- sort(unique(data$Cluster_32PCs[!is.na(data$Cluster_32PCs)]))

# Make a full grid of time × species
full_grid <- expand.grid(
  time_padded = all_times,
  Cluster_32PCs = all_species
)

# Sum detections
detection_15min <- data %>%
  filter(!is.na(Cluster_32PCs)) %>%
  group_by(time_padded, Cluster_32PCs) %>%
  summarise(detection_count = sum(Presence == 1, na.rm = TRUE), .groups = "drop")

# Merge with full grid and fill in zeros
detection_15min_full <- full_grid %>%
  left_join(detection_15min, by = c("time_padded", "Cluster_32PCs")) %>%
  mutate(detection_count = ifelse(is.na(detection_count), 0, detection_count))

# Plot
ggplot(detection_15min_full, aes(x = time_padded, y = detection_count, fill = as.factor(Cluster_32PCs))) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Total Detections at 15-Minute Intervals by Sonotype",
    x = "Time of Day (HHMMSS)",
    y = "Total Detection Count",
    fill = "Sonotype"
  ) +
  theme(axis.text.x = element_text(angle = 90, size = 6))

#######################################################################################
#Mapping other covariates#
#######################################################################################
covariates <- read.csv("AM_covs.csv")
covs_numeric<-covariates[ , c(5, 14, 19, 34, 35, 36)]
covs_numeric$cattle_30min_event_rate[is.na(covs_numeric$cattle_30min_event_rate)] <- 0
covs_numeric$shoat_30min_event_rate[is.na(covs_numeric$shoat_30min_event_rate)] <- 0
covs_numeric <- as.data.frame(scale(covs_numeric))
col_1 <- covariates[, 1, drop = FALSE]
col_2 <- covariates[, 2, drop = FALSE]
col_3 <- covariates[, 3, drop = FALSE]
covs_numeric <- cbind(col_1, col_2, col_3, covs_numeric)
head(covs_numeric)
#filter covs_numeric by cluster1_sites, cluster1_sites not found, should I use AM covs here then



# Read shapefile with consistent CRS
map_data <- st_read("study_area_shapefile_OMCcombined.shp")





#Preliminary plot

# Create base map with specified coordinate limits
base_map <- ggplot() +
  geom_sf(data = map_data, fill = "lightgrey", color = "black", alpha = 0.3) +
  coord_sf(xlim = c(34.75, 35.5), ylim = c(-1.6, -1.0)) +
  theme_minimal() +
  labs(
    title = "Audio Moth Stations",
    x = "Longitude",
    y = "Latitude"
  ) +
  geom_point(data = covs_numeric, 
             aes(x = longitude, y = latitude), 
             color = "black",  # Set fixed color outside aes()
             size = 0.1,
             alpha = 1
  )

base_map

#Need to add AM stations

####
#water

# Add point heatmap layer
water_heatmap <- base_map +
  geom_point(
    data = covs_numeric, 
    aes(x = longitude, y = latitude, color = waterdist_short), 
    size = 5,  # Point size
    alpha = 0.5  # Transparency
  ) +
  scale_color_gradient(
    low = "darkblue", 
    high = "white",     
    name = "Distance from water"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )

print(water_heatmap)


####
#open

# Add point heatmap layer
open_heatmap <- base_map +
  geom_point(
    data = covs_numeric, 
    aes(x = longitude, y = latitude, color = propopen500m), 
    size = 5,  # Point size
    alpha = 0.5  # Transparency
  ) +
  scale_color_viridis_c(
    option = "inferno",  # or "viridis", "inferno", "magma", "cividis"
    name = "Proportion of openness",
    direction = -1       # Optional: reverses the scale
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )

print(open_heatmap)

####
#humstructures

# Add point heatmap layer
hum_heatmap <- base_map +
  geom_point(
    data = covs_numeric, 
    aes(x = longitude, y = latitude, color = humdist_short), 
    size = 5,  # Point size
    alpha = 0.5  # Transparency
  ) +
  scale_color_viridis_c(
    option = "inferno",  # or "viridis", "inferno", "magma", "cividis"
    name = "Distance to humans",
    direction = -1       # Optional: reverses the scale
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )

print(hum_heatmap)


####
#Grazing - cattle

# Add point heatmap layer
cattle_heatmap <- base_map +
  geom_point(
    data = covs_numeric, 
    aes(x = longitude, y = latitude, color = cattle_30min_event_rate), 
    size = 5,  # Point size
    alpha = 0.5  # Transparency
  ) +
  scale_color_gradient(
    low = "white", 
    high = "purple",     
    name = "Cattle Grazing"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )

print(cattle_heatmap)

####
#Grazing - shoat

# Add point heatmap layer
shoat_heatmap <- base_map +
  geom_point(
    data = covs_numeric, 
    aes(x = longitude, y = latitude, color = shoat_30min_event_rate), 
    size = 5,  # Point size
    alpha = 0.5  # Transparency
  ) +
  scale_color_gradient(
    low = "white", 
    high = "forestgreen",     
    name = "Shoat Grazing"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )

print(shoat_heatmap)
