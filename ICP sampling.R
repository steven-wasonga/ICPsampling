# Geographically regulated sampling (ICP design)
# Based on: https://gh.bmj.com/content/7/10/e009500 

setwd("C:/Users/HP/Documents/GNSAMPLING")

#load("ICPsampling.RData")

# Load libraries
library(geosphere)
library(ggplot2)
library(sp)
library(sf)
library(dplyr)
library(parallel)
library(data.table)

# For reproducibility
set.seed(9365)


# 1. LOAD AND PREPARE DATA


# Load village data
kenya_villagesraw <- st_read("kenya_villages.shp")
kenya_villagesraw <- kenya_villagesraw |>
  st_drop_geometry() |>
  as.data.frame()

# rename columns to lowercase
colnames(kenya_villagesraw) <- tolower(colnames(kenya_villagesraw))
kenya_villagesraw <- kenya_villagesraw |>
  dplyr::select(name, lon = longdd, lat = latdd) |>
  dplyr::mutate(id = row_number(), .before = 1)

head(kenya_villagesraw)
cat("Total villages loaded:", nrow(kenya_villagesraw), "\n")

# Data validation and cleaning
#remove duplicated villages
kenya_villages_clean <- kenya_villagesraw[
  !duplicated(kenya_villagesraw[, c("name", "lon", "lat")]),
]

coordinates <- data.frame(
  lon = kenya_villages_clean$lon,
  lat = kenya_villages_clean$lat
)
# Check coordinates within Kenya boundaries (approx: lon 33.9-41.9, lat -4.7-5.5)
cat("Lon range:", range(coordinates$lon), "\n")
cat("Lat range:", range(coordinates$lat), "\n")


# Remove NAs
coordinates <- coordinates[complete.cases(coordinates), ]
# Remove duplicated coordinates
duplicates <- duplicated(coordinates[, c("lon", "lat")])
coordinates <- coordinates[!duplicates, ]
cat("Villages after removing NAs and duplicates:", nrow(coordinates), "\n")

# Create spatial object for villages
villages_sf <- st_as_sf(coordinates, coords = c("lon", "lat"), crs = 4326)


# 2. LOAD AND PROCESS WARDS/SHAPEFILES


# Load wards shapefile
wards <- st_read("wards.shp")
wards <- st_make_valid(wards)

# Clean county names
wards$county <- tolower(wards$county)
wards$county <- gsub("[-_]", " ", wards$county)
wards$county <- trimws(wards$county)

# Assign each village to a ward/county
villages_with_ward <- st_join(villages_sf, wards, join = st_within)

# Handle villages not assigned to any ward
missing <- is.na(villages_with_ward$pop2009)
if (any(missing)) {
  missing_villages <- villages_sf[missing, ]
  nearest_wards <- st_nearest_feature(missing_villages, wards)
  villages_with_ward$county[missing] <- wards$county[nearest_wards]
  villages_with_ward$pop2009[missing] <- wards$pop2009[nearest_wards]
}
cat("Villages successfully assigned to wards:", sum(!is.na(villages_with_ward$pop2009)), "\n")


# 3. FILTER OUT COUNTIES COVERED IN THE CLASSICAL SAMPLING TECHNIQUE +NAIROBI

exclude_counties <- c("turkana", "kilifi", "garissa", "tharaka nithi", 
                      "baringo", "kitui", "siaya", "nairobi")

# Filter out excluded counties
villages_with_ward <- villages_with_ward[!villages_with_ward$county %in% exclude_counties, ]
cat("Villages after excluding the 7 counties:", nrow(villages_with_ward), "\n")

# Load counties shapefile for visualization
counties <- st_read("County.shp")
counties <- st_make_valid(counties)


# 4.ICP SAMPLING starts here
# Extract coordinates for sampling
villages_with_ward <- villages_with_ward %>%
  mutate(indices = row_number())
all_indices <- villages_with_ward$indices  #unique identifications


# ------------------------------------------------------------
# Function 1: Select primary points with minimum distance (δ)
# ------------------------------------------------------------
select_primary_points <- function(coords, n_primary, delta_m, max_attempts = 1000) {
  # n_primary: Number of primary points to select
  # delta_m: Minimum distance between primary points (in meters)
  # max_attempts: Maximum random attempts before using fallback strategy
  
  # Get total number of available villages/points
  coords <- st_coordinates(villages_with_ward)
  n_total <- nrow(coords)
  
  # Create vector to store indices of selected primary points
  # These indices refer to unique IDs in the 'villages_with_ward' data
  selected_indices <- integer(n_primary)
  
  # Create matrix to store coordinates of selected points
  selected_coords <- matrix(NA, nrow = n_primary, ncol = 2)
  colnames(selected_coords) <- c("X", "Y")
  
  # Define all possible indices (1 through total number of villages)
  all_indices <- 1:n_total
  
  cat("Selecting", n_primary, "primary villages with minimum spacing", 
      delta_m/1000, "km...\n")
  
  #progress bar
  pb <- txtProgressBar(min = 0, max = n_primary, style = 3)
  
  # MAIN LOOP: Select each primary point one by one
  for (i in 1:n_primary) {
    attempts <- 0          # Counter for random attempts
    point_found <- FALSE   # Flag to track if suitable point found
    
    # Keep trying to find a point until max attempts reached or point found
    while (attempts < max_attempts && !point_found) {
      
      # STEP 1: Randomly select a candidate point from remaining villages
      # setdiff() removes already selected indices from all available indices
      candidate_idx <- sample(setdiff(all_indices, selected_indices[1:(i-1)]), 1)
      candidate_coord <- coords[candidate_idx, ]  # Get coordinates of candidate
      
      # STEP 2: Check if candidate meets distance requirements
      if (i == 1) {
        # First point has no distance constraints - accept immediately
        selected_indices[i] <- candidate_idx
        selected_coords[i, ] <- candidate_coord
        point_found <- TRUE
      } else {
        # For points 2 and onward, check distance to all previously selected points
        
        # Calculate distances using Haversine formula (accounts for Earth's curvature)
        distances <- distHaversine(
          candidate_coord,                           # New candidate point
          selected_coords[1:(i-1), , drop = FALSE]   # All previously selected points
        )
        
        # Check if minimum distance to any existing point is >= delta_m
        if (min(distances) >= delta_m) {
          # Candidate is sufficiently far from all existing points - accept it
          selected_indices[i] <- candidate_idx
          selected_coords[i, ] <- candidate_coord
          point_found <- TRUE
        }
      }
      attempts <- attempts + 1  # Increment attempt counter
    }
    
    # STEP 3: Fallback if no suitable point found after max attempts
    if (!point_found) {
      # Get all villages that haven't been selected yet
      available <- setdiff(all_indices, selected_indices[1:(i-1)])
      
      if (length(available) > 0) {
        # For each available village, calculate its minimum distance to existing points
        min_dists <- sapply(available, function(idx) {
          if (i == 1) return(Inf)  # If first point, return infinite distance
          
          # Calculate distances from this village to all selected points
          dists <- distHaversine(
            coords[idx, ],                            # Current available village
            selected_coords[1:(i-1), , drop = FALSE]  # All selected points
          )
          min(dists)  # Return minimum distance to any selected point
        })
        
        # Select the village with the MAXIMUM of these minimum distances
        # This picks the point that's "most isolated" from existing points
        best_idx <- available[which.max(min_dists)]
        selected_indices[i] <- best_idx
        selected_coords[i, ] <- coords[best_idx, ]
        cat("\nWarning: Used fallback for point", i, "\n")
      } else {
        stop("No available villages left - cannot select more primary points")
      }
    }
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close progress bar
  close(pb)
  
  # Return results as a list
  return(list(
    indices = selected_indices,   # Indices of selected primary points
    coords = selected_coords,     # Coordinates of selected points
    n_selected = n_primary        # Number of points selected (should equal n_primary)
  ))
}

# ------------------------------------------------------------
# Function 2: Add close pairs within radius ζ
# ------------------------------------------------------------
add_close_pairs <- function(primary_indices, primary_coords, all_coords, 
                            all_indices, zeta_m) {
  # primary_indices: Indices of primary points (from select_primary_points)
  # primary_coords: Coordinates of primary points
  # all_coords: Coordinates of ALL villages/points in dataset
  # all_indices: Indices of ALL villages (usually 1:nrow(all_coords))
  # zeta_m: Maximum radius for close pairs (in meters)
  
  n_primary <- length(primary_indices)  # Number of primary points
  
  # Create storage vectors:
  close_pair_indices <- integer(n_primary)  # Indices of close pair villages
  assignments <- integer(n_primary)         # Which primary each pair is assigned to
  pair_distances <- numeric(n_primary)      # Distance between primary and its pair
  
  cat("Selecting", n_primary, "close pairs within", zeta_m/1000, "km radius...\n")
  pb <- txtProgressBar(min = 0, max = n_primary, style = 3)  # Progress bar
  
  # MAIN LOOP: Find close pair for each primary point
  for (i in 1:n_primary) {
    # Get current primary point's index and coordinates
    primary_idx <- primary_indices[i]
    primary_coord <- primary_coords[i, ]
    
    # STEP 1: Create list of excluded villages (cannot be chosen as close pair)
    # Exclude: 1) All primary villages, 2) Close pairs already selected for previous primaries
    excluded <- c(primary_indices, close_pair_indices[1:(i-1)])
    available <- setdiff(all_indices, excluded)  # Remaining eligible villages
    
    # Check if any villages are available
    if (length(available) == 0) {
      warning("No available villages for close pair ", i)
      next  # Skip to next iteration
    }
    
    # STEP 2: Calculate distances from current primary to all available villages
    distances <- distHaversine(
      primary_coord,                       # Current primary point
      all_coords[available, , drop = FALSE]  # All available villages
    )
    
    # Identify villages within the specified radius (zeta_m)
    within_radius <- which(distances <= zeta_m)
    
    # STEP 3: Select a close pair
    if (length(within_radius) > 0) {
      # CASE A: One or more villages within radius - randomly select one
      selected <- if (length(within_radius) == 1) {
        within_radius  # Only one choice
      } else {
        sample(within_radius, 1)  # Randomly pick one from within radius
      }
      
      # Record the selected village and its distance
      pair_idx <- available[selected]
      pair_dist <- distances[selected]
    } else {
      # CASE B: No villages within radius - use closest available village
      closest <- which.min(distances)  # Index of closest village
      pair_idx <- available[closest]
      pair_dist <- distances[closest]
      
      # Warning message with distance in km
      cat("\nNo village within radius for primary", i, 
          "- using closest at", round(pair_dist/1000, 1), "km\n")
    }
    
    # STEP 4: Store results
    close_pair_indices[i] <- pair_idx      # Index of selected close pair
    assignments[i] <- primary_idx          # Which primary this pair belongs to
    pair_distances[i] <- pair_dist         # Distance between primary and pair
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  close(pb)  # Close progress bar
  
  # Remove any zeros (if loop broke early due to no available villages)
  valid <- close_pair_indices > 0
  
  # Return results as a list
  return(list(
    indices = close_pair_indices[valid],      # Indices of close pair villages
    assignments = assignments[valid],         # Corresponding primary indices
    distances = pair_distances[valid],        # Distances between pairs
    n_pairs = sum(valid)                      # Number of valid pairs found
  ))
}



# 5. RUN ICP SAMPLING
# Step 1: set ICP parameters
n_primary <- 150     # Number of primary sample points
k <- 150             # Number of close pairs 
delta_km <- 40       # Minimum distance between primary points (km)
delta_m <- delta_km * 1000 #since distances are measured in m
zeta_km <- 15        # Radius for close pair selection (km)
zeta_m <- zeta_km * 1000
n_total <- n_primary + k

# Step 2: Select primary points
primary_result <- select_primary_points(
  coords <- st_coordinates(villages_with_ward),
  n_primary = n_primary,
  delta_m = delta_m
)

cat("Selected", primary_result$n_selected, "primary villages\n")

# Step 3: Add close pairs
pair_result <- add_close_pairs(
  primary_indices = primary_result$indices,
  primary_coords = primary_result$coords,
  all_coords = st_coordinates(villages_with_ward),
  all_indices = villages_with_ward$indices,
  zeta_m = zeta_m
)

cat("Selected", length(pair_result$indices), "close pairs\n")



# Create results data frame
primary_points <- villages_with_ward[primary_result$indices, ]
primary_points$type <- "Primary"
primary_points$pair_id <- 1:nrow(primary_points)

close_pairs <- villages_with_ward[pair_result$indices, ]
close_pairs$type <- "Close Pair"
close_pairs$pair_id <- 1:nrow(close_pairs)

# Combine all sampled points
all_sampled <- rbind(primary_points, close_pairs)
head(all_sampled)

# 6. VISUALIZATION
# Load counties shapefile 
counties <- st_read("County.shp")
counties <- st_make_valid(counties)

#plot counties
p<-ggplot() +
  geom_sf(data = counties, 
          fill = "white", 
          color = "gray40", 
          linewidth = 0.3) +
  
  # Plot all sampled points
  geom_sf(data = all_sampled, 
          aes(color = type),
          size = 1) +         
  scale_color_manual(values = c("Primary" = "blue", 
                                "Close Pair" = "orange")) +
  theme_minimal() +
  labs(title = "Sampled Villages",
       color = "Type")  # Legend title

#view plot
p


# 6. SAVE RESULTS
#plot
ggsave("sampling_map.pdf", plot = p, width = 10, height = 8)

#villages per county
county_table <- as.data.frame(table(all_sampled$county))
colnames(county_table) <- c("County", "Number_of_Villages_Sampled")
head(county_table)
write.csv(county_table, "county_village_counts.csv", row.names = FALSE)

#add village names from the original dataset
kenya_villages_sf <- st_as_sf(
  kenya_villagesraw,
  coords = c("lon", "lat"),
  crs = 4326
)
all_sampled <- st_transform(all_sampled, 4326)

all_sampled_named <- st_join(
  all_sampled,
  kenya_villages_sf["name"],
  join = st_contains
)
all_sampled_named <- all_sampled_named[, c(
  "county", "subcounty", "ward", "name", "type", "geometry"
)]
head(all_sampled_named)
all_sampled_named_df <- st_drop_geometry(all_sampled_named)
write.csv(
  all_sampled_named_df,
  "ICP_sampled_villages.csv",
  row.names = FALSE
)

#######=================================================
#summary statistics
cat(
  "\n=== ICP SAMPLING SUMMARY ===\n",
  "Total villages in the raw data: ", nrow(kenya_villagesraw), "\n",
  "After cleaning: ", nrow(coordinates), "\n",
  "After excluding the 8 counties: ", nrow(villages_with_ward), "\n",
  "Primary points selected: ", n_primary, "\n",
  "Close pairs selected: ", length(pair_result$indices), "\n",
  "Total number of villages selected: ", n_primary + length(pair_result$indices), "\n",
  "Minimum distance between primaries (δ): ", delta_km, " km\n",
  "Close pair radius (ζ): ", zeta_km, " km\n",
  "\n=== SAMPLING COMPLETE ===\n",
  sep = ""
)


#=============
#save.image("ICPsampling.RData")





