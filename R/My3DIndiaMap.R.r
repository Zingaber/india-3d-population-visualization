# ============================================================================
# ULTIMATE INDIA 3D POPULATION VISUALIZATION
# Combining Local XYZ Data + HDX Integration + Advanced 3D Techniques
# ============================================================================

# INSTALLATION GUIDE - Run once
# install.packages(c("tidyverse", "sf", "stars", "httr", "jsonlite", "R.utils", 
#                   "rayshader", "rayrender", "raster", "elevatr", "giscoR", 
#                   "terra", "magick", "rgl", "pracma", "scales"))

# ============================================================================
# LOAD REQUIRED PACKAGES
# ============================================================================

library(tidyverse)      # Data manipulation
library(sf)             # Spatial data handling  
library(stars)          # Raster processing
library(raster)         # Raster operations
library(terra)          # Modern raster handling
library(rayshader)      # 3D visualization
library(rayrender)      # Advanced rendering
library(elevatr)        # Elevation data
library(giscoR)         # Country boundaries
library(magick)         # Image processing
library(rgl)            # 3D graphics
library(pracma)         # Mathematical functions
library(scales)         # Scale formatting

# ============================================================================
# CONFIGURATION AND SETUP
# ============================================================================

# Directory setup
DATA_DIR <- "C:/Users/josze/MyRWorkspace/Download"
OUTPUT_DIR <- "C:/Users/josze/MyRWorkspace/Download/ultimate_india_maps"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Your existing data file path
LOCAL_DATA_PATH <- file.path(DATA_DIR, "ind_pd_2020_1km_ASCII_XYZ/ind_pd_2020_1km_ASCII_XYZ.csv")

# High-quality rendering settings
MAP_RESOLUTION <- 6000    # Very high resolution
EXPORT_WIDTH <- 4000      # 4K width
EXPORT_HEIGHT <- 3000     # 3K height
WINDOW_SIZE <- 1600       # Large window

# India-specific projection (optimized for the subcontinent)
INDIA_CRS <- "+proj=laea +lat_0=20 +lon_0=78 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# ============================================================================
# ENHANCED DATA LOADING AND PROCESSING
# ============================================================================

load_and_process_india_data <- function() {
  cat("=== LOADING AND PROCESSING INDIA POPULATION DATA ===\n")
  
  # Load your existing XYZ data
  cat("Loading local XYZ data...\n")
  pop_xyz <- read_csv(
    LOCAL_DATA_PATH,
    col_types = cols(x = col_double(), y = col_double(), z = col_double()),
    skip = 1
  ) %>%
    drop_na() %>%
    filter(z > 0)  # Remove zero population areas
  
  cat("Loaded", scales::comma(nrow(pop_xyz)), "population points\n")
  cat("Population range:", scales::comma(min(pop_xyz$z)), "to", scales::comma(max(pop_xyz$z)), "\n")
  
  # Get India boundaries for proper masking
  cat("Fetching India boundaries...\n")
  india_borders <- giscoR::gisco_get_countries(country = "IN", resolution = "3")
  india_bbox <- st_bbox(india_borders)
  
  cat("India bounding box:\n")
  print(india_bbox)
  
  # Filter data to India extent with buffer
  pop_filtered <- pop_xyz %>%
    filter(
      x >= india_bbox["xmin"] - 0.5,
      x <= india_bbox["xmax"] + 0.5,
      y >= india_bbox["ymin"] - 0.5,
      y <= india_bbox["ymax"] + 0.5
    )
  
  cat("Filtered to", scales::comma(nrow(pop_filtered)), "points within India bounds\n")
  
  # Create high-resolution raster
  cat("Creating high-resolution raster...\n")
  pop_raster <- rasterFromXYZ(pop_filtered, crs = "+proj=longlat +datum=WGS84")
  
  # Crop and mask to India boundaries
  pop_raster_cropped <- crop(pop_raster, india_borders)
  pop_raster_masked <- mask(pop_raster_cropped, india_borders)
  
  # Get elevation data for terrain context
  cat("Adding elevation context...\n")
  elevation_data <- get_elev_raster(india_borders, z = 7)  # High resolution
  elevation_cropped <- crop(elevation_data, extent(pop_raster_masked))
  elevation_resampled <- resample(elevation_cropped, pop_raster_masked)
  
  return(list(
    population_raster = pop_raster_masked,
    elevation_raster = elevation_resampled,
    boundaries = india_borders,
    raw_data = pop_filtered
  ))
}

# ============================================================================
# ADVANCED DATA PROCESSING FOR SPIKE REDUCTION
# ============================================================================

process_for_3d_visualization <- function(pop_raster, elev_raster, method = "hybrid") {
  cat("=== PROCESSING DATA FOR 3D VISUALIZATION ===\n")
  
  # Convert to matrices
  pop_matrix <- raster_to_matrix(pop_raster)
  elev_matrix <- raster_to_matrix(elev_raster)
  
  # Handle missing values
  pop_matrix[is.na(pop_matrix)] <- 0
  elev_matrix[is.na(elev_matrix)] <- 0
  
  cat("Matrix dimensions:", dim(pop_matrix)[1], "x", dim(pop_matrix)[2], "\n")
  cat("Population range:", scales::comma(min(pop_matrix)), "to", scales::comma(max(pop_matrix)), "\n")
  
  # Advanced smoothing to reduce spikes
  cat("Applying advanced smoothing...\n")
  
  if (method == "gaussian") {
    # Gaussian smoothing
    pop_smooth <- apply_gaussian_smooth(pop_matrix, sigma = 2)
    
  } else if (method == "adaptive") {
    # Adaptive smoothing based on local density
    pop_smooth <- apply_adaptive_smooth(pop_matrix)
    
  } else if (method == "hybrid") {
    # Hybrid approach: different processing for different density levels
    pop_smooth <- apply_hybrid_processing(pop_matrix)
    
  } else {
    # Simple moving average
    pop_smooth <- apply_simple_smooth(pop_matrix, window = 5)
  }
  
  # Advanced scaling options
  cat("Applying intelligent scaling...\n")
  pop_scaled <- apply_intelligent_scaling(pop_smooth, method = "cube_root_log")
  
  # Normalize elevation for subtle base
  elev_normalized <- (elev_matrix - min(elev_matrix)) / (max(elev_matrix) - min(elev_matrix))
  elev_base <- elev_normalized * 3  # Subtle elevation base
  
  # Combine population spikes with terrain base
  combined_matrix <- elev_base + pop_scaled
  
  cat("Final matrix range:", round(min(combined_matrix), 2), "to", round(max(combined_matrix), 2), "\n")
  
  return(list(
    population_matrix = pop_matrix,
    elevation_matrix = elev_matrix,
    processed_matrix = pop_scaled,
    combined_matrix = combined_matrix,
    smooth_matrix = pop_smooth
  ))
}

# ============================================================================
# SMOOTHING ALGORITHMS
# ============================================================================

apply_gaussian_smooth <- function(matrix, sigma = 1.5) {
  # Create Gaussian kernel
  kernel_size <- ceiling(sigma * 6)
  if (kernel_size %% 2 == 0) kernel_size <- kernel_size + 1
  
  x <- seq(-(kernel_size-1)/2, (kernel_size-1)/2, 1)
  kernel <- outer(x, x, function(x, y) exp(-(x^2 + y^2) / (2 * sigma^2)))
  kernel <- kernel / sum(kernel)
  
  # Apply convolution
  if (require(pracma, quietly = TRUE)) {
    return(conv2(matrix, kernel, "same"))
  } else {
    return(apply_simple_smooth(matrix, window = kernel_size))
  }
}

apply_adaptive_smooth <- function(matrix, base_window = 3, max_window = 9) {
  result <- matrix
  rows <- nrow(matrix)
  cols <- ncol(matrix)
  
  for (i in 2:(rows-1)) {
    for (j in 2:(cols-1)) {
      # Calculate local variance to determine smoothing level
      local_values <- as.vector(matrix[max(1, i-1):min(rows, i+1), 
                                      max(1, j-1):min(cols, j+1)])
      local_var <- var(local_values, na.rm = TRUE)
      
      # Adaptive window size based on local variance
      if (is.na(local_var) || local_var < quantile(matrix, 0.7, na.rm = TRUE)) {
        window <- base_window
      } else {
        window <- max_window
      }
      
      # Apply smoothing
      half_window <- floor(window / 2)
      i_range <- max(1, i - half_window):min(rows, i + half_window)
      j_range <- max(1, j - half_window):min(cols, j + half_window)
      
      result[i, j] <- mean(matrix[i_range, j_range], na.rm = TRUE)
    }
  }
  
  return(result)
}

apply_hybrid_processing <- function(matrix) {
  # Identify high, medium, and low density areas
  q75 <- quantile(matrix, 0.75, na.rm = TRUE)
  q95 <- quantile(matrix, 0.95, na.rm = TRUE)
  
  result <- matrix
  
  # Heavy smoothing for very high density areas (reduce extreme spikes)
  high_density_mask <- matrix > q95
  if (sum(high_density_mask, na.rm = TRUE) > 0) {
    high_smooth <- apply_gaussian_smooth(matrix, sigma = 3)
    result[high_density_mask] <- high_smooth[high_density_mask]
  }
  
  # Moderate smoothing for medium density areas
  med_density_mask <- matrix > q75 & matrix <= q95
  if (sum(med_density_mask, na.rm = TRUE) > 0) {
    med_smooth <- apply_gaussian_smooth(matrix, sigma = 1.5)
    result[med_density_mask] <- med_smooth[med_density_mask]
  }
  
  # Light smoothing for low density areas (preserve detail)
  low_density_mask <- matrix <= q75
  if (sum(low_density_mask, na.rm = TRUE) > 0) {
    low_smooth <- apply_simple_smooth(matrix, window = 3)
    result[low_density_mask] <- low_smooth[low_density_mask]
  }
  
  return(result)
}

apply_simple_smooth <- function(matrix, window = 3) {
  result <- matrix
  half_window <- floor(window / 2)
  rows <- nrow(matrix)
  cols <- ncol(matrix)
  
  for (i in (half_window + 1):(rows - half_window)) {
    for (j in (half_window + 1):(cols - half_window)) {
      i_range <- (i - half_window):(i + half_window)
      j_range <- (j - half_window):(j + half_window)
      result[i, j] <- mean(matrix[i_range, j_range], na.rm = TRUE)
    }
  }
  
  return(result)
}

# ============================================================================
# INTELLIGENT SCALING METHODS
# ============================================================================

apply_intelligent_scaling <- function(matrix, method = "cube_root_log") {
  if (method == "cube_root") {
    # Cube root scaling (gentler than square root)
    scaled <- matrix^(1/3)
    
  } else if (method == "fourth_root") {
    # Even gentler fourth root scaling
    scaled <- matrix^(1/4)
    
  } else if (method == "log_plus") {
    # Log scaling with offset
    scaled <- log10(matrix + 1)
    
  } else if (method == "cube_root_log") {
    # Hybrid: cube root for low values, log for high values
    threshold <- quantile(matrix, 0.8, na.rm = TRUE)
    scaled <- ifelse(matrix <= threshold, 
                    matrix^(1/3), 
                    log10(matrix + 1) * (threshold^(1/3) / log10(threshold + 1)))
    
  } else if (method == "quantile_based") {
    # Quantile-based scaling to reduce outlier impact
    q99 <- quantile(matrix, 0.99, na.rm = TRUE)
    scaled <- pmin(matrix, q99)  # Cap at 99th percentile
    scaled <- scaled^(1/3)
    
  } else {
    scaled <- matrix
  }
  
  # Normalize to 0-1
  scaled <- (scaled - min(scaled, na.rm = TRUE)) / (max(scaled, na.rm = TRUE) - min(scaled, na.rm = TRUE))
  
  # Scale for 3D height (much more conservative)
  scaled <- scaled * 15  # Reduced from original 100+ to 15
  
  return(scaled)
}

# ============================================================================
# ENHANCED COLOR PALETTES
# ============================================================================

create_ultimate_color_palettes <- function() {
  palettes <- list(
    # India flag inspired with gradient
    tricolor_gradient = colorRampPalette(c(
      "#000080", "#4169E1", "#87CEEB", "#F0F8FF", 
      "#FFE4B5", "#FFA500", "#FF6347", "#FF4500"
    ))(256),
    
    # Population density spectrum
    population_spectrum = colorRampPalette(c(
      "#F7FCF0", "#E0F3DB", "#CCEBC5", "#A8DDB5", 
      "#7BCCC4", "#4EB3D3", "#2B8CBE", "#0868AC", "#084081"
    ))(256),
    
    # Terrain + Population
    terrain_pop = colorRampPalette(c(
      "#2E7D32", "#4CAF50", "#8BC34A", "#CDDC39",
      "#FFEB3B", "#FFC107", "#FF9800", "#FF5722", "#D32F2F"
    ))(256),
    
    # Sunset India
    sunset_india = colorRampPalette(c(
      "#1A237E", "#3F51B5", "#2196F3", "#03DAC6",
      "#4CAF50", "#8BC34A", "#CDDC39", "#FFC107", 
      "#FF9800", "#FF5722", "#E91E63"
    ))(256),
    
    # Monochrome elegant
    mono_elegant = colorRampPalette(c(
      "#F8F9FA", "#E9ECEF", "#DEE2E6", "#CED4DA",
      "#ADB5BD", "#6C757D", "#495057", "#343A40", "#212529"
    ))(256),
    
    # Vibrant contrast
    vibrant_contrast = colorRampPalette(c(
      "#E8F5E8", "#A5D6A7", "#66BB6A", "#4CAF50",
      "#388E3C", "#2E7D32", "#1B5E20", "#0D4619"
    ))(256)
  )
  
  return(palettes)
}

# ============================================================================
# ULTIMATE 3D RENDERING FUNCTION
# ============================================================================

create_ultimate_3d_map <- function(data_list, palette_name = "population_spectrum", 
                                  quality = "high", style = "realistic") {
  
  cat("=== CREATING ULTIMATE 3D VISUALIZATION ===\n")
  cat("Palette:", palette_name, "| Quality:", quality, "| Style:", style, "\n")
  
  # Get processed matrices
  combined_matrix <- data_list$combined_matrix
  pop_matrix <- data_list$population_matrix
  
  # Get color palette
  palettes <- create_ultimate_color_palettes()
  colors <- palettes[[palette_name]]
  
  # Create base hillshade
  cat("Generating base hillshade...\n")
  hillshade <- height_shade(combined_matrix, texture = colors)
  
  # Add multiple shadow layers for realism
  if (style == "realistic") {
    cat("Adding realistic lighting...\n")
    
    # Ambient shade
    ambient <- ambient_shade(combined_matrix, zscale = 0.03, maxsearch = 200)
    hillshade <- add_shadow(hillshade, ambient, 0.4)
    
    # Ray-traced shadows
    ray_shadow <- ray_shade(combined_matrix, zscale = 0.03, lambert = TRUE)
    hillshade <- add_shadow(hillshade, ray_shadow, 0.8)
    
    # Lamb shading for soft lighting
    lamb_shade <- lamb_shade(combined_matrix, zscale = 0.03)
    hillshade <- add_shadow(hillshade, lamb_shade, 0.3)
    
  } else if (style == "artistic") {
    # Artistic style with enhanced contrast
    ambient <- ambient_shade(combined_matrix, zscale = 0.05)
    hillshade <- add_shadow(hillshade, ambient, 0.6)
  }
  
  # Close existing windows
  try(rgl::close3d(), silent = TRUE)
  
  # Render settings based on quality
  if (quality == "ultra") {
    window_size <- 2000
    samples <- 500
    zscale <- 0.015
  } else if (quality == "high") {
    window_size <- 1600
    samples <- 300
    zscale <- 0.02
  } else {
    window_size <- 1200
    samples <- 150
    zscale <- 0.025
  }
  
  # Create 3D plot
  cat("Rendering 3D plot...\n")
  plot_3d(
    hillshade,
    heightmap = combined_matrix,
    zscale = zscale,
    solid = TRUE,
    shadow = TRUE,
    shadowdepth = -50,
    windowsize = c(window_size, window_size),
    phi = 30,                    # Optimal viewing angle for India
    theta = 45,                  # Good rotation for coastlines
    zoom = 0.65,                 # Balanced zoom
    background = "#87CEEB",      # Sky blue
    shadowcolor = "#2F4F4F"      # Dark slate gray
  )
  
  # Create output filename
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  base_filename <- sprintf("india_ultimate_%s_%s_%s", palette_name, quality, timestamp)
  
  # Multiple viewing angles
  views <- list(
    perspective = list(phi = 30, theta = 45, zoom = 0.65, suffix = "perspective"),
    aerial = list(phi = 85, theta = 0, zoom = 0.8, suffix = "aerial"),
    coastline = list(phi = 25, theta = 15, zoom = 0.7, suffix = "coastline"),
    diagonal = list(phi = 40, theta = 60, zoom = 0.6, suffix = "diagonal")
  )
  
  rendered_files <- c()
  
  for (view_name in names(views)) {
    view <- views[[view_name]]
    
    cat("Creating", view_name, "view...\n")
    
    # Set camera angle
    render_camera(phi = view$phi, theta = view$theta, zoom = view$zoom)
    
    # Quick snapshot
    snapshot_file <- file.path(OUTPUT_DIR, paste0(base_filename, "_", view$suffix, ".png"))
    render_snapshot(snapshot_file, width = 1600, height = 1200)
    rendered_files <- c(rendered_files, snapshot_file)
    
    # High-quality render for main view
    if (view_name == "perspective") {
      hq_file <- file.path(OUTPUT_DIR, paste0(base_filename, "_", view$suffix, "_HQ.png"))
      
      render_highquality(
        filename = hq_file,
        lightdirection = c(225, 315, 135),    # Multiple light sources
        lightaltitude = c(45, 80, 25),
        lightintensity = c(300, 500, 200),
        lightcolor = c("white", "#FFFFCC", "#E6E6FA"),
        width = EXPORT_WIDTH,
        height = EXPORT_HEIGHT,
        samples = samples,
        interactive = FALSE,
        backgroundhigh = "#87CEEB",
        backgroundlow = "#F0F8FF"
      )
      
      rendered_files <- c(rendered_files, hq_file)
    }
  }
  
  cat("Rendered files:\n")
  for (file in rendered_files) {
    cat("-", basename(file), "\n")
  }
  
  return(rendered_files)
}

# ============================================================================
# COMPARISON AND ANALYSIS FUNCTIONS
# ============================================================================

create_processing_comparison <- function(data_list) {
  cat("=== CREATING PROCESSING METHOD COMPARISON ===\n")
  
  pop_raster <- data_list$population_raster
  elev_raster <- data_list$elevation_raster
  
  methods <- c("gaussian", "adaptive", "hybrid")
  scaling_methods <- c("cube_root", "fourth_root", "cube_root_log")
  
  comparison_files <- c()
  
  for (method in methods) {
    for (scaling in scaling_methods) {
      cat("Processing with", method, "smoothing and", scaling, "scaling...\n")
      
      processed <- process_for_3d_visualization(pop_raster, elev_raster, method = method)
      processed$combined_matrix <- apply_intelligent_scaling(processed$smooth_matrix, method = scaling) + 
                                  (processed$elevation_matrix - min(processed$elevation_matrix)) / 
                                  (max(processed$elevation_matrix) - min(processed$elevation_matrix)) * 3
      
      # Quick visualization
      colors <- create_ultimate_color_palettes()$population_spectrum
      hillshade <- height_shade(processed$combined_matrix, texture = colors)
      
      try(rgl::close3d(), silent = TRUE)
      plot_3d(hillshade, heightmap = processed$combined_matrix, zscale = 0.02, 
              windowsize = c(800, 800), phi = 30, theta = 45, zoom = 0.65)
      
      filename <- file.path(OUTPUT_DIR, sprintf("comparison_%s_%s.png", method, scaling))
      render_snapshot(filename, width = 800, height = 600)
      comparison_files <- c(comparison_files, filename)
    }
  }
  
  cat("Created", length(comparison_files), "comparison images\n")
  return(comparison_files)
}

generate_statistics_report <- function(data_list) {
  cat("=== GENERATING COMPREHENSIVE STATISTICS REPORT ===\n")
  
  pop_matrix <- data_list$population_matrix
  raw_data <- data_list$raw_data
  
  stats <- list(
    total_features = nrow(raw_data),
    total_population = sum(raw_data$z, na.rm = TRUE),
    max_population = max(raw_data$z, na.rm = TRUE),
    mean_population = mean(raw_data$z, na.rm = TRUE),
    median_population = median(raw_data$z, na.rm = TRUE),
    std_population = sd(raw_data$z, na.rm = TRUE),
    matrix_dimensions = paste(dim(pop_matrix)[1], "x", dim(pop_matrix)[2]),
    geographic_extent = list(
      min_lon = min(raw_data$x),
      max_lon = max(raw_data$x),
      min_lat = min(raw_data$y),
      max_lat = max(raw_data$y)
    )
  )
  
  # Quantile analysis
  quantiles <- quantile(raw_data$z, c(0.5, 0.75, 0.9, 0.95, 0.99), na.rm = TRUE)
  
  # Create report dataframe
  report_df <- data.frame(
    Metric = c("Total Features", "Total Population", "Max Population", "Mean Population",
               "Median Population", "Std Deviation", "Matrix Size", "Min Longitude",
               "Max Longitude", "Min Latitude", "Max Latitude", "50th Percentile",
               "75th Percentile", "90th Percentile", "95th Percentile", "99th Percentile"),
    Value = c(
      scales::comma(stats$total_features),
      scales::comma(round(stats$total_population)),
      scales::comma(round(stats$max_population)),
      scales::comma(round(stats$mean_population, 2)),
      scales::comma(round(stats$median_population, 2)),
      scales::comma(round(stats$std_population, 2)),
      stats$matrix_dimensions,
      round(stats$geographic_extent$min_lon, 3),
      round(stats$geographic_extent$max_lon, 3),
      round(stats$geographic_extent$min_lat, 3),
      round(stats$geographic_extent$max_lat, 3),
      scales::comma(round(quantiles[1])),
      scales::comma(round(quantiles[2])),
      scales::comma(round(quantiles[3])),
      scales::comma(round(quantiles[4])),
      scales::comma(round(quantiles[5]))
    )
  )
  
  # Save report
  report_file <- file.path(OUTPUT_DIR, "india_population_statistics.csv")
  write_csv(report_df, report_file)
  
  cat("Statistics report saved to:", report_file, "\n")
  
  # Print summary
  cat("\nKEY STATISTICS:\n")
  cat("Total Population:", scales::comma(round(stats$total_population)), "\n")
  cat("Max Population Density:", scales::comma(round(stats$max_population)), "\n")
  cat("Geographic Coverage:", round(stats$geographic_extent$max_lon - stats$geographic_extent$min_lon, 1), 
      "° longitude x", round(stats$geographic_extent$max_lat - stats$geographic_extent$min_lat, 1), "° latitude\n")
  cat("Matrix Resolution:", stats$matrix_dimensions, "\n")
  
  return(report_df)
}

# ============================================================================
# MAIN EXECUTION PIPELINE
# ============================================================================

run_ultimate_india_pipeline <- function(create_comparisons = FALSE, quality = "high") {
  cat("\n")
  cat(rep("=", 70), "\n")
  cat("ULTIMATE INDIA 3D POPULATION VISUALIZATION PIPELINE\n")
  cat("Combining Local Data + Advanced Processing + Stunning 3D Rendering\n")
  cat(rep("=", 70), "\n\n")
  
  # Step 1: Load and process data
  cat("STEP 1: Loading and processing data...\n")
  data_list <- load_and_process_india_data()
  
  # Step 2: Advanced processing
  cat("\nSTEP 2: Advanced 3D processing...\n")
  processed_data <- process_for_3d_visualization(
    data_list$population_raster, 
    data_list$elevation_raster, 
    method = "hybrid"
  )
  
  # Combine processed data
  final_data <- c(data_list, processed_data)
  
  # Step 3: Generate statistics
  cat("\nSTEP 3: Generating statistics...\n")
  stats_report <- generate_statistics_report(final_data)
  
  # Step 4: Create visualizations
  cat("\nSTEP 4: Creating ultimate 3D visualizations...\n")
  
  palettes <- names(create_ultimate_color_palettes())
  all_rendered_files <- c()
  
  # Create visualizations with different palettes
  for (palette in palettes[1:3]) {  # Limit to top 3 palettes
    cat("\nCreating visualization with", palette, "palette...\n")
    rendered_files <- create_ultimate_3d_map(
      final_data, 
      palette_name = palette, 
      quality = quality,
      style = "realistic"
    )
    all_rendered_files <- c(all_rendered_files, rendered_files)
  }
  
  # Step 5: Create comparisons (optional)
  if (create_comparisons) {
    cat("\nSTEP 5: Creating processing comparisons...\n")
    comparison_files <- create_processing_comparison(data_list)
    all_rendered_files <- c(all_rendered_files, comparison_files)
  }
  
  # Step 6: Final summary
  cat("\n")
  cat(rep("=", 70), "\n")
  cat("ULTIMATE INDIA 3D VISUALIZATION COMPLETED!\n")
  cat(rep("=", 70), "\n")
  cat("Output directory:", OUTPUT_DIR, "\n")
  cat("Total files created:", length(all_rendered_files), "\n")
  cat("\nFiles created:\n")
  for (file in all_rendered_files) {
    cat("✓", basename(file), "\n")
  }
  
  cat("\nRecommended viewing order:\n")
  cat("1. *_perspective_HQ.png - Main high-quality visualization\n")
  cat("2. *_aerial.png - Top-down view for boundaries\n")
  cat("3. *_coastline.png - Coastal geography emphasis\n")
  cat("4. india_population_statistics.csv - Data analysis\n")
  
  return(all_rendered_files)
}

# ============================================================================
# QUICK START FUNCTIONS
# ============================================================================

# Quick function for immediate results
quick_india_3d <- function() {
  cat("Creating quick India 3D visualization...\n")
  data_list <- load_and_process_india_data()
  processed_data <- process_for_3d_visualization(data_list$population_raster, data_list$elevation_raster)
  final_data <- c(data_list, processed_data)
  
  rendered_files <- create_ultimate_3d_map(final_data, quality = "high", style = "realistic")
  
  cat("Quick visualization complete! Check:", basename(rendered_files[1]), "\n")
  return(rendered_files)
}

# Function to test different scaling methods
test_scaling_methods <- function() {
  cat("Testing different scaling methods...\n")
  data_list <- load_and_process_india_data()
  
  scaling_methods <- c("cube_root", "fourth_root", "cube_root_log", "quantile_based")
  test_files <- c()
  
  for (method in scaling_methods) {
    cat("Testing", method, "scaling...\n")
    
    # Process with current scaling method
    pop_matrix <- raster_to_matrix(data_list$population_raster)
    pop_matrix[is.na(pop_matrix)] <- 0
    pop_smooth <- apply_hybrid_processing(pop_matrix)
    pop_scaled <- apply_intelligent_scaling(pop_smooth, method = method)
    
    # Quick visualization
    colors <- create_ultimate_color_palettes()$population_spectrum
    hillshade <- height_shade(pop_scaled, texture = colors)
    
    try(rgl::close3d(), silent = TRUE)
    plot_3d(hillshade, heightmap = pop_scaled, zscale = 0.03, 
            windowsize = c(800, 800), phi = 30, theta = 45, zoom = 0.65)
    
    filename <- file.path(OUTPUT_DIR, paste0("scaling_test_", method, ".png"))
    render_snapshot(filename, width = 800, height = 600)
    test_files <- c(test_files, filename)
  }
  
  cat("Scaling tests complete! Check files:", paste(basename(test_files), collapse = ", "), "\n")
  return(test_files)
}

# ============================================================================
# INTERACTIVE FUNCTIONS
# ============================================================================

# Interactive palette selector
select_palette_interactive <- function() {
  palettes <- create_ultimate_color_palettes()
  cat("Available color palettes:\n")
  for (i in seq_along(palettes)) {
    cat(sprintf("%d. %s\n", i, names(palettes)[i]))
  }
  
  cat("\nEnter palette number (1-", length(palettes), "): ")
  choice <- as.numeric(readline())
  
  if (choice >= 1 && choice <= length(palettes)) {
    selected_palette <- names(palettes)[choice]
    cat("Selected palette:", selected_palette, "\n")
    return(selected_palette)
  } else {
    cat("Invalid choice. Using default palette.\n")
    return("population_spectrum")
  }
}

# Create custom palette
create_custom_palette <- function(colors) {
  cat("Creating custom palette with", length(colors), "colors\n")
  custom_palette <- colorRampPalette(colors)(256)
  return(custom_palette)
}

# ============================================================================
# BATCH PROCESSING AND AUTOMATION
# ============================================================================

batch_create_all_styles <- function() {
  cat("=== BATCH CREATING ALL VISUALIZATION STYLES ===\n")
  
  # Load data once
  data_list <- load_and_process_india_data()
  processed_data <- process_for_3d_visualization(data_list$population_raster, data_list$elevation_raster)
  final_data <- c(data_list, processed_data)
  
  # All palettes
  palettes <- names(create_ultimate_color_palettes())
  
  # All quality settings
  qualities <- c("high", "ultra")
  
  # All styles
  styles <- c("realistic", "artistic")
  
  batch_files <- c()
  total_combinations <- length(palettes) * length(qualities) * length(styles)
  current_combination <- 0
  
  cat("Creating", total_combinations, "visualization combinations...\n")
  
  for (palette in palettes) {
    for (quality in qualities) {
      for (style in styles) {
        current_combination <- current_combination + 1
        cat(sprintf("\n[%d/%d] Creating %s-%s-%s...\n", 
                   current_combination, total_combinations, palette, quality, style))
        
        tryCatch({
          files <- create_ultimate_3d_map(final_data, palette_name = palette, 
                                        quality = quality, style = style)
          batch_files <- c(batch_files, files)
        }, error = function(e) {
          cat("Error creating", palette, quality, style, ":", e$message, "\n")
        })
      }
    }
  }
  
  cat("\nBatch processing complete!", length(batch_files), "files created\n")
  return(batch_files)
}

# ============================================================================
# EXPORT AND SHARING FUNCTIONS
# ============================================================================

create_web_gallery <- function(image_files) {
  cat("=== CREATING WEB GALLERY ===\n")
  
  html_content <- paste0(
    "<!DOCTYPE html>\n",
    "<html>\n<head>\n",
    "<title>Ultimate India 3D Population Maps</title>\n",
    "<style>\n",
    "body { font-family: Arial, sans-serif; margin: 20px; background: #f5f5f5; }\n",
    ".gallery { display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 20px; }\n",
    ".image-card { background: white; padding: 15px; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }\n",
    ".image-card img { width: 100%; height: auto; border-radius: 4px; }\n",
    ".image-title { font-weight: bold; margin: 10px 0 5px 0; }\n",
    ".image-desc { color: #666; font-size: 0.9em; }\n",
    "h1 { text-align: center; color: #333; }\n",
    "</style>\n</head>\n<body>\n",
    "<h1>Ultimate India 3D Population Density Maps</h1>\n",
    "<p style='text-align: center; color: #666;'>Generated on ", Sys.Date(), "</p>\n",
    "<div class='gallery'>\n"
  )
  
  for (file in image_files) {
    filename <- basename(file)
    title <- gsub("india_ultimate_|_\\d{8}_\\d{6}", "", filename)
    title <- gsub("_", " ", title)
    title <- tools::toTitleCase(title)
    
    html_content <- paste0(html_content,
      "<div class='image-card'>\n",
      "<img src='", filename, "' alt='", title, "'>\n",
      "<div class='image-title'>", title, "</div>\n",
      "<div class='image-desc'>High-resolution 3D population visualization</div>\n",
      "</div>\n"
    )
  }
  
  html_content <- paste0(html_content, "</div>\n</body>\n</html>")
  
  gallery_file <- file.path(OUTPUT_DIR, "india_3d_gallery.html")
  writeLines(html_content, gallery_file)
  
  cat("Web gallery created:", gallery_file, "\n")
  return(gallery_file)
}

create_image_montage <- function(image_files, output_name = "india_montage.png") {
  cat("=== CREATING IMAGE MONTAGE ===\n")
  
  if (!require(magick, quietly = TRUE)) {
    cat("magick package required for montage creation\n")
    return(NULL)
  }
  
  # Read images
  images <- image_read(image_files[1:min(4, length(image_files))])  # Limit to 4 images
  
  # Resize to consistent size
  images <- image_resize(images, "800x600!")
  
  # Create montage
  montage <- image_montage(images, tile = "2x2", geometry = "800x600+10+10")
  
  # Add title
  montage <- image_annotate(montage, "Ultimate India 3D Population Maps", 
                           size = 30, gravity = "north", color = "white")
  
  # Save montage
  montage_file <- file.path(OUTPUT_DIR, output_name)
  image_write(montage, montage_file)
  
  cat("Montage created:", montage_file, "\n")
  return(montage_file)
}

# ============================================================================
# DOCUMENTATION AND METADATA
# ============================================================================

generate_comprehensive_documentation <- function() {
  cat("=== GENERATING COMPREHENSIVE DOCUMENTATION ===\n")
  
  doc_content <- paste0(
    "# Ultimate India 3D Population Visualization Project\n\n",
    "## Overview\n",
    "This project creates stunning 3D population density maps of India using advanced data processing and visualization techniques.\n\n",
    "## Features\n",
    "- **Multiple smoothing algorithms**: Gaussian, adaptive, and hybrid processing\n",
    "- **Intelligent scaling**: Cube root, fourth root, and hybrid scaling methods\n",
    "- **Advanced 3D rendering**: Multiple lighting, shadows, and quality settings\n",
    "- **6 artistic color palettes**: From traditional to modern styles\n",
    "- **Multiple viewing angles**: Perspective, aerial, coastline, and diagonal views\n",
    "- **High-quality exports**: Up to 4K resolution with anti-aliasing\n\n",
    "## Data Source\n",
    "- **File**: `ind_pd_2020_1km_ASCII_XYZ.csv`\n",
    "- **Resolution**: 1km grid cells\n",
    "- **Format**: XYZ coordinates with population density\n",
    "- **Coverage**: Complete India with boundaries\n\n",
    "## Processing Pipeline\n",
    "1. **Data Loading**: Read XYZ data and filter to India boundaries\n",
    "2. **Boundary Masking**: Use official India boundaries from giscoR\n",
    "3. **Elevation Context**: Add terrain data for geographic realism\n",
    "4. **Advanced Smoothing**: Reduce spikes while preserving detail\n",
    "5. **Intelligent Scaling**: Apply appropriate scaling for 3D visualization\n",
    "6. **Multi-layer Rendering**: Combine population, terrain, and lighting\n",
    "7. **Quality Export**: Generate multiple views and resolutions\n\n",
    "## Key Functions\n",
    "\n",
    "### Main Execution\n",
    "- `run_ultimate_india_pipeline()` - Complete pipeline with all features\n",
    "- `quick_india_3d()` - Fast single visualization\n",
    "- `batch_create_all_styles()` - Generate all combinations\n\n",
    "### Processing Options\n",
    "- `apply_gaussian_smooth()` - Gaussian kernel smoothing\n",
    "- `apply_adaptive_smooth()` - Variance-based adaptive smoothing\n",
    "- `apply_hybrid_processing()` - Multi-level processing\n",
    "- `apply_intelligent_scaling()` - Various scaling methods\n\n",
    "### Visualization\n",
    "- `create_ultimate_3d_map()` - Main 3D rendering function\n",
    "- `create_ultimate_color_palettes()` - Color scheme definitions\n",
    "- `create_processing_comparison()` - Method comparison\n\n",
    "### Analysis\n",
    "- `generate_statistics_report()` - Comprehensive data analysis\n",
    "- `test_scaling_methods()` - Compare scaling methods\n\n",
    "## Color Palettes\n",
    "1. **tricolor_gradient** - India flag inspired\n",
    "2. **population_spectrum** - Scientific population density\n",
    "3. **terrain_pop** - Terrain with population overlay\n",
    "4. **sunset_india** - Warm sunset colors\n",
    "5. **mono_elegant** - Monochromatic elegance\n",
    "6. **vibrant_contrast** - High contrast vibrant\n\n",
    "## Quality Settings\n",
    "- **high**: 1600x1600 window, 300 samples, fast render\n",
    "- **ultra**: 2000x2000 window, 500 samples, maximum quality\n\n",
    "## Output Files\n",
    "- `*_perspective_HQ.png` - Main high-quality view\n",
    "- `*_aerial.png` - Top-down boundary view\n",
    "- `*_coastline.png` - Coastal emphasis\n",
    "- `*_diagonal.png` - Diagonal artistic view\n",
    "- `india_population_statistics.csv` - Data analysis\n",
    "- `india_3d_gallery.html` - Web gallery\n\n",
    "## Technical Specifications\n",
    "- **Coordinate System**: Lambert Azimuthal Equal Area (India-centered)\n",
    "- **Smoothing Kernels**: Gaussian, adaptive, hybrid\n",
    "- **Scaling Methods**: Cube root, fourth root, logarithmic, hybrid\n",
    "- **Rendering Engine**: rayshader with rayrender\n",
    "- **Export Resolution**: Up to 4000x3000 pixels\n",
    "- **Lighting**: Multi-source realistic lighting\n\n",
    "## Performance Notes\n",
    "- High-quality renders may take 5-15 minutes\n",
    "- Ultra-quality renders may take 15-30 minutes\n",
    "- Memory usage: 2-8GB depending on resolution\n",
    "- Recommended: 16GB RAM for ultra-quality\n\n",
    "## Usage Examples\n",
    "\n",
    "```r\n",
    "# Quick start - single high-quality map\n",
    "files <- quick_india_3d()\n\n",
    "# Complete pipeline with comparisons\n",
    "all_files <- run_ultimate_india_pipeline(create_comparisons = TRUE, quality = \"ultra\")\n\n",
    "# Test different scaling methods\n",
    "test_files <- test_scaling_methods()\n\n",
    "# Create custom visualization\n",
    "data_list <- load_and_process_india_data()\n",
    "processed <- process_for_3d_visualization(data_list$population_raster, data_list$elevation_raster)\n",
    "final_data <- c(data_list, processed)\n",
    "custom_files <- create_ultimate_3d_map(final_data, palette_name = \"sunset_india\", quality = \"ultra\")\n```\n\n",
    "## Troubleshooting\n",
    "- **Out of memory**: Reduce quality setting or close other applications\n",
    "- **Slow rendering**: Use 'high' instead of 'ultra' quality\n",
    "- **RGL issues**: Update rgl package and graphics drivers\n",
    "- **File not found**: Check DATA_DIR path configuration\n\n",
    "Generated on: ", Sys.Date(), "\n",
    "R Version: ", R.version.string, "\n"
  )
  
  # Save documentation
  doc_file <- file.path(OUTPUT_DIR, "ULTIMATE_INDIA_3D_DOCUMENTATION.md")
  writeLines(doc_content, doc_file)
  
  cat("Comprehensive documentation saved:", doc_file, "\n")
  return(doc_file)
}

# ============================================================================
# FINAL EXECUTION INSTRUCTIONS
# ============================================================================

cat("\n")
cat(rep("=", 80), "\n")
cat("ULTIMATE INDIA 3D POPULATION VISUALIZATION - READY TO RUN\n")
cat(rep("=", 80), "\n")

cat("\n🚀 QUICK START OPTIONS:\n")
cat("1. quick_india_3d()                           # Fast single visualization\n")
cat("2. run_ultimate_india_pipeline()              # Complete pipeline\n")
cat("3. run_ultimate_india_pipeline(quality='ultra') # Maximum quality\n")
cat("4. test_scaling_methods()                     # Compare scaling approaches\n")

cat("\n🎨 CUSTOMIZATION OPTIONS:\n")
cat("1. select_palette_interactive()               # Choose color palette\n")
cat("2. batch_create_all_styles()                  # Generate all combinations\n")
cat("3. create_processing_comparison()             # Compare methods\n")

cat("\n📊 ANALYSIS AND EXPORT:\n")
cat("1. generate_statistics_report()               # Data analysis\n")
cat("2. create_web_gallery(files)                  # Web gallery\n")
cat("3. create_image_montage(files)               # Image montage\n")
cat("4. generate_comprehensive_documentation()     # Full documentation\n")

cat("\n⚙️ CONFIGURATION:\n")
cat("- Output directory: ", OUTPUT_DIR, "\n")
cat("- Data file: ", LOCAL_DATA_PATH, "\n")
cat("- Export resolution: ", EXPORT_WIDTH, "x", EXPORT_HEIGHT, "\n")

cat("\n🎯 RECOMMENDED FIRST RUN:\n")
cat("run_ultimate_india_pipeline(create_comparisons = FALSE, quality = 'high')\n")

cat("\n", rep("=", 80), "\n")
cat("ULTIMATE INDIA 3D VISUALIZATION SCRIPT LOADED - EXECUTE WHEN READY!\n")
cat(rep("=", 80), "\n\n")

# ============================================================================
# END OF ULTIMATE INDIA 3D VISUALIZATION SCRIPT
# ============================================================================