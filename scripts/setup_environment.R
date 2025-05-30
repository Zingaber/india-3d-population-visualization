# ============================================================================
# INDIA 3D POPULATION VISUALIZATION - ENVIRONMENT SETUP
# Generated on: 2025-05-30
# ============================================================================

cat('Setting up India 3D Population Visualization environment...\n')

# Core packages (required)
core_packages <- c(
  "tidyverse",
  "dplyr",
  "tidyr",
  "readr",
  "sf",
  "stars",
  "raster",
  "terra",
  "elevatr",
  "giscoR",
  "rayshader",
  "rayrender",
  "rgl",
  "magick",
  "imager",
  "pracma",
  "scales"
)

# Optional packages
optional_packages <- c(
  "httr",
  "jsonlite",
  "R.utils"
)

# Function to install packages safely
install_if_missing <- function(packages, type = 'core') {
  new_packages <- packages[!(packages %in% installed.packages()[,'Package'])]
  if(length(new_packages)) {
    cat('Installing', type, 'packages:', paste(new_packages, collapse = ', '), '\n')
    install.packages(new_packages, dependencies = TRUE)
  } else {
    cat('All', type, 'packages already installed\n')
  }
}

# Install packages
install_if_missing(core_packages, 'core')
install_if_missing(optional_packages, 'optional')

# Load core packages and check
cat('\nLoading packages...\n')
success <- sapply(core_packages, function(pkg) {
  tryCatch({
    library(pkg, character.only = TRUE)
    TRUE
  }, error = function(e) {
    cat('Failed to load:', pkg, '\n')
    FALSE
  })
})

# Summary
if(all(success)) {
  cat('\n✓ All packages installed and loaded successfully!\n')
  cat('✓ Environment ready for India 3D Population Visualization\n')
} else {
  cat('\n⚠ Some packages failed to load. Please check installation.\n')
}

# Display system info
cat('\n=== SYSTEM INFORMATION ===\n')
cat('R Version:', R.Version()$version.string, '\n')
cat('Platform:', R.Version()$platform, '\n')
cat('RAM:', round(as.numeric(system('wmic computersystem get TotalPhysicalMemory /value', intern=TRUE)[2]) / 1024^3), 'GB (approximate)\n')

