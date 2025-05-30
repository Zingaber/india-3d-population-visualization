# map3D.R
# =============================
# India Population Density 3D Map (Smoothed + Color Gradient)
# =============================

# Load required libraries
library(readr)
library(dplyr)
library(tidyr)       # For drop_na()
library(raster)
library(rayshader)
library(rgl)
library(imager)      # For image smoothing

# Step 1: Load XYZ population data
file_path <- "C:/Users/josze/MyRWorkspace/Download/ind_pd_2020_1km_ASCII_XYZ/ind_pd_2020_1km_ASCII_XYZ.csv"

pop_xyz <- read_csv(
  file = file_path,
  col_types = cols(
    x = col_double(),
    y = col_double(),
    z = col_double()
  ),
  skip = 1   # skip header if needed
) %>%
  drop_na()

# Step 2: Convert to raster and then to matrix
pop_raster <- rasterFromXYZ(pop_xyz)
pop_raster[is.na(pop_raster)] <- 0
pop_matrix <- raster_to_matrix(pop_raster)

cat("Dimensions of matrix are:", dim(pop_matrix)[1], "x", dim(pop_matrix)[2], "\n")

# Step 3: Normalize, log-scale, smooth
pop_matrix <- log10(pop_matrix + 1)                          # Log scaling to reduce outliers
pop_matrix <- (pop_matrix - min(pop_matrix)) / 
              (max(pop_matrix) - min(pop_matrix))           # Normalize to 0-1
pop_matrix <- as.matrix(isoblur(as.cimg(pop_matrix), 1))    # Smooth with Gaussian blur
pop_matrix <- pop_matrix * 100                              # Optional exaggeration

# Step 4: Define color palette
color_palette <- colorRampPalette(c("#F7FCF0", "#41B6C4", "#253494"))(256)

# Step 5: Generate hillshade (colored overlay)
hillshade <- height_shade(pop_matrix, texture = color_palette)

# Step 6: Render 3D plot
rgl::open3d()

plot_3d(
  hillshade,
  heightmap = pop_matrix,
  zscale = 0.08,                      # 🔽 Lower value = less vertical exaggeration
  solid = FALSE,
  windowsize = c(1400, 1400),
  phi = 40,
  theta = 30,
  zoom = 0.6,
  background = "white"
)

# Step 7: Export high-quality image
render_highquality(
  filename = "india_population_3d_spike.png",
  lightdirection = 225,
  lightaltitude = 60,
  lightintensity = 450,
  interactive = FALSE,
  width = 2000,
  height = 2000
)

# Optional quick snapshot
render_snapshot("india_population_3d_snapshot.png")
