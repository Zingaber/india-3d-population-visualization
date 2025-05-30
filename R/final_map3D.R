
# Load required libraries
library(readr)
library(dplyr)
library(tidyr)     # For drop_na()
library(raster)
library(rayshader)
library(rgl)
library(imager)    # For smoothing

# Step 1: Load XYZ population data
file_path <- "C:/Users/josze/MyRWorkspace/Download/ind_pd_2020_1km_ASCII_XYZ/ind_pd_2020_1km_ASCII_XYZ.csv"

pop_xyz <- read_csv(file_path, col_names = c("x", "y", "z")) %>%
  mutate(across(everything(), as.numeric)) %>%
  drop_na()

# Step 2: Convert to raster and then to matrix
pop_raster <- rasterFromXYZ(pop_xyz)
pop_raster[is.na(pop_raster)] <- 0
pop_matrix <- raster_to_matrix(pop_raster)

# Step 3: Normalize, smooth, and scale
pop_matrix <- log10(pop_matrix + 1)
pop_matrix <- (pop_matrix - min(pop_matrix)) / (max(pop_matrix) - min(pop_matrix))

# Apply smoothing
pop_matrix <- as.matrix(isoblur(as.cimg(pop_matrix), 2))

# Optional exaggeration for visualization
pop_matrix <- pop_matrix * 100

# Step 4: Define improved color palette (greenish contrast)
color_palette <- colorRampPalette(c("#EDF8FB", "#B2E2E2", "#66C2A4", "#238B45"))(256)

# Step 5: Open 3D device
rgl::open3d()

# Step 6: Render 3D plot
hillshade <- height_shade(pop_matrix, texture = color_palette)

plot_3d(
  hillshade,
  heightmap = pop_matrix,
  zscale = 0.03,
  solid = FALSE,
  windowsize = c(1400, 1400),
  phi = 40,
  theta = 30,
  zoom = 0.6,
  background = "white"
)

# Step 7: Export a high-quality PNG
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
