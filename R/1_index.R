
# -------------------------------------------------------------------------
# Properties of the Normalised Spatial Dispersion Index 
# -------------------------------------------------------------------------

library(pbapply)
library(dplyr)

source("R/utils.R")


# Illustration of the index -----------------------------------------------

# The total number of birds within a square.
n <- 20

# Low spatial dispersion: all the birds in a single segment.
x1 <- c(0, 0, n, 0, 0, 0, 0, 0, 0, 0)
xy <- inv_obs(x1)
plot_obs(xy)
nsdi(x1, 1000)

# Spatial randomness
x2 <- csr(n)
xy <- inv_obs(x2)
plot_obs(xy)
nsdi(x2, 1000)

# Uniform distribution
x3 <- rep(2, 10)
xy <- inv_obs(x3)
plot_obs(xy)
nsdi(x3, 1000)


# Test the independence of the index with respect to total abundance --------

# Generate a spatial pattern for each abundance.
# The returned value is the average of 10,000 repetitions.
# Two boundry cases are considered: complete randomness and complete clustering.

# Range of abundances tested
n <- seq(3, 100, 1)

# 1. Under the assumption of complete spatial randomness ----

# Setting the seed for reproducibility
set.seed(1)

ir <- pbsapply(n, function(ni) {
  xi <- csr(ni)
  nsdi(xi, n_sim = 10000)
})

op <- par(mar = c(5, 5, 2, 2))
plot(n, ir, xlab = "Total abundance", ylab = "Spatial dispersion", cex.lab = 1.5)
lines(n, ir, col = adjustcolor(1, alpha = 0.2))
abline(h = 1)
par(op)


# 2. Under the assumption of complete clustering ----

# Setting the seed for reproducibility
set.seed(1)

ic <- pbsapply(n, function(ni) {
  xi <- sample(c(ni, rep(0, 9)))
  nsdi(xi, n_sim = 10000)
})

op <- par(mar = c(5, 5, 2, 2))
plot(n, ic, xlab = "Total abundance", ylab = "Spatial dispersion", cex.lab = 1.5)
lines(n, ic, col = adjustcolor(1, alpha = 0.2))
abline(h = mean(ic))
par(op)


# Generating the dummy data -----------------------------------------------

# The raw bird counts are provided by BirdLife Poland and the Chief Inspectorate 
# for Environmental Protection (GIOŚ). The GIOŚ data policy does not permit
# the release of raw data. 
# For this reason, we are only making the processed data available 
# (in the 'data/data.RDS' and 'data/data.csv' files), which includes densities estimated using the 
# distance sampling method and the Normalised Spatial Dispersion Index (NSDI) calculated.  The original observers' and plots' IDs were pseudonymised.

# However, to enable the replication of the entire data processing chain, 
# we generate sample data below, which is an exact replica of the original data. 
# We then use this data to demonstrate how to calculate the NSDI.

n_years <- 10
n_plots <- 5
n_segments <- 10
years <- seq(n_years) + 1999

base_grid <- expand.grid(year = years, id = seq(n_plots), segment = seq(n_segments))
base_grid$id <- factor(base_grid$id)
base_grid <- transform(base_grid, id_year = interaction(id, year, sep = "_"))
base_grid <- transform(base_grid, id_year_seg = interaction(id_year, segment , sep = "_"))

# Empirical frequency distribution of observed counts
count <- c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
           22, 23, 24)
prob <- c(0.3435, 0.2716, 0.1262, 0.0927, 0.0543, 0.024, 0.0272, 0.0208, 0.0064,
          0.0096, 0.0016, 0.0064, 0.0016, 0.0032, 0.0032, 0.0016, 0.0016, 0.0016,
          0, 0, 0.0016, 0.0016)

# Setting the seed for reproducibility
set.seed(1)

# Generate random counts using the observed frequency
base_grid$n <- sample(count, size = nrow(base_grid), replace = TRUE, prob = prob)

# Generate random observer's ID per plot
obs_id <- factor(sample(1:10000, n_plots))
base_grid$observer_id <- obs_id[match(base_grid$id, unique(base_grid$id))]

# Final check
head(base_grid[base_grid$id == 1, ])
tail(base_grid[base_grid$id == 1, ])

# Add random structure ----

ddat <- base_grid |>
  dplyr::group_by(id) |>
  dplyr::arrange(year , .by_group = TRUE)  |>
  dplyr::group_modify(~ {
    size <- sample(nrow(.x) , 1)
    dplyr::slice_sample(.x , n = size)
  }) |>
  dplyr::ungroup()


# Calculate the index -----------------------------------------------------

ddat_splitted <- split(ddat, ddat$id, drop = TRUE)

set.seed(1)
ddat_splitted <- pblapply(ddat_splitted, disp, n_sim = 1000)
# ddat_splitted

ddat <- do.call(rbind.data.frame, ddat_splitted)

summary(ddat)
