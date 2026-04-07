
# -------------------------------------------------------------------------
# Fit a Generalised Additive Mixed Model
# -------------------------------------------------------------------------

library(mgcv)

# Data --------------------------------------------------------------------

data <- readRDS("data/data.RDS")

data <- droplevels(data)
summary(data)

summary(data$dispersion)

hist(data$dispersion, breaks = 30)
hist(log(data$dispersion), breaks = 30)


# GAMM --------------------------------------------------------------------

nthreads <- 8L  # no. of cores to be used

ctrl <- list(ncv.threads = nthreads)
# ctrl <- list(maxit = 1000, epsilon = 1e-05)


# List of predictors:
dput(names(data))

v <- c("log_Lapwing", "log_Crow", "Arable_land", "Grasslands_Pastures", "Precip_winter", "Tmin_winter")

# f <- formula(paste0("dispersion ~ ", paste0("s(", v, ", k = k)", collapse = " + ")))

f <- dispersion ~ s(log_Lapwing, k = k) + s(log_Crow, k = k) + ti(log_Crow, log_Lapwing, k = k) + s(Arable_land, k = k) + s(Grasslands_Pastures, k = k) + s(Precip_winter, k = k) + s(Tmin_winter, k = k)

# f <- dispersion ~ s(log_Crow, log_Lapwing, k = k) + s(Arable_land, k = k) + s(Grasslands_Pastures, k = k) + s(Precip_winter, k = k) + s(Tmin_winter, k = k)


# f <- dispersion ~ s(log_Lapwing, k = k) + s(log_Crow, by = n_cat, k = k) + s(Arable_land, k = k) + s(Grasslands_Pastures, k = k) + s(Precip_winter, k = k) + s(Tmin_winter, k = k)


# f01 <- update(f, ". ~ . + s(observer_id, bs = 're') + s(plot_id, bs = 're') + s(year, k = 10, bs = 'gp', m = c(3, 1)) + s(x, y, k = 20, bs = 'gp')")

f01 <- update(f, ". ~ . + s(observer_id, bs = 're') + s(plot_id, bs = 're')")

f01

k <- 6

# Full model
m <- gam(f01, data, family = Gamma("log"), control = ctrl)

summary(m)

op <- par(mfrow = c(2, 2)); gam.check(m, rep = 100); par(op)


# Response curves
op <- par(mar = c(5, 4, 1, 1))
plot.gam(m, scale = 0, scheme = 2, pages = 1, seWithMean = TRUE, residuals = TRUE, trans = I, pch = 21, cex = 0.3, by.resids = TRUE, all.terms = TRUE, too.far = 0.2)
par(op)

# Response curves - no points
op <- par(mar = c(5, 4, 1, 1))
plot.gam(m, scale = -1, scheme = 2, pages = 1, seWithMean = TRUE, residuals = FALSE, trans = I, pch = 21, cex = 0.3, all.terms = TRUE, too.far = 0.2)
par(op)

saveRDS(m, file = ("data/gamm_model.RDS"))


