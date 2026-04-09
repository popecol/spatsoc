
# -------------------------------------------------------------------------
# Visualisation of the results
# -------------------------------------------------------------------------

source("R/utils.R")

dir <- dir()
if(sum(dir == "Figures") == 0) {
  dir.create("Figures")
} 

# Data & Setup ------------------------------------------------------------

data <- readRDS("data/data.RDS")

m <- readRDS("data/gamm_model.RDS")

model_terms <- terms(m)
is_re <- vapply(m$smooth, inherits, logical(1), what = "random.effect")
v <- all.vars(model_terms)[!is_re][-1]

v_name <- c("Lapwing density", "Crow density", "Arable land", "Grassland", "Winter precipitation", "Winter temperature")
n_var <- length(v)
s_table <- data.frame(summary(m)[["s.table"]])
nam <- rownames(s_table)
p <- s_table[grep("^s", nam), ][["p.value"]]
p_val <- format_pval(p, show_stars = TRUE)


# Figure 1 ----------------------------------------------------------------
# Frequency distribution of the Normalised Spatial Dispersion Index (NSDI)

col <- adjustcolor("black", alpha = 0.5)

ylab <- "Normalised Spatial Dispersion Index"

mult <- 1.5

# Dimensions in millimetres:
w <- 85; h <- 80

# Dimensions in inches * multiplier:
wi <- round(mult * w / 25.4, 1)
hi <- round(mult * h / 25.4, 1)

cairo_pdf("Figures/Figure_1.pdf", width = wi, height = hi, symbolfamily = "OpenSymbol")

op <- par(mar = c(5, 5, 4, 2))

hist(data$dispersion, breaks = "Scott", main = "",
     col = col, xlab = ylab, cex.lab = 1.5, cex.axis = 1.2, ylab = "Frequency")
box(col = "black", lty = "solid")

par(op)

dev.off()


# Figure 2 ----------------------------------------------------------------
# Marginal main effects of socio-ecological predictors 

col <- adjustcolor("black", alpha = 0.8)

mult <- 1.3

# Dimensions in millimetres:
w <- 170; h <- 150

# Dimensions in inches * multiplier:
wi <- round(mult * w / 25.4, 1)
hi <- round(mult * h / 25.4, 1)

cairo_pdf("Figures/Figure_2.pdf", width = wi, height = hi, symbolfamily = "OpenSymbol")

op <- par(mfrow = c(2, 3), oma = c(0, 5, 0, 0), mar = c(5, 2, 2, 2)) 

for (i in seq_along(v)) {
  term_plot(m, v[i], col_main = col, xlab = v_name[i], ylab = "", cex.lab = 2, ylim = c(0.2, 1.2), cex.axis = 1.5) # ylim = c(0.28, 1.15)
  legend("topleft", p_val[i], bty = "n", cex = 1.5)
  mtext(letters[i], side = 3, line = 0.5, cex = 1.2, adj = 0, font = 2) 
}
mtext(ylab, side = 2, line = 2, outer = TRUE, cex = 1.7)

par(op)

dev.off()


# Figure 3 ----------------------------------------------------------------
# Interaction effects between local conspecific density and predation risk

# Levels of crow abundance
quantile(data$Crow, probs = c(0, 0.5, 0.95), na.rm = TRUE)
crow_lev <- c(0, 4)
log1p(crow_lev)

# Levels of lapwing abundance
quantile(data$Lapwing, probs = c(0, 0.5, 0.95), na.rm = TRUE)
lap_lev <- c(3, 9)
log1p(lap_lev)

cols_lap <- c("#56B4E9", adjustcolor("#56B4E9", alpha = 0.5))

cols_crow <- c("#D55E00",adjustcolor("#D55E00", alpha = 0.5))

mult <- 1.5

# Dimensions in millimetres:
w <- 110; h <- 75

# Dimensions in inches * multiplier:
wi <- round(mult * w / 25.4, 1)
hi <- round(mult * h / 25.4, 1)

# cex.lab <- 1.5

cairo_pdf("Figures/Figure_3.pdf", width = wi, height = hi, symbolfamily = "OpenSymbol")

op <- par(mfrow = c(1, 2), mar = c(5, 5, 2, 2))

# Panel (a)
term_plot(m, ~ log_Lapwing + log_Crow, group_levels = log1p(crow_lev), xlab = "Lapwing density", ylab = ylab, cex.lab = 1.2, col_pal = cols_crow, ylim = c(0.2, 1.2))

mtext(letters[1], side = 3, line = 0.5, cex = 1.2, adj = 0, font = 2)

legend("topleft", title = "Crow abundance", legend = crow_lev, col = cols_crow, lty = 1:2, lwd = 2, bty = "n")

# Panel (b)
term_plot(model = m, formula = ~ log_Crow + log_Lapwing, group_levels = log1p(lap_lev), xlab = "Crow density", ylab = "", cex.lab = 1.2, col_pal = cols_lap, ylim = c(0.2, 1.2))

mtext(letters[2], side = 3, line = 0.5, cex = 1, adj = 0, font = 2)

legend("topleft", title = "Lapwing abundance", legend = lap_lev, col = cols_lap, lty = 1:2, lwd = 2, bty = "n")

par(op)

dev.off()
