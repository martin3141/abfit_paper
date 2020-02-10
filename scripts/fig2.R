source("~/abfit/abfit.R")

library(ggplot2)
library(cowplot)

theme_set(theme_cowplot(font_size = 10))

# change the working directory to the source file location
# when "Sourcing" from the RStudio GUI
if (Sys.getenv("RSTUDIO") == "1" & !is.null(parent.frame(2)$ofile)) {
  this.dir <- dirname(parent.frame(2)$ofile)
  setwd(this.dir)
}

set.seed(1)

calc_yhat <- function(y, spline_basis, deriv_mat, lambda) {
  inv_mat <- MASS::ginv(rbind(spline_basis, lambda ^ 0.5 * deriv_mat))
  alpha   <- inv_mat %*% c(y, rep(0, dim(deriv_mat)[1]))
  yhat    <- spline_basis %*% alpha
  return(yhat)
}

# generate a smooth function
x      <- seq(0, 4.6 * pi, length.out = 200)
ytrue  <- sin(x)
y      <- ytrue + rnorm(200, 0, 0.3)  # add noise

comps        <- 50                    # number of spline functions in the basis
spline_basis <- bbase(length(x), comps - 3)
deriv_mat    <- diff(diag(comps), differences = 2)

auto_smo_res <- auto_pspline_smoother(y, spline_basis, deriv_mat)

y_true_lab <- "true : y = sin(x)"
y1 <- calc_yhat(y, spline_basis, deriv_mat, 0.0)
y1_lab <- "λ=0, ED=50"
y2 <- calc_yhat(y, spline_basis, deriv_mat, 0.1)
y2_lab <- "λ=0.1, ED=33"
y3 <- calc_yhat(y, spline_basis, deriv_mat, 20)
y3_lab <- "λ=20, ED=12"
y4 <- calc_yhat(y, spline_basis, deriv_mat, 500)
y4_lab <- "λ=500, ED=6"
y5 <- calc_yhat(y, spline_basis, deriv_mat, 1e8)
y5_lab <- "λ=1e8, ED=2"

df_est <- rbind(data.frame(x = x, y = y,     c1 = "data",     c2 = y_true_lab),
                data.frame(x = x, y = ytrue, c1 = "estimate", c2 = y_true_lab),
                data.frame(x = x, y = y,     c1 = "data",     c2 = y1_lab),
                data.frame(x = x, y = y1,    c1 = "estimate", c2 = y1_lab),
                data.frame(x = x, y = y,     c1 = "data",     c2 = y2_lab),
                data.frame(x = x, y = y2,    c1 = "estimate", c2 = y2_lab),
                data.frame(x = x, y = y,     c1 = "data",     c2 = y3_lab),
                data.frame(x = x, y = y3,    c1 = "estimate", c2 = y3_lab),
                data.frame(x = x, y = y,     c1 = "data",     c2 = y4_lab),
                data.frame(x = x, y = y4,    c1 = "estimate", c2 = y4_lab),
                data.frame(x = x, y = y,     c1 = "data",     c2 = y5_lab),
                data.frame(x = x, y = y5,    c1 = "estimate", c2 = y5_lab))
    
p1 <- ggplot(data = df_est, aes(x = x, y = y, col = c1)) + geom_line() + 
             facet_wrap(~ c2, ncol = 2) +
             scale_colour_manual(values = c("black", "red")) + 
             theme(legend.title = element_blank(),
                   legend.position = c(0.44, 0.75))

lambda_start <- 1e-5
lambda_end   <- 1e8
lambda_n     <- 100
lambda_vec   <- 10 ^ (seq(log10(lambda_start), log10(lambda_end),
                              length.out = lambda_n))
ed_vec    <- rep(NA, lambda_n)
resid_vec <- rep(NA, lambda_n)
error_vec <- rep(NA, lambda_n)
aic_vec   <- rep(NA, lambda_n)

for (n in 1:lambda_n) {
  ed_vec[n]    <- calc_ed_from_lambda(spline_basis, deriv_mat, lambda_vec[n])
  yhat         <- calc_yhat(y, spline_basis, deriv_mat, lambda_vec[n])
  resid_vec[n] <- sum((y - yhat) ^ 2)
  error_vec[n] <- sum((ytrue - yhat) ^ 2)
  aic_vec[n]   <- log(resid_vec[n]) + 2 * ed_vec[n] / 200
}

error_vec[60:lambda_n] <- NA # helps with plot clarity

df_diag <- rbind(data.frame(x = lambda_vec, y = resid_vec, c1 = "residual"),
                 data.frame(x = lambda_vec, y = ed_vec,    c1 = "ED"),
                 data.frame(x = lambda_vec, y = error_vec, c1 = "error"),
                 data.frame(x = lambda_vec, y = aic_vec,   c1 = "AIC"))

breaks <- c(1e-5, 1e-2, 1, 20, 1000, 1e6, 1e8)
labs   <- c("1e-5", "1e-2", "1", "20", "1000", "1e6", "1e8")
p2 <- ggplot(data = df_diag, aes(x = x, y = y)) + geom_line() + 
       facet_wrap(~ c1, nrow = 4, scales = "free_y") + 
       scale_x_continuous(trans='log10', breaks = breaks, labels = labs) + 
       xlab("λ") + ylab("Value")


full_plot <- plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12,
                       rel_widths = c(3,2))

# print(full_plot)

cairo_pdf("../figures/fig2.pdf", width = 6.92, height = 5.5)
print(full_plot)
dev.off()

tiff("../figures/fig2.tiff", width = 300 * 6.92, height = 300 * 5.5,
     pointsize = 10, res = 300)
print(full_plot)
dev.off()
