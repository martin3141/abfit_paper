library(cowplot)
source("~/abfit/abfit.R")

# change the working directory to the source file location
# when "Sourcing" from the RStudio GUI
if (Sys.getenv("RSTUDIO") == "1" & !is.null(parent.frame(2)$ofile)) {
  print(parent.frame(2)$ofile)
  this.dir <- dirname(parent.frame(2)$ofile)
  setwd(this.dir)
}

set.seed(10)

calc_yhat_sc_basis <- function(y, spline_basis, deriv_mat, lambda) {
  inv_mat   <- MASS::ginv(rbind(spline_basis, lambda ^ 0.5 * deriv_mat))
  alpha     <- inv_mat %*% c(y, rep(0, dim(deriv_mat)[1]))
  yhat      <- spline_basis %*% alpha
  alpha_mat <- matrix(rep(alpha, dim(spline_basis)[1]), dim(spline_basis)[1],
                      length(alpha), byrow = TRUE)
  sc_basis  <- spline_basis * alpha_mat
  return(list(yhat = yhat, sc_basis = sc_basis))
}

# generate a smooth function
x      <- seq(0, pi, length.out = 200)
ytrue  <- sin(x) + 1.0

y      <- ytrue + rnorm(200, 0, 0.2)  # add noise

comps        <- 15                    # number of spline functions in the basis
spline_basis <- bbase(length(x), comps - 3)
deriv_mat    <- diff(diag(comps), differences = 2)


res   <- calc_yhat_sc_basis(y, spline_basis, deriv_mat, 0.0)
y_est <- res$yhat


p1 <- function() {
  par(bty = "L", mar = c(3.5, 2, 1.5, 0.5), mgp = c(2.0, 0.7, 0), cex = 0.7)
  offset <- matrix(rep(0:14, 200), nrow(spline_basis), ncol(spline_basis),
                   byrow = TRUE)
  
  matplot(x, spline_basis + offset / 1.2, type = "l", lty = 1, col = c("blue"),
          yaxt = "n", ylab = "")
}

p2 <- function() {
  par(bty = "L", mar = c(3.5, 4, 1.5, 0.5), mgp = c(2.0, 0.7, 0), cex = 0.7)
  matplot(x, res$sc_basis, type = "l", lty = 1, col = c("blue"),
          ylim = c(0, 2.4), ylab = "y")
  
  lines(x, y)
  lines(x, y_est, col = "red")
  legend(2.1, 2.6, c("data", "spline fit", "spline basis"),
         col = c("black", "red", "blue"), lty = c(1,1,1),
         xpd = TRUE, box.lty = 0, y.intersp = 1.2)
}

full_plot <- plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12,
                       rel_widths = c(0.92, 1))

# print(full_plot)

cairo_pdf("../figures/fig1.pdf", width = 6.92, height = 3.4, pointsize = 10)
print(full_plot)
dev.off()
