library(spant)   # mrs processing
library(ggplot2) # fancy plots
library(cowplot) # fancier plots

# change the working directory to the source file location
# when "Sourcing" from the RStudio GUI
if (Sys.getenv("RSTUDIO") == "1" & !is.null(parent.frame(2)$ofile)) {
  this.dir <- dirname(parent.frame(2)$ofile)
  setwd(this.dir)
}

# mrs data file
mrs_f    <- "../data/2D_MRSI.rds"

# read the mrs data
mrs_data <- readRDS(mrs_f)

# extract a voxel with asymmetric lineshape
asy_vox <- get_voxel(mrs_data, 9, 12)

# simulate the basis set
basis <- sim_basis_1h_brain(seq_slaser_ideal,
                            acq_paras = get_acq_paras(mrs_data), TE1 = 0.012,
                            TE2 = 0.012, TE3 = 0.016, xlim = c(0.5, 4.1))

# default fitting options
opts <- abfit_opts()
res_asym <- fit_mrs(asy_vox, opts = opts, basis = basis)

# restrict the asym parameter to force a symmetric lineshape
opts <- abfit_opts(max_asym = 1e-4)
res_sym <- fit_mrs(asy_vox, opts = opts, basis = basis)

sym_plot_fn  <- function() {
  par(cex = 0.75)
  plot(res_sym)
}

asym_plot_fn <- function() {
  par(cex = 0.75)
  plot(res_asym)
}

full_plot <- plot_grid(asym_plot_fn, sym_plot_fn, labels = c('A', 'B'),
                       label_size = 12, rel_widths = c(1, 1), ncol = 2)

cairo_pdf("../figures/fig10.pdf", width = 6.92, height = 3.5)
print(full_plot)
dev.off()

tiff("../figures/fig10.tiff", width = 300 * 6.92, height = 300 * 3.5,
     pointsize = 10, res = 300)
print(full_plot)
dev.off()