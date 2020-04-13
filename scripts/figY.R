library(spant)   # mrs processing
library(ggplot2) # fancy plots
library(cowplot) # fancier plots
library(doParallel)

# change the working directory to the source file location
# when "Sourcing" from the RStudio GUI
if (Sys.getenv("RSTUDIO") == "1" & !is.null(parent.frame(2)$ofile)) {
  this.dir <- dirname(parent.frame(2)$ofile)
  setwd(this.dir)
}

# are we going to run the fitting in parallel?
parallel_fits <- TRUE

# use the same number of parallel jobs as we have available cores
jobs <- detectCores()

# mrs data file
mrs_f     <- "../data/2D_MRSI_artefact.rds"

# read the mrs data
mrs_data <- read_mrs(mrs_f, format = "rds")

# extract a row of voxels
mrs_data_cropped <- mrs_data %>% get_subset(x_set = 5:12, y_set = 7)

# simulate the basis set
basis <- sim_basis_1h_brain(seq_slaser_ideal,
                            acq_paras = get_acq_paras(mrs_data), TE1 = 0.012,
                            TE2 = 0.012, TE3 = 0.016, xlim = c(0.5, 4.1))

# all default fitting options
opts <- abfit_opts()

if (parallel_fits) {
  # clusters are platform dependant
  clust_type <- ifelse(.Platform$OS.type == "unix", "FORK", "PSOCK") 
  cl <- makeCluster(jobs, type = clust_type)
  registerDoParallel(cl)
}

res <- fit_mrs(mrs_data_cropped, opts = opts, basis = basis,
               parallel = parallel_fits, time = FALSE)

if (parallel_fits) stopCluster(cl)

plot_cex <- 0.6

plot_1 <- function() {
  par(cex = plot_cex)
  plot(res, x_pos = 1, restore_def_par = FALSE)
  ed_pppm <- sprintf('%.1f', res$res_tab$bl_ed_pppm[1])
  text(3.6, -7000, adj = 0, paste("ED per\nppm =", ed_pppm))
}

plot_2 <- function() {
  par(cex = plot_cex)
  plot(res, x_pos = 2, restore_def_par = FALSE)
  ed_pppm <- sprintf('%.1f', res$res_tab$bl_ed_pppm[2])
  text(1.6, 12000, adj = 0, paste("ED per\nppm =", ed_pppm))
}

plot_3 <- function() {
  par(cex = plot_cex)
  plot(res, x_pos = 3, restore_def_par = FALSE)
  ed_pppm <- sprintf('%.1f', res$res_tab$bl_ed_pppm[3])
  text(1.6, 8000, adj = 0, paste("ED per\nppm =", ed_pppm))
}

plot_4 <- function() {
  par(cex = plot_cex)
  plot(res, x_pos = 4, restore_def_par = FALSE)
  ed_pppm <- sprintf('%.1f', res$res_tab$bl_ed_pppm[4])
  text(1.6, 10200, adj = 0, paste("ED per\nppm =", ed_pppm))
}

plot_5 <- function() {
  par(cex = plot_cex)
  plot(res, x_pos = 5, restore_def_par = FALSE)
  ed_pppm <- sprintf('%.1f', res$res_tab$bl_ed_pppm[5])
  text(1.6, 10400, adj = 0, paste("ED per\nppm =", ed_pppm))
}

plot_6 <- function() {
  par(cex = plot_cex)
  plot(res, x_pos = 6, restore_def_par = FALSE)
  ed_pppm <- sprintf('%.1f', res$res_tab$bl_ed_pppm[6])
  text(1.6, 10000, adj = 0, paste("ED per\nppm =", ed_pppm))
}

plot_7 <- function() {
  par(cex = plot_cex)
  plot(res, x_pos = 7, restore_def_par = FALSE)
  ed_pppm <- sprintf('%.1f', res$res_tab$bl_ed_pppm[7])
  text(1.6, 10000, adj = 0, paste("ED per\nppm =", ed_pppm))
}

plot_8 <- function() {
  par(cex = plot_cex)
  plot(res, x_pos = 8, restore_def_par = FALSE)
  ed_pppm <- sprintf('%.1f', res$res_tab$bl_ed_pppm[8])
  text(1.6, 8800, adj = 0, paste("ED per\nppm =", ed_pppm))
}

labs <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H') 
full_plot <- plot_grid(plot_1, plot_2, plot_3, plot_4, plot_5, plot_6, plot_7,
                       plot_8, labels = labs, label_size = 10, ncol = 4,
                       label_x = -0.02)

cairo_pdf("../figures/figY.pdf", width = 6.92, height = 3.5)
print(full_plot)
dev.off()

tiff("../figures/figY.tiff", width = 300 * 6.92, height = 300 * 3.5, res = 300)
print(full_plot)
dev.off()