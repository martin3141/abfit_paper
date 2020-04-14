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
mrs_f     <- "../data/2D_MRSI.rds"

# mri data file
mri_f     <- "../data/T1_vol_deface.nii.gz"

# segmented mri data file
mri_seg_f <- "../data/brain_seg.nii.gz"

# read the mrs data
mrs_data <- read_mrs(mrs_f, "rds")

# extract the central 8x8 matrix of spectra
mrs_data_cropped <- mrs_data %>% crop_xy(8, 8)

# crop the corner voxels as they are close to the diagonal saturation bands
mask <- matrix(FALSE, 8, 8)
mask[1, 1] <- mask[1, 8] <- mask[8, 1] <- mask[8, 8] <- TRUE
mrs_data_cropped <- mrs_data_cropped %>% mask_xy_mat(mask)

# simulate the basis set
basis <- sim_basis_1h_brain(seq_slaser_ideal,
                            acq_paras = get_acq_paras(mrs_data), TE1 = 0.012,
                            TE2 = 0.012, TE3 = 0.016, xlim = c(0.5, 4.1))

# all default fitting options
opts <- abfit_opts()

fname <- "../data/fig9.rds"

if (file.exists(fname)) {  # don't recalc unless we have to
  cat("Reading precomputed results :", fname, "\n")
  res <- readRDS(fname) 
} else {
  if (parallel_fits) {
    # clusters are platform dependant
    clust_type <- ifelse(.Platform$OS.type == "unix", "FORK", "PSOCK") 
    cl <- makeCluster(jobs, type = clust_type)
    registerDoParallel(cl)
  }
  
  res <- fit_mrs(mrs_data_cropped, opts = opts, basis = basis,
                 parallel = parallel_fits, time = FALSE)
  
  if (parallel_fits) stopCluster(cl)
  cat("Saving precomputed results :", fname, "\n")
  saveRDS(res, fname)
}

# print data quality summary
summary(res)

# extract some metabolite maps
tnaa_tcr_metab_map  <- get_fit_map(res, "tNAA") / get_fit_map(res, "tCr")
tcho_tcr_metab_map  <- get_fit_map(res, "tCho") / get_fit_map(res, "tCr")
glx_tcr_metab_map   <- get_fit_map(res, "Glx")  / get_fit_map(res, "tCr")

# read the segmented T1 MRI
seg_data <- readNifti(mri_seg_f)

# get the approximate MRSI excitation region
mrsi_mask     <- get_mrsi_voi(mrs_data %>% crop_xy(10, 10), seg_data)

# mask the segmented data by the excitation region
seg_data_crop <- seg_data * mrsi_mask

# calculate the PSF for convolution with the segementation data
mat_size <- Nx(mrs_data)
psf_ker  <- Re(get_2d_psf(FOV = mrs_data$resolution[2] * mat_size + 1,
               mat_size = mat_size))

# calculate the CSF, WM, GM contribution to each voxel
seg_res <- get_mrsi2d_seg(mrs_data_cropped, seg_data_crop, psf_ker)

# make a data frame combining metabolite levels and GM fraction
metab_gmf <- data.frame(gmf = seg_res$GMF,
                        tnaa_tcr  = as.vector(tnaa_tcr_metab_map),
                        tcho_tcr  = as.vector(tcho_tcr_metab_map),
                        glx_tcr   = as.vector(glx_tcr_metab_map))

# remove any rows with missing values (eg masked voxels at the corners)
metab_gmf <- na.omit(metab_gmf)

# calculate some linear models
tnaa_tcr_res <- lm(tnaa_tcr ~ gmf, metab_gmf)
tcho_tcr_res <- lm(tcho_tcr ~ gmf, metab_gmf)
glx_tcr_res  <- lm(glx_tcr  ~ gmf, metab_gmf)

# extract some R2 and p vals for the linear model
tnaa_p_val <- summary(tnaa_tcr_res)$coefficients[2,4]
tnaa_R2    <- summary(tnaa_tcr_res)$r.squared

tcho_p_val <- summary(tcho_tcr_res)$coefficients[2,4]
tcho_R2    <- summary(tcho_tcr_res)$r.squared

glx_p_val  <- summary(glx_tcr_res)$coefficients[2,4]
glx_R2     <- summary(glx_tcr_res)$r.squared

# set ggplot2 theme
theme_set(theme_cowplot(font_size = 10))

tnaa_lab <- c(paste("italic(R)^2 ==", round(tnaa_R2, 2)),
              paste("italic(p)*'-value' ==", sprintf("%.1e", tnaa_p_val)))
            
tnaa_plt <- ggplot(metab_gmf, aes(x = gmf, y = tnaa_tcr)) + geom_point() +
            geom_smooth(method = "lm") +
            annotate("text", x = 43, y = c(1.98, 1.9), label = tnaa_lab,
                     parse = T, hjust = 0) + xlab("Gray matter fraction (%)") +
            ylab("tNAA / tCr")

tcho_lab <- c(paste("italic(R)^2 ==", round(tcho_R2, 2)),
              paste("italic(p)*'-value' ==", sprintf("%.1e", tcho_p_val)))
            
tcho_plt <- ggplot(metab_gmf, aes(x = gmf, y = tcho_tcr)) + geom_point() +
            geom_smooth(method = "lm") +
            annotate("text", x = 43, y = c(0.37, 0.35), label = tcho_lab,
                     parse = T, hjust = 0) + xlab("Gray matter fraction (%)") +
            ylab("tCho / tCr")

glx_lab <- c(paste("italic(R)^2 ==", round(glx_R2, 2)),
             paste("italic(p)*'-value' ==", sprintf("%.1e", glx_p_val)))
            
glx_plt <- ggplot(metab_gmf, aes(x = gmf, y = glx_tcr)) + geom_point() +
           geom_smooth(method = "lm") +
           annotate("text", x = 43, y = c(1.5, 1.4), label = glx_lab,
                    parse = T, hjust = 0) + xlab("Gray matter fraction (%)") +
           ylab("Glx / tCr")

# resice MRI to match MRSI
mri_data <- readNifti(mri_f)
mri_data_resliced <- reslice_to_mrs(mri_data, mrs_data)

# generate 8x8 metabolite map
mrsi_volume_map <- get_mrsi_voi(mrs_data_cropped, mri_data_resliced,
                                tnaa_tcr_metab_map)

im_plot_fn <- function() {
  ortho3(mri_data_resliced[,,90:256], mrsi_volume_map[,,90:256], 
         xyz = c(105, 112, 102), alpha = 0.6, zlim = c(300, 2000),
         zlim_ol = c(0.6, 1.8), bg = "white", mar = rep(0.5, 4), ch_lwd = 0.5)
}

full_plot <- plot_grid(im_plot_fn, tnaa_plt, tcho_plt, glx_plt,
                       labels = c('A', 'B', 'C', 'D'), label_size = 12,
                       rel_widths = c(1,1), ncol = 2)

cairo_pdf("../figures/figS7.pdf", width = 6.92, height = 5.5, pointsize = 10)
print(full_plot)
dev.off()

tiff("../figures/figS7.tiff", width = 300 * 6.92, height = 300 * 5.5, res = 300,
     pointsize = 10)
print(full_plot)
dev.off()
