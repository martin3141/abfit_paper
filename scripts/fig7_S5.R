library(spant)
library(doParallel)
library(ggplot2)
library(cowplot)

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

theme_set(theme_cowplot(font_size = 10))

ft <- def_acq_paras()$ft

ala    <- get_mol_paras("ala")
asp    <- get_mol_paras("asp")
cr     <- get_mol_paras("cr")
gaba   <- get_mol_paras("gaba")
glc    <- get_mol_paras("glc")
gln    <- get_mol_paras("gln")
glu    <- get_mol_paras("glu")
gpc    <- get_mol_paras("gpc")
gsh    <- get_mol_paras("gsh")
ins    <- get_mol_paras("ins")
lac    <- get_mol_paras("lac")
naa    <- get_mol_paras("naa")
naag   <- get_mol_paras("naag")
pch    <- get_mol_paras("pch")
pcr    <- get_mol_paras("pcr")
sins   <- get_mol_paras("sins")
tau    <- get_mol_paras("tau")
mm_exp <- get_mol_paras("mm_3t", ft)

metab_basis_list <- list(ala, asp, cr, gaba, glc, gln, glu, gpc, gsh, ins, lac,
                         naa, naag, pch, pcr, sins, tau)

full_basis_list  <- append(metab_basis_list, list(mm_exp))

full_basis <- sim_basis(full_basis_list, pul_seq = seq_slaser_ideal,
                        xlim = c(0.5, 4.2))

# metab values from de Graff book

amps <- c( 0.80,  # 1  Ala
           1.00,  # 2  Asp
           7.50,  # 3  Cr
           1.50,  # 4  GABA
           1.50,  # 5  Glc
           4.50,  # 6  Gln
           9.25,  # 7  Glu
           1.00,  # 8  GPC
           2.25,  # 9  GSH
           6.50,  # 10 Ins
           0.60,  # 11 Lac
          12.25,  # 13 NAA
           1.50,  # 14 NAAG
           0.60,  # 15 PCh
           4.25,  # 16 PCr
           0.35,  # 17 sIns
           4.00,  # 18 Tau
          30.00)  # 12 MMexp

set.seed(1)

# simulate mrs data
lb_para    <- 6
noise_N    <- 32
metab_mm   <- basis2mrs_data(full_basis, sum_elements = TRUE, amp = amps)
amps_no_mm <- c(amps[1:(length(amps) - 1)], 0)
metab      <- basis2mrs_data(full_basis, sum_elements = TRUE, amp = amps_no_mm)
metab      <- lb(metab, lb_para)
broad_sig  <- sim_resonances(freq = 1.3, amp = 0, lw = 100, lg = 1) # zero

mrs_data_nn    <- lb(metab_mm, lb_para)     # no noise data
mrs_data_noise <- sim_noise(sd = 2.0, fd = FALSE, dyns = noise_N)
mrs_data       <- rep_dyn(mrs_data_nn, noise_N) + mrs_data_noise

ed_pppm_start <- 2.01 / (4 - 0.2)
ed_pppm_end   <- 15
ed_pppm_N     <- 15
ed_pppm_vec   <- 10 ^ (seq(log10(ed_pppm_start), log10(ed_pppm_end),
                           length.out = ed_pppm_N))

lip09  <- get_mol_paras("lip09",  ft)
lip13a <- get_mol_paras("lip13a", ft)
lip13b <- get_mol_paras("lip13b", ft)
lip20  <- get_mol_paras("lip20",  ft)
mm09   <- get_mol_paras("mm09",   ft)
mm12   <- get_mol_paras("mm12",   ft)
mm14   <- get_mol_paras("mm14",   ft)
mm17   <- get_mol_paras("mm17",   ft)
mm20   <- get_mol_paras("mm20",   ft)

lip_mm_basis_list <- list(lip09, lip13a, lip13b, lip20, mm09, mm12, mm14, mm17,
                          mm20)

approx_basis_list <- append(metab_basis_list, lip_mm_basis_list)
approx_basis <- sim_basis(approx_basis_list, pul_seq = seq_slaser_ideal,
                        xlim = c(0.5, 4.2))

fname <- "../data/fig8_S5.rds"        # precomputed results

if (file.exists(fname)) {  # don't recalc unless we have to
  cat("Reading precomputed results :", fname, "\n")
  res_list <- readRDS(fname) 
} else {
  if (parallel_fits) {
    # clusters are platform dependant
    clust_type <- ifelse(.Platform$OS.type == "unix", "FORK", "PSOCK") 
    cl <- makeCluster(jobs, type = clust_type)
    registerDoParallel(cl)
  }
  
  res_list <- vector(mode = "list", length = (ed_pppm_N + 1))
  for (n in 1:ed_pppm_N) {
    
    opts  <- abfit_opts(auto_bl_flex = FALSE, bl_ed_pppm = ed_pppm_vec[n],
                        bl_comps_pppm = 25, pre_fit_bl_ed_pppm = 10)
    
    res_list[[n]] <- fit_mrs(mrs_data, method = "abfit", opts = opts,
                             basis = approx_basis, parallel = parallel_fits,
                             time = FALSE)
  }

  opts  <- abfit_opts(max_bl_ed_pppm = 15, bl_comps_pppm = 25,
                      pre_fit_bl_ed_pppm = 10)
  
  res_list[[n + 1]] <- fit_mrs(mrs_data, method = "abfit", opts = opts,
                               basis = approx_basis, parallel = parallel_fits,
                               time = FALSE)
  
  if (parallel_fits) stopCluster(cl)
  cat("Saving precomputed results :", fname, "\n")
  saveRDS(res_list, fname)
}

error_vec        <- rep(NA, ed_pppm_N)
sd_error_vec     <- rep(NA, ed_pppm_N)

# calc amp est errors
for (n in 1:ed_pppm_N) {
  amp_inds <- c(6:22)
  fit_amp_mat  <- res_list[[n]]$res_tab[amp_inds]
  true_amp_mat <- matrix(amps[-18], nrow(fit_amp_mat), ncol(fit_amp_mat),
                         byrow = TRUE)
  
  error           <- (true_amp_mat - fit_amp_mat) ^ 2
  mean_error      <- mean(rowSums(error))
  sd_error        <- sd(rowSums(error))
  error_vec[n]    <- mean_error
  sd_error_vec[n] <- sd_error
}

df <- data.frame(ed_pppm_vec, error_vec, sd_error_vec)

mean_ed_pppm <- mean(res_list[[16]]$res_tab$bl_ed_pppm)

breaks <- c(0.5, 1, 2, 3, 5, 10, 25)

p1 <- ggplot(data = df, aes(x = ed_pppm_vec)) + 
             geom_line(aes(y = error_vec)) + geom_point(aes(y = error_vec)) +
             geom_errorbar(aes(ymin = error_vec - sd_error_vec,
                           ymax = error_vec + sd_error_vec), width = .02) +
             scale_x_continuous(trans = "log10", breaks = breaks) + 
             xlab("Baseline ED per ppm") + 
             ylab("Metabolite estimate error") +
             geom_vline(xintercept = mean_ed_pppm, linetype = "dashed")

p2 <- function() {
  par(cex = 0.75)
  plot(res_list[[1]], restore_def_par = FALSE)
}

# find a fit closest to the automically determined mean ED pppm
closest_dyn <- which.min((mean_ed_pppm - res_list[[16]]$res_tab$bl_ed_pppm) ^ 2)

p3 <- function() {
  par(cex = 0.75)
  plot(res_list[[16]], restore_def_par = FALSE, dyn = closest_dyn)
}

p4 <- function() {
  par(cex = 0.75)
  plot(res_list[[15]], restore_def_par = FALSE)
}

full_plot <- plot_grid(p1, p2, p3, p4, labels = c('A', 'B', 'C', 'D'),
                       label_size = 12, rel_widths = c(1,1,1,1), ncol = 2,
                       label_x = -0.01)

# print(full_plot)

cairo_pdf("../figures/fig8.pdf", width = 6.92, height = 5.5)
print(full_plot)
dev.off()

tiff("../figures/fig8.tiff", width = 300 * 6.92, height = 300 * 5.5, res = 300)
print(full_plot)
dev.off()

# supp. figure

true_full  <- crop_spec(td2fd(zf(mrs_data_nn)), xlim = c(4, 0.2))
true_metab <- crop_spec(td2fd(zf(metab)), xlim = c(4, 0.2))
true_bl    <- crop_spec(td2fd(zf(broad_sig)), xlim = c(4, 0.2))
true_mm    <- crop_spec(td2fd(zf(mrs_data_nn - metab - broad_sig)),
                        xlim = c(4, 0.2))
true_noise <- crop_spec(td2fd(zf(get_dyns(mrs_data_noise, closest_dyn))),
                        xlim = c(4, 0.2))

dummy <- true_full
fit_tab  <- res_list[[16]]$fits[[closest_dyn]]

est_full <- dummy
est_full$data[1,1,1,1,1,1,] <- fit_tab$Fit + fit_tab$Baseline

est_bl <- dummy
est_bl$data[1,1,1,1,1,1,] <- fit_tab$Baseline

est_mm <- dummy
est_mm$data[1,1,1,1,1,1,] <- rowSums(fit_tab[22:30])

est_metab <- dummy
est_metab$data[1,1,1,1,1,1,] <- rowSums(fit_tab[5:21])

labs <- c("true", "est.", "resid.", "noise")
# full
sp1 <- function() {
  stacked_data <- append_dyns(true_full, est_full, true_full - est_full,
                              true_noise)
  stackplot(stacked_data, xlim = c(4, 0.2), restore_def_par = FALSE,
            y_offset = 10, bl_lty = 2, labels = labs, right_marg = 3)
}

# metab
sp2 <- function() {
  stacked_data <- append_dyns(true_metab, est_metab, true_metab - est_metab,
                              true_noise)  
  stackplot(stacked_data, xlim = c(4, 0.2), restore_def_par = FALSE,
            y_offset = 20, bl_lty = 2, labels = labs, right_marg = 3)
}

# bl
sp3 <- function() {
  stacked_data <- append_dyns(true_bl, est_bl, true_bl - est_bl, true_noise)  
  stackplot(stacked_data, xlim = c(4, 0.2), restore_def_par = FALSE,
            y_offset = 300, bl_lty = 2, labels = labs, right_marg = 3)
}

# mm
sp4 <- function() {
  stacked_data <- append_dyns(true_mm, est_mm, true_mm - est_mm, true_noise)  
  stackplot(stacked_data, xlim = c(4, 0.2), restore_def_par = FALSE,
            y_offset = 200, bl_lty = 2, labels = labs, right_marg = 3)
}

full_plot_supp <- plot_grid(sp1, sp2, sp3, sp4, labels = c('A', 'B', 'C', 'D'),
                            label_size = 12, rel_widths = c(1,1,1,1), ncol = 2)

print(full_plot_supp)

cairo_pdf("../figures/figS5.pdf", width = 6.92, height = 5.5, pointsize = 10)
print(full_plot_supp)
dev.off()

tiff("../figures/figS5.tiff", width = 300 * 6.92, height = 300 * 5.5,
     res = 300, pointsize = 10)
print(full_plot_supp)
dev.off()