library(spant)
library(ggplot2)
library(cowplot)

theme_set(theme_cowplot(font_size = 10))

# change the working directory to the source file location
# when "Sourcing" from the RStudio GUI
if (Sys.getenv("RSTUDIO") == "1" & !is.null(parent.frame(2)$ofile)) {
  this.dir <- dirname(parent.frame(2)$ofile)
  setwd(this.dir)
}

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
lb_para   <- 6
noise_N   <- 32
metab_mm  <- basis2mrs_data(full_basis, sum_elements = TRUE, amp = amps)
broad_sig <- sim_resonances(freq = 1.3, amp = 150, lw = 100, lg = 1)

mrs_data_nn    <- lb(metab_mm, lb_para) + broad_sig    # no noise data
mrs_data_noise <- sim_noise(sd = 2.0, fd = FALSE, dyns = noise_N)
mrs_data       <- rep_dyn(mrs_data_nn, noise_N) + mrs_data_noise

ed_pppm_start <- 1
ed_pppm_end   <- 25
ed_pppm_N     <- 10
ed_pppm_vec   <- 10 ^ (seq(log10(ed_pppm_start), log10(ed_pppm_end),
                       length.out = ed_pppm_N))

fname <- "../data/fig3.rds"        # precomputed results

if (file.exists(fname)) {  # don't recalc unless we have to
  cat("Reading precomputed results :", fname, "\n")
  res_list <- readRDS(fname) 
} else {
  res_list <- vector(mode = "list", length = ed_pppm_N)
  for (n in 1:ed_pppm_N) {
    
    
    opts  <- abfit_opts(auto_bl_flex = FALSE, bl_ed_pppm = ed_pppm_vec[n],
                        maxiters = 0, maxiters_pre = 0, init_damping = lb_para,
                        bl_comps_pppm = 40, pre_align = FALSE)
    
    res_list[[n]] <- fit_mrs(mrs_data, method = "abfit", opts = opts,
                             basis = full_basis, time = FALSE)
  }
  cat("Saving precomputed results :", fname, "\n")
  saveRDS(res_list, fname)
}

error_vec        <- rep(NA, ed_pppm_N)
sd_error_vec     <- rep(NA, ed_pppm_N)
aic_vec          <- rep(NA, ed_pppm_N)
sd_aic_vec       <- rep(NA, ed_pppm_N)
aic_mod_vec      <- rep(NA, ed_pppm_N)
sd_aic_mod_vec   <- rep(NA, ed_pppm_N)
smo_res_fit_list <- vector(mode = "list", length = ed_pppm_N)
resid_list       <- vector(mode = "list", length = ed_pppm_N)

# calc amp est errors and smooth residuals to check for modelling quality
for (n in 1:ed_pppm_N) {
  amp_inds <- c(6:22)
  temp_tab     <- res_list[[n]]$res_tab
  fit_amp_mat  <- temp_tab[amp_inds]
  true_amp_mat <- matrix(amps[-18], nrow(fit_amp_mat), ncol(fit_amp_mat),
                         byrow = TRUE)
  
  error           <- (true_amp_mat - fit_amp_mat) ^ 2
  mean_error      <- mean(rowSums(error))
  sd_error        <- sd(rowSums(error)) * 5 # inflate sds by x5  for clarity
  error_vec[n]    <- mean_error
  sd_error_vec[n] <- sd_error
  
  aic             <- log(temp_tab$res.deviance) + 2 * temp_tab$bl_ed_pppm * 
                     temp_tab$ppm_range / temp_tab$fit_pts
  
  aic_vec[n]      <- mean(aic)
  sd_aic_vec[n]   <- sd(aic)
  
  aic_mod         <- log(temp_tab$res.deviance) + 10 * temp_tab$bl_ed_pppm * 
                     temp_tab$ppm_range / temp_tab$fit_pts
  
  aic_mod_vec[n]    <- mean(aic_mod)
  sd_aic_mod_vec[n] <- sd(aic_mod)
}

df <- data.frame(ed_pppm_vec, error_vec, sd_error_vec, aic_vec, sd_aic_vec,
                 aic_mod_vec, sd_aic_mod_vec)

breaks <- c(1, 2, 3, 5, 10, 25)

p1 <- ggplot(data = df, aes(x = ed_pppm_vec)) + 
             geom_line(aes(y = error_vec)) + geom_point(aes(y = error_vec)) +
             geom_errorbar(aes(ymin = error_vec - sd_error_vec,
                           ymax = error_vec + sd_error_vec), width = .02) +
             scale_x_continuous(trans='log10', breaks = breaks) +
             xlab("Baseline ED per ppm") +  ylab("Metabolite estimate error")

aic_df <- rbind(data.frame(ed_pppm_vec, Value = aic_vec, group = "AIC"),
                data.frame(ed_pppm_vec, Value = aic_mod_vec,
                           group = "mAIC"))

p2 <- ggplot(data = aic_df, aes(x = ed_pppm_vec, y = Value, col = group)) +
      geom_line() + geom_point() + 
      scale_x_continuous(trans='log10', breaks = breaks) + 
      theme(legend.title = element_blank(), legend.position = c(0.2, 0.9)) + 
      xlab("Baseline ED per ppm") 

p3 <- function() {
  par(cex = 0.75)
  plot(res_list[[1]], restore_def_par = FALSE)
}

p4 <- function() {
  par(cex = 0.75)
  plot(res_list[[5]], restore_def_par = FALSE)
}

p5 <- function() {
  par(cex = 0.75)
  plot(res_list[[6]], restore_def_par = FALSE)
}

p6 <- function() {
  par(cex = 0.75)
  plot(res_list[[10]], restore_def_par = FALSE)
}

full_plot <- plot_grid(p1, p2, p3, p4, p5, p6,
                       labels = c('A', 'B', 'C', 'D', 'E', 'F'),
                       label_size = 12, rel_widths = c(1,1,1),
                       nrow = 2)

# print(full_plot)

cairo_pdf("../figures/fig3.pdf", width = 6.92, height = 5.0)
print(full_plot)
dev.off()

tiff("../figures/fig3.tiff", width = 300 * 6.92, height = 300 * 5.0,
     pointsize = 10, res = 300)
print(full_plot)
dev.off()