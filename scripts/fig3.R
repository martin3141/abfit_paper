library(spant)
library(ggplot2)
library(cowplot)

theme_set(theme_cowplot(font_size = 10))

# change the working directory to the source file location
# when "Sourcing" from the RStudio GUI
if (Sys.getenv("RSTUDIO") == "1") {
  this.dir <- dirname(parent.frame(2)$ofile)
  setwd(this.dir)
}

source("~/abfit/abfit.R")

ft <- def_acq_paras()$ft

ala  <- get_mol_paras("ala")
asp  <- get_mol_paras("asp")
cr   <- get_mol_paras("cr")
gaba <- get_mol_paras("gaba")
glc  <- get_mol_paras("glc")
gln  <- get_mol_paras("gln")
glu  <- get_mol_paras("glu")
gpc  <- get_mol_paras("gpc")
gsh  <- get_mol_paras("gsh")
ins  <- get_mol_paras("ins")
lac  <- get_mol_paras("lac")
mm   <- get_mol_paras("mm_3t", ft)
naa  <- get_mol_paras("naa")
naag <- get_mol_paras("naag")
pch  <- get_mol_paras("pch")
pcr  <- get_mol_paras("pcr")
sins <- get_mol_paras("sins")
tau  <- get_mol_paras("tau")

basis_list <- list(ala, asp, cr, gaba, glc, gln, glu, gpc, gsh, ins, lac, mm,
                   naa, naag, pch, pcr, sins, tau)

full_basis <- sim_basis(basis_list, pul_seq = seq_slaser_ideal,
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
          30.00,  # 12 MM
          12.25,  # 13 NAA
           1.50,  # 14 NAAG
           0.60,  # 15 PCh
           4.25,  # 16 PCr
           0.35,  # 17 sIns
           4.00)  # 18 Tau

set.seed(1)

# simulate mrs data
lb_para   <- 6
noise_N   <- 32
metab_mm  <- basis2mrs_data(full_basis, sum_elements = TRUE, amp = amps)
broad_sig <- sim_resonances(freq = 1.3, amp = 150, lw = 100, lg = 1,
                                phase = 0)
mrs_data_nn    <- lb(metab_mm, lb_para) + broad_sig    # no noise data
mrs_data_noise <- sim_noise(sd = 2.0, fd = FALSE, dyns = noise_N)
mrs_data       <-  rep_dyn(mrs_data_nn, noise_N) + mrs_data_noise

ed_pppm_start <- 2
ed_pppm_end   <- 25
ed_pppm_N     <- 10
ed_pppm_vec   <- 10 ^ (seq(log10(ed_pppm_start), log10(ed_pppm_end),
                          length.out = ed_pppm_N))

fname <- "fig3.rds"        # precomputed results

if (file.exists(fname)) {  # don't recalc unless we have to
  cat("Reading precomputed results :", fname, "\n")
  res_list <- readRDS(fname) 
} else {
  res_list <- vector(mode = "list", length = ed_pppm_N)
  for (n in 1:ed_pppm_N) {
    
    
    opts  <- abfit_opts(auto_bl_flex = FALSE, ed_pppm = ed_pppm_vec[n],
                        maxiters = 0, maxiters_pre = 0, init_damping = lb_para,
                        bl_comps_pppm = 40)
    
    res_list[[n]] <- fit_mrs(mrs_data, method = "abfit", opts = opts,
                             basis = full_basis)
  }
  cat("Saving precomputed results :", fname, "\n")
  saveRDS(res_list, fname)
}

error_vec        <- rep(NA, ed_pppm_N)
sd_error_vec     <- rep(NA, ed_pppm_N)
smo_res_fit_list <- vector(mode = "list", length = ed_pppm_N)
resid_list       <- vector(mode = "list", length = ed_pppm_N)

comps     <- round(40 * (4 - 0.5))
bl_bas    <- bbase(458, comps)
bl_comps  <- comps + 3
deriv_mat <- diff(diag(bl_comps), lag = 1, differences = 2)

# calc amp est errors and smooth residuals to check for modelling quality
for (n in 1:ed_pppm_N) {
  amp_inds <- c(6:23)
  amp_inds <- amp_inds[-12]               # remove MM
  fit_amp_mat  <- res_list[[n]]$res_tab[amp_inds]
  true_amp_mat <- matrix(amps[-12], nrow(fit_amp_mat), ncol(fit_amp_mat),
                         byrow = TRUE)
  
  error           <- (true_amp_mat - fit_amp_mat) ^ 2
  mean_error      <- mean(rowSums(error))
  sd_error        <- sd(rowSums(error))
  error_vec[n]    <- mean_error
  sd_error_vec[n] <- sd_error
  
 
  # auto-smooth residual to find a good baseline freedom estimate
  fit_tab               <- res_list[[n]]$fits$`1.fit`
  resid_list[[n]]       <- fit_tab$Data - fit_tab$Fit - fit_tab$Baseline
  smo_res_fit_list[[n]] <- auto_pspline_smoother(resid_list[[n]], bl_bas,
                                                 deriv_mat)
}

df <- data.frame(ed_pppm_vec, error_vec, sd_error_vec)

breaks <- c(2, 3, 5, 10, 25)

p1 <- ggplot(data = df, aes(x = ed_pppm_vec)) + 
             geom_line(aes(y = error_vec)) + geom_point(aes(y = error_vec)) +
             geom_errorbar(aes(ymin = error_vec - sd_error_vec,
                           ymax = error_vec + sd_error_vec), width = .02) +
             scale_x_continuous(trans='log10', breaks = breaks) + xlab("ED per PPM") + 
             ylab("Metabolite estimate error")

resid1   <- resid_list[[1]]
smo1     <- smo_res_fit_list[[1]]$yhat
resid2   <- resid_list[[2]]
smo2     <- smo_res_fit_list[[2]]$yhat
resid3   <- resid_list[[4]]
smo3     <- smo_res_fit_list[[4]]$yhat

lab1 <- "fit ED = 7, fit residual ED = 25"
lab2 <- "fit ED = 9, fit residual ED = 18"
lab3 <- "fit ED = 16, fit residual ED = 2"

ppm_sc <- res_list[[1]]$fits$`1.fit`$PPMScale
resid_smo_df <- rbind(data.frame(x = ppm_sc, y = resid1, c1 = "fit residual",
                                 ED = lab1),
                      data.frame(x = ppm_sc, y = smo1,   c1 = "smoothed",
                                 ED = lab1),
                      data.frame(x = ppm_sc, y = resid2, c1 = "fit residual",
                                 ED = lab2),
                      data.frame(x = ppm_sc, y = smo2,   c1 = "smoothed",
                                 ED = lab2),
                      data.frame(x = ppm_sc, y = resid3, c1 = "fit residual",
                                 ED = lab3),
                      data.frame(x = ppm_sc, y = smo3,   c1 = "smoothed",
                                 ED = lab3))
                     
p3 <- ggplot(resid_smo_df, aes(x, y, col = c1)) + geom_line() +
             facet_wrap(~ED, nrow = 3) +
             scale_x_reverse(breaks = c(1, 1.5, 2, 2.5, 3, 3.5),
                             expand = c(0,0)) + 
             scale_colour_manual(values = c("black", "red")) +
             xlab("Chemical shift (ppm)") + ylab("Intensity (au)") +
             theme(legend.title = element_blank(), legend.position="top",
                   axis.text.y = element_blank(),
                   axis.title.y = element_blank(),
                   legend.box.margin = margin(-2,0,-4,10),
                   legend.margin = margin(-2,0,-4,10)) 

res_ed_vec <- as.numeric(lapply(smo_res_fit_list,`[[`,"ED"))
df_ed      <- data.frame(ed_pppm_vec, res_ed_vec)
p2 <- ggplot(df_ed, aes(x = ed_pppm_vec, y = res_ed_vec)) + geom_line() +
      geom_point() + scale_x_continuous(trans='log10', breaks = breaks) +
      xlab("ED per PPM") + ylab("Fit residual ED")

p4 <- function() {
  par(cex = 0.75)
  plot(res_list[[1]], restore_def_par = FALSE)
}

p5 <- function() {
  par(cex = 0.75)
  plot(res_list[[2]], restore_def_par = FALSE)
}

p6 <- function() {
  par(cex = 0.75)
  plot(res_list[[4]], restore_def_par = FALSE)
}

full_plot <- plot_grid(p1, p2, p3, p4, p5, p6,
                       labels = c('A', 'B', 'C', 'D', 'E', 'F'),
                       label_size = 12, rel_widths = c(0.9,0.9,1,1,1,1), ncol = 3)

print(full_plot)

cairo_pdf("fig3.pdf", width = 6.92, height = 5.5)
print(full_plot)
dev.off()