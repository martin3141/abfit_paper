# change the working directory to the source file location
# when "Sourcing" from the RStudio GUI
if (Sys.getenv("RSTUDIO") == "1" & !is.null(parent.frame(2)$ofile)) {
  this.dir <- dirname(parent.frame(2)$ofile)
  setwd(this.dir)
}

# install any necessary packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(spant, ggplot2, cowplot, doParallel, gridGraphics)

system("Rscript fig1.R")
system("Rscript fig2.R")
system("Rscript fig3.R")
system("Rscript fig4_S1.R")
system("Rscript fig5_S2.R")
system("Rscript fig6_S3.R")
system("Rscript fig7_S4.R")
system("Rscript fig8_S5.R")
system("Rscript fig9_S6.R")
system("Rscript fig10.R")
system("Rscript figS7.R")
system("Rscript figS8.R")
