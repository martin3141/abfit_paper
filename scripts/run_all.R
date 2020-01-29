# change the working directory to the source file location
# when "Sourcing" from the RStudio GUI
if (Sys.getenv("RSTUDIO") == "1" & !is.null(parent.frame(2)$ofile)) {
  this.dir <- dirname(parent.frame(2)$ofile)
  setwd(this.dir)
}

system("Rscript fig1.R")
system("Rscript fig2.R")
system("Rscript fig3.R")
system("Rscript fig4.R")
system("Rscript fig5.R")
system("Rscript fig6.R")
system("Rscript fig7.R")
system("Rscript fig8.R")