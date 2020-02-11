# Script and data to reproduce the results presented in the ABfit paper

## Local installation instructions

### Dependencies
- R (all required packages will be automatically installed)
- git (optional but recommended)

### Install on Linux or OSX

```
git clone https://github.com/martin3141/abfit_paper.git
cd abfit_paper/scripts
Rscript run_all.R
```

### Install on other (or without git)

Download the repository zip using the green "Clone or download" button above. Extract the archive and run run_all.R (located in the scripts directory) with the R interpreter.

### Notes

All paper figures will be generated in the "figures" directory once rendered.

Analysis will take some time, however fitting can be performed in parallel to speed things up. Edit the "jobs" variable in the figX.R scripts to match the number of available CPU cores (default = 4). All analyses can be run under 15 mins on a machine with 32 cores.
