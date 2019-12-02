#!/bin/bash

#SBATCH --mail-type BEGIN,END
#SBATCH --ntasks 32
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --time 00:30:00
#SBATCH --account wilsonmp-mrs-analysis
#SBATCH --qos castles

# mail-type options : NONE, BEGIN, END, FAIL
# qos options       : bbshort, bbdefault, castles

module purge
module load bluebear
module load bear-apps/2019a
module load spant/0.18.0-foss-2019a-R-3.6.0

Rscript fig4.R
Rscript fig5.R
Rscript fig6.R
Rscript fig7.R
Rscript fig8.R
