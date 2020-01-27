#!/bin/bash

# batch script to run on bluebear cluster

#SBATCH --mail-type BEGIN,END
#SBATCH --ntasks 32
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --time 01:00:00
#SBATCH --account wilsonmp-mrs-analysis
#SBATCH --qos castles

# mail-type options : NONE, BEGIN, END, FAIL
# qos options       : bbshort, bbdefault, castles

module purge
module load bluebear
module load bear-apps/2019a
module load spant/0.19.0-foss-2019a-R-3.6.0
module load gridGraphics/0.4-1-foss-2019a-R-3.6.0

Rscript run_all.R
