 library(TADCompare)
 library(data.table)
# Convert to sparse 3-column matrix using cooler2sparse from HiCcompare
dir = "/mnt/isilon/sfgi/pahlm/projects/Tcell_activation/data/tad_coords/"

cd4_unstim = fread(paste(dir, "CD4_unstim_10K.txt", sep="/"))
cd4_8hr = fread(paste(dir, "CD4_8hr_10K.txt", sep="/"))
cd4_24hr = fread(paste(dir, "CD4_24hr_10K.txt", sep="/"))

cd4_unstim_sparse_mat <- HiCcompare::cooler2sparse(cd4_unstim)
cd4_8hr_sparse_mat <- HiCcompare::cooler2sparse(cd4_8hr)
cd4_24hr_sparse_mat <- HiCcompare::cooler2sparse(cd4_24hr)


TD_Compare_8hr_unstim <- TADCompare(cd4_unstim_sparse_mat, cd4_8hr_sparse_mat, resolution = 10000)
TD_Compare_2 <- TADCompare(cd4_8hr_sparse_mat, cd4_24hr_sparse_mat, resolution = 10000)

save(TD_Compare_8hr_unstim, file="TD_Compare_8hr_unstim.RData")

#!/bin/bash
source ~/.bashrc
conda activate py37
 cooler dump --join /mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_unstimulated_2reps/cool/merge/naiveT_unstimulated_2reps.10K.cool > CD4_unstim_10K.txt
 cooler dump --join /mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_8hr_2reps/cool/merge/naiveT_8hr_2reps.10K.cool > CD4_8hr_10K.txt
 cooler dump --join /mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_24hr_2reps/cool/merge/naiveT_24hr_2reps.10K.cool  > CD4_24hr_10K.txt

# Read in data
cool_mat1 <- read.table("Zuin.HEK293.50kb.Control.txt")
cool_mat2 <- read.table("Zuin.HEK293.50kb.Depleted.txt")

# Convert to sparse 3-column matrix using cooler2sparse from HiCcompare
sparse_mat1 <- HiCcompare::cooler2sparse(cool_mat1)
sparse_mat2 <- HiCcompare::cooler2sparse(cool_mat2)
