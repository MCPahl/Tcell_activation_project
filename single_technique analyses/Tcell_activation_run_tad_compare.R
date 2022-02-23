 library(TADCompare)
 library(data.table)
# Convert to sparse 3-column matrix using cooler2sparse from HiCcompare
dir = "/mnt/isilon/sfgi/pahlm/projects/Tcell_activation/data/tad_coords/"

cd4_unstim = fread(paste(dir, "CD4_unstim_10K.txt", sep="/"))
cd4_8hr = fread(paste(dir, "CD4_8hr_10K.txt", sep="/"))
cd4_24hr = fread(paste(dir, "CD4_24hr_10K.txt", sep="/"))

cd4_unstim_sparse_mat <- HiCcompare::cooler2sparse(cd4_unstim)
cd4_8hr_sparse_mat <- HiCcompare::cooler2sparse(cd4_8hr)
TD_Compare_8hr_unstim <- TADCompare(cd4_unstim_sparse_mat, cd4_8hr_sparse_mat, resolution = 10000)
save(TD_Compare_8hr_unstim, file="TD_Compare_8hr_unstim.RData")
