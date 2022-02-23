library("parseIbed")
library("DFbedtools")
library("data.table")

#Read Promoter annotation
data(prom_bed)
names(prom_bed) = c("pro_chr", "pro_start","pro_end", "pro_anno", "strand", "gene_id")
prom_bed$pro_anno = paste(prom_bed$pro_anno, prom_bed$gene_id, sep="+")


load("data/norm_counts/atac_fpkm.Rdata")
load("data/loop_calls/fithic_loops_3resolutions.Rdata")

ocrs = fpkm[,c(1:4)] %>% unique()
names(fithic_loops) = c("chr_a", "start_a", "end_a","chr_b", "start_b", "end_b", "n_reads", "fdr", "condition")

prom_ocr_hic_annotated = annotate_bedpe2geneOCR(fithic_loops, ocrs, prom_bed)

save(prom_ocr_hic_annotated, file="intersection/promoter_ocr_hic_annotated.RData")

