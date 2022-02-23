######################################################################################################################################
library(tidyverse)
load("intersection/promoter_ocr_hic_annotated.RData")

#Get promoter-ocr connections
gene_ocr_linakges = unique(data.frame(
	anno = c(gsub("\\+.*\\+", "+", prom_ocr_hic_annotated$anno_a),gsub("\\+.*\\+", "+", prom_ocr_hic_annotated$anno_b)), 
	ocr = c(prom_ocr_hic_annotated$ocr_b,prom_ocr_hic_annotated$ocr_a),
	condition =  prom_ocr_hic_annotated$condition
))
gene_ocr_linakges = gene_ocr_linakges[complete.cases(gene_ocr_linakges),]

gene_ocr_linkages = gene_ocr_linakges %>% mutate(
	ocr_id = gsub(".*:", "", ocr), 
	chr = gsub(":.*", "", ocr), 
	start = gsub(":.*", "", gsub(":OCR.*", "", gsub("chr\\d*:", "", gsub("chrX:", "", ocr)))),
	end = gsub(".*:", "", gsub(":OCR.*", "", gsub("chr\\d*:", "", gsub("chrX:", "", ocr)))),
	gene_name = gsub("\\+.*", "", anno),
	gene_id = gsub(".*\\+", "", anno)
) %>%
select(-ocr, -anno) %>%
group_by(ocr_id, chr, start, end, gene_name, gene_id) %>%
summarize(hic_contact_condition = paste(condition, collapse=":"))

save(gene_ocr_linkages, file="intersection/gene_ocr_annotation.RData")
######################################################################################################################################

library(tidyverse)
load("intersection/gene_ocr_annotation.RData")
load("data/comparisions/rna_differential_analysis.Rdata")

rna_DE_wide = DE_df %>% select(gene_name, gene_id, logFC, FDR, comp) %>% 
	nest(logFC, FDR, .key='value_col') %>%
	spread(key = comp, value = value_col) %>% 
 	unnest(CD4_8hr_vs_unstim) %>% 
 	rename(RNAseq_CD4_8hr_vs_unstim_logFC = logFC, RNAseq_CD4_8hr_vs_unstim_FDR = FDR) %>% 
 	unnest(CD4_24hr_vs_8hr) %>% 
 	rename(RNAseq_CD4_24hr_vs_8hr_logFC = logFC, RNAseq_CD4_24hr_vs_8hr_FDR = FDR) %>%
 	unnest(CD4_24hr_vs_unstim) %>% 
 	rename(RNAseq_CD4_24hr_vs_unstim_logFC = logFC, RNAseq_CD4_24hr_vs_unstim_FDR = FDR) %>%


 gene_ocr_linkages_rna = left_join(gene_ocr_linkages, DE_wide)


 load("data/comparisions/atac_differential_analysis.Rdata")

atac_DE_wide = DE_df %>% select(id, logFC, FDR, comp) %>% 
	nest(logFC, FDR, .key='value_col') %>%
	spread(key = comp, value = value_col) %>% 
 	unnest(CD4_8hr_vs_unstim) %>% 
 	rename(ATACseq_CD4_8hr_vs_unstim_logFC = logFC, ATACseq_CD4_8hr_vs_unstim_FDR = FDR) %>% 
 	unnest(CD4_24hr_vs_8hr) %>% 
 	rename(ATACseq_CD4_24hr_vs_8hr_logFC = logFC, ATACseq_CD4_24hr_vs_8hr_FDR = FDR) %>%
 	unnest(CD4_24hr_vs_unstim) %>% 
 	rename(ATACseq_CD4_24hr_vs_unstim_logFC = logFC, ATACseq_CD4_24hr_vs_unstim_FDR = FDR) %>%
 	rename(ocr_id = id)

gene_ocr_linkages_rna_atac = left_join(gene_ocr_linkages_rna, atac_DE_wide)

gene_ocr_linkages_rna_atac = gene_ocr_linkages_rna_atac[is.na(gene_ocr_linkages_rna_atac$RNAseq_CD4_24hr_vs_8hr_logFC)==F,]
save(gene_ocr_linkages_rna_atac, file = "intersection/gene_ocr_linkages_differential.RData")
write.csv(gene_ocr_linkages_rna_atac, file = "tables/gene_ocr_linkages_differential.csv", quote=F, row.names=F)


rna.kcluster  = read.csv("tables/Tcell_activation_kmeansCluster_geneList.csv")

cluster_ocr_linkages = left_join(rna.kcluster[,1:3], gene_ocr_linkages_rna_atac )

save(cluster_ocr_linkages, file = "intersection/gene_ocr_linkages_differential_filtered_rnaseq_kmeans_clusters.RData")
write.csv(cluster_ocr_linkages, file = "tables/gene_ocr_linkages_differential_filtered_rnaseq_kmeans_clusters.csv", quote=F, row.names=F)