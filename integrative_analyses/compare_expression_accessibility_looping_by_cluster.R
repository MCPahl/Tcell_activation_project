#Retreive differential data from each, then get cluster means and plot.

library(tidyverse)

#load looping data
load("intersection/promoter_ocr_hic_annotated.RData")
prom_ocr_hic_annotated = prom_ocr_hic_annotated %>% filter(grepl("4kb", condition))

#Reformat to OCR-gene
anno_key = bind_rows(
	annoA_ocrB = prom_ocr_hic_annotated %>%
	 	select(chr_a, start_a, end_a, chr_b, start_b, end_b, anno_a, ocr_b) %>%
	 	dplyr::rename(anno = anno_a, ocr = ocr_b),
	annoA_ocrB = prom_ocr_hic_annotated %>%
	 	select(chr_a, start_a, end_a, chr_b, start_b, end_b, anno_b, ocr_a) %>%
	 	dplyr::rename(anno = anno_b, ocr = ocr_a)) %>% 
	mutate(gene_name = gsub("\\+.*", "", anno), gene_id = gsub(".*\\+", "", anno), ocr = gsub(".*:", "", ocr)) %>% select(-anno) %>% filter(complete.cases(.)) %>%
	unique()

#plot mean number of genes connected per ocr, plot mean ocr connected per gene
ocr_to_gene_count = anno_key%>% group_by(ocr) %>% summarize(n_gene = length(unique(gene_id)))
gene_to_ocr_count = anno_key%>% group_by(gene_id) %>% summarize(n_ocr = length(unique(ocr)))

ocr_to_gene_med = median(ocr_to_gene_count$n_gene)
gene_to_ocr_med = median(gene_to_ocr_count$n_ocr)


pdf("plots/allOCR_gene_pairs_median_number_genes_per_ocr.pdf", width = 5, height = 2)
ggplot(ocr_to_gene_count, aes(x = n_gene))+
	geom_histogram(binwidth=1, fill = "cyan")+
	theme_bw()+ 
	geom_vline(xintercept = ocr_to_gene_med, color = "red")
dev.off()

pdf("plots/allOCR_gene_pairs_median_number_ocr_per_gene.pdf", width = 5, height = 2)
ggplot(gene_to_ocr_count, aes(x = n_ocr))+
	geom_histogram(binwidth=1, fill = "orange")+
	theme_bw()+ 
	geom_vline(xintercept = gene_to_ocr_med, color = "red")
dev.off()



#Add RNA-seq differential
load("data/comparisions/rna_differential_analysis.Rdata")

DE_df = DE_df %>% select(gene_name, gene_id, logFC, FDR, comp) %>% filter(comp != "CD4_24hr_vs_unstim") %>% dplyr::rename(RNA_logFC = logFC, RNA_FDR = FDR)

anno_key = left_join(anno_key, DE_df) %>% filter(complete.cases(.))

#Add ATACSeq Differential
load("data/comparisions/atac_differential_analysis.Rdata")

DE_df = DE_df %>% select(id, logFC, FDR, comp) %>% filter(comp != "CD4_24hr_vs_unstim") %>% dplyr::rename(ATAC_logFC = logFC, ATAC_FDR = FDR) %>% dplyr::rename(ocr = id)

anno_key = left_join(anno_key, DE_df) %>% filter(complete.cases(.))

#Add HiC differential 
load("/mnt/isilon/sfgi/suc1/analyses/wells/hiC/compare/naiveT_stim/data/4000/hr24_vs_hr8.results.Rdata")

hic_result_1 = result %>% select(chr, region1, region2, logFC, p.adj) %>% 
	mutate(chr = paste0("chr", chr), comp = "CD4_24hr_vs_8hr") %>%
	dplyr::rename(chr_a = chr, start_a = region1, start_b = region2, HiC_logFC = logFC, HiC_FDR = p.adj)

load("/mnt/isilon/sfgi/suc1/analyses/wells/hiC/compare/naiveT_stim/data/4000/hr8_vs_unstim.results.Rdata")

hic_result_2 = result %>% select(chr, region1, region2, logFC, p.adj) %>% 
	mutate(chr = paste0("chr", chr), comp = "CD4_8hr_vs_unstim") %>%
	dplyr::rename(chr_a = chr, start_a = region1, start_b = region2, HiC_logFC = logFC, HiC_FDR = p.adj)

hic_result = bind_rows(hic_result_1, hic_result_2)

anno_key = left_join(anno_key, hic_result)  %>% filter(complete.cases(.))

kclusters = read.csv("tables/Tcell_activation_kmeansCluster_geneList.csv") %>% select(gene_name, gene_id, cluster)


anno_key_diff = anno_key %>% filter(ATAC_FDR < 0.05 & abs(ATAC_logFC) > 1 & HiC_FDR < 0.05 & abs(HiC_logFC) > 1)

save(anno_key, anno_key_diff, file = "intersection/gene_ocr_linkages_differental_annotated.RData")

kclusters_diff = inner_join(kclusters, anno_key) %>% filter(complete.cases(.))

x = kclusters_diff %>% 
	group_by(gene_name, gene_id, comp, cluster) %>% 
	summarize(RNA_logFC = median( RNA_logFC), ATAC_logFC = median(ATAC_logFC), HiC_logFC = median(HiC_logFC)) 

x$comp = factor(x$comp, levels =c("CD4_8hr_vs_unstim", "CD4_24hr_vs_8hr"))
pdf("plots/Tcell_compare_promoterConnectedOCR_ATAC_RNA_seq_logFC.pdf")
ggplot(x, aes(x = ATAC_logFC, y = RNA_logFC))+
	geom_point()+
	geom_smooth()+
	theme_bw()+
	facet_grid(~comp)
dev.off()

x$class = ifelse(x$RNA_logFC>0, "up", "down")



xt = cbind(x, t)

pdf("plots/Tcell_compare_promoterConnectedOCR_ATAC_RNA_seq_logFC.pdf")
ggplot(x, aes(x = ATAC_logFC, y = RNA_logFC, color=class, group = class))+
	#geom_point()+
	#geom_smooth()+
	theme_bw()+
	facet_grid(~comp+cluster)
dev.off()

x.t = x %>% group_by(comp) %>%
	 mutate(ATAC_logFC_group = factor(cut(ATAC_logFC,breaks=5, labels=FALSE)), HiC_logFC_group = factor(cut(HiC_logFC,breaks=5, labels=FALSE)))


pdf("plots/Tcell_compare_promoterConnectedOCR_ATAC_RNA_seq_logFC.pdf")
ggplot(x.t, aes(x =  HiC_logFC_group , y = RNA_logFC))+
	geom_boxplot()+
	#geom_point()+
	theme_bw()+
	facet_grid(vars(comp), vars(cluster), scales = "free")+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



pdf("plots/Tcell_compare_promoterConnectedOCR_RNA_HiC_seq_logFC.pdf")
ggplot(x.t, aes(x = HiC_logFC_group, y = RNA_logFC))+
	geom_boxplot()+
	#geom_point()+
	theme_bw()+
	facet_grid(~comp+cluster, scales = "free")+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("plots/Tcell_compare_promoterConnectedOCR_ATAC_HiC_seq_logFC.pdf")
ggplot(x.t, aes(x = ATAC_logFC_group, y = HiC_logFC))+
	geom_boxplot()+
	#geom_point()+
	theme_bw()+
	facet_grid(vars(comp), vars(cluster), scales = "free")+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()




pdf("plots/Tcell_compare_promoterConnectedOCR_ATAC_HiC_seq_logFC.pdf")
ggplot(x.t, aes(x = ATAC_logFC_group, y = HiC_logFC))+
	geom_boxplot()+
	#geom_point()+
	theme_bw()+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()




anno_key_diff_long = anno_key_diff %>% 
	select(gene_name, RNA_logFC, ATAC_logFC, HiC_logFC, comp) %>%
	gather(key = experiment, value = logFC, RNA_logFC, ATAC_logFC, HiC_logFC) %>%
	mutate(class = paste(comp, experiment, sep="_"))	

pdf("plots/Differential_loop_chromatin_contacts_expression_summary_heatmap.pdf", width=20, height=20)
ggplot(anno_key_diff_long, aes(x=class, y = gene_name, fill = logFC))+
	geom_tile()+
	scale_fill_gradient2(low="navy", mid="white", high="red")+
	theme_bw()
dev.off()




