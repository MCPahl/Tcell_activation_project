library(tidyverse)
library(viridis)
load("intersection/promoter_ocr_hic_annotated.RData")

load("data/comparisions/atac_differential_analysis.Rdata")
atac_de = DE_df
load("data/comparisions/hic_differential_analysis.Rdata")
load("data/comparisions/rna_differential_analysis.Rdata")
rna_de = DE_df



prom_ocr_hic_annotated1 = prom_ocr_hic_annotated %>%
mutate(gene_name_a = gsub("\\+.*", "", anno_a), gene_id_a = gsub(".*\\+", "", anno_a), gene_name_b = gsub("\\+.*", "", anno_b), gene_id_b = gsub(".*\\+", "", anno_b)) %>%
select(-c(anno_a, anno_b)) %>%
unique()

sigDE_consensus_loop = sigDE_consensus_loop %>% 
select(chr_a, start_a, end_a, chr_b, start_b, end_b, resol, comp, logFC, p.adj) %>%
rename(hic_comp = comp, hic_logFC = logFC, hic_p.adj = p.adj)

atac_de = atac_de %>% 
select(id, chr, start, end, logFC, FDR, comp) %>% 
mutate(ocr_id = paste(chr, start, end, id, sep=":")) %>% 
select(ocr_id, logFC, FDR, comp)

rna_de = rna_de %>% 
select(gene_name, gene_id,logFC, FDR, comp)

prom_ocr_hic_annotated2 = left_join(prom_ocr_hic_annotated1, sigDE_consensus_loop )

atac_de_a = atac_de %>%
rename(ocr_a = ocr_id, ocr_a_logFC = logFC, ocr_a_FDR = FDR, ocr_comp_a =comp) 

atac_de_b = atac_de %>%
rename(ocr_b = ocr_id, ocr_b_logFC = logFC, ocr_b_FDR = FDR, ocr_comp_b=comp) 

prom_ocr_hic_annotated3 = prom_ocr_hic_annotated2 %>% left_join(atac_de_a) %>% left_join(atac_de_b)

rna_de_a = rna_de %>%
	rename(gene_name_a = gene_name, gene_id_a = gene_id, rna_logFC_a = logFC, rna_FDR_a = FDR, rna_comp_a = comp)

rna_de_b = rna_de %>%
	rename(gene_name_b = gene_name, gene_id_b = gene_id, rna_logFC_b = logFC, rna_FDR_b = FDR, rna_comp_b = comp)

prom_ocr_hic_annotated4 = prom_ocr_hic_annotated3 %>% left_join(rna_de_a) %>% left_join(rna_de_b)

prom_ocr_hic_annotated5 = prom_ocr_hic_annotated4 %>% 
	filter(
		(rna_comp_a == rna_comp_b | is.na(rna_comp_a) | is.na(rna_comp_b) ) &
		(ocr_comp_a == ocr_comp_b | is.na(ocr_comp_a) | is.na(ocr_comp_b) ) & 
		(ocr_comp_b == rna_comp_a | ocr_comp_a == rna_comp_b)) %>%
	unique()

prom_ocr_hic_annotated6 = data.frame(
	gene_name = c(prom_ocr_hic_annotated5$gene_name_a, prom_ocr_hic_annotated5$gene_name_b), 
	gene_id = c(prom_ocr_hic_annotated5$gene_id_a, prom_ocr_hic_annotated5$gene_id_b), 
	ocr_id = c(prom_ocr_hic_annotated5$ocr_a, prom_ocr_hic_annotated5$ocr_b), 
	rna_comp = c(prom_ocr_hic_annotated5$rna_comp_a, prom_ocr_hic_annotated5$rna_comp_b), 
	ocr_comp = c(prom_ocr_hic_annotated5$ocr_comp_a, prom_ocr_hic_annotated5$ocr_comp_b), 
	rna_logFC = c(prom_ocr_hic_annotated5$rna_logFC_a, prom_ocr_hic_annotated5$rna_logFC_b),
	rna_FDR = c(prom_ocr_hic_annotated5$na_FDR_a, prom_ocr_hic_annotated5$rna_FDR_b),
	ocr_logFC = c(prom_ocr_hic_annotated5$ocr_a_logFC, prom_ocr_hic_annotated5$ocr_b_logFC),
	ocr_FDR = c(prom_ocr_hic_annotated5$ocr_a_FDR, prom_ocr_hic_annotated5$ocr_b_FDR), 
	hic_comp = c(prom_ocr_hic_annotated5$hic_comp, prom_ocr_hic_annotated5$hic_comp),
	hic_logFC = c(prom_ocr_hic_annotated5$hic_logFC, prom_ocr_hic_annotated5$hic_logFC),
	hic_p.adj = c(prom_ocr_hic_annotated5$hic_p.adj, prom_ocr_hic_annotated5$hic_p.adj)
) 

prom_ocr_hic_annotated7 = prom_ocr_hic_annotated6 %>% filter( (is.na(gene_id) | is.na(ocr_id)) == F) %>% unique() %>% mutate(hic_sig = ifelse(is.na(hic_p.adj)==F & hic_logFC >0, "increased_contact", ifelse(is.na(hic_p.adj)==F &hic_logFC<0, "decreased_contact", "stable_contact")))
prom_ocr_hic_annotated7 = prom_ocr_hic_annotated7[is.na(prom_ocr_hic_annotated7$ocr_comp)==F,]
prom_ocr_hic_annotated7 = prom_ocr_hic_annotated7[is.na(prom_ocr_hic_annotated7$rna_comp)==F,]
prom_ocr_hic_annotated7 = prom_ocr_hic_annotated7 %>% filter(ocr_comp %in% c("CD4_8hr_vs_unstim", "CD4_24hr_vs_8hr")) %>% mutate(ocr_comp = factor(ocr_comp, levels=c("CD4_8hr_vs_unstim", "CD4_24hr_vs_8hr") ))



pdf("plots/Tcell_atac_hic_boxplot.pdf")
ggplot(prom_ocr_hic_annotated7, aes(x = ocr_comp, y= ocr_logFC, fill = hic_sig))+
geom_boxplot()
dev.off()



pdf("plots/Tcell_atac_rna_seq_corr.pdf")
ggplot(prom_ocr_hic_annotated7, aes(x= ocr_logFC, y = rna_logFC))+
geom_density_2d_filled(contour_var = "ndensity") +
#scale_fill_viridis()+
facet_grid(hic_sig ~ ocr_comp, scales="free")+
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
xlim(-7.5,7.5)+
ylim(-7.5,7.5)+
geom_vline(xintercept=0, color = "red",linetype=3)+
geom_hline(yintercept=0, color = "red", linetype=3)
dev.off()


write.table(prom_ocr_hic_annotated7, file = "tables/Tcell_atac_rna_seq_corr_hic.txt", quote=F, row.names=F)





