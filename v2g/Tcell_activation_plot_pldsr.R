library("tidyverse")
library("viridis")

clusters = read_csv("tables/Tcell_activation_kmeansCluster_geneList.csv")

cluster_top_genes = clusters %>% 
	group_by(cluster) %>%
	mutate(rank = rank(-score)) %>% 
	filter(rank <= 100) %>% 
	ungroup() %>%
	gather(key = sample, value = tpm, CD4_unstim_ND543, CD4_unstim_ND581, CD4_unstim_ND589, CD4_8hr_ND543, CD4_8hr_ND581, CD4_8hr_ND589, CD4_24hr_ND543, CD4_24hr_ND581, CD4_24hr_ND589) %>%
	mutate(timepoint = ifelse(grepl("unstim", sample), "unstim", ifelse(grepl("8hr",sample), "8hr","24hr"))) %>% 
	mutate(class = paste0(as.character(cluster),"_", rank)) %>% 
	group_by(rank, cluster, gene_name, gene_id, timepoint, class) %>%
	summarize(tpm = mean(tpm)) %>% 
	ungroup() %>% 
	group_by(gene_name, gene_id) %>%
	mutate(z = scale(tpm))


cluster_top_genes$gene_name = factor(cluster_top_genes$gene_name, levels= rev(unique(cluster_top_genes$gene_name[order(cluster_top_genes$class)])))
cluster_top_genes$timepoint = factor(cluster_top_genes$timepoint, levels= c("unstim", "8hr", "24hr"))

genes = c("RELB", "IL21R", "NFKBIA", "PER2", "NFKB1", "IL20RA", "BMPR1B", "RELN", "NPAS3", "CTNNBL", "MTOR", "HNRNPM","ZNF75D", "KLF2", "KLF7", "FOXO4", "CDK2", "TRAIP1", "AKT1", "FBXO4", "GAMT")


cluster_top_genes.selected = cluster_top_genes[cluster_top_genes$gene_name %in% genes,]
pdf("plots/core_cluster_genes_subset.pdf", useDingbats=FALSE)
ggplot(cluster_top_genes.selected , aes(x = timepoint, y = gene_name, fill = z))+
	geom_tile()+
	scale_fill_viridis(discrete=FALSE)+
	theme_minimal()
dev.off()
