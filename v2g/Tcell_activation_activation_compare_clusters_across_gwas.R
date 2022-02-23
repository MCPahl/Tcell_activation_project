library(tidyverse)
v2g = read.delim("tables/variant_to_gene_mapping_allTimepoints.txt")
traits = c("ALG", "CEL", "ECZ", "JIA", "PSO", "RA", "SLE", "T1D", "UC")

v2g_count_class = v2g %>% 
	select(gene_name, gene_id, rsid, gwas, cluster, hla_region) %>%
	filter(hla_region == FALSE) %>%
	filter(gwas %in% traits) %>% 
	unique() %>%
	group_by(gwas, cluster) %>%
	tally() %>%
	data.frame()

v2g_count_class$gwas = factor(v2g_count_class$gwas, levels = rev(unique(v2g_count_class$gwas)))

pdf("plots/v2g_counts.pdf")
ggplot(v2g_count_class, aes(x=gwas, y = n, fill = cluster))+
geom_bar(position="stack", stat="identity")+
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
coord_flip()
dev.off()









library(ComplexHeatmap)

v2g_count_genes = v2g %>% 
	select(gene_name, gene_id, rsid, gwas, hla_region) %>%
	filter(hla_region == FALSE) %>%
	filter(gwas %in% traits) %>% 
	unique()

v2g_count_genes_groups = split(v2g_count_genes, v2g_count_genes$gwas)

v2g_count_genes_groups_upset = lapply(v2g_count_genes_groups , function(x){
	unique(paste(x$gene_name, x$gene_id, sep="|")
	})

names(v2g_count_genes_groups_upset) = names(v2g_count_genes_groups)

gwas.matrix = list_to_matrix(v2g_count_genes_groups_upset)
 gwas.matrix1 = as.data.frame(gwas.matrix)
 gwas.matrix1$sum = unname(c(apply(gwas.matrix1, 1, sum)))

write.csv(gwas.matrix1, file = "tables/v2g_implicated_by_gwas.csv", quote=F, row.names=F)
gwas.matrix_comb = make_comb_mat(gwas.matrix)

pdf("plots/Upset_GWAS_genes.pdf", width=12, height=4,)
UpSet(gwas.matrix_comb,
	#set_order = c("AaM","Bipolar"BMI", "HEIGHT", "MDD", "SLEEP"),
	pt_size = unit(5, "mm"), lwd = 1,
	comb_col = c("red", "blue", "black", "darkgreen", "orange", "purple")[comb_degree(gwas.matrix)],
	comb_order = order(comb_size(gwas.matrix)),
	#set_on_rows
	)
dev.off()
