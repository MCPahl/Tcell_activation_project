library(tidyverse)

load("/mnt/isilon/sfgi/pahlm/annotationFiles/gtex/v8/gtex_median_tissue_expression_tpm.RData")
druggable_genes = read.delim("tables/drugs/tcell_activation_v2g_all_drug_gene_results_table.txt")

druggable_genes = druggable_genes[druggable_genes$DistinctDrugCount > 0,]

load("data/norm_counts/rnaseq_tpm.Rdata")

tpm_mean = tpm %>% group_by(gene_name, gene_id, condition) %>% 
summarize(tpm_mean = mean(tpm), tpm_sd = sd(tpm))


 gtex.tissue.med[gtex.tissue.med$gene_name %in% druggable_genes$Gene,]

 gtex.filtered =  gtex.tissue.med[gtex.tissue.med$gene_name %in% druggable_genes$Gene,]

gtex.filtered$score = apply(gtex.filtered[,c(50,56)], 1, sum)/ apply(gtex.filtered[,-c(1,2,50,56)], 1, sum)
gtex.filtered = gtex.filtered[,c(1:2,57)]

library("RISmed")

genes_druggable = read.delim("tables/drugs/tcell_activation_v2g_all_drug_gene_results_table.txt")
genes_druggable = genes_druggable[genes_druggable$DistinctDrugCount>0,]

genes = genes_druggable$Gene

terms = c("autoimmune", "allergy", "celiac", "eczema", "ibd", "irritable bowel", "juvenile idiopathic arthritis", "psoriasis", "rheumatoid arthritis",
"sle", "lupus", "t1d", "Type 1 Diabetes", "uc", "ulcerative colitis")

terms = paste(terms, collapse=" OR ")
terms = gsub("$", ")", gsub("^", "(", terms))

x=list()
for(i in seq_along(genes)){
search_topic <- paste(genes[[i]], terms, sep = " AND ")
search_query <- EUtilsSummary(search_topic, retmax=36000)
x[[i]] <- summary(search_query)
Sys.sleep(5)
}

lit.df = data.frame(gene_name = genes, pubmed_count = unlist(lapply(x,length)))


comp_gtex = inner_join(lit.df, gtex.filtered )

write.csv(comp_gtex, file ="tables/drugs/drug_list_gtex_pubmed.csv", quote=F, row.names=F)


pdf("plots/drug_targets_literature_query.pdf", useDingbats=FALSE, width=20, height = 20)
ggplot(comp_gtex, aes(x = log2(pubmed_count+1), y=score, label=gene_name))+
geom_point()+
geom_text()
dev.off()


comp_ourdat = left_join(lit.df, tpm_mean)
comp_ourdat = comp_ourdat[complete.cases(comp_ourdat),]

comp_ourdat$condition = factor(comp_ourdat$condition, levels = c("CD4_unstim", "CD4_8hr", "CD4_24hr"))
write.csv(comp_ourdat, file ="tables/drugs/drug_list_rnaseq_pubmed.csv", quote=F, row.names=F)

pdf("plots/drug_targets_literature_query_tpmOurData.pdf", useDingbats=FALSE, width=20, height = 20)
ggplot(comp_ourdat, aes(x = log2(pubmed_count+1), y=log2(tpm_mean+1), label=gene_name))+
#geom_point()+
geom_text()+
facet_grid(~condition)+
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

pubmed = x 
save(pubmed, file = "pubmedlookup.Rdata")


load("data/comparisions/rna_differential_analysis.Rdata")


sigDE_df_wide = sigDE_df %>% select(gene_name, gene_id, logFC, comp) %>% spread(key = comp, value = logFC)

lit.df_de = left_join(lit.df, sigDE_df_wide )


lit.df_de$CD4_24hr_vs_8hr[is.na(lit.df_de$CD4_24hr_vs_8hr)] = 0
lit.df_de$CD4_24hr_vs_unstim[is.na(lit.df_de$CD4_24hr_vs_unstim)] = 0
lit.df_de$CD4_8hr_vs_unstim[is.na(lit.df_de$CD4_8hr_vs_unstim)] = 0

lit.df_de = lit.df_de %>% gather(key=comp, value=logFC, CD4_8hr_vs_unstim, CD4_24hr_vs_8hr, CD4_24hr_vs_unstim)
pdf("plots/drug_targets_literature_query_differential_expression_OurData.pdf", useDingbats=FALSE, width=20, height = 20)
ggplot(lit.df_de, aes(x = log2(pubmed_count+1), y=logFC, label=gene_name))+
#geom_point()+
geom_text()+
facet_grid(~comp)+
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()


