x = read.csv("tables/drugs/drug_list_rnaseq_pubmed.csv")


z = z[complete.cases(z),]
p = z %>% select(gene_name, pubmed_count, gene_id, BF,gwas) %>% group_by(gene_name, pubmed_count, gene_id, gwas) %>% summarize(BF= max(BF)) %>%unique()

pdf("plots/drug_targets_literature_query_vs_BF.pdf", useDingbats=FALSE, width=20, height = 20)
ggplot(p, aes(x = log2(pubmed_count+1), y=log10(BF), label=gene_name))+
#geom_point()+
geom_text()+
facet_wrap(~gwas)+
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

q = z %>% select(gene_name, pubmed_count, gene_id, pp, gwas) %>% group_by(gene_name, pubmed_count, gene_id, gwas)%>% summarize(pp= max(pp)) %>% unique()

pdf("plots/drug_targets_literature_query_vs_pp.pdf", useDingbats=FALSE, width=20, height = 20)
ggplot(q, aes(x = log2(pubmed_count+1), y=log10(pp), label=gene_name))+
#geom_point()+
geom_text()+
facet_wrap(~gwas)+
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
