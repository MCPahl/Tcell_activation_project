
library(tidyverse)

v2g = read.delim("tables/variant_to_gene_mapping_allTimepoints.txt")
traits = c("ALG", "CEL", "ECZ", "IBD", "JIA", "PSO", "RA", "SLE", "T1D", "UC")

v2g = v2g[v2g$gwas %in% traits & v2g$hla_region==FALSE,]

tally = v2g %>% 
	group_by(rsid, gwas) %>%
	summarize(BF_max = max(BF), gene_tally = length(unique(gene_id))) %>%
	unique()

pdf("plots/gwas_geneNumber_vs_BayesFactor.pdf", useDingbats=FALSE)
ggplot(tally, aes(y = log10(BF_max),x=gene_tally))+
 geom_point()+
 xlim(0,100)+
 facet_wrap(vars(gwas), nrow=2)+ 
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()



proxy_tally = v2g %>% 
	group_by(rsid, gwas) %>%
	summarize(pp = pp, gene_tally = length(unique(gene_id)))%>%
	unique()

pdf("plots/proxy_vs_pp.pdf", useDingbats=FALSE)
ggplot(proxy_tally, aes(x = gene_tally , y=log10(pp)))+
 geom_point()+
 #xlim(0,100)+
 facet_wrap(vars(gwas), nrow=2)+ 
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()



s = v2g %>% 
	group_by(rsid, gwas) %>%
	summarize(BFxpp= BF*pp, gene_tally = length(unique(gene_id)))%>%
	unique()

pdf("plots/proxy_vs_ppBF.pdf", useDingbats=FALSE)
ggplot(proxy_tally, aes(x = gene_tally , y=log10(BFxpp)))+
 geom_point()+
 #xlim(0,100)+
 facet_wrap(vars(gwas), nrow=2)+ 
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()


#!/bin/bash
source ~/.bashrc
cd /mnt/isilon/sfgi/pahlm/projects/Tcell_activation/tf/hint_atac/hint_out
rgt-motifanalysis matching --organism=hg19 \
	--input-files /mnt/isilon/sfgi/pahlm/projects/Tcell_activation/tf/hint_atac/hint_out/CD4_24hr.bed \
	--motif-dbs ~/rgtdata/motifs/jasparcore_vert2020/ \
	--pseudocounts 0.8 
rgt-motifanalysis matching --organism=hg19 \
	--input-files /mnt/isilon/sfgi/pahlm/projects/Tcell_activation/tf/hint_atac/hint_out/CD4_Unstim.bed \
	--motif-dbs ~/rgtdata/motifs/jasparcore_vert2020/ \
	--pseudocounts 0.8 
rgt-motifanalysis matching --organism=hg19 \
	--input-files /mnt/isilon/sfgi/pahlm/projects/Tcell_activation/tf/hint_atac/hint_out/CD4_8hr.bed \
	--motif-dbs ~/rgtdata/motifs/jasparcore_vert2020/ \
	--pseudocounts 0.8 
