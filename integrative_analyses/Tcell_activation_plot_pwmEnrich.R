library(tidyverse)
library("PWMEnrich")
library(viridis)
library(ComplexHeatmap)

load("pwm_enrichment.RData")

anno = read.delim("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/jasper2020//TF_motif.anno.txt")
anno$id = NULL

anno = anno %>%
	dplyr::rename(id = motif_id)

out = lapply(pwm_enrichment_cluster_linked_ocrs, function(cluster_enrichment){
	groupSig <- as.data.frame(groupReport(cluster_enrichment))
	groupSig$fdr = p.adjust(groupSig$p.value,method="BH")
	groupSig = left_join(groupSig, anno)
})

names(out) = names(pwm_enrichment_cluster_linked_ocrs)

for(i in seq_along(out)){
	out[[i]]$pwm_enriched_cluster = names(out)[i]
}

out = do.call("rbind", out)
row.names(out) = NULL


clusters = read.csv("../../tables/Tcell_activation_kmeansCluster_geneList.csv")
clusters = clusters[,c(1:3)]
clusters$gene_id = gsub("\\..*", "", clusters$gene_id)

out = left_join(out, clusters)

out$cluster[is.na(out$cluster)] = "not_differential"


matr = out %>% select(TF, pwm_enriched_cluster, id, fdr) %>% mutate(fdr_t = -log10(fdr+2.225074e-308)) %>% select(-fdr) %>% unique() %>%
	spread(key = pwm_enriched_cluster, value = fdr_t) 

matr.m = matr[,-c(1:2)]
row.names(matr.m) = matr$TF



pdf("../../plots/tf_heatmap.pdf", useDingbats=F, width=50, height = 150)
Heatmap(as.matrix(matr.m),
  show_row_names=TRUE,
  col=viridis(10000),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  )
 dev.off()

#filtered

geneList = unique(out[out$fdr<0.01,]$TF)

matr = out %>% select(TF, pwm_enriched_cluster, id, fdr) %>% filter(TF %in% geneList) %>% mutate(fdr_t = -log10(fdr+2.225074e-308)) %>% select(-fdr) %>% unique() %>%
	spread(key = pwm_enriched_cluster, value = fdr_t) 

matr.m = matr[,-c(1:2)]
row.names(matr.m) = matr$TF


pdf("../../plots/tf_filtered.pdf", useDingbats=F, width=50, height = 150)
Heatmap(as.matrix(matr.m),
  show_row_names=TRUE,
  col=viridis(10000),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  )
 dev.off()

