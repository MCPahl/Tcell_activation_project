library(tidyverse)

load("/mnt/isilon/sfgi/pahlm/projects/Tcell_activation/intersection/gene_ocr_linkages_differential_filtered_rnaseq_kmeans_clusters.RData")

load("intersection/summarized_differential_conditionMean_normalized_zscores.Rdata")
cluster_ocr_linkages = cluster_ocr_linkages[complete.cases(cluster_ocr_linkages),]

head(anno_IF)

cluster_ocr_linkages = cluster_ocr_linkages %>%
	 select(gene_name, gene_id, cluster, ocr_id) %>%
	 rename(ocr = ocr_id) %>% 
	 unique()


dat = left_join(cluster_ocr_linkages, anno_IF) 



anno_IF_filt = dat %>% 
	filter(complete.cases(.)) %>%
    select(cluster, chr_a, start_a, end_a, chr_b, start_b, end_b, ocr, gene_name, gene_id, HiC_CD4_unstim,  HiC_CD4_8hr, HiC_CD4_24hr, ATAC_CD4_unstim, ATAC_CD4_8hr, ATAC_CD4_24hr  ,RNA_CD4_unstim,  RNA_CD4_8hr,  RNA_CD4_24hr) %>%
    mutate(anno = paste0(gene_name, ":", ocr)) %>% 
    select(anno, cluster, HiC_CD4_unstim, HiC_CD4_8hr, HiC_CD4_24hr, ATAC_CD4_unstim,ATAC_CD4_8hr, ATAC_CD4_24hr, RNA_CD4_unstim, RNA_CD4_8hr, RNA_CD4_24hr) %>% 
    gather(key = class, value = zscore,-anno, -cluster) %>% 
    group_by(anno, class, cluster) %>%
    summarize(zscore = min(zscore)) %>%
    ungroup() %>% 
    mutate(value_type = gsub("_.*", "", class)) %>% 
    mutate(class = factor(class, levels = c("RNA_CD4_unstim", "RNA_CD4_8hr", "RNA_CD4_24hr", 
        "ATAC_CD4_unstim", "ATAC_CD4_8hr", "ATAC_CD4_24hr",
        "HiC_CD4_unstim", "HiC_CD4_8hr", "HiC_CD4_24hr")))%>%
    mutate(value_type = factor(value_type, levels = c("RNA", "ATAC", "HiC"))) %>%
    group_by(anno, cluster, value_type) %>% 
    mutate(zscore_perc = (zscore - min(zscore))/(max(zscore)-min(zscore))) %>%
    ungroup()


library(viridis)

pdf("plots/top_differential_zscore.pdf", width = 15, height = 15)
ggplot(anno_IF_filt, aes(x = class, y =anno, fill = zscore_perc))+
geom_tile()+
scale_fill_viridis(discrete=FALSE)+
theme_bw()+
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
 facet_wrap(~cluster+value_type, scale="free", nrow=5)+
 theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
dev.off()



