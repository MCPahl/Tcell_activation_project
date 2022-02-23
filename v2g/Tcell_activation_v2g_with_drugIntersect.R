library(tidyverse)
library(GenomicRanges)

load("intersection/gene_ocr_annotation.RData") #gene_ocr_linkages

snp_dir = "data/gwas/credsets/"
snp_files = list.files(snp_dir)

autoimmune_credset = lapply(snp_files, function(snp_file){
	load(paste0(snp_dir, snp_file))
	for(i in seq_along(sub_sumstat)){
		sentinel_name = names(sub_sumstat)[i]
		sub_sumstat[[i]]$sentinel = sentinel_name 
		sub_sumstat[[i]]$chr = paste("chr", sub_sumstat[[i]]$chr, sep="")
		sub_sumstat[[i]] = sub_sumstat[[i]][order(sub_sumstat[[i]]$pp, decreasing=TRUE),]
		sub_sumstat[[i]]$pp_cumsum = cumsum(sub_sumstat[[i]]$pp)
		under_thresh=sub_sumstat[[i]]$pp_cumsum < 0.95
		cred_set_index=(length(under_thresh[under_thresh==TRUE])+1)
		sub_sumstat[[i]] = sub_sumstat[[i]][c(1:cred_set_index),]
	}
	credset_out_table = do.call("rbind", sub_sumstat)
	credset_out_table$gwas = gsub("\\..*", "", snp_file)
	credset_out_table
})
autoimmune_credset = do.call("rbind", autoimmune_credset)


row.names(autoimmune_credset)=NULL


autoimmune_credset$loc = paste0(autoimmune_credset$chr, ":", autoimmune_credset$pos, "-", autoimmune_credset$pos)

autoimmune_loci = GRanges(autoimmune_credset$loc)

gene_connected_ocrs = GRanges(gene_ocr_linkages)

index = as.data.frame(findOverlaps(gene_connected_ocrs, autoimmune_loci))


v2g = data.frame(gene_ocr_linkages[index[,1],], autoimmune_credset[index[,2],])
 kmean_rna = read.csv("tables/Tcell_activation_kmeansCluster_geneList.csv")
v2g = left_join(v2g, kmean_rna[,c(1:3)])
v2g$cluster[is.na(v2g$cluster)] = "not_differential"

#annotate hla

v2g$hla_region = (v2g$chr == "chr6") & (v2g$pos > 28477797) & (v2g$pos <33448354)

write.table(v2g, file = "tables/variant_to_gene_mapping_allTimepoints.txt", quote=F, row.names=F, sep="\t")

#Check for drugs that target genes


library("rDGIdb")
genes = unique(v2g$gene_name)
result <- queryDGIdb(genes)

results_table = resultSummary(result)
results_table_detailed = detailedResults(result)

gene_results_table =  byGene(result)

write.table(results_table, file = "tables/drugs/tcell_activation_v2g_all_drug_info_results_table.txt", sep="\t", quote=F, row.names=F)
write.table(results_table_detailed ,file = "tables/drugs/tcell_activation_v2g_all_drug_info_results_table_details.txt", sep="\t", quote=F, row.names=F)
write.table(gene_results_table, file = "tables/drugs/tcell_activation_v2g_all_drug_gene_results_table.txt", sep="\t", quote=F, row.names=F)

