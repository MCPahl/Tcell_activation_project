 library(tidyverse)
 load("data/norm_counts/atac_fpkm.Rdata")
 load("data/norm_counts/rnaseq_tpm.Rdata")
 load("intersection/gene_ocr_linkages_differental_annotated.RData")
 kclusters = read.csv("tables/Tcell_activation_kmeansCluster_geneList.csv") %>% select(gene_name, gene_id, cluster)

data.dir = "/mnt/isilon/sfgi/suc1/analyses/wells/hiC/compare/naiveT_stim/data/4000/"
files = list.files(data.dir)
files = files[grepl("chr", files)]

paths = paste0(data.dir, files, "/normIF.Rdata")

anno_key_diff = anno_key %>% filter(ATAC_FDR < 0.05  & HiC_FDR < 0.05)


chroms = paste0("chr", c(1:22, "X"))
anno_IF = lapply(chroms, function(chrom){
    path = paths[grepl(paste0("/", chrom, "/"), paths)]
    print(path)

    load(path)
        normIF = normIF %>% 
        rename(chr_a = chr, start_a = region1, end_a = region2) %>% 
        select(-D) %>%
        mutate(chr_a = paste0("chr", chr_a))

        anno_key_diff_tmp = anno_key_diff %>% 
        inner_join(normIF) %>% 
        rowwise %>%
        mutate(naiveT_unstimulated_mean = sum(c_across(c(starts_with('naiveT_unstimulated'))))/2, 
        naiveT_8hr_mean = sum(c_across(c(starts_with('naiveT_8hr'))))/2, 
        naiveT_24hr_mean = sum(c_across(c(starts_with('naiveT_24hr'))))/2) %>%
        ungroup() %>% 
        select(-naiveT_24hr_ND517, -naiveT_24hr_TMP442, -naiveT_8hr_ND517, -naiveT_8hr_TMP442, -naiveT_unstimulated_ND517, -naiveT_unstimulated_TMP442)
})

anno_IF = do.call("rbind", anno_IF )

anno_IF = anno_IF %>%
   gather(key = HiC_norm, value = mean_IF, naiveT_unstimulated_mean, naiveT_8hr_mean, naiveT_24hr_mean) %>% 
   group_by(chr_a, start_a, end_a, chr_b, start_b, end_b, ocr, gene_name, gene_id,comp) %>%
   mutate(mean_IF = scale(mean_IF)) %>% 
   ungroup() %>%
   spread(key = HiC_norm, value = mean_IF) %>% 
   rename(HiC_CD4_unstim = naiveT_unstimulated_mean, HiC_CD4_8hr = naiveT_8hr_mean, HiC_CD4_24hr = naiveT_24hr_mean)

mean_fpkm = fpkm %>%
    group_by(condition, id) %>%
    summarize(fpkm = mean(fpkm)) %>%
    mutate(fpkm = scale(fpkm)) %>%
    spread(key = condition, value = fpkm) %>%
    rename(ocr = id, ATAC_CD4_24hr= CD4_24hr,ATAC_CD4_8hr = CD4_8hr, ATAC_CD4_unstim = CD4_unstim)


mean_tpm = tpm %>%
    group_by(gene_name, gene_id, condition) %>%
    summarize(tpm= mean(tpm)) %>%
    mutate(tpm = scale(tpm)) %>%
    spread(key = condition, value = tpm) %>%
    rename(RNA_CD4_24hr= CD4_24hr,RNA_CD4_8hr = CD4_8hr, RNA_CD4_unstim = CD4_unstim)


anno_IF = anno_IF %>% left_join(mean_fpkm) %>% left_join(mean_tpm)

save(anno_IF, file = "intersection/summarized_differential_conditionMean_normalized_zscores.Rdata") %>% 

v2g = read.table("tables/variant_to_gene_mapping_allTimepoints.txt",header=TRUE) %>% rename(ocr = ocr_id)

#Dynamic genes & snp pairs
 gene = c(
"SIK1",
"IL2",
"PARK7",
"DUSP5",
"PARK7",
"CLEC2D",
"TRIP10",
"GPR108",
"CLEC2D"
)
var = c(
"rs73380290",
"rs62322750",
"rs35554366",
"rs180826373",
"rs12121408",
"rs10844503",
"rs1077667",
"rs1077667", 
"chr1:7962954",
"rs10844503")


v2g_filt = v2g %>%
    filter(gene_name %in% gene & rsid %in% var) %>%
    unique()


anno_IF_filt = anno_IF %>% 
    select(chr_a, start_a, end_a, chr_b, start_b, end_b, ocr, gene_name, gene_id, HiC_CD4_unstim,  HiC_CD4_8hr, HiC_CD4_24hr, ATAC_CD4_unstim, ATAC_CD4_8hr, ATAC_CD4_24hr  ,RNA_CD4_unstim,  RNA_CD4_8hr,  RNA_CD4_24hr) %>%
    inner_join(v2g_filt) %>% 
    mutate(pp_z = scale(pp), anno = paste0(rsid, ":", gene_name)) %>% 
    select(anno, pp_z, HiC_CD4_unstim, HiC_CD4_8hr, HiC_CD4_24hr, ATAC_CD4_unstim,ATAC_CD4_8hr, ATAC_CD4_24hr, RNA_CD4_unstim, RNA_CD4_8hr, RNA_CD4_24hr) %>% 
    gather(key = class, value = zscore,-anno) %>% 
    group_by(anno, class) %>%
    summarize(zscore = min(zscore)) %>%
    ungroup() %>% 
    mutate(value_type = gsub("_.*", "", class)) %>% 
    mutate(class = factor(class, levels = c("pp_z", "RNA_CD4_unstim", "RNA_CD4_8hr", "RNA_CD4_24hr", 
        "ATAC_CD4_unstim", "ATAC_CD4_8hr", "ATAC_CD4_24hr",
        "HiC_CD4_unstim", "HiC_CD4_8hr", "HiC_CD4_24hr")))%>%
    mutate(value_type = factor(value_type, levels = c("pp", "RNA", "ATAC", "HiC"))) 




library(viridis)

pdf("plots/top_differential_zscore.pdf", width = 15)
ggplot(anno_IF_filt, aes(x = class, y =anno, fill = zscore))+
geom_tile()+
scale_fill_viridis(discrete=FALSE)+
theme_bw()+
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
 facet_wrap(~value_type, scale="free", nrow=1)
dev.off()



