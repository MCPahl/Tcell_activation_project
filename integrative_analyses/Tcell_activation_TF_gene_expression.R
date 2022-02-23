#unstim and 8hr
library(tidyverse)

tf_dat = read.delim("tf/hint_atac/hint_out/differential/CD4_8hr_vs_Unstim/differential_statistics.txt")


load("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/jasper2020/motif_anno.Rdata")
motif_anno = motif_anno %>%
	mutate(Motif = paste(motif_id, TF_old, sep=".")) %>%
	select(Motif, motif_id, gene_name, gene_id, TF_family)
tf_dat_anno = inner_join(tf_dat, motif_anno)


load("data/comparisions/rna_differential_analysis.Rdata")
DE_df_comp = DE_df %>% select(gene_name, gene_id, logFC, FDR, comp) %>% filter(comp == "CD4_8hr_vs_unstim") %>% mutate(gene_id = gsub("\\.\\d*", "", gene_id))

tf_dat_anno_comp = inner_join(tf_dat_anno, DE_df_comp)
tf_dat_anno_comp = tf_dat_anno_comp %>%
	mutate(class = ifelse(FDR <0.05 & P_values< 0.05, "both", 
	 	ifelse(FDR <0.05, "expression",
	 		ifelse(P_values< 0.05, "footprint", "NS"))))


write.table(tf_dat_anno_comp, file = "tables/TF_footprinting_differential_comparison_8hr_vs_unstim.tsv", sep="\t", quote=F, row.names=F)


highlight = c("MA0472.2.EGR2", "MA1633.1.BACH1", "MA1101.2.BACH2", "MA1138.1.FOSL2::JUNB", "MA0099.3.FOS::JUN", "MA1634.1.BATF", "MA0089.2.NFE2L1",
	"MA0732.1.EGR3", "MA1112.2.NR4A1", "MA1419.1.IRF4", "MA0517.1.STAT1::STAT2", "MA0050.2.IRF1", "MA0059.1.MAX::MYC", "MA0519.1.Stat5a::Stat5b", "MA1515.1.KLF2",
	"MA0508.3.PRDM1", "MA0506.1.NRF1", "MA0769.2.TCF7", "MA0098.3.ETS1", "MA0603.1.Arntl", "MA1564.1.SP9",
	"MA1564.1.SP9", "MA1622.1.Smad2::Smad3", "MA0157.2.FOXO3", "MA0850.1.FOXP3", "MA1522.1.MAZ")


pdf("plots/TF_footprinting_differential_comparison_8hr_vs_unstim.pdf")
ggplot(tf_dat_anno_comp, aes(x = logFC , y =  -1 * TF_Activity, color = class))+
geom_point(aes(size = Num))+
scale_color_manual(values=c("red",  "orange", "yellow4", "grey" ))+
geom_text(aes(label=ifelse(Motif %in% highlight, Motif, '')),hjust=0,vjust=0, color = 'black', size = 4)+
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

####
#24 and 8hr


library(tidyverse)

tf_dat = read.delim("tf/hint_atac/hint_out/differential/CD4_24hr_vs_8hr/differential_statistics.txt")


load("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/jasper2020/motif_anno.Rdata")
motif_anno = motif_anno %>%
	mutate(Motif = paste(motif_id, TF_old, sep=".")) %>%
	select(Motif, motif_id, gene_name, gene_id, TF_family)
tf_dat_anno = inner_join(tf_dat, motif_anno)


load("data/comparisions/rna_differential_analysis.Rdata")
DE_df_comp = DE_df %>% select(gene_name, gene_id, logFC, FDR, comp) %>% filter(comp == "CD4_24hr_vs_8hr") %>% mutate(gene_id = gsub("\\.\\d*", "", gene_id))

tf_dat_anno_comp = inner_join(tf_dat_anno, DE_df_comp)
tf_dat_anno_comp = tf_dat_anno_comp %>%
	mutate(class = ifelse(FDR <0.05 & P_values< 0.05, "both", 
	 	ifelse(FDR <0.05, "expression",
	 		ifelse(P_values< 0.05, "footprint", "NS"))))


write.table(tf_dat_anno_comp, file = "tables/TF_footprinting_differential_comparison_24hr_vs_8hr.tsv", sep="\t", quote=F, row.names=F)


highlight = c("MA0472.2.EGR2", "MA0637.1.CENPB", "MA0871.2.TFEC", "MA0732.1.EGR3", "MA0506.1.NRF1", "MA1515.1.KLF2", "MA1522.1.MAZ", "MA0079.4.SP1", "MA0516.2.SP2",
	"MA1513.1.KLF15", "MA0685.1.SP4", "MA1564.1.SP9", "MA1653.1.ZNF148", "MA0753.2.ZNF740", "MA0741.1.KLF16",
	"MA0099.3.FOS::JUN", "MA0024.3.E2F1", "MA0796.1.TGIF1", "MA0693.2.VDR", "MA0071.1.RORA", "MA0478.1.FOSL2", "MA0508.3.PRDM1",
	"MA1418.1.IRF3", "MA0463.2.BCL6", "MA0151.1.Arid3a","MA0523.1.TCF7L2", "MA1101.2.BACH2" 
	"MA0474.2.ERG", "MA0888.1.EVX2")

pdf("plots/TF_footprinting_differential_comparison_24hr_vs_8hr.pdf")
ggplot(tf_dat_anno_comp, aes(x = logFC , y =  -1 * TF_Activity, color = class))+
geom_point(aes(size = Num))+
scale_color_manual(values=c("red",  "orange", "yellow4", "grey" ))+
geom_text(aes(label=ifelse(Motif %in% highlight, Motif, '')),hjust=0,vjust=0, color = 'black', size = 4)+
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

