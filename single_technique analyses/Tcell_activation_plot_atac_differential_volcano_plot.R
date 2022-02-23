
library(tidyverse)
load("data/comparisions/atac_differential_analysis.Rdata")
DE_df = DE_df %>% filter(comp == "CD4_8hr_vs_unstim" | comp == "CD4_24hr_vs_8hr" ) %>% 
	mutate(sig = ifelse(FDR <0.05, "sig", "ns")) %>%
	mutate(comp = factor(comp, levels = c("CD4_8hr_vs_unstim", "CD4_24hr_vs_8hr")))


pdf("plots/atac_ocr_differential_comparison_volcanoPlot.pdf", useDingbats=FALSE)
ggplot(DE_df, aes(x = logFC, y = -log10(FDR), color = sig))+
	geom_point()+
	scale_color_manual(values = c("grey", "red"))+
	facet_grid(~comp)+
	theme_minimal()+
	theme(aspect.ratio=1)+
	theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
#As png
png("plots/atac_ocr_differential_comparison_volcanoPlot.png", width =6, height = 4, units="in", res=300)
ggplot(DE_df, aes(x = logFC, y = -log10(FDR), color = sig))+
	geom_point()+
	scale_color_manual(values = c("grey", "red"))+
	facet_grid(~comp)+
	theme_minimal()+
	theme(aspect.ratio=1)+
	theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+

dev.off()

