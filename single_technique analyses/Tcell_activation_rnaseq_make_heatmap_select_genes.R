#make plot
library(tidyverse)

dir = "tables"
files = list.files(dir)
files = files[grepl("Reactome_Pathway_Tcell_activation_cluster", files)]


enrichment = lapply(files, function(file){
	path = paste(dir, file, sep="/")
	dat = read.csv(path)
	dat$cluster = gsub("Reactome_Pathway_Tcell_activation_", "", gsub(".csv", "", file))
	dat
})

enrichment = do.call("rbind", enrichment)

enrichment = enrichment %>% mutate(log2_enrichment=log2(observed/expected+1), Score= -log10(p.adjust+1e-100)) %>% group_by(cluster) %>%
mutate(rank = rank(-Score,ties.method = "first"))

highlighted_pathways= c(
"REACTOME_METABOLISM_OF_RNA", 
"REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR",
"REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
"REACTOME_TRANSLATION",
"REACTOME_DEATH_RECEPTOR_SIGNALLING",
"REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
"REACTOME_MET_ACTIVATES_PTK2_SIGNALING",
 "REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES",
 "REACTOME_INFECTIOUS_DISEASE",
 "REACTOME_REGULATION_OF_MRNA_STABILITY_BY_PROTEINS_THAT_BIND_AU_RICH_ELEMENTS",
 "REACTOME_REGULATION_OF_APOPTOSIS",
 "REACTOME_CELLULAR_RESPONSES_TO_STRESS",
 "REACTOME_DNA_REPAIR",
 "REACTOME_UNWINDING_OF_DNA",
 "REACTOME_BETA_OXIDATION_OF_OCTANOYL_COA_TO_HEXANOYL_COA",
 "REACTOME_CELL_CYCLE",
 "REACTOME_CHROMOSOME_MAINTENANCE"
)

pdf("plots/Reactome_Pathway_Tcell_activation_cluster_scatterplot.pdf", useDingbats=FALSE)
ggplot(enrichment, aes(x = log2_enrichment, y = Score, color = cluster))+
	geom_point(aes(size = observed))+
	geom_text(aes(label=ifelse(Pathway %in% highlighted_pathways,gsub("_", " ", gsub("REACTOME_", "", as.character(Pathway))),'')),hjust=0,vjust=0, color = 'black', size = 2)+
	theme_minimal()+
	geom_hline(yintercept= 1.30103, linetype="dashed", color = "red")+
	facet_wrap(~cluster)+
	theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

