#Gather pldsr results
library(tidyverse)
library("viridis")

dir = "tables/ldsr_ocrs/out"

files = list.files(dir)
files = files[grepl("results", files)]

pldsr = lapply(files, function(file){
	dat = read.delim(paste0(dir, "/", file))
	dat = dat[1,]
	dat$Category = gsub("_hic_linked_ocrs.results", "", file)
	dat$disease = gsub("_.*", "", dat$Category)
	dat$cluster = gsub(".*_", "", dat$Category)
	dat
})

 pldsr = do.call("rbind",  pldsr)
 pldsr$fdr = p.adjust(pldsr$Enrichment_p, method="BH")
 pldsr$z = -1 * qnorm(pldsr$Enrichment_p)
 pldsr$sig = ifelse(pldsr$fdr >0.01, "", ifelse(pldsr$fdr > 0.001, "**", "***"))

 pldsr$disease = factor( pldsr$disease, levels = rev(unique( pldsr$disease)))

pdf("plots/pldsr_output.pdf", useDingbats=FALSE, height=5, width=5)
 ggplot( pldsr, aes(x= cluster, y = disease, fill = z))+
 	geom_tile()+
	scale_fill_viridis(discrete=FALSE)+
	geom_text(aes(label = sig), size = 12, color = "white") +
	theme_minimal()
dev.off()


