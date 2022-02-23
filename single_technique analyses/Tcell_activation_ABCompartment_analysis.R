library(tidyverse)
library(GenomicRanges)
dir = "/mnt/isilon/sfgi/pahlm/projects/Tcell_activation/data/ab_compartment_coords"
files = list.files(dir)

ab_compart = lapply(files, function(f){
	path = paste(dir, f, sep="/")
	dat = read.delim(path, header=TRUE)
	n = gsub( "_2reps.eigen.40K.cis.vecs.tsv", "", f)
	dat$class = n
	dat
})

ab_compart = do.call("rbind", ab_compart )
ab_compart = ab_compart[complete.cases(ab_compart),]
ab_compart = ab_compart %>% mutate(call = ifelse(E1 > 0, "A", "B"))

ab_compart = ab_compart %>% mutate(chrom = factor(chrom, 
	levels = c(paste("chr", 1:22, sep=""), "chrX", "chrY")))

ab_class = ab_compart %>% select(chrom, start, end, call, class) %>%
mutate(class = gsub("unstimulated", "unstim", gsub("naiveT", "CD4", class)), start = start + 1) %>% 
dplyr::rename(chr =  chrom)


#save pdf for figure
pdf("plots/hic_abcompartments.pdf", useDingbats=FALSE, width =10 ,height=20)
ggplot(ab_compart, aes(x = (start+end)/2, y = E1))+
	geom_col(position="identity", aes(fill= call))+
	scale_fill_manual(values=c("indianred2", "steelblue"))+
facet_grid(chrom ~ class)+
theme_minimal()+
theme(
    panel.border=element_rect(colour="black",size=1, fill=NA),
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())
dev.off()

#save png for presentations
png("plots/hic_abcompartments.png",  width =10 ,height=20, units="in", res =300)
ggplot(ab_compart, aes(x = (start+end)/2, y = E1))+
	geom_col(position="identity", aes(fill= call))+
	scale_fill_manual(values=c("indianred2", "steelblue"))+
facet_grid(chrom ~ class)+
theme_minimal()+
theme(
    panel.border=element_rect(colour="black",size=1, fill=NA),
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())
dev.off()

######
#Compare to expression data

gene_ref = read.delim("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.gene.bed",header=F, skip=5)
gene_ref = gene_ref %>% 
	dplyr::rename(chr= V1, start=V2, end=V3, gene=V4, score= V5, strand = V6) %>% 
	mutate(start = start+1, gene_name = gsub("\\+.*", "", gene), gene_id = gsub(".*\\+", "", gene)) %>% 
	select(-gene) %>% 
	mutate(tss = ifelse(strand == "+", start, end)) %>%
	mutate(start = tss, end = tss)

load("data/norm_counts/rnaseq_tpm.Rdata")

tpm_mean = tpm %>% 
	group_by(condition,gene_id, gene_name, gene_type) %>%
	 summarize(tpm_mean = mean(tpm)) %>%
	 ungroup()

tpm_position_appended = left_join(tpm_mean, gene_ref)

index = as.data.frame(findOverlaps(GRanges(tpm_position_appended), GRanges(ab_class)))

ab_tpm_merged = data.frame(tpm_position_appended[index[,1],], ab_class[index[,2], ])
ab_tpm_merged = ab_tpm_merged[ab_tpm_merged$class==ab_tpm_merged$condition,]


ab_tpm_merged$class = factor(ab_tpm_merged$class, levels= c("CD4_unstim", "CD4_8hr", "CD4_24hr"))

pdf("plots/hic_abcompartments_expression_trend.pdf", useDingbats=FALSE)
ggplot(ab_tpm_merged, aes(x = call, y = log2(tpm_mean+1)))+
geom_violin(aes(fill=call))+
geom_boxplot(width = 0.1, outlier.shape = NA)+
	scale_fill_manual(values=c("indianred2", "steelblue"))+
facet_grid(~class)+
theme_minimal()+
theme(panel.border=element_rect(colour="black",size=1, fill=NA))
dev.off()

######
#Compare to atac-seq data

load("data/norm_counts/atac_fpkm.Rdata")

fpkm_mean = fpkm %>% 
	group_by(condition,chr, start, end, id) %>%
	 summarize(fpkm_mean = mean(fpkm)) %>%
	 ungroup()

fpkm_mean = fpkm_mean %>% mutate(midpoint = (start+end)/2) %>%
	mutate(start = midpoint, end = midpoint)

index = as.data.frame(findOverlaps(GRanges(fpkm_mean), GRanges(ab_class)))

fpkm_mean_ab = data.frame(fpkm_mean[index[,1],], ab_class[index[,2],])
fpkm_mean_ab = fpkm_mean_ab[fpkm_mean_ab$condition==fpkm_mean_ab$class,]
fpkm_mean_ab$class = factor(fpkm_mean_ab$class, levels= c("CD4_unstim", "CD4_8hr", "CD4_24hr"))

pdf("plots/hic_abcompartments_accessibility_trend.pdf", useDingbats=FALSE)
ggplot(fpkm_mean_ab, aes(x = call, y = log2(fpkm_mean+1)))+
geom_violin(aes(fill=call))+
geom_boxplot(width = 0.1, outlier.shape = NA)+
	scale_fill_manual(values=c("indianred2", "steelblue"))+
facet_grid(~class)+
theme_minimal()+
theme(panel.border=element_rect(colour="black",size=1, fill=NA))
dev.off()

######
#Identify compartments swapping sign or staying consistant
ab_class_timecourse = ab_class %>%
	group_by(chr, start, end) %>%
	summarize(call_timecourse =paste( call[order(factor(class, levels= c("CD4_unstim", "CD4_8hr", "CD4_24hr")))], collapse=">")) %>% 
	filter(grepl(".>.>.", call_timecourse)) 

ab_class_timecourse_tally = ab_class_timecourse %>%
	group_by(call_timecourse) %>%
	tally() %>%
	mutate(percent = n/sum(n) * 100) 

pdf("plots/hic_abcompartments_timecourse_changes.pdf")
ggplot(ab_class_timecourse_tally, aes(x= "", y= n, fill = call_timecourse))+
	geom_bar(stat="identity")+
  	coord_polar("y", start=0)+
  	theme_minimal()+
theme(
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())
dev.off()

#Identify genes in expression changing regions

index = as.data.frame(findOverlaps(GRanges(ab_class_timecourse), GRanges(tpm_position_appended )))

ab_class_tpm = data.frame(ab_class_timecourse[index[,1],], tpm_position_appended[index[,2],])

dynamic = c("A>B>A","A>B>B", "B>B>A", "B>A>A", "A>A>B", "B>A>B")

dab_class_tpm = ab_class_tpm[ab_class_tpm$call_timecourse %in% dynamic,]

ggplot(dab_class_tpm, aes(x = condition, y = tpm_mean))






library(ggplot2)

HG.test = function(grp,anno){
  universe = unique(unlist(anno)) #Set of genes that have been annotated, number of balls 
  gset = grp
  gset = gset[gset %in% universe]
  tmp = lapply(anno,function(j){
      p = length(j[j %in% gset]) #white balls drawn from an urn
      m = length(j) #white balls in the urn
      n = length(universe) - m #Number of black balls in the urn
      k = length(gset) #Number of balls drawn from the urn
      expected = ceiling(length(j)/length(universe)*length(gset)) #Number of white balls expected to be drawn, rounding up (replace ceiling with floor to round down, or remove to leave the decimal)
      hg.pval = phyper(p,m,n,k,lower.tail=FALSE) #Calculate p value
      genes = paste(j[j %in% gset], collapse=":") #Combine the list of genes
      data.frame(p.val = hg.pval,genes= genes, observed=p, expected=expected) #Put the data together in a data.frame
    })
   out = do.call("rbind",tmp) #Merge the list of dataframes into a single dataframe
   out$Pathway = names(tmp) #Convert row.names to a column
   row.names(out) = NULL #Remove the row.names
   out$p.adjust = p.adjust(out$p.val,method="BH") #Calculate p.value after FDR correction
   out #Print results
   out = out[out$p.adjust<0.05 & out$observed>0,]
   out = out[order(out$p.val,decreasing=FALSE), c(5,3,4,1,6,2)]
}

#Reactome
load("/mnt/isilon/sfgi/pahlm/annotationFiles/msigdb_v7.0_GMTs/c2.cp.reactome.v7.0.symbols.Rdata")

c2.cp.reactome.nohistone = lapply(c2.cp.reactome, function(x){
    x[grepl("HIST", x)==FALSE]
})


geneset_enrichment_out = lapply(dynamic, function(d){
	x = unique(ab_class_tpm[ab_class_tpm$call_timecourse == d & ab_class_tpm$tpm_mean >1,]$gene_name)
	geneset.hg <- HG.test(x, c2.cp.reactome.nohistone) 
	
	if(nrow(geneset.hg)> 0){
	geneset.hg$class = d
	geneset.hg 
}
})

 geneset_enrichment_out = do.call("rbind",  geneset_enrichment_out) 
  geneset_enrichment_out  = geneset_enrichment_out   %>% mutate(log2_enrichment=log2(observed/expected+1), Score= -log10(p.adjust+1e-100)) %>% group_by(class) %>%
mutate(rank = rank(-Score,ties.method = "first"))

highlighted_pathways = unique(geneset_enrichment_out$Pathway)

pdf("plots/Reactome_Pathway_Tcell_AB_compartment.pdf", useDingbats=FALSE)
ggplot(geneset_enrichment_out, aes(x = log2_enrichment, y = Score, color = class))+
	geom_point(aes(size = observed))+
	geom_text(aes(label=ifelse(Pathway %in% highlighted_pathways,gsub("_", " ", gsub("REACTOME_", "", as.character(Pathway))),'')),hjust=0,vjust=0, color = 'black', size = 2)+
	theme_minimal()+
	#geom_hline(yintercept= 1.30103, linetype="dashed", color = "red")+
	facet_wrap(~class)+
	theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
