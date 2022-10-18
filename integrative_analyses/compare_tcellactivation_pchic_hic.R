library(tidyverse)
library("parseIbed")
library("DFbedtools")
library(GenomicRanges)

#Load format HiC
load("data/loop_calls/fithic_loops_3resolutions.Rdata")
fithic_loops = fithic_loops %>% dplyr::rename(chr_a = chrA, start_a = startA, end_a = endA, chr_b = chrB, start_b = startB, end_b = endB)

#Load format pcHiC
pchic <- read.delim("external_data/Burren_2017_GB_PCHiC_ActivatedTCells/13059_2017_1285_MOESM5_ESM")
pchic_reformat = pchic %>% select(baitChr, baitStart, baitLength, oeChr, oeStart, oeLength,Total_CD4_Activated, Total_CD4_NonActivated) %>%
	mutate(baitEnd = baitStart+baitLength, oeEnd = oeStart+ oeLength) %>%
	filter(baitChr == oeChr) %>% 
	select(-c(baitLength, oeLength)) %>%
	dplyr::rename(chr_a = baitChr, start_a = baitStart, end_a = baitEnd, chr_b = oeChr, start_b = oeStart, end_b = oeEnd) %>%
	mutate(loop_stage = ifelse(Total_CD4_Activated > 5 & Total_CD4_NonActivated > 5 ,"both", ifelse(Total_CD4_Activated > 5, "activated", "nonactivated")),
		dist = (start_b+ end_b)/2 - (start_a+ end_a)/2)

#Loop through each HiC condition (unstim, 8hr, 24hr at 1k,2k,4k loop calls), and PCHIC condition (activated, nonactivated, both)

#Find overlap of both ends (GRanges each side to get each side's index, then take intersect of both indexes)
#Then count number of shared loops, number of 

hic_conds = unique(fithic_loops$condition)
pchic_lstages = unique(pchic_reformat$loop_stage)

loop_comparison = lapply(hic_conds, function(hic_cond){
	hic_gene_loop = unique(fithic_loops[fithic_loops$condition == hic_cond,])
	x = lapply(pchic_lstages, function(pchic_lstage){
		pchic_gene_loops = unique(pchic_reformat[pchic_reformat$loop_stage == pchic_lstage,])

		chr_a_overlap = data.frame(findOverlaps(
			GRanges(hic_gene_loop$chr_a, ranges = IRanges(hic_gene_loop$start_a, hic_gene_loop$end_a)),
			GRanges(pchic_gene_loops$chr_a, ranges = IRanges(pchic_gene_loops$start_a, pchic_gene_loops$end_a))
		))
		
		chr_b_overlap = data.frame(findOverlaps(
			GRanges(hic_gene_loop$chr_b, ranges = IRanges(hic_gene_loop$start_b, hic_gene_loop$end_b)),
			GRanges(pchic_gene_loops$chr_b, ranges = IRanges(pchic_gene_loops$start_b, pchic_gene_loops$end_b))
		))
		shared_loops = inner_join(chr_a_overlap, chr_b_overlap)
		data.frame(hic_condition = hic_cond,
			pchic_conditon = pchic_lstage,
			hic_loops_n = nrow(hic_gene_loop), 
			pchic_loop_n = nrow(pchic_gene_loops), 
			pchic_shared_sideA = length(unique(chr_a_overlap$subjectHits)), 
			pchic_shared_sideB = length(unique(chr_b_overlap$subjectHits)), 
			shared_loops = nrow(shared_loops)
		)
	})
	x = do.call("rbind", x)
})
loop_comparison = do.call("rbind", loop_comparison)

#Filter to only 4k resolution and get % shared
loop_comparison = loop_comparison %>%
	 mutate(hic_res = gsub(".*_", "",  hic_condition), hic_timepoint = gsub("_\\dkb", "", hic_condition)) %>%
	 filter(hic_res == "4kb") %>%
	 mutate(pchic_loop_shared = shared_loops/pchic_loop_n * 100) %>%
	 
	 loop_comparison$hic_timepoint = factor(loop_comparison$hic_timepoint, levels = c("HiC_24hr", "HiC_8hr","HiC_0hr" ))

#Plot as heatmap
pdf("plots/Tcell_activation_pchic_tcell_activation.pdf")
ggplot(loop_comparison, aes(x = pchic_conditon, y = hic_timepoint, fill = pchic_loop_shared))+
	geom_tile()+
	scale_fill_gradient2(low="darkblue", high="darkgreen")+
	theme_bw()
dev.off()

#Write table
 write.table(loop_comparison, file = "Burren_loop_comparison.txt", quote=F, row.names=F, sep="\t")


####Compare both specific

fithic_loops_4kb = fithic_loops %>% 
	filter(grepl("4kb", condition)) %>%
	mutate(loop = paste0(chr_a, start_a, end_a, "_", chr_b, start_b, end_b))

unstim_loops = setdiff(fithic_loops_4kb[fithic_loops_4kb$condition== "HiC_0hr_4kb",]$loop,
	union(fithic_loops_4kb[fithic_loops_4kb$condition== "HiC_8hr_4kb",]$loop,
	fithic_loops_4kb[fithic_loops_4kb$condition== "HiC_24hr_4kb",]$loop))

hr8_loops = setdiff(fithic_loops_4kb[fithic_loops_4kb$condition== "HiC_8hr_4kb",]$loop,
	union(fithic_loops_4kb[fithic_loops_4kb$condition== "HiC_0hr_4kb",]$loop,
	fithic_loops_4kb[fithic_loops_4kb$condition== "HiC_24hr_4kb",]$loop))

hr24_loops = setdiff(fithic_loops_4kb[fithic_loops_4kb$condition== "HiC_24hr_4kb",]$loop,
	union(fithic_loops_4kb[fithic_loops_4kb$condition== "HiC_0hr_4kb",]$loop,
	fithic_loops_4kb[fithic_loops_4kb$condition== "HiC_8hr_4kb",]$loop))

shared_loops = intersect(fithic_loops_4kb[fithic_loops_4kb$condition== "HiC_24hr_4kb",]$loop,
	intersect(fithic_loops_4kb[fithic_loops_4kb$condition== "HiC_0hr_4kb",]$loop,
	fithic_loops_4kb[fithic_loops_4kb$condition== "HiC_8hr_4kb",]$loop))


specific = bind_rows(
fithic_loops_4kb %>% filter(condition== "HiC_0hr_4kb" & loop %in% unstim_loops ) %>% mutate(condition = "unstim_specific"),
fithic_loops_4kb %>% filter(condition== "HiC_8hr_4kb" & loop %in% hr8_loops ) %>% mutate(condition = "8hr_specific"),
fithic_loops_4kb %>% filter(condition== "HiC_24hr_4kb" & loop %in% hr24_loops ) %>% mutate(condition = "24hr_specific"),
fithic_loops_4kb %>% filter(condition== "HiC_0hr_4kb" & loop %in% shared_loops ) %>% mutate(condition = "shared")
)

hic_conds_specific = unique(specific $condition)


loop_comparison = lapply(hic_conds_specific, function(hic_cond){
	hic_gene_loop = unique(specific[specific$condition == hic_cond,])
	x = lapply(pchic_lstages, function(pchic_lstage){
		pchic_gene_loops = unique(pchic_reformat[pchic_reformat$loop_stage == pchic_lstage,])

		chr_a_overlap = data.frame(findOverlaps(
			GRanges(hic_gene_loop$chr_a, ranges = IRanges(hic_gene_loop$start_a, hic_gene_loop$end_a)),
			GRanges(pchic_gene_loops$chr_a, ranges = IRanges(pchic_gene_loops$start_a, pchic_gene_loops$end_a))
		))
		
		chr_b_overlap = data.frame(findOverlaps(
			GRanges(hic_gene_loop$chr_b, ranges = IRanges(hic_gene_loop$start_b, hic_gene_loop$end_b)),
			GRanges(pchic_gene_loops$chr_b, ranges = IRanges(pchic_gene_loops$start_b, pchic_gene_loops$end_b))
		))
		shared_loops = inner_join(chr_a_overlap, chr_b_overlap)
		data.frame(hic_condition = hic_cond,
			pchic_conditon = pchic_lstage,
			hic_loops_n = nrow(hic_gene_loop), 
			pchic_loop_n = nrow(pchic_gene_loops), 
			pchic_shared_sideA = length(unique(chr_a_overlap$subjectHits)), 
			pchic_shared_sideB = length(unique(chr_b_overlap$subjectHits)), 
			shared_loops = nrow(shared_loops)
		)
	})
	x = do.call("rbind", x)
})
loop_comparison = do.call("rbind", loop_comparison)


#Filter to only 4k resolution and get % shared
loop_comparison = loop_comparison %>%
	 mutate(hic_res = gsub(".*_", "",  hic_condition), hic_timepoint = gsub("_\\dkb", "", hic_condition)) %>%
	 mutate(pchic_loop_shared = shared_loops/pchic_loop_n * 100) 
	 
	 loop_comparison$hic_timepoint = factor(loop_comparison$hic_timepoint, levels = rev(c("shared","unstim_specific", "8hr_specific", "24hr_specific")))

#Plot as heatmap
pdf("plots/Tcell_activation_pchic_tcell_activation_specific.pdf")
ggplot(loop_comparison, aes(x = pchic_conditon, y = hic_timepoint, fill = pchic_loop_shared))+
	geom_tile()+
	scale_fill_gradient2(low="darkblue", high="darkgreen")+
	theme_bw()
dev.off()

#Write table
 write.table(loop_comparison, file = "Burren_loop_comparison_specific.txt", quote=F, row.names=F, sep="\t")



 Cornelis Blauwendraat
