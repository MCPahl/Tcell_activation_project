library(tidyverse)

load("data/comparisions/hic_differential_analysis.Rdata")

comp_direction = sigDE_consensus_loop %>% 
select(chr_a, start_a, end_a, chr_b, start_b, end_b, logFC, p.adj, comp) %>%
filter(comp == "hr8_vs_unstim" | comp == "hr24_vs_hr8") %>% 
unique() %>%
mutate(direction = ifelse(logFC > 0, "up", "down")) %>%
group_by(comp, direction) %>% 
tally() %>%
mutate(n = ifelse(direction == "down", n*-1, n)) %>%
mutate(comp = factor(comp, levels = c("hr8_vs_unstim", "hr24_vs_hr8")))

pdf("plots/hic_loop_sig_comparison.pdf", width = 4, height = 2)
ggplot(comp_direction, aes(x = comp, y=n, fill = direction))+
geom_bar(stat = "identity", position="dodge", width = 0.5)+
scale_fill_manual(values = c("blue", "red"))+
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
	

load("data/loop_calls/fithic_loops_3resolutions.Rdata")
fithic_loops_tally = fithic_loops %>% group_by(condition) %>% tally() %>% 
mutate(timepoint = ifelse(grepl("0hr", condition), "CD4_unstim", ifelse(grepl("8hr", condition), "CD4_8hr", "CD4_24hr")), resolution = gsub(".*_", "", condition)) %>%
mutate(timepoint = factor(timepoint, levels = c("CD4_unstim", "CD4_8hr", "CD4_24hr")))

pdf("plots/hic_loop_sig_tally.pdf", width = 4, height = 2)
ggplot(fithic_loops_tally, aes(x = timepoint, y=n, fill = resolution))+
geom_bar(stat = "identity", position="dodge", width = 0.5)+
scale_fill_manual(values = c("grey", "darkgrey", "black"))+
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()