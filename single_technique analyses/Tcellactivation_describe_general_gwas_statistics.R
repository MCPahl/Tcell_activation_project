library(tidyverse)

gwas = read_delim("tables/gwas_summary.txt") %>% 
	rename(sentinels = "sumStat snp_n", N = "EUR N", Ncase = "EUR case", Ncontrol = "EUR control" )

gwas$Trait = factor(gwas$Trait, levels = rev(unique(gwas$Trait)))

pdf("plots/gwas_signal_count.pdf")
ggplot(gwas, aes(x = Trait, y = sentinels))+
geom_bar(stat="identity")+
theme_bw()+
coord_flip()
dev.off()


gwas_long = gwas %>%
	gather(key = group, value = n_stratified, Ncase, Ncontrol ) %>% 
	select(Trait, group, n_stratified)

pdf("plots/gwas_n.pdf")
ggplot(gwas_long, aes(x = Trait, y = n_stratified, fill = group))+
geom_bar(stat="identity", position= "dodge")+
theme_bw()+
coord_flip()
dev.off()