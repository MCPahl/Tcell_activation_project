#R version 4.1.1
setwd("/mnt/isilon/sfgi/pahlm/projects/Tcell_activation")
library(tidyverse)
library(cluster)
library(vegan)


#Set seed 
set.seed(777)

#Load data
load("data/norm_counts/rnaseq_tpm.Rdata")
load("data/comparisions/rna_differential_analysis.Rdata")

#Get differential genes
de.genes = unique(sigDE_df$gene_id)
z = tpm %>% filter((gene_type %in% "noncoding_RNA_short")==F) %>% select(gene_name, gene_id, sample, tpm) %>% distinct() %>% spread(key= sample, value = tpm)
tpm_var <- apply(z[,3:length(z)], 1, var)
tpm_mean <- apply(z[,3:length(z)], 1, mean)
diag = data.frame(gene_id = z$gene_id, gene_name=z$gene_name, tpm_mean= tpm_mean, tpm_var = tpm_var)
diag$diff = diag$gene_id %in% de.genes

#Plot the mean and variance across samples
pdf("plots/diagnosis/rnaseq_tpm_mean_vs_variance.pdf", useDingbats=FALSE)
ggplot(diag, aes(x=log2(tpm_mean), y=log2(tpm_var), colour = diff))+
geom_point(alpha=0.25)+
scale_color_manual(values = c('#999999','#E69F00'))+
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()


#filter to differentially expressed genes (any timepoint)
z.filt = z %>% filter(gene_id %in% de.genes)

#Change from tibble to matrix (sample x gene)
z.filt.matrix = z.filt %>% mutate(id = paste(gene_name, gene_id, sep="|"))
z.filt.matrix =as.data.frame(z.filt.matrix )
row.names(z.filt.matrix)= z.filt.matrix$id
z.filt.matrix$gene_name=NULL
z.filt.matrix$gene_id=NULL
z.filt.matrix$id=NULL
z.filt.matrix = as.matrix(z.filt.matrix)

#Center and scale
z.scaled.matrix = t(scale(t(z.filt.matrix)))

#Identify clusters
z.hc = hclust(as.dist(1-cor(z.scaled.matrix, method="spearman")), method="complete")
z.tree = as.dendrogram(z.hc, method="average")

#Plot relationship between samples
pdf("plots/diagnosis/rnaseq_tpm_diff_clustering.pdf", useDingbats=FALSE, height=16, width=12)
plot(z.tree,
     main = "Sample Clustering",
     ylab = "Height")
dev.off()

#Testing metrics for determining K

#SSE: sum of squared error, sum of squared distance between each cluster member and it's centroid, vary the number of centers (clusters)
#Adding new clusters should reduce the SSE since the distance between the center and any point in the cluster will be smaller
#Look for inflection point where you don't improve this metric, find when it's not worth adding additional clusters
#Lower values are better

wss <- (nrow(z.scaled.matrix -1)*sum(apply(z.scaled.matrix ,2,var)))
for (i in 2:25) wss[i] <- sum(kmeans(z.scaled.matrix,centers=i)$withinss)
done

wss = data.frame(centers = 1:25, within_group_sum_of_squares = wss)

pdf("plots/diagnosis/rnaseq_wss.pdf", useDingbats=FALSE)
ggplot(wss, aes(x=centers, y=within_group_sum_of_squares))+
geom_line()+
geom_point()+
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

#Silhouette width
#Another metric for cluster quality
#Sil describes how similar a gene is to own cluster compared to other clusters. Higher values are better
sil <- rep(0,25)

#Takes around ~30 min
for(i in 2:25){
  k1to20 <- kmeans(z.scaled.matrix, centers = i, nstart = 25, iter.max = 20)
  ss <- silhouette(k1to20$cluster, dist(z.scaled.matrix))
  sil[i] <- mean(ss[, 3])
}
sil = data.frame(centers = 1:25, sil_width = sil)

pdf("plots/diagnosis/hic_sil.pdf", useDingbats=FALSE)
ggplot(sil, aes(x=centers, y=sil_width))+
geom_line()+
geom_point()+
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

#Calinsky-Harabasz index
#Based on a intra and inter cluster sum of squares, maximizing the index indicates how well seperated the clusters are
fit <- cascadeKM(z.scaled.matrix, 1, 25, iter = 100)

pdf("plots/diagnosis/hic_calinsky.pdf", useDingbats=FALSE)
plot(fit, sortg = TRUE, grpmts.plot = TRUE)
dev.off()

#Gap statistic
#Compares the log SSE with it's expectation under a null distribution, identifying where the number of clusters where the gap between the log(wss) and null ref is highest

gap <- clusGap(z.scaled.matrix, kmeans, 25, B = 100, verbose = interactive())

plot(gap, main = "Gap statistic")
abline(v=which.max(gap$Tab[,3]), lty = 2)

#Perform K means clustering
kClust = kmeans(z.scaled.matrix, centers = 10, nstart = 1000, iter.max=100)
kClusters <- kClust$cluster

#Plot number of genes per cluster
cluster_counts = summary(as.factor(kClusters))
cluster_counts = data.frame(cluster = names(cluster_counts), counts=cluster_counts)
pdf("plots/Cluster_counts.pdf", useDingbats=FALSE)
plot = ggplot(cluster_counts, aes(x= cluster, y = counts, fill=cluster))+
  geom_bar(stat='identity')+
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
print(plot)
dev.off()

#Identify the centroids of clusters per timepoint
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}
kClustcentroids <- sapply(levels(factor(kClusters)), clust.centroid, z.scaled.matrix, kClusters)

#Check how the centroids correlate with one another to determine if the number of clusters should be reduced
kClustcentroids_cor = cor(kClustcentroids) %>% melt()

#Plot the clusters
pdf("plots/diagnosis/rnaseq_kmeansClusters_correlation.pdf", useDingbats=FALSE)
ggplot(kClustcentroids_cor, aes(x=X1, y=X2, fill= value))+
geom_tile()+    
geom_text(aes(label = round(value, 2)), color="white")+
scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  xlab("Cluster") + ylab("Cluster")
dev.off()


#Prepare graph of the pattern of centroids overtime (representing mean and sd of replicates)
library(reshape)
#get in long form for plotting
Kmolten <- melt(kClustcentroids)
colnames(Kmolten) <- c('sample','cluster','value')
Kmolten$sample = factor(Kmolten$sample, levels=c("CD4_unstim_ND543", "CD4_unstim_ND581", "CD4_unstim_ND589", "CD4_8hr_ND543",  "CD4_8hr_ND581", "CD4_8hr_ND589", "CD4_24hr_ND543", "CD4_24hr_ND581", "CD4_24hr_ND589"))
Kmolten$class = gsub("_ND.*", "", Kmolten$sample)
Kmolten$class = factor(Kmolten$class, levels=c("CD4_unstim", "CD4_8hr", "CD4_24hr"))
Kmolten_rep_summary = Kmolten %>% group_by(class,cluster) %>% summarize(mean=mean(value), sd=sd(value))

#plot centroid pattern
pdf("plots/diagnosis/rnaseq_kmeansClusters_centroid_behavior.pdf", useDingbats=FALSE)
ggplot(Kmolten_rep_summary, aes(x=class,y=mean, group=cluster, colour=as.factor(cluster))) + 
  geom_point() + 
  geom_line() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster Expression by Time",color = "Cluster")+
  theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
facet_grid(rows = vars(cluster))
dev.off()

#Identify correlation between each gene with the centroid (identify the core genes vs poorly matching genes)
k.index = 1:5
kMeans_cluster_score <- lapply(k.index, function(index){
  core<- Kmolten[Kmolten$cluster==index,]
  K <- (z.scaled.matrix[kClusters==index,])
  corscore <- function(x){
    cor(x,core$value)
  }
  score <- apply(K, 1, corscore)
  Kmolten <- melt(K)
  colnames(Kmolten) <- c('gene','sample','value')
  Kmolten <- merge(Kmolten,score, by.x='gene',by.y='row.names', all.x=T)
  colnames(Kmolten) <- c('gene','sample','value','score')
  Kmolten$order_factor <- 1:length(Kmolten$gene)
  Kmolten <- Kmolten[order(Kmolten$score),]
  Kmolten$cluster = index
  Kmolten
  })
kMeans_cluster_score = do.call("rbind", kMeans_cluster_score)
kMeans_cluster_score$class= gsub("_ND.*", "", kMeans_cluster_score$sample)

#Take the mean of biological replicates
kMeans_cluster_score$class = factor(kMeans_cluster_score$class, levels=c("CD4_unstim", "CD4_8hr", "CD4_24hr"))
kMeans_cluster_score_mean = kMeans_cluster_score %>% group_by(gene, cluster, class) %>% summarize(mean=mean(value), score = unique(score)) %>% ungroup()

#Prepare and write table with clusters and genes 
kMeans_cluster_score_mean_unique = kMeans_cluster_score_mean %>% mutate(gene_id = gsub(".*\\|", "", gene), gene_name= gsub("\\|.*", "", gene)) %>% select(gene_name, gene_id, cluster, score) %>% distinct()
z.filt.k_anno = left_join(z.filt, kMeans_cluster_score_mean_unique)
z.filt.k_anno = z.filt.k_anno  %>% relocate(gene_name, gene_id, cluster, score, CD4_unstim_ND543, CD4_unstim_ND581, CD4_unstim_ND589, CD4_8hr_ND543, CD4_8hr_ND581, CD4_8hr_ND589, CD4_24hr_ND543, CD4_24hr_ND581, CD4_24hr_ND589)
write.csv(z.filt.k_anno, file = "tables/Tcell_activation_kmeansCluster_geneList.csv", quote=F, row.names=F)

# Plot all genes in cluster, black line for centroid, color indicates how close to centroid
p2 <- ggplot(kMeans_cluster_score_mean , aes(x=class,y=mean)) + 
  geom_line(aes(color=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  #this adds the core 
  geom_line(data=Kmolten_rep_summary, aes(class,mean, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster Expression by Time",color = "Score")+  
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  facet_grid(rows = vars(cluster))
pdf("plots/diagnosis/rnaseq_kmeansClusters_cluster4_genes.pdf", useDingbats=FALSE)
p2
dev.off()

#Hierachal Clustering to check heatmap and similarity with kmeans
#make the matrix
dist <- cor(t(z.scaled.matrix), method="pearson")
#make the tree
hr <- hclust(as.dist(1-dist), method="complete") # Cluster rows by Pearson correlation.

#draw heatmap
png("plots/diagnosis/rnaseq_heatmap_hr_clustering.png" ,width=3.25,height=3.25,units="in",res=1200)
heatmap.2(dist,
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hr),
          scale="row",
          margins = c(2, 2),
          cexCol = 0.7,
          labRow = F,
          labCol = F,
          main = "Heatmap",
          trace = "none"
)
dev.off()

#Cut the tree to make 5 clusters
hclustk4 = cutree(hr, k=5)
library('dendextend')
Tree_hr = as.dendrogram(hr, method="complete")

#Plot the comparison between k-means and Hierachal Clustering 
pdf("plots/diagnosis/rnaseq_compare_kmean_hier_clustering.pdf", useDingbats=FALSE)
plot(Tree_hr,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")
bars <- cbind(hclustk4, kClusters)
colored_bars(bars, Tree_hr, sort_by_labels_order = T, y_shift=-0.1, rowLabels = c("Treecut",'K-means'),cex.rowLabels=0.7)
dev.off()


##Functionally annotate clusters with REACTOME enrichment

#Load the annotation files, each located in the directory
load("/mnt/isilon/sfgi/pahlm/annotationFiles/msigdb_v7.0_GMTs/c2.cp.reactome.v7.0.symbols.Rdata")

#Function for performing hypergeometric enrichment 
HG.test = function(grp,anno){
  universe = unique(unlist(anno)) #Set of genes that have been annotated, number of balls 
  gset = grp
  gset = gset[gset %in% universe]
  tmp = lapply(anno,function(j){
      p = length(j[j %in% gset]) #white balls drawn from an urn
      m = length(j) #white balls in the urn
      n = length(universe) - m #Number of black balls in the urn
      k = length(gset) #Number of balls drawn from the urn
     # expected = ceiling(length(j)/length(universe)*length(gset)) #Number of white balls expected to be drawn, rounding up (replace ceiling with floor to round down, or remove to leave the decimal)
      expected = ceiling(m/length(universe)*k)
      hg.pval = phyper(p,m,n,k,lower.tail=FALSE) #Calculate p value
      genes = paste(j[j %in% gset], collapse=":") #Combine the list of genes
      pathway_size = length(j)
      data.frame(p.val = hg.pval,genes= genes, observed=p, expected=expected, pathway_size) #Put the data together in a data.frame
    })
   out = do.call("rbind",tmp) #Merge the list of dataframes into a single dataframe
   out$Pathway = names(tmp) #Convert row.names to a column
   row.names(out) = NULL #Remove the row.names
   out$p.adjust = p.adjust(out$p.val,method="BH") #Calculate p.value after FDR correction
   out #Print results
}

#Run through clusters
for(i in 1:5){
cluster_number = i
genelist= unique(gsub("\\|ENSG.*", "", names(kClusters[kClusters==cluster_number])))\
#Run the test on my gene list and the reactome dataset
geneset.hg <- HG.test(genelist, c2.cp.reactome)
geneset.hg$Score =  -log10(geneset.hg$p.adjust)
#Sort the table by score
geneset.hg = geneset.hg[order(geneset.hg$Score,decreasing=TRUE),]
#This line makes sure they come up in the right order
geneset.hg$Pathway = factor(geneset.hg$Pathway, levels= geneset.hg$Pathway)
#Reorder columns
geneset.hg = geneset.hg[,c(6,3,4,5,1,7,8,2)]
write.csv(geneset.hg, file = paste0("tables/Reactome_Pathway_Tcell_activation_cluster_", cluster_number, ".csv"), quote=F, row.names=F)

#Plot reactome results
pdf(paste0("plots/Reactome_Pathway_Tcell_activation_cluster_", cluster_number, ".pdf"), useDingbats=FALSE)
plot = ggplot(geneset.hg[c(1:10),], aes(x= Pathway, y = Score))+
  geom_bar(stat='identity')+
  coord_flip()+
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_x_discrete(limits=rev)
print(plot)
dev.off()
}