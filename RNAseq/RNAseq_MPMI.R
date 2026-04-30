# Load the required libraries
library(limma)
library(edgeR)
library(qvalue)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# Load the input read count data "01_counttable_SUPP1331.txt" can be downloaded from the GitHub page.
data = read.delim(file.choose(), header = T, row.names = 1) #read "01_counttable_SUPP1331.txt"
data = data[1:32548, ]#Remove mitochondrial and chloroplast genes since poly A enrichment was performed during library preparation.
tail(data)

# TMM normalization
deg <- DGEList(data)
nf <- calcNormFactors(deg)  

# Exclude lowly expressed genes
c.sum = rowSums(data)    
s= c.sum[c.sum>=ncol(data)*10]
dat.mat = data[names(s),]
nrow(dat.mat)
#[1] 16899 

#Create design matrix
sample <- colnames(dat.mat) 
geno <- as.factor(paste(unlist(lapply(sample, function(x){strsplit(x, "_")[[1]][1]})))) 
design = model.matrix(~0 + geno) 
rownames(design) = colnames(dat.mat)
colnames(design) = levels(geno)
design

# voom transformation 
v <- voomWithQualityWeights(dat.mat, design, lib.size=colSums(data)*nf$samples$norm.factors, plot = T)
LogReadCounts = v$E
write.table(LogReadCounts, file = "LogReadCounts.txt", row.names = T, col.names = NA, quote = F, sep = '\t')

# fitting the linear model
fit = lmFit(v, design) 

# Write out a table containing the mean expression levels and standard errors for each gene
fitted_mean = fit$coef[,1:ncol(fit)]
write.table(fitted_mean, file="fitted_mean.txt", row.names=T, col.names=NA, sep="\t", quote=F)
fitted_se = fit$stdev.unscaled * fit$sigma 
write.table(fitted_se, file='fitted_se.txt', row.names=T, col.names=NA, sep="\t", quote=F)

# Specify comparisons and calculate q values
contrast.matrix1 <- makeContrasts(DC3000 - mock, SUPP1331 - mock, delta - mock, levels=design)
fit1 = contrasts.fit(fit, contrast.matrix1) 
eb.fit1 = eBayes(fit1)
write.table(eb.fit1$coef, file='ExpChange.txt', row.names=T, col.names=NA, sep="\t", quote=F)
q.val = qvalue(eb.fit1$p.val)$qvalues
write.table(q.val, file='q_val.txt', row.names=T, col.names=NA, sep="\t", quote=F)

# Extract differentially expressed genes based on > a 2-fold change in mean expression levels and q-value < 0.05
DEG_q0.05_2fold = subset(eb.fit1$coef, (q.val[, 1] < 0.05 & abs(eb.fit1$coef[, 1]) > 1)|(q.val[, 2] < 0.05 & abs(eb.fit1$coef[, 2]) > 1)|(q.val[, 3] < 0.05 & abs(eb.fit1$coef[, 3]) > 1)) 
nrow(DEG_q0.05_2fold)
write.table(DEG_q0.05_2fold,  file="Fig5A.txt", row.names = T, col.names = NA, sep = "\t", quote = F)
#Please note that this output table was used for generating heatmaps using cluster 3.0 (http://bonsai.hgc.jp/~mdehoon/software/cluster/software.htm)


# Select upregulated genes
DC3000_up = subset(eb.fit1$coef, q.val[, 1] < 0.05 & eb.fit1$coef[, 1] > 1)
SUPP_up = subset(eb.fit1$coef, q.val[, 2] < 0.05 & eb.fit1$coef[, 2] > 1)
delta_up = subset(eb.fit1$coef, q.val[, 3] < 0.05 & eb.fit1$coef[, 3] > 1)
nrow(DC3000_up)
#[1] 3077
nrow(SUPP_up)
#[1] 2364
nrow(delta_up)
#[1] 2568

# Select downregulated genes
DC3000_down = subset(eb.fit1$coef, q.val[, 1] < 0.05 & eb.fit1$coef[, 1] < -1)
SUPP_down = subset(eb.fit1$coef, q.val[, 2] < 0.05 & eb.fit1$coef[, 2] < -1)
delta_down = subset(eb.fit1$coef, q.val[, 3] < 0.05 & eb.fit1$coef[, 3] < -1)
nrow(DC3000_down)
#[1] 3236
nrow(SUPP_down)
#[1] 2052
nrow(delta_down)
#[1] 2360

#GO enrichment analysis for Fig. 5B
library(clusterProfiler)
library(org.At.tair.db)  
library(ggplot2)
library(dplyr)

#These clusters of genes were selected based on the expression patterns and are available from the GitHub page
cluster1 <- read.table(file.choose(), header = FALSE, stringsAsFactors = FALSE)[,1]
cluster2 <- read.table(file.choose(), header = FALSE, stringsAsFactors = FALSE)[,1]
cluster3 <- read.table(file.choose(), header = FALSE, stringsAsFactors = FALSE)[,1]
cluster4 <- read.table(file.choose(), header = FALSE, stringsAsFactors = FALSE)[,1]
gene_lists <- list(clusterI = cluster1, clusterII = cluster2, clusterIII = cluster3, clusterIV = cluster4)

cc <- compareCluster(geneCluster = gene_lists,
                     fun = "enrichGO",
                     OrgDb = org.At.tair.db,
                     keyType = "TAIR",
                     ont = "BP",
                     pvalueCutoff=0.05)

# simplify and extract top 5
cc_s <- simplify(cc, cutoff = 0.7, by = "p.adjust", select_fun = min)
cc_df <- as.data.frame(cc_s)
cc_top5 <- cc_df %>%
  group_by(Cluster) %>%
  arrange(p.adjust, pvalue, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

Fig5B = ggplot(cc_top5, aes(x = Cluster, 
                    y = Description,
                    color = -log10(p.adjust),
                    size = Count)) +
  geom_point() +  theme_bw() +  labs(x = NULL, y = NULL, size = "Gene count", color = "log.p" ) +
  theme(axis.text.y = element_text(size = 8))

ggsave(file = "Fig5B.png", plot = Fig5B, dpi = 300, width = 6, height = 5)


##SUPP1331 vs delta, more than 2 fold change
SUPP_vs_delta = subset(DEG_q0.05_2fold, abs(DEG_q0.05_2fold[ ,2] - DEG_q0.05_2fold[ ,3]) > 1)
nrow(SUPP_vs_delta)
#[1] 624
write.table(SUPP_vs_delta,  file="Fig5C.txt", row.names = T, col.names = NA, sep = "\t", quote = F)

#Gene set enrichment analysis
GSEA_SUPP_vs_delta = eb.fit1$coefficients[ ,2] - eb.fit1$coefficients[ ,3]
geneList <- sort(GSEA_SUPP_vs_delta, decreasing=TRUE)
gsea_go <- gseGO(
  geneList     = geneList,
  OrgDb        = org.At.tair.db,
  keyType      = "TAIR",
  ont          = "BP",
  pvalueCutoff = 0.1
)
#Please note that GSEA did not identify any statistically significantly enriched GO terms

#GO enrichment analysis for Fig. 5D
clusterA <- read.table(file.choose(), header = FALSE, stringsAsFactors = FALSE)[,1]
clusterB <- read.table(file.choose(), header = FALSE, stringsAsFactors = FALSE)[,1]
clusterC <- read.table(file.choose(), header = FALSE, stringsAsFactors = FALSE)[,1]
clusterD <- read.table(file.choose(), header = FALSE, stringsAsFactors = FALSE)[,1]
gene_lists <- list(clusterA = clusterA, clusterB = clusterB, clusterC = clusterC, clusterD = clusterD)

# GO enrichment against each gene list
cc <- compareCluster(geneCluster = gene_lists,
                     fun = "enrichGO",
                     OrgDb = org.At.tair.db,
                     keyType = "TAIR",
                     ont = "BP",
                     pvalueCutoff=0.1)

# simplify and extract top 5
cc_s <- simplify(cc, cutoff = 0.7, by = "p.adjust", select_fun = min)
cc_df <- as.data.frame(cc_s)

cc_top5 <- cc_df %>%
  group_by(Cluster) %>%
  arrange(p.adjust, pvalue, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

#dotplot
Fig5D = ggplot(cc_top5, aes(x = Cluster, 
                    y = Description,
                    size = Count,
                    color = -log10(p.adjust))) +
  geom_point() +  theme_bw() +  labs(x = NULL, y = NULL, size = "Gene count", color = "log.p") + theme(axis.text.y = element_text(size = 8))

ggsave(file = "Fig5D.png", plot = Fig5D, dpi = 300, width = 5, height = 5)
