rm(list = ls())

## load data 
counts <- readRDS("TCGA-SKCM_rawCounts.RDS")
metadata <- readRDS("TCGA-SKCM_metadata.RDS")
table(metadata$pharmaceutical_therapy_type)
table(metadata$treatment_best_response)
identical(rownames(metadata), colnames(counts))


## ditribution of a sample
# base plot - histogram 
hist(counts[, 1])
hist(counts[, 1], xlim = c(0, 25000), ylim = c(0, 2500), breaks = 800)
# base plot - probability density plot 
plot(density(counts[, 1]),  xlim = c(0, 25000), main="sample density plot")

# ggplot2
library(ggplot2)
ggplot(counts) +
  geom_histogram(aes(x = counts[, 1]), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Frequency") 
  # xlim(0, 25000) +
  # ylim(0, 2000)
# ggplot2 - probability density plot
library(ggpubr)
ggdensity(counts[,1]) +
  xlim(0, 2500)
# 
# Cosa possiamo notare da questo plot:
#   La maggior parte dei geni ha un numero di conte basso 
# La distribuzione ha una lunga coda perchè non esiste nessun limite superiore
# La dinamica della distribuzione ha un range ampio 


## ditribution of a gene
# base plot
hist(as.numeric(counts[1,]), breaks = 50)
plot(density(as.numeric(counts[1,])), main="gene density plot")

# ggplot2
counts_2 <- as.data.frame(t(counts))
ggplot(counts_2) +
  geom_histogram(aes(x = as.numeric(counts_2[,1])), stat = "bin", bins = 50)
ggdensity(as.numeric(counts_2[,1])) +
  xlim(0, 2500)

## check for homoscedasticity
# plotting mean vs variance
mean_counts <- apply(counts, 1, mean)
variance_counts <- apply(counts, 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) 
  # geom_abline(intercept = 0, slope = 1, color="red")

# Con dati di RNAseq è vero che abbiamo un numero di geni molto alto e 
# la probabilità di prenderne uno in particolare è bassa, però succederà 
# sempre che i dati sono sovradispersi 
# questo accade perchè c'è comunque una bella variabilità tra un campione e 
# l'altro, anche per campioni che appartengono allo stesso fenotipo 


# boxplot 
my.gene.expression <- counts[2,]
my.gene.expression <- as.data.frame(t(my.gene.expression))

my.gene.expression$Response <- metadata$treatment_best_response

ggplot(my.gene.expression, aes(x = my.gene.expression$Response,
                               y = my.gene.expression[,1],
                               color = Response)) + 
  geom_point() +
  theme_bw() +
  xlab("Response") +
  ylab("Expression Level")
# idealmente dovremmo cercare di capire da dove viene questa variabilità 
# e aggiunsstare in base ad essa 


## normalise data 
# Il fatto di lavorare con una distribuzione non normale comunque per noi sarà
# sempre un problema perché molte analisi statistiche assunomo normalità
# e omoschedasticità. 
# Ci sono situazioni in cui non possiamo usare una distribuzione negativa 
# binomiale, per esempio se vogliamo fare una correllazione di pearson o una 
# PCA abbiamo bisogno di una distribuzione normale. 

sizeFactors <- colSums(counts)  
sizeFactors <- sizeFactors / mean(sizeFactors)  # normalise
normCounts <- counts / sizeFactors
counts_log2 <- log2(normCounts + 1)
library(DESeq2)
counts_m <- as.matrix(counts)
# counts_rlog <- rlog(counts_m)
counts_rlog <- readRDS("TCGA-SKCM_counts_rlog.RDS")
counts_vst <- vst(counts_m)

# op <- par(mfrow=c(1,4))
plot(density(as.numeric(counts[1,])), main="raw counts", cex.main=2)
plot(density(as.numeric(counts_log2[1,])), main="log2(x+1)", cex.main=2)
plot(density(as.numeric(counts_rlog[1,])), main="rlog", cex.main=2)
plot(density(as.numeric(counts_vst[1,])), main="vst", cex.main=2)
# par(op)

# library(ggpubr)
# ggdensity(counts[,1])  
# ggdensity(counts_vst[1,])  
# ggdensity(counts_log2[,1])  

## check for homoscedasticity
# plotting mean vs variance
# mean_counts_log2 <- apply(counts_log2, 1, mean)
# variance_counts_log2 <- apply(counts_log2, 1, var)
# df_log2 <- data.frame(mean_counts_log2, variance_counts_log2)
# 
# mean_counts_rlog <- apply(counts_rlog, 1, mean)
# variance_counts_rlog <- apply(counts_rlog, 1, var)
# df_rlog <- data.frame(mean_counts_rlog, variance_counts_rlog)
# 
# mean_counts_vst <- apply(counts_vst, 1, mean)
# variance_counts_vst <- apply(counts_vst, 1, var)
# df_vst <- data.frame(mean_counts_vst, variance_counts_vst)
# 
# 
# ggplot(df_log2) +
#   geom_point(aes(x=mean_counts_log2, y=variance_counts_log2)) + 
# geom_abline(intercept = 0, slope = 1, color="red")
# 
# ggplot(df_rlog) +
#   geom_point(aes(x=mean_counts_rlog, y=variance_counts_rlog)) + 
#   geom_abline(intercept = 0, slope = 1, color="red")
# 
# ggplot(df_vst) +
#   geom_point(aes(x=mean_counts_vst, y=variance_counts_vst)) + 
#   geom_abline(intercept = 0, slope = 1, color="red")
# 

## PCA
counts_vst <- counts_vst[which(apply(counts_vst, 1, var) != 0), ]

p <- PCAtools::pca(counts_vst, metadata = metadata)

biplot(p, x="PC1", y="PC2", lab = NULL,
       # colkey = c(Male='gold2', Female='purple'),
       colby = 'treatment_best_response')



## DEGs
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ treatment_best_response)
dds <- DESeq(dds)
deseq2_res <- results(dds)
deseq2_res <- as.data.frame(deseq2_res)

library(easylabel)
easyVolcano(deseq2_res)

## pathway analysis 
library(enrichR)
downregulated_genes <- subset(deseq2_res, deseq2_res$padj < 0.05)
downregulated_genes <- subset(deseq2_res, downregulated_genes$log2FoldChange < 0)
downregulated_genes <- rownames(downregulated_genes)
head(downregulated_genes)

dbs <- listEnrichrDbs()

downregulated_pathways <- enrichr(genes = downregulated_genes, 
                                  databases = c("Reactome_2016",
                                                 "GO_Biological_Process_2021",
                                                 "GO_Molecular_Function_2021",
                                                 "GO_Cellular_Component_2021",
                                                 "OMIM_Disease"))


## DEG-Gs
library(DEGGs)
subnetworks_object <- generate_subnetworks(normalised_counts = as.data.frame(counts_vst), 
                                           metadata = metadata,
                                           subgroup_variable = "treatment_best_response",
                                           use_qvalues = FALSE) 
View_interactive_subnetwork(subnetworks_object)


## cemitool 
library(CEMiTool)
sample_annot <- data.frame("Class" = metadata$treatment_best_response,
                           "SampleName" = rownames(metadata))

cem <- cemitool(as.data.frame(counts_vst), sample_annot, verbose=TRUE)
cem <- mod_gsea(cem)
cem <- plot_gsea(cem)
cem <- plot_profile(cem)

gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)
cem <- mod_ora(cem, gmt_in)
cem <- plot_ora(cem)

int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
interactions_data(cem) <- int_df # add interactions
cem <- plot_interactions(cem) 


generate_report(cem, title = "cemitool_report", directory="./", force = TRUE)
