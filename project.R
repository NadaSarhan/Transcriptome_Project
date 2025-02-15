# download libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

################################################

#import library
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(biomaRt)

################################################

trans_samples = read.delim("~/NGS/project/all_samples/alignment_samples/count_matrix.tsv")

# Remove version numbers if they exist (e.g., ENST00000367770.8 → ENST00000367770)
trans_samples$target_id <- sub("\\..*", "", trans_samples$target_id)

# Connect to Ensembl database
# For human genes
mart <- useMart("ensembl", 
                dataset="hsapiens_gene_ensembl", 
                host="https://useast.ensembl.org")  # US East mirror

# Check whether target_id values are transcript or gene IDs
head(trans_samples$target_id)  # Print a few IDs to check

# Convert transcript ID to gene name
conversion <- getBM(attributes=c("ensembl_transcript_id", "external_gene_name"),
                    filters="ensembl_transcript_id",
                    values=trans_samples$target_id,
                    mart=mart)

# Merge gene names with the count matrix
new_matrix <- merge(conversion, 
                      trans_samples, 
                      by.x="ensembl_transcript_id", 
                      by.y="target_id", 
                     all.y=TRUE)

# Save updated count matrix
write.table(new_matrix, 
            "~/NGS/project/all_samples/alignment_samples/count_matrix_with_gene_names.tsv", 
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE)

################################################

# Remove rows with any NA values
count_matrix <- na.omit(new_matrix)

# Remove duplicate rows (genes_id)
count_matrix = count_matrix[!duplicated(count_matrix$external_gene_name), ]# Now set row names safely

rownames(count_matrix) = count_matrix$external_gene_name
dim(trans_samples)

count_matrix$ensembl_transcript_id = NULL
count_matrix$external_gene_name = NULL
dim(count_matrix)

# Save updated count matrix
write.table(count_matrix, 
            "~/NGS/project/all_samples/alignment_samples/count_matrix_with_gene_names.tsv", 
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE)

ncbi_matrix = read.delim("~/NGS/project/all_samples/alignment_samples/count_matrix_ncbi_geo/GSE223598_raw_counts_GRCh38.p13_NCBI.tsv", 
                         row.names = 1)
dim(ncbi_matrix)
dim(count_matrix)
################################################

# obtain data into matrix
count_matrix = as.matrix(count_matrix)
genes = row.names(count_matrix)

# convert each value to integer number
count_matrix = apply(count_matrix, 2, as.integer)
row.names(count_matrix) = genes

dim(count_matrix)
str(count_matrix)  # Check structure
class(count_matrix)  # Should return "matrix"
head(count_matrix)  # Preview first few rows

# load metadata
trans_meta = read.csv("~/NGS/project/all_samples/normalization_pathway/SraRunTable.csv",
                row.names = 1)
dim(trans_meta)

#Filter_other_virus_Cases
trans_meta = trans_meta %>% select(-Assay.Type, -AvgSpotLen,
                                   -Bases, -BioProject,
                                   -BioSample, -Bytes,
                                   -cell_line, -cell_type, 
                                   -Center.Name, -Collection_Date,
                                   -Consent, -DATASTORE.filetype,
                                   -DATASTORE.provider, -DATASTORE.region,
                                   -Experiment, -geo_loc_name_country, 
                                   -geo_loc_name_country_continent, -geo_loc_name,
                                   -Instrument, -Library.Name, 
                                   -LibraryLayout, -LibrarySelection,
                                   -LibrarySource, -Organism, 
                                   -Platform, -ReleaseDate, 
                                   -create_date, -version, 
                                   -Sample.Name, -source_name, 
                                   -SRA.Study)

data= subset(as.matrix(count_matrix), 
             select= as.character(rownames(trans_meta)))

dim(count_matrix)
dim(trans_meta)

str(count_matrix)
str(trans_meta)

trans_meta$genotype <- gsub(" ", "_", trans_meta$genotype)  # Replace spaces with "_"
trans_meta$genotype <- as.factor(trans_meta$genotype)  # Ensure it's a factor
# Check if levels are now clean
levels(trans_meta$genotype)

trans_meta$treatment <- gsub(" ", "_", trans_meta$treatment)  # Replace spaces with "_"
trans_meta$treatment <- as.factor(trans_meta$treatment)  # Ensure it's a factor
# Check if levels are now clean
levels(trans_meta$treatment)

# create DESeq2 object (dds_object)
deseq2Object = DESeqDataSetFromMatrix(countData = count_matrix, 
                                      colData = trans_meta, 
                                      design = ~ genotype + treatment)
# run DESeq2 (dds)
deseq2Run = DESeq(deseq2Object)

# differentially expressed genes (DEGs) and how many genes were 
# filtered in the DESeq2 analysis.
summary(results(deseq2Run))

# Before running results(), check that your conditions are properly recognized:
levels(trans_meta$genotype)   # Check genotype levels
levels(trans_meta$treatment)   # Check treatment levels

# differential analysis statistics
res = results(deseq2Run, 
              contrast = c("genotype", "Parental", "Venetoclax_resistant"))
res = results(deseq2Run, 
              contrast = c("treatment", "Decitabine", "DMSO"))

# remove NaN
res = as.data.frame(res[complete.cases(res), ])

# check the overall results
summary(res)
nrow(res)
nrow(res[res$padj < 0.05, ])  # Number of significant DEG

# add Thresholds:
# This filters the results from res (DESeq2 output) to 
# extract significant differentially expressed genes (DEGs) based on 
# statistical and biological significance.
deseq_deg = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
dim(deseq_deg)
summary(deseq_deg)

# save my results for further analysis:
write.csv(res, "~/NGS/project/all_samples/normalization_pathway/DESeq2_results.csv")
View(read.csv("~/NGS/project/all_samples/normalization_pathway/DESeq2_results.csv"))
write.csv(deseq_deg, "~/NGS/project/all_samples/normalization_pathway/DESeq2_significant_DEGs.csv")
View(read.csv("~/NGS/project/all_samples/normalization_pathway/DESeq2_significant_DEGs.csv"))

################################################

# create normal transformation
ntd = normTransform(deseq2Run)
# Normalize counts using variance stabilizing transformation (recommended)
vsd = vst(deseq2Run, blind = FALSE)

# extrace sigDif genes
significantGeneNames = rownames(deseq_deg)

################################################

# extract test matrix to plot heatmap
testMatrix = as.matrix(assay(ntd))[significantGeneNames[1:200], ]

# Make sure that significantGeneNames is not empty and contains at least 15 genes:
length(significantGeneNames)
head(significantGeneNames)
# Ensure that ntd (normalized transformation) exists and contains data:
dim(ntd)
head(assay(ntd))
# check testMatrix:
dim(testMatrix)
testMatrix
# If testMatrix contains NA values, it can cause pheatmap() to fail.
# Check for NA:
sum(is.na(testMatrix))

dev.off()  # Reset the plotting device

# scale = "row" normalizes each row independently (Z-score transformation).
# This helps in comparing gene expression patterns across samples by standardizing values.
# Without scale = "row", values remain in their original scale.
# This makes it easier to compare up/down-regulated genes, regardless of absolute expression levels.
# If you want to keep raw expression values, remove scale = "row" (like in the second version).
pheatmap(testMatrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         annotation_col = trans_meta, 
         scale = "row",
         plot = TRUE)

################################################

# Extract PCA data
pca_data = plotPCA(ntd, 
                   intgroup = c("treatment", "genotype"), 
                   ntop=500, 
                   returnData = TRUE)

# Create the PCA plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = treatment, shape = genotype)) +
  geom_point(size = 4) +  # Adjust point size as needed
  theme_minimal() +
  labs(
    x = paste0("PC1: ", round(attr(pca_data, "percentVar")[1] * 100, 1), "% variance"),
    y = paste0("PC2: ", round(attr(pca_data, "percentVar")[2] * 100, 1), "% variance")
  ) +
  scale_color_manual(values = c("red", "blue")) +  # Customize colors for treatment
  scale_shape_manual(values = c(16, 17))           # Customize shapes for genotype

################################################

# Convert DESeq2 results to data frame
res$gene = rownames(res)

# Volcano plot
ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme(legend.position = "none")

################################################

# next step is to Perform pathway analysis (KEGG, Reactome, GO) to find affected pathways.

# load packages
library(clusterProfiler)
library(org.Hs.eg.db)  # Database for human genes
library(igraph)
library(qvalue)
library(GO.db)
library(GOSemSim)
library(gson)
library(AnnotationDbi)

# 1- Step 1: Prepare Your Gene List for Pathway Analysis

# Convert Gene Symbols to Entrez IDs
# Extract gene symbols from DESeq2 results
conversion <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                    filters = "hgnc_symbol",
                    values = rownames(deseq_deg),  # Significant DEGs
                    mart = mart)

# Merge Entrez IDs with DESeq2 results
deseq_deg <- merge(deseq_deg, conversion, by.x = "row.names", by.y = "hgnc_symbol", all.x = TRUE)
colnames(deseq_deg)[1] <- "Gene_Symbol"

dim(deseq_deg)

##########
# 2- Perform KEGG Pathway Enrichment Analysis
# KEGG (Kyoto Encyclopedia of Genes and Genomes) helps identify which molecular pathways are significantly affected in your data.

# Extract Entrez IDs
entrez_ids <- deseq_deg$entrezgene_id

# Run KEGG Pathway Analysis
kegg_results <- enrichKEGG(gene = entrez_ids, 
                           organism = "hsa",  # "hsa" is for human pathways
                           keyType = "ncbi-geneid",
                           pvalueCutoff = 0.05)  # Adjusted p-value threshold
# This identifies KEGG pathways significantly affected in your dataset.
# View results
head(kegg_results)

# 3- Visualize KEGG Pathways
# This helps you identify affected pathways such as cancer-related signaling pathways.
# Dot plot of KEGG results
dotplot(kegg_results, showCategory = 20, title = "Top 20 Enriched KEGG Pathways")
# Bar plot of KEGG results
barplot(kegg_results, showCategory = 20, title = "KEGG Pathway Enrichment")
# View a Specific KEGG Pathway in a Web Browser
# This will open the KEGG pathway diagram showing how genes are involved.
browseKEGG(kegg_results, "hsa05200")  # Example: KEGG pathway for cancer

#######################################################################

# Perform Reactome Pathway Analysis
# Functional Enrichment Using ReactomePA (Pathway-Specific)
# ReactomePA is best for pathway enrichment analysis, especially for cancer-related pathways.
# Install & Load ReactomePA
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ReactomePA")
library(ReactomePA)
library(aplot)

# Run Reactome Pathway Analysis
reactome_results <- enrichPathway(gene = entrez_ids, 
                                  organism = "human",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05)
# View results
head(reactome_results)

# Visualize Reactome Pathway Results
# Reactome helps in understanding disease-specific biological interactions.
# Dot Plot
dotplot(reactome_results, showCategory = 20, title = "Top 20 Reactome Pathways")
# Bar Plot
barplot(reactome_results, showCategory = 20, title = "Reactome Pathway Enrichment")

#######################################################################

# Perform Gene Ontology (GO) Analysis
# Gene Ontology (GO) analysis helps understand biological processes (BP), 
# molecular functions (MF), and cellular components (CC) affected in your data.

# Run GO Enrichment Analysis
go_results <- enrichGO(gene = entrez_ids,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",  # "BP" = Biological Process, "MF" = Molecular Function, "CC" = Cellular Component
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05)
# View results
head(go_results)
# Visualize GO Results
# GO analysis tells you which biological functions are disrupted in your data.
dotplot(go_results, showCategory = 20, title = "Top 20 Enriched GO Biological Processes")

################################
