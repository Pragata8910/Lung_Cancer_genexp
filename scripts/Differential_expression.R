##Data processing and model fitting
#*Do note that sig_genes_0.05 df has genes with strict padj cutoff of <0.05
#*res_filtered_.2 contains genes with lenient cutoff of padj <0.2
#Note that all the data extracted or generated are stored in data folder

##Required libraries
library(DESeq2)
library(dplyr)
library(biomaRt)
library(reshape2)
#creating a Deseq data-frame
dds <- DESeqDataSetFromMatrix(countData = count_raw,
                              colData = metadata_new,
                              design = ~ condition)
#preferred pre_filtering with a sum threshold of 10
keep <- rowSums(counts(dds))>= 10
dds <- dds [keep,]
#setting the factor level
dds$condition <- relevel(dds$condition, ref = "normal")

##estimating size factor and fitting the best model
dds <-  DESeq(dds)
result <- results(dds)
summary(result)

#Annotating the ensembl id of the genes using Ensembl mirror(Note that the database is according to ch19(ch37))
mart <- useMart("ENSEMBL_MART_ENSEMBL", 
                host = "https://grch37.ensembl.org",  
                dataset = "hsapiens_gene_ensembl")
annotations <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
                     filters = "ensembl_gene_id",
                     values = rownames(result),
                     mart = mart)
#merging
res_annotated <- merge(as.data.frame(result), annotations, by.x = "row.names", by.y = "ensembl_gene_id")
head(res_annotated)

# Filter for genes with padj < 0.05
sig_genes_0.05 <- res_annotated[which(res_annotated$padj < 0.05 & !is.na(result$padj)),]
#saving the results
write.csv(as.data.frame(sig_genes_0.05),"data/significant_genes_.05.csv", row.names = FALSE)

#Filtering with padj<0.2 and extracting significant genes
filtered_.2 <- res_annotated[!is.na(res_annotated$padj) & res_annotated$padj < 0.2,]
res_filtered_.2 <- filtered_.2$hgnc_symbol
#Applying additional filter for lof2FC >0.5 or < -0.5
res_significant <- filtered_.2[!is.na(filtered_.2$log2FoldChange) &(filtered_.2$log2FoldChange > 0.5 | filtered_.2$log2FoldChange < -0.5),]
#saving the results
write.csv(as.data.frame(res_significant),"data/sig.results_log2fc&pval_.2.csv", row.names = FALSE)


#Checking for the expression of key genes of interest
genes_of_interest <- c("EGFR","RB1","ALK","ROS1","BRAF","STK11","NF1","BRCA1","BRCA2","MET","KRAS","MYC","TP53","ATM","SOX2","NOTCH1","CDKN2A","PIK3CA","CD274","MMP9","CYP1B1","FOXP3","PDCD1","CTLA4","RAD51")
key_results <- res_annotated[res_annotated$hgnc_symbol %in% genes_of_interest,]
#Generating the normalized count matrix
normalized_counts <- counts(dds, normalized =TRUE)
#placing gene names in counts matrix
matched_genes <-  match(rownames(normalized_counts), res_annotated$Row.names)
valid_indices <-  !is.na(matched_genes)
rownames(normalized_counts)[valid_indices] <- res_annotated$hgnc_symbol[matched_genes[valid_indices]]
#extracting normalized counts for these genes
specific_gene_counts <-  normalized_counts[genes_of_interest,]

