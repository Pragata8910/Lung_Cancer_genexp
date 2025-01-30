###Visualizing the processed result data
##We also visualise some previously mentioned special genes in various conditions.
#Please refer to README.md for running instructions and modification of the functions and plots for different cases.
#All the plots generated here are saved in the results folder.
#Loading the libraries.
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(pheatmap)

#Generating the MA plot & saving it
pdf("MA_plot.pdf", width = 8, height = 6)
plotMA(result, ylim= c(-5,5), main= "main MA Plot")  
abline(h= c(-1,1), col = "blue", lty= 2)
dev.off()

#Generating the Volcano Plot
volcano_data <- res_annotated
volcano_data$significant <- ifelse(volcano_data$padj <0.2, "Significant","Not Significant")

volcano <- ggplot(volcano_data, aes(x = log2FoldChange, y= -log10(padj), color= significant)) +
  geom_point(alpha= 0.7, size= 1.5)+
  scale_color_manual(values = c("Significant"= "red", "Not Significant"= "grey"))+
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y= "-Log10(padj)")+
  theme_minimal()+
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"), 
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90")
  )
#Saving the volcano plot
ggsave("Volcano_plot_.2.png", volcano, width = 12, height = 12,dpi = 300)

#Generating the heatmap
#Taking normalised count for significant genes
significant_genes <- res_significant$Row.names
norm_counts <- counts(dds, normalized = TRUE)[significant_genes, ]

heatmap <- pheatmap(norm_counts,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_rownames = TRUE,
         main = "Heatmap of Significant Genes")
#Saving the heatmap
ggsave("Heatmap_significant_genes.pdf", heatmap, height = 10, width =12)

##Code to make a function to visualize some previously specified genes in various conditions
box_gene_expression <- function(normalized_counts, genes, condition){
  sig_counts <- normalized_counts[genes, ]
  sig_counts_df <- as.data.frame(sig_counts)
  sig_counts_df$Gene <- rownames(sig_counts_df)
  sig_counts_long <- melt(sig_counts_df, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
  colnames(sig_counts_long) <- c("Gene","Sample","Expression")
  metadata_new$Sample <- rownames(metadata_new)
  metadata_subset <- metadata_new[,c("Sample",condition)]
  colnames(metadata_subset) <- c("Sample","Condition")
  
  sig_counts_long <- merge(sig_counts_long, metadata_subset, by = "Sample")
  
  ggplot(sig_counts_long, aes(x = condition, y = Expression, fill = condition))+
    geom_boxplot(alpha = 0.7)+
    facet_wrap(~ Gene, scales = "free_y")+
    scale_y_log10()+
    theme_minimal()+
    labs(title = paste("Expression of Significant Genes by:", condition),
         x= condition, y= "Normalised Expression (log scale)")+
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"), 
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_line(color = "gray90")
    )
}
#saving the boxplot
sig_boxplots <- box_gene_expression(normalized_counts, genes_of_interest, "Smoking")
ggsave("results/special_genes_box_smoking.png", sig_boxplots, height = 10, width = 12, dpi = 500)

#barplot
bar_gene_expression <- function(normalized_counts, genes, condition){
    sig_counts <- normalized_counts[genes, ]
    sig_counts_df <- as.data.frame(sig_counts)
    sig_counts_df$Gene <- rownames(sig_counts_df)
    sig_counts_long <- melt(sig_counts_df, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
    colnames(sig_counts_long) <- c("Gene","Sample","Expression")
    metadata_new$Sample <- rownames(metadata_new)
    metadata_subset <- metadata_new[,c("Sample",condition)]
    colnames(metadata_subset) <- c("Sample","Condition")
    
    sig_counts_long <- merge(sig_counts_long, metadata_subset, by = "Sample")
    #mean and standard error
    sig_counts_summary <- sig_counts_long %>%
      group_by(Gene, Condition) %>%
      summarise(
        mean_expr = mean(Expression),
        se = sd(Expression)/sqrt(n()),
        .groups = 'drop'
      )
    #error bars
    ggplot(sig_counts_summary, aes(x = Condition, y = mean_expr, fill = Condition)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
      geom_errorbar(aes(ymin = mean_expr - se, ymax = mean_expr + se),
                    width = 0.2, position = position_dodge(0.9)) +
      facet_wrap(~ Gene, scales = "free_y") +
      scale_y_log10() +
      theme_minimal() +
      labs(title = paste("Mean Expression of Significant Genes by:", condition),
           x = condition, 
           y = "Mean Normalised Expression (log scale)")+
      theme(
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"), 
        panel.grid.major = element_line(color = "gray80"),
        panel.grid.minor = element_line(color = "gray90")
      )
  }
#saving
sig_barplot <- bar_gene_expression(normalized_counts, genes_of_interest, "histology")
ggsave("results/special_genes_bar_histo.png", sig_barplot, height = 10, width = 12, dpi = 500)
