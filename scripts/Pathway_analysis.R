##We would be performing pathway enrichment analysis for the res_filtered_.2 dataset made previously
#We would be performing two pathway analysis, namely GO and KEGG enrichment analysis
#Please note that the res_filtered_.2 dataset we are using for this has been made with a lenient padj cutoff of 0.2
#All the plots generated in this code are saved in the results folder
#Please refer to README.md for instructions to this code

#Loading the required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
#converting to Entrez ID
gene_entrez <- bitr(res_filtered_.2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

#GO analysis(Note that here we have p value cut off as 0.1 and padj as 0.2 to check for biologically signifact enriched pathways)
go_result <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.2,
  readable = TRUE
)
#Visualization of GO analysis
dotplot(go_result, showCategory = 35) + ggtitle("GO Enrichment Analysis")
barplot(go_result, showCategory = 35)+ ggtitle ("Top GO Terms")
#Saving
write.csv(as.data.frame(go_result), file = "data/GO_results.csv", row.names = TRUE)

#KEGG enrichment analysis
kegg_result <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.2
)
#Visualization
dotplot(kegg_result, showCategory = 20) + ggtitle("KEGG Pathway Analysis")
barplot(kegg_result, showCategory = 20) + ggtitle("Top KEGG Pathways")
#Saving
write.csv(as.data.frame(kegg_result), file = "data/KEGG_results.csv", row.names = TRUE)
#saving the visualizations
ggsave("results/GO_barplot.png", (barplot(go_result, showCategory = 35)+ ggtitle ("Top GO Terms")), width= 15, height = 15, dpi =300)
ggsave("results/KEGG_barplot.png",(barplot(kegg_result, showCategory = 20) + ggtitle("Top KEGG Pathways")), width = 12, height =12, dpi = 300)
