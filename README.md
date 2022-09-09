# MetaPhOR
Metabolic Pathway Analysis of RNA

MetaPhOR was developed to enable users to assess metabolic dysregulation using transcriptomic-level data (RNA-sequencing and Microarray data) and produce publication-quality figures. A list of differentially expressed genes (DEGs), which includes fold change and p value, from DESeq2 [1] or limma [2], can be used as input, with sample size for MetaPhOR, and will produce a data frame of scores for each KEGG pathway. These scores represent the magnitude and direction of transcriptional change within the pathway, along with estimated p-values [3].MetaPhOR then uses these scores to visualize metabolic profiles within and between samples through a variety of mechanisms, including: bubble plots, heatmaps, and pathway models.  
