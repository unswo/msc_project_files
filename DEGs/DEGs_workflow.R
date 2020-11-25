#!/usr/bin/Rscript
# Differential expression analysis
# Adrienne Unsworth 190033848, 2020
#
# Perform standard analysis from transcript abundance estimates obtained from Stringtie
# Could easily be adapted for use with salmon by altering tximport command or raw counts from HTseq/Featurecounts.
#
# Import libraries
library(biomaRt)
library(GenomicFeatures)
library(tximport)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(fgsea)
library(genefilter)
library(pheatmap)
library(apeglm)
library(gprofiler2)
library(EnhancedVolcano)
# Define useful function(s)
`%notin%` <- Negate(`%in%`)
# Set working directory
setwd("path/to/transcript/files")

# Define biomaRt marts that will be used later on
# Full mart, other marts will be drived from this
ensembl_hs_mart <-
  useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_df <-
  getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
        mart = ensembl_hs_mart)

# data fram of ncRNAs, mtRNAs, rRNAs and pseudogenes
ensembl_ncrna <-
  getBM(
    attributes = c("transcript_biotype", "external_gene_name"),
    mart = ensembl_hs_mart,
    filters = 'biotype',
    values = c(
      'lncRNA',
      'rRNA',
      'Mt_rRNA',
      'Mt_tRNA',
      'miRNA',
      'rRNA_pseudogene',
      'pseudogene'
    )
  )

# define protein coding genes
ensembl_coding <- getBM(
  attributes = c("transcript_biotype", "external_gene_name"),
  mart = ensembl_hs_mart,
  filters = 'biotype',
  values = c('protein_coding')
)

# example of filtering mart for specific function, ie as defined by a GO term
# ensembl_fusion <-
#   getBM(
#     attributes = c("ensembl_gene_id", "go_id", 'hgnc_symbol'),
#     filters = 'go',
#       values = 'GO:0003700',
#     mart = ensembl_hs_mart
#   )
# TF <- unique(ensembl_fusion$hgnc_symbol)

# another example of creating a dataframe of ensembl IDs from hgnc symbols
# ensembl_splice <-
#   getBM(
#     attributes = c("ensembl_gene_id", 'hgnc_symbol'),
#     filters = 'hgnc_symbol',
#     values = c(
#       'VEGFA',
#       'KLF6',
#       'BCL2L2',
#       'FGFR2',
#       'TMPRSS2',
#       'ERG',
#       'AR',
#       'PTEN',
#       'TP53',
#       'RB1',
#       'BRCA1',
#       'BRCA2'
#     ),
#     mart = ensembl_hs_mart
#   )


# Import samples; _t_.ctab is specific to stringtie files
samples_stringtie <-
  list.files(pattern = '_t_data.ctab')
# project data acquired from ENA website, uncomment to download
#download.file('https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJNA411786&result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,submitted_ftp,submitted_galaxy,sra_ftp,sra_galaxy,cram_index_ftp,cram_index_galaxy&download=txt', destfile = 'PRJNA411786.txt')
#download.file('https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJNA552058&result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,submitted_ftp,submitted_galaxy,sra_ftp,sra_galaxy,cram_index_ftp,cram_index_galaxy&download=txt', destfile = 'PRJNA552058.txt')
samples_PRJNA411786 <-
  read.delim(('PRJNA411786.txt'))
samples_PRJNA552058 <-
  read.delim(('PRJNA552058.txt'))
samples_info <-
  rbind(samples_PRJNA411786, samples_PRJNA552058)

# #remove extra info - ftp links etc so we're only interested in col 1-9
samples_info <- samples_info[, 1:9]
names(samples_stringtie) <-
  samples_info$run_accession

# pull annotation from ctab files
tmp <-
  read_tsv(samples_stringtie[1], col_types = cols())
tx2gene <- tmp[, c("t_name", "gene_name")]
txi.genes.stringtie <-
  tximport(samples_stringtie, type = 'stringtie', tx2gene = tx2gene)

# does.. something? seems superfluous
sampleTable <- read.csv('proj_design.csv')
sampleTable <- sampleTable[,-1]
sampleTable$Patient <-
  as.factor(sampleTable$Patient)
rownames(sampleTable) <-
  colnames(txi.genes.stringtie$counts)

# DESeq2 - more or less standard DESeq workflow
# ~ BioProject + Condition is the overall design of the experiment
# The first parameter is the batch effect
dds <-
  DESeqDataSetFromTximport(txi.genes.stringtie, sampleTable, ~ BioProject + Condition)

dds_2 <- DESeq(dds)
res <- results(dds_2)
# log2fold shrink
res <-
  lfcShrink(dds_2,
            coef = "Condition_T_vs_N",
            res = res,
            type = 'apeglm')

# estimate dispersion
plotDispEsts(dds_2,
             ylim = c(1e-10, 1e2),
             xlab = 'Mean of normalised counts',
             ylab = 'Dispersion')


res_df <- as.data.frame(res)

res_df$ensembl_gene_id <- rownames(res_df)

# Not needed unless need to manually annotate external gene names ie counts/abundances from another tool
# res_df <-
#   merge.data.frame(
#     res_df,
#     ensembl_df,
#     by = intersect(colnames(res_df), colnames(ensembl_df)),
#     all.x = TRUE,
#     all.y = FALSE
#   )

res_df2 <-
  filter(res_df,!is.na(res_df$padj))
gseaInput <-
  filter(res_df2,!is.na(rownames(res_df2)))
ranks <- gseaInput$log2FoldChange
names(ranks) <- gseaInput$ensembl_gene_id
ranks <- sort(ranks, decreasing = TRUE)

# Range of gmt files for different enrichment analyses
#
# gmt_list <-
#   gmtPathways('h.all.v7.1.symbols.gmt')
# kegg <-
#   gmtPathways('c2.cp.kegg.v7.1.symbols.gmt')
# GO_bp <-
#   gmtPathways('c5.bp.v7.1.symbols.gmt')
gseaRes <-
  fgsea(gmt_list, ranks, nperm = 1000)
gseaResTidy <-
  gseaRes %>% as_tibble() %>% arrange(desc(NES))
topPathwaysUp <-
  gseaRes[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown <-
  gseaRes[ES < 0][head(order(pval), n = 10), pathway]
topPathways <-
  c(topPathwaysUp, rev(topPathwaysDown))

# Different option for GSEA visuals
#plot.new()
#plotGseaTable(gmt_list[topPathways], ranks, gseaRes,
#              gseaParam = 0.5)

# Define top DE genes for heatmap
rownames(res_df) <-
  res_df$external_gene_name
topGenes <- arrange(res_df, padj)
topGenes <-
  filter(topGenes, padj <= 0.01, abs(log2FoldChange) >= 2)
rownames(topGenes) <-
  topGenes$ensembl_gene_id
topGenes <- rownames(topGenes)

#GO overrepresentation analysis
gostres <-
  gost(
    topGenes[1:5000],
    organism = 'hsapiens',
    sources = c('GO', 'KEGG'),
    evcodes = TRUE
  )
gostplot(gostres)

counts_from_dds <- counts(dds_2)

# Pinched from:
# DESeq results to pathways in 60 Seconds with the fgsea package
# Stephen Turner
ggplot(gseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.05)) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       title = "Hallmark pathways NES from GSEA") +
  theme_minimal()

# Filter for protein coding only
dds_2 <-
  dds_2[rownames(dds_2) %in% ensembl_coding$external_gene_name, ]
res <- results(dds_2)
res <-
  lfcShrink(dds_2,
            coef = "Condition_T_vs_N",
            res = res,
            type = 'apeglm')
plotMA(res, alpha = 0.01, ylim = c(-5, 5))
res_df <- as.data.frame(res)


DE_genes <- subset(res, res$padj < 0.01)
DE_genes <-
  DE_genes[order(abs(DE_genes$log2FoldChange), decreasing = TRUE), ]

DE_genes <-
  DE_genes[rownames(DE_genes) %in% ensembl_coding$external_gene_name, ]
DE_genes <- as.data.frame(DE_genes)

# write to csv
#write.csv2(DE_genes, file = 'DE_genes.csv')

res_df <- as.data.frame(DE_genes)
res_df$ensembl_gene_id <- rownames(res_df)
res_df2 <-
  filter(res_df,!is.na(res_df$padj))

#plot gene counts for gene with smallest padj
rownames(counts_interest) <- NULL
p <- ggplot(counts_interest,
            aes(
              x = Patient,
              y = count,
              color = Condition,
              group = Condition
            )) +
  geom_point() + stat_summary(fun = mean, geom = "line") +
  scale_y_log10()
p <-
  p + facet_wrap(~ gene, scales = 'fixed') + labs(x = 'Gene', y = 'Normalised counts')
p

# QC data
rlog_data_blind <-
  rlogTransformation(dds_2, blind = TRUE)

# Euclidean distances of samples
dist_rl <- dist(t(assay(rlog_data_blind)))
dist_rl <- as.matrix(dist_rl)
heatmap(dist_rl)

# PCA analysis
plotPCA(rlog_data_blind, intgroup = c('Condition', 'BioProject'))

top_genes <- head(order(DE_genes$padj), 20)

# gene counts
rownames(res_df) <-
  res_df$external_gene_name
topGenes <- arrange(res_df, padj)
topGenes <-
  filter(topGenes, padj <= 0.01, abs(log2FoldChange) >= 2)
rownames(topGenes) <-
  topGenes$ensembl_gene_id
topGenes <- rownames(topGenes)
plotcounts <- NULL
counts_df <- NULL
plottedGenes <- topGenes[1:20]
for (x in plottedGenes) {
  plotcounts[[x]] <-
    plotCounts(
      dds,
      gene = x,
      intgroup = c("Condition"),
      returnData = TRUE
    )
  plotcounts[[x]][, 'gene'] <-
    rep(x, length(plotcounts[[x]][, 'count']))
  counts_df <-
    rbind(counts_df, as.data.frame(plotcounts[[x]]))
  rownames(counts_df) <- NULL
}


p <-
  ggplot(data = counts_df, aes(x = gene, y = count)) + geom_boxplot(aes(fill =
                                                                          Condition))
p + facet_wrap(~ gene, scales = 'free') + labs(x = 'Gene', y = 'Normalised counts')

## counts for specific genes of interest
# counts_df2 <- NULL
# gene_interest <- ensembl_splice$hgnc_symbol
# for (x in gene_interest) {
#   plotcounts[[x]] <-
#     plotCounts(
#       dds,
#       gene = x,
#       intgroup = c("Condition"),
#       returnData = TRUE
#     )
#   plotcounts[[x]][, 'gene'] <-
#     rep(x, length(plotcounts[[x]][, 'count']))
#   # y <-
#   #   ensembl_df$external_gene_name[ensembl_df$ensembl_gene_id == x]
#   # plotcounts[[x]][, 'gene_name'] <-
#   #   rep(y, length(plotcounts[[x]][, 'count']))
#   counts_df2 <-
#     rbind(counts_df2, as.data.frame(plotcounts[[x]]))
#   rownames(counts_df2) <- NULL
# }
#
# p <-
#   ggplot(data = counts_df2, aes(x = gene, y = count)) + geom_boxplot(aes(fill =
#                                                                            Condition))
# p + facet_wrap( ~ gene, scales = 'free') + labs(x = 'Gene', y = 'Normalised counts')

# Heatmap of most variable genes
rld <- vst(dds_2)
# topVarGenes <- head(order(-rowVars(assay(rld))), 20)
# mat <- assay(rld)[topVarGenes,]
# mat <- mat - rowMeans(mat)
# df <-
#   as.data.frame(colData(rld)[, c("BioProject", "Condition")])
# annotation_colours <-
#   list(Condition = c(N = "lightblue", T = "red"))
# pheatmap(mat, annotation_col = df,)#annotation_colors = annotation_colours)

# Heatmap of top genes
mat <- assay(rld)[topGenes[1:30], ]
mat <- mat - rowMeans(mat)
df <-
  as.data.frame(colData(rld)[, c("BioProject", "Condition")])
annotation_colours <-
  list(Condition = c(N = "lightblue", T = "red"))
pheatmap(mat, annotation_col = df, )

#VolcanoPlot
# Better visualisation option than MA plot
EnhancedVolcano(
  res,
  lab = rownames(res),
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-6, 6),
  title = 'Tumour versus Normal tissues',
  pCutoff = 0.01,
  FCcutoff = 2,
  legendPosition = 'bottom',
  caption = 'Fold change cut off = 2.0, adjusted p-value = 0.01',
  subtitle = 'Differential expression'
)

#STRING list
# List of DEGs submitted to STRING
cat(rownames(DE_genes[abs(DE_genes$log2FoldChange) >= 2,]), sep = '\n')
