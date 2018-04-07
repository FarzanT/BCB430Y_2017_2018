# KM_Plot_mapping.R
download.file("http://kmplot.com/analysis/studies/probe%20sets%20to%20genes.txt", "KMplot_probedata.txt")
kmp <- fread("KMplot_probedata.txt", header = F)

brca_parad <- readRDS("Paradoxical_Genes/DESeq2_paradoxical_TCGA-BRCA.rds")
brca_mirna <- readRDS("miRNA_Differential_Expression/DESeq2_miRNA_TCGA-BRCA.rds")

# Convert ENSG to HGNC
gene_dict <- fread("gene_dictionary.txt")
cur_data <- gene_dict[ensembl_gene_id %in% brca_parad$over_exp_parad,
              c("ensembl_gene_id", "external_gene_name")]

# Print (UP) paradoxical genes in BRCA
cat(cur_data$external_gene_name, sep = ",")
# Print deregulated miRNAs in BRCA
cur_results <-
    results(
        object = brca_mirna,
        lfcThreshold = 2,
        alpha = 0.01,
        pAdjustMethod = "fdr",
        independentFiltering = T
    )
diff_mirna <- data.table(miRNA = rownames(cur_results), as.data.table(cur_results))
dereg_mirna <- diff_mirna[cur_results$padj <= 0.01, ]$miRNA
cat(dereg_mirna, sep = ",")


# Convert to probe ids
cur_sub <- kmp[V2 %in% cur_data$external_gene_name]
# Select the first 65
cat(cur_sub$V1[1:65], sep = "\n")