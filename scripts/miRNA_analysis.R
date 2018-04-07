# miRNA_analysis.R

if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}
if (!require(qdapTools)) {
    install.packages("qdapTools")
    library(qdapTools)
}
if (!require(clusterProfiler)) {
    BiocInstaller::biocLite("clusterProfiler")
    library(clusterProfiler)
}
if (!require(biomaRt)) {
    BiocInstaller::biocLite("biomaRt")
    library(biomaRt)
}
if (!require(targetscan.Hs.eg.db)) {
    BiocInstaller::biocLite("targetscan.Hs.eg.db")
    library(targetscan.Hs.eg.db)
}
if (!require(stringr)) {
    install.packages("stringr")
    library(stringr)
}
if (!require(SummarizedExperiment)) {
    install.packages("SummarizedExperiment")
    library(SummarizedExperiment)
}
if (!require(DESeq2)) {
    install.packages("DESeq2")
    library(DESeq2)
}
if (!require(parallel)) {
    install.packages("parallel")
    library(parallel)
}
if (!require(compiler)) {
    install.packages("compiler")
    library(compiler)
}
if (!require(ppcor)) {
    install.packages("ppcor")
    library(ppcor)
}
if (!require(org.Hs.eg.db)) {
    BiocInstaller::biocLite("org.Hs.eg.db")
    library(org.Hs.eg.db)
}



# Find the targets of the _differentially expressed_ miRNA's in each cancer type
mirna_files <- list.files("SummarizedExperiments/miRNA/", full.names = T)
# annotated_exp <- list.files("AnnotatedExperiments/", full.names = T)
mRNA_files <- list.files("SummarizedExperiments/mRNA", full.names = T)
dir.create("miRNA_Differential_Expression")
dir.create("miRNA_Tidy_Expression_Files")

idx <- 4
mirna_diffex <- function(idx) {
    proj <- gsub(pattern = ".*//(TCGA-\\w{1,4})_.*",
                 replacement = "\\1",
                 x = mirna_files[idx])
    # Load annotated experiment file containing samples with both mRNA and CNA
    # data
    # if (identical(annotated_exp[grep(pattern = proj, x = annotated_exp)], character(0))) {
    #     return(paste0("Project ", proj, " does not have matching mRNA/CNA data"))
    # }
    cur_mrna <- get(load(mRNA_files[grepl(pattern = proj, mRNA_files)]))
    # annot_file <- readRDS(annotated_exp[grep(pattern = proj, x = annotated_exp)])
    if (!("TP" %in% cur_mrna$shortLetterCode & "NT" %in% cur_mrna$shortLetterCode)) {
        return(paste0("Project ", proj, " does not have both normal and tumor samples"))
    }
    # Subset
    cur_mrna <- cur_mrna[, cur_mrna$shortLetterCode %in% c("NT", "TP")]
    # Load miRNA file
    cur_file <- as.data.table(get(load(mirna_files[idx])))
    # Find the type of samples
    # annot_file$shortLetterCode
    
    # Extract sample barcodes from current miRNA file and the annotation file
    annot_barcodes <- colnames(cur_mrna)
    mirna_barcodes <- unique(gsub(pattern = ".*(TCGA.*)", replacement = "\\1",
                              x = colnames(cur_file))[-1])
    
    # Split barcodes into sub-identifiers
    mirna_split <- as.data.table(str_split(mirna_barcodes, pattern = "-",
                                              simplify = T))
    mirna_tissue_split = str_split(mirna_split$V4, "", simplify = T)
    
    mirna_split$participant <- paste(mirna_split$V1,
                                     mirna_split$V2,
                                     mirna_split$V3,
                                     paste0(mirna_tissue_split[, 1],
                                            mirna_tissue_split[, 2]), sep = "-")
    annot_split <- as.data.table(str_split(annot_barcodes, pattern = "-",
                                                simplify = T))
    annot_tissue_split = str_split(annot_split$V4, "", simplify = T)
    
    annot_split$participant <- paste(annot_split$V1,
                                     annot_split$V2,
                                     annot_split$V3,
                                    paste0(annot_tissue_split[, 1],
                                           annot_tissue_split[, 2]), sep = "-")
    # Find the overlapping participant IDs, extract matches from the original
    # (indices haven't changed)
    match_idx = mirna_split$participant %in% annot_split$participant
    matches <- mirna_barcodes[match_idx]
    # Find tissue type from the last numbers of the barcode
    matches_tissue_type <-
        ifelse(test = as.integer(paste0(
            mirna_tissue_split[match_idx, 1],
            mirna_tissue_split[match_idx, 2])) > 10,
        yes = "NT",
        no = "TP")
    if (length(table(matches_tissue_type)) == 1) {
        return(paste0("The miRNA file for project ", proj, " is only of one type"))
    }
    # Extract the read counts from the miRNA file
    matching_counts <- cur_file[, c("miRNA_ID", paste0("read_count_", matches)), with = F]
    sample_data <- data.table(miRNA = matches, tissue_type = matches_tissue_type)
    # Createa a DESeqDataSet for DESeq2, exclude the first column (miRNA IDs)
    my_deseq <- DESeqDataSetFromMatrix(
        countData = matching_counts[,-1],
        colData = sample_data,
        design = ~ tissue_type
    )
    # Add miRNA IDs and meta data
    rownames(my_deseq) <- matching_counts$miRNA_ID
    metadata(my_deseq)$title <- paste0("miRNA count data for ", proj)
    
    # Save
    saveRDS(my_deseq, file = paste0("miRNA_Tidy_Expression_Files/miRNA_Exp_",
                                    proj, ".rds"))
    
    # Run differential miRNA expression analysis:
    # Set number of workers
    myBPPARAM <- BiocParallel::MulticoreParam(workers = detectCores())

    cur_difEx <- DESeq(object = my_deseq, parallel = T, BPPARAM = myBPPARAM)

    # Save results
    saveRDS(cur_difEx, file = paste0("miRNA_Differential_Expression/DESeq2_miRNA_", proj, ".rds"))

    # Stop the cluster and run a garbage collection
    BiocParallel::bpstop(myBPPARAM)
    gc()
}
mirna_diffex <- compiler::cmpfun(mirna_diffex)
# mirna_diffex(idx = 3)
my_results <- lapply(X = 1:length(mirna_files), FUN = mirna_diffex)

# ==== Find consistently deregulated miRNAs and their targets ====
mirna_difex_files <- list.files("miRNA_Differential_Expression/", full.names = T, pattern = "DESeq")
dir.create("miRNA_DEA_Plots")
dir.create("Signif_miRNA")

idx <- 2
min_lfc <- 2
p_value <- 0.05
mirna_analysis <- function(idx, min_lfc, p_value) {
    # Load file
    cur_file <- readRDS(mirna_difex_files[idx])
    # Extract project name
    proj <- gsub(pattern = ".*miRNA_(TCGA-\\w{1,4})\\.rds", replacement = "\\1",
                 x = mirna_difex_files[idx])
    num_NT <- sum(cur_file$tissue_type == "NT")
    num_TP <- sum(cur_file$tissue_type == "TP")
    # Find deregulated miRNAs
    cur_results <-
        results(
            object = cur_file,
            lfcThreshold = min_lfc,
            alpha = p_value,
            pAdjustMethod = "fdr",
            independentFiltering = T
        )
    # Create an MA plot
    png(filename = paste0("miRNA_DEA_Plots/MA_miRNA_DESeq2_", proj, ".png"))
    DESeq2::plotMA(object = cur_results, ylim = c(-10, 10),
                   main = paste0("MA plot of miRNA differential expression in ", proj),
                   sub = paste0("# normal: ", num_NT, " | # tumor: ", num_TP,
                                " | min lfc: ", min_lfc, " | max p-value: ", p_value))
    # Save plot
    dev.off()
    
    # Extract deregulated miRNAs
    cur_sig_mirna <- cur_results[(abs(cur_results$log2FoldChange) >= min_lfc &
                                     cur_results$padj <= p_value) %in% TRUE, ]
    cur_sig_signs <- sign(cur_sig_mirna$log2FoldChange)
    
    cur_over <- rownames(cur_sig_mirna)[cur_sig_signs == 1]
    cur_under <- rownames(cur_sig_mirna)[cur_sig_signs == -1]
    
    # Save over and under expressed miRNAs
    saveRDS(object = list(under_exp_mirna = cur_under,
                          over_exp_mirna = cur_over),
            file = paste0("Signif_miRNA/Signif_miRNA_DESeq2_", proj, ".rds"))
    
    # Create a volcano plot
    cur_volcano_data <- data.table(mirna_names = rownames(cur_results),
                                   logFC = cur_results$log2FoldChange, 
                                   negLog10 = -log10(cur_results$padj))
    cur_volcano_data$color <- "black"
    cur_volcano_data[cur_volcano_data$mirna_names %in% c(cur_over, cur_under)]$color <- "red"
    
    png(filename = paste0("miRNA_DEA_Plots/Volcano_miRNA_DESeq2_", proj, ".png"))
    plot(
        x = cur_volcano_data$logFC,
        y = cur_volcano_data$negLog10,
        pch = 20,
        cex = 0.5,
        col = cur_volcano_data$color,
        xlab = "Log2 Fold Change",
        ylab = "- log10 p-value",
        main = paste0("Volcano plot of miRNAs in ", proj),
        sub = paste0("# normal: ", num_NT, " | # tumor: ", num_TP,
                     " | min lfc: ", min_lfc, " | max p-value: ", p_value)
    )
    dev.off()
}
mirna_analysis <- compiler::cmpfun(mirna_analysis)

mc.results <-
    mclapply(
        X = 1:length(mirna_difex_files),
        FUN = mirna_analysis,
        mc.preschedule = T,
        mc.cores = detectCores(),
        min_lfc = 2,
        p_value = 0.01
    )

# ==== Use a hypergeometric test to assess the significance of deregulated miRNA
# and paradoxical genes for each cancer project ====
mirdip <- fread("Fa_Ta_mirDIP_v_3_1520371014/Unique_mirdip.txt")

# Convert the target entrez IDs to ensembl IDs
gene_dict <- fread("gene_dictionary.txt")
gene_dict[entrezgene %in% humanMirnaTargets[[1]]]$ensembl_gene_id
humanMirnaTargets <- lapply(
    X = humanMirnaTargets,
    FUN = function(x) {gene_dict[entrezgene %in% x]$ensembl_gene_id}
)

# ==== Check whether paradoxical genes are targets of deregulated miRNAs ====
dereg_mirna <- list.files("Signif_miRNA/", full.names = T)
paradoxical_genes <- list.files("Paradoxical_Genes/", full.names = T)
mirna_exp_files <- list.files("miRNA_Tidy_Expression_Files", full.names = T)
deseq_files <- list.files("DESeq_Results/", full.names = T)

idx <- 16
min_lfc <- 2
p_value <- 0.01
paradoxical_analysis <- function(idx, min_lfc, p_value) {
    cur_paradox <- readRDS(paradoxical_genes[idx])
    # Extract project name
    proj <- gsub(pattern = ".*paradoxical_(TCGA-\\w{1,4})\\.rds", replacement = "\\1",
                 x = paradoxical_genes[idx])
    mirna_match <- grepl(pattern = proj, x = dereg_mirna, fixed = T)
    if (sum(mirna_match) == 0) {
        return(paste0("Project ", proj, " does not have relevant miRNA data"))
    }
    # Load deregulated miRNAs
    cur_mirna <- readRDS(dereg_mirna[mirna_match])
    # Change mir to miR, discard numbers that indicate miRNA with identical mature sequences
    # TODO: It's the mature sequences that matter during mRNA regulation, correct?
    # under_exp_mirna <- gsub(pattern = "mir", replacement = "miR", x = cur_mirna$under_exp_mirna, fixed = T)
    # over_exp_mirna <- gsub(pattern = "mir", replacement = "miR", x = cur_mirna$over_exp_mirna, fixed = T)
    under_exp_mirna <-
        sapply(X = cur_mirna$under_exp_mirna, function(x)
            gsub(pattern = "hsa-mir-(\\d{1,}\\w*)-*.*", replacement = "hsa-miR-\\1", x = x), USE.NAMES = F)
    over_exp_mirna <-
        sapply(X = cur_mirna$over_exp_mirna, function(x)
            gsub(pattern = "hsa-mir-(\\d{1,}\\w*)-*.*", replacement = "hsa-miR-\\1", x = x), USE.NAMES = F)

    # Convert to data.table
    under_exp_mirna <- unique(data.table(miRNA = under_exp_mirna))
    over_exp_mirna <- unique(data.table(miRNA = over_exp_mirna))
    
    # Find their targets (all matches)
    if (nrow(under_exp_mirna) != 0) {
        under_targets <- merge(under_exp_mirna, mirdip, by = "miRNA")
    } else {
        under_targets <- data.table(miRNA = character())
    }
    if (nrow(over_exp_mirna) != 0) {
        over_targets <- merge(over_exp_mirna, mirdip, by = "miRNA")
    } else {
        over_targets <- data.table(miRNA = character())
    }
    
    # for (cur in under_exp_mirna) {
    #     mirdip[miRNA == cur]
    # }
    # under_match_idx <- fmatch(x = under_exp_mirna, table = mirdip$miRNA)
    # over_match_idx <- fmatch(x = over_exp_mirna, table = mirdip$miRNA)
    # under_targets <- humanMirnaTargets[under_exp_mirna]
    # over_targets <- humanMirnaTargets[over_exp_mirna]
    
    # Get the original expression levels for miRNAs and genes to calculate correlations
    cur_mirna_exp <- readRDS(paste0("miRNA_Tidy_Expression_Files/miRNA_Exp_", proj, ".rds"))
    cur_gene_exp <- readRDS(deseq_files[grepl(pattern = proj, x = deseq_files, ignore.case = T)])
    
    # Get the deregulated genes
    dereg <- results(object = cur_gene_exp, lfcThreshold = min_lfc, independentFiltering = T,
                     pAdjustMethod = "fdr", alpha = p_value)
    dereg_genes <- rownames(dereg[(dereg$padj <= p_value & abs(dereg$log2FoldChange) >= min_lfc) %in% T,])
    
    # Calculate ratio of deregulated genes that are targets of deregulated miRNAs
    dereg_ratio <- sum(dereg_genes %in% c(unlist(under_targets, use.names = F),
                           unlist(over_targets, use.names = F))) / length(dereg_genes)
    
    # Number of paradoxically over/under expressed genes that are targets of
    # under/over expressed miRNAs
    over_target_size <- sum(cur_paradox$under_exp_parad %in% unique(over_targets$ENSG))
    under_target_size <- sum(cur_paradox$over_exp_parad %in% unique(under_targets$ENSG))
    
    # Hypergeometric test of significance
    # x, q - vector of quantiles representing the number of white balls drawn
    # without replacement from an urn which contains both black and white balls.
    # m	- the number of white balls in the urn.
    # n	- the number of black balls in the urn.
    # k	- the number of balls drawn from the urn.
    p_hyp_geom_test <- 1 - phyper(q = over_target_size + under_target_size,
                                  k = length(unlist(cur_paradox)),
                                  m = length(dereg_genes),
                                  n = nrow(cur_gene_exp) - length(dereg_genes))
    
    return(data.table(proj = proj, dereg_ratio = dereg_ratio, over_target_size = over_target_size,
               under_target_size = under_target_size,
               p_hyp_geom_test = p_hyp_geom_test))
}

paradoxical_analysis <- compiler::cmpfun(paradoxical_analysis)

mc.results <- mclapply(
    X = 1:length(paradoxical_genes),
    FUN = paradoxical_analysis,
    mc.preschedule = F,
    mc.cores = 16,
    min_lfc = 2,
    p_value = 0.01
)

mc.results <- mc.results[-6]
all_results <- rbindlist(mc.results)

fwrite(all_results, "HyperGeom_Ratio_miRNA_targets.csv")