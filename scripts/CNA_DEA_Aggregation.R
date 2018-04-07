# CNA_DEA_Aggregation.R
if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}
setDTthreads(threads = detectCores())
if (!require(SummarizedExperiment)) {
    BiocInstaller::biocLite("SummarizedExperiment")
    library(SummarizedExperiment)
}
if (!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
}
if (!require(compiler)) {
    install.packages("compiler")
    library(compiler)
}
if (!require(parallel)) {
    install.packages("parallel")
    library(parallel)
}
if (!require(stringr)) {
    install.packages("stringr")
    library(stringr)
}
if (!require(genefilter)) {
    BiocInstaller::biocLite("genefilter")
    library(genefilter)
}
if (!require(WGCNA)) {
    BiocInstaller::biocLite(c("GO.db", "preprocessCore", "impute")) 
    install.packages("WGCNA")
    library(WGCNA)
    allowWGCNAThreads(nThreads = 8)
}
if (!require(IdeoViz)) {
    BiocInstaller::biocLite("IdeoViz")
    library(IdeoViz)
}

# ideo_table <- getIdeo("hg38")

# ==== Categorize segment mean values to loss, gain and normal ====
# Note that the functions in DEA_plots.R must be run before using the following
# script

# List all files
full_cna_files <- list.files("CNA_Gene_Pairs/Full/", full.names = T)
partial_cna_files <- list.files("CNA_Gene_Pairs/Partial/", full.names = T)
limma_dea_files <- list.files("SigGenes/", pattern = "limma", full.names = T)
deseq_dea_files <- list.files("SigGenes/", pattern = "DESeq", full.names = T)
dir.create("CNA_Plots")
dir.create("Paradoxical_Genes")
dir.create("CNA_Ideograms")
dir.create("Paradoxical_Gene_Plots/")
gene_dict <- fread("gene_dictionary.txt")

idx <- 1
upper_margin <- 0.2
lower_margin <- 0.2

tumor_percentile <- 0.4
normal_percentile <- 0.9
p_thresh <- 0.05
cna_analysis <- function(idx, upper_margin, lower_margin, p_thresh,
                         tumor_percentile, normal_percentile) {
    # Load a SummarizedExperiment file, retrieve the segment means and convert to DT
    curFile <- readRDS(full_cna_files[idx])
    # Get the project name
    proj <- gsub(pattern = ".*//(TCGA-\\w{1,4})_.*",
                 replacement = "\\1",
                 x = full_cna_files[idx])
    
    # Find the associated file containing the differentially expressed genes for 
    # this project
    # limma_dea <- grepl(pattern = paste0(".*", proj, ".*"), x = limma_dea_files)
    deseq_dea <- grepl(pattern = paste0(".*", proj, ".*"), x = deseq_dea_files)
    
    if (sum(deseq_dea) == 0) {
        # The associated CNA data has no matching DEA results
        return(paste0(proj, " has no matching DEA files"))
    } else {
        # limma_dea <- readRDS(limma_dea_files[limma_dea])
        deseq_dea <- readRDS(deseq_dea_files[deseq_dea])
    }
    
    # The current segment means from the CNA data
    curSegMean <- SummarizedExperiment::assay(curFile, 1)
    
    # Find the average segment means, by tissue type
    sample_names <- colnames(curSegMean)
    sample_split <- str_split_fixed(sample_names, "-", n = 7)
    # Get tissue type
    sample_type <- as.integer(gsub(pattern = "(\\d\\d)\\w",
                                   replacement = "\\1",
                                   x = sample_split[, 4]
    ))
    
    # Annotate samples based on type
    # TP: Solid tumor 1
    # TM: Metastatic 6
    # NB: Blood derived normal 10
    # NT: Solid Tissue normal 11
    sample_type[sample_type == 1] <- "TP"
    sample_type[sample_type == 6] <- "TM"
    sample_type[sample_type == 10] <- "NB"
    sample_type[sample_type == 11] <- "NT"
    
    # Create a 'dictionary'
    all_samples <- data.table(aliquot_id = sample_names, type = sample_type)
    
    NT_ids <- all_samples[type == "NT"]$aliquot_id
    NT_idx <- which(colnames(curSegMean) %in% NT_ids)
    
    TP_ids <- all_samples[type == "TP"]$aliquot_id
    TP_idx <- which(colnames(curSegMean) %in% TP_ids)
    
    # Subset only TP and NT for now
    NT_TP_segmeans <- curSegMean[, c(NT_idx, TP_idx)]
    NT_idx = 1:length(NT_idx)
    TP_idx = (length(NT_idx) + 1):(length(NT_idx) + length(TP_idx))
    
    # TODO: Should outliers be removed? Along which axis?
    # TODO: Why some genes have NA segment means? No coverage in CNA data?
    # TODO: Replace the NA values with the average of segment means in similar
    # tissue types; this allows us to use the 'genefilter' package to calculate
    # two-group t-tests for specified columns using efficient C code
    
    # Get na positions by tissue type
    # all_NTTP_na <- which(is.na(NT_TP_segmeans), arr.ind = T)
    # all_NT_na <- all_NTTP_na[all_NTTP_na[, 2] %in% NT_idx, ]
    # all_TP_na <- all_NTTP_na[all_NTTP_na[, 2] %in% TP_idx, ]
    # 
    # NT_TP_segmeans[all_NT_na] <- base::rowMeans(curSegMean[, NT_idx], na.rm = T)[all_NT_na[, 1]]
    # NT_TP_segmeans[all_TP_na] <- base::rowMeans(curSegMean[, TP_idx], na.rm = T)[all_TP_na[, 1]]
    
    NT_TP_facs <- factor(c(rep(x = "NT", length(NT_idx)), rep(x = "TP", length(TP_idx))))
    # Run t-tests by row (gene) and tissue type
    # NT_TP_ttest <- genefilter::rowttests(NT_TP_segmeans, fac = NT_TP_facs)
    
    NT_TP_segmeans <- as.data.table(NT_TP_segmeans, keep.rownames = T)
    
    # Calculate t-test p-values in parallel:
    # Define function
    # t_func <- function(row) {
    #     t.test(x = NT_TP_segmeans[row, NT_idx + 1, with = F],
    #            y = NT_TP_segmeans[row, TP_idx + 1, with = F],
    #            na.action = na.omit,
    #            alternative = "two.sided", mu = 0,
    #            var.equal = F)$p.value
    # }
    # # Compile and execute
    # t_func <- compiler::cmpfun(t_func)
    # t_results <- mclapply(
    #     X = 1:nrow(NT_TP_segmeans),
    #     FUN = t_func,
    #     mc.preschedule = T,
    #     mc.cores = detectCores(),
    #     mc.cleanup = T
    # )
    # # Append back to data.table
    # NT_TP_segmeans$t_pvalue <- unlist(t_results)
    # # Do a multiple comparison adjustment
    # # TODO: is this necessary given each gene is considered separately?
    # NT_TP_segmeans$t_pvalue_adj <- p.adjust(p = NT_TP_segmeans$t_pvalue,
    #                                     method = "fdr")
    
    # Calculate KS-test p-values
    ks_func <- function(row) {
        ks.test(x = unlist(NT_TP_segmeans[row, NT_idx + 1, with = F]),
               y = unlist(NT_TP_segmeans[row, TP_idx + 1, with = F]),
               alternative = "two.sided")$p.value
    }
    # Compile and execute
    ks_func <- compiler::cmpfun(ks_func)
    ks_results <- mclapply(
        X = 1:nrow(NT_TP_segmeans),
        FUN = ks_func,
        mc.preschedule = T,
        mc.cores = detectCores(),
        mc.cleanup = T
    )
    # Append back to data.table
    NT_TP_segmeans$ks_pvalue <- unlist(ks_results)
    # Do a multiple comparison adjustment
    # TODO: is this necessary given each gene is considered separately?
    NT_TP_segmeans$ks_pvalue_adj <- p.adjust(p = NT_TP_segmeans$ks_pvalue,
                                        method = "fdr")
    
    # Find means by type (there shouldn't be any NA values)
    # Add 1 to account for 'rn' column
    NT_TP_segmeans$NT_mean <- base::rowMeans(NT_TP_segmeans[, NT_idx + 1, with = F],
                                             na.rm = T)
    NT_TP_segmeans$TP_mean <- base::rowMeans(NT_TP_segmeans[, TP_idx + 1, with = F],
                                             na.rm = T)
    
    # Find standard deviation by type
    NT_TP_segmeans[, NT_sd := sd(.SD, na.rm = T), by = "rn", .SDcols = NT_idx + 1]
    NT_TP_segmeans[, TP_sd := sd(.SD, na.rm = T), by = "rn", .SDcols = TP_idx + 1]
    
    # Find median by type
    NT_TP_segmeans[, NT_median := median(unlist(.SD), na.rm = T),
                   by = "rn", .SDcols = NT_idx + 1]
    NT_TP_segmeans[, TP_median := median(unlist(.SD), na.rm = T),
                   by = "rn", .SDcols = TP_idx + 1]
    
    # Add actual quantiles to the data.table
    # NT_TP_segmeans$NT_actual_pcntl_up <- NT_row_quantiles_upper
    # NT_TP_segmeans$NT_actual_pcntl_lo <- NT_row_quantiles_lower
    # 
    # NT_TP_segmeans$TP_actual_pcntl_up <- TP_row_quantiles_upper
    # NT_TP_segmeans$TP_actual_pcntl_lo <- TP_row_quantiles_lower
    
    # Calculate and add empirical quantiles to the data.table
    NT_TP_segmeans[, NT_empirical_pcntl_up :=
                       qnorm(p = normal_percentile,
                             mean = .SD$NT_mean,
                             sd = .SD$NT_sd, lower.tail = F), by = "rn"]
    NT_TP_segmeans[, NT_empirical_pcntl_lo :=
                       qnorm(p = normal_percentile,
                             mean = .SD$NT_mean,
                             sd = .SD$NT_sd, lower.tail = T), by = "rn"]
    NT_TP_segmeans[, TP_empirical_pcntl_up :=
                       qnorm(p = tumor_percentile,
                             mean = .SD$TP_mean,
                             sd = .SD$TP_sd, lower.tail = F), by = "rn"] 
    NT_TP_segmeans[, TP_empirical_pcntl_lo :=
                       qnorm(p = tumor_percentile,
                             mean = .SD$TP_mean,
                             sd = .SD$TP_sd, lower.tail = T), by = "rn"] 
    # Subset for faster access
    NT_TP_segmeans <- NT_TP_segmeans[, c("rn", "NT_mean", "TP_mean",
                                         "NT_median", "TP_median",
                                         "NT_sd", "TP_sd",
                                         #"NT_actual_pcntl_up",
                                         #"NT_actual_pcntl_lo",
                                         "NT_empirical_pcntl_up",
                                         "NT_empirical_pcntl_lo",
                                         #"TP_actual_pcntl_up",
                                         #"TP_actual_pcntl_lo",
                                         "TP_empirical_pcntl_up",
                                         "TP_empirical_pcntl_lo",
                                         "t_pvalue",
                                         "ks_pvalue",
                                         "t_pvalue_adj",
                                         "ks_pvalue_adj")]
    colnames(NT_TP_segmeans)[1] <- "ensembl_gene_id"
    
    # Add chromosomal location to each gene
    NT_TP_segmeans <- merge(NT_TP_segmeans,
                            gene_dict[, c("ensembl_gene_id", "chromosome_name",
                                          "start_position", "end_position",
                                          "strand")])
    NT_TP_segmeans$strand <- ifelse(NT_TP_segmeans$strand > 0, yes = "+", no = "-")
    NT_TP_segmeans$chromosome_name <- paste0("chr", NT_TP_segmeans$chromosome_name)
    
    # ==== Ideogram plotting using IdeoViz ====
    # Create GRanges for plotting with the IdeoViz package
    # NT_TP_segmeans_gr <- makeGRangesFromDataFrame(
    #     df = NT_TP_segmeans,
    #     keep.extra.columns = T,
    #     ignore.strand = F,
    #     seqnames.field = "chromosome_name",
    #     start.field = "start_position",
    #     end.field = "end_position",
    #     strand.field = "strand"
    # )
    # 
    # # Create chromosome name vector
    # chromosome_names <- c(paste0("chr", seq(1, 22)), "chrX")
    # # Loop indices
    # i <- 1
    # image_idx <- 1
    # while (i < 23) {
    #     png(filename = paste0("CNA_Ideograms/CNA_Ideogram_", proj, "_",
    #                           image_idx, ".png"),
    #         res = 320, width = 3000, height = 3000)
    #     if (i == 19) {
    #         max_chrom = 4
    #     } else{
    #         max_chrom = 6
    #     }
    #     IdeoViz::plotOnIdeo(
    #         chrom = chromosome_names[i:(i + max_chrom)],
    #         ideoTable = ideo_table,
    #         values_GR = NT_TP_segmeans_gr,
    #         val_range = c(min(c(
    #             NT_TP_segmeans$NT_mean,
    #             NT_TP_segmeans$TP_mean)),
    #         max(c(
    #             NT_TP_segmeans$NT_mean,
    #             NT_TP_segmeans$TP_mean))),
    #         value_cols = c("NT_mean", "TP_mean"),
    #         plotType = "lines",
    #         col = c("green", "red"),
    #         addOnetoStart = T,
    #         vertical = F,
    #         addScale = F,
    #         plot_title = paste0("Medians of Segment Means by Chromosome for ", proj),
    #         ylab = "Averaged Segment Mean, across all samples by tissue type",
    #         xlab = "Position on chromosome"
    #     )
    #     dev.off()
    #     
    #     i <- i + max_chrom
    #     image_idx <- image_idx + 1
    # }
    # graphics.off()
    # ==== Segment mean plot preparation ====
    # Convert chromosome names to numbers to be able to sort by them
    NT_TP_segmeans$chromosome_num <- gsub(pattern = "chr(\\d{1,2})",
                                          replacement = "\\1",
                                          x = NT_TP_segmeans$chromosome_name)
    NT_TP_segmeans$chromosome_num <- gsub(pattern = "chrX", replacement = "23",
                                          x = NT_TP_segmeans$chromosome_num)
    NT_TP_segmeans$chromosome_num <- as.integer(NT_TP_segmeans$chromosome_num)
    # Order by chromosome number and start position
    data.table::setorder(x = NT_TP_segmeans, chromosome_num, start_position)
    # Add index
    NT_TP_segmeans$index <- 1:nrow(NT_TP_segmeans)
    # Find the indices of the last gene on each chromosome
    chrom_lines <- vector(length = 23, mode = "integer")
    for (i in 1:23) {
        chrom_lines[i] <- which(NT_TP_segmeans$chromosome_num == i)[1]
    }
    chrom_lines <- c(chrom_lines, nrow(NT_TP_segmeans))
    
    # graphics.off()
    # ==== Plot of segment means with no thresholds ====
    # png(filename = paste0("CNA_Plots/SegMean_Dist_No_Thresh_", proj, ".png"),
    #     res = 220, height = 1000, width = 2000)
    # plot(x = NT_TP_segmeans$index, y = NT_TP_segmeans$NT_mean, col = "blue",
    #                pch = 20, cex = 0.01, ylim = c(-.5, .5), xaxt = "n",
    #      xlab = "", ylab = "")
    # points(x = NT_TP_segmeans$index, y = NT_TP_segmeans$TP_mean,
    #        col = "black", xaxt = "n", ylab = "",
    #      xlab = "", ylim = c(-1, 1), cex = 0.01)
    # for (i in 1:length(chrom_lines)) {
    #     abline(v = chrom_lines[i], col = "black", lwd = 0.1)
    # }
    # text(labels = chromosome_names,
    #      x = ((chrom_lines + chrom_lines[-1])/2)[-24],
    #      y = 0.55, srt = 45, xpd = T, adj = 0, cex = 0.6)
    # abline(h = log2(2/2), col = "red", lwd = 1)
    # 
    # title(main = paste0("Gene-wise average segment mean distribution in ", proj),
    #       ylab = "Averaged Segment Mean", mgp = c(2,0,1),
    #       xlab = paste0("Genes ordered by chromosome location. No thresholds.\nBlue: average normal (solid) | Black: average tumor (solid" ))
    # dev.off()
    
    # ==== Selection of aberrated genes ====
    # upper_margin = 0.2
    # lower_margin = 0.2
    NT_TP_segmeans[, color := "black"]
    # NT_TP_segmeans[(TP_empirical_pcntl_up >
    #                    (NT_empirical_pcntl_lo + upper_margin)) &
    #                    (t_pvalue_adj < p_thresh),
    #                color := "red", by = "ensembl_gene_id"]
    NT_TP_segmeans[(TP_empirical_pcntl_up >
                        (NT_empirical_pcntl_lo + upper_margin)) &
                       (t_pvalue_adj < p_thresh),
                   color := "red", by = "ensembl_gene_id"]
    
    # NT_TP_segmeans[(TP_empirical_pcntl_lo <
    #                     (NT_empirical_pcntl_up - lower_margin)) &
    #                    (t_pvalue < p_thresh),
    #                color := "red", by = "ensembl_gene_id"]
    NT_TP_segmeans[(TP_empirical_pcntl_lo <
                        (NT_empirical_pcntl_up - lower_margin)) &
                       (t_pvalue_adj < p_thresh),
                   color := "red", by = "ensembl_gene_id"]
    NT_TP_segmeans[color == "red"]
    # NT_TP_segmeans$normal_mean
    # NT_TP_segmeans$tumor_actual_percentile
    # NT_TP_segmeans$normal_empirical_percentile[1:10]
    # NT_TP_segmeans$tumor_empirical_percentile[1:10]
    # NT_TP_segmeans[chromosome_num == 1]
    # TODO: Ensure the following respects the two-tailed quantile calculation
    # i.e. the addition/subtraction of upper and lower margins actually distinguishes
    # the genes of interes
    
    # ==== Plot of segment means with thresholds ====
    png(filename = paste0("CNA_Plots/SegMean_Dist_With_Thresh_", proj, ".png"),
        res = 220, height = 1000, width = 2000)
    plot(x = NT_TP_segmeans$index, y = NT_TP_segmeans$NT_median, col = "blue",
         pch = 20, cex = 0.01, ylim = c(-.5, .5), xaxt = "n",
         xlab = "", ylab = "")
    points(x = NT_TP_segmeans$index, y = NT_TP_segmeans$TP_median,
           col = NT_TP_segmeans$color, xaxt = "n", ylab = "",
           xlab = "", ylim = c(-1, 1), cex = 0.01)
    
    for (i in 1:length(chrom_lines)) {
        abline(v = chrom_lines[i], col = "black", lwd = 0.1)
    }
    text(labels = chromosome_names,
         x = ((chrom_lines + chrom_lines[-1])/2)[-24],
         y = 0.55, srt = 45, xpd = T, adj = 0, cex = 0.6)
    abline(h = log2(2/2), col = "red", lwd = 0.5)
    
    title(main = paste0("Gene-wise median of segment mean distribution in ", proj),
          ylab = "Medians of Segment Means", mgp = c(2,0,1),
          xlab = paste0("Genes ordered by chromosome location. Applied upper and lower thresholds.\nBlue: average normal | Red: Consistently aberrated\n t-test p-value < ", p_thresh))
    dev.off()
    
    # ==== Extraction of genes that are consistently in regions of CNAs ====
    cna_genes <- NT_TP_segmeans[color == "red",
                                           c("ensembl_gene_id",
                                             "NT_empirical_pcntl_lo",
                                             "NT_empirical_pcntl_up",
                                             "TP_empirical_pcntl_lo",
                                             "TP_empirical_pcntl_up")]
    # Convert the seg-means that fall outside the cutoffs to zero
    # NT_TP_segmeans$tumor_mean[!(NT_TP_segmeans$tumor_mean > upper_margin |
    #                       NT_TP_segmeans$tumor_mean < lower_margin)] <- 0
    # # And those that are not significant
    # NT_TP_segmeans$tumor_mean[NT_TP_segmeans$t_pvalue < 0.05] <- 0
    # upper_margin <- 0.5
    # lower_margin <- -0.5
    # sum((NT_TP_segmeans$tumor_mean > upper_margin |
    #           NT_TP_segmeans$tumor_mean < lower_margin) & NT_TP_segmeans$tumor_mean > 0.05)
    # sum(!((NT_TP_segmeans$tumor_mean > upper_margin |
    #       NT_TP_segmeans$tumor_mean < lower_margin) & NT_TP_segmeans$t_pvalue < 0.05))
    
    # Add the sign of the aberrations ()
    cna_genes$gain <- cna_genes$TP_empirical_pcntl_up > cna_genes$NT_empirical_pcntl_lo
    cna_genes$loss <- cna_genes$NT_empirical_pcntl_up > cna_genes$TP_empirical_pcntl_lo
    # table(cna_genes$gain)
    # table(cna_genes$loss)

    # cna_genes$gain <- ifelse(test = all_up_signs == 1, yes = 1, no = -1)
    # cna_genes$loss <- ifelse(test = all_lo_signs == -1, yes = 1, no = -1)
    
    # Indicate whether each gene is part of limma's under or over expressed set
    # NT_TP_segmeans$limma_over = ifelse(test = NT_TP_segmeans$Gene %in% limma_dea$over_exp,
    #                                yes = 1, no = -1)
    # NT_TP_segmeans$limma_under = ifelse(test = NT_TP_segmeans$Gene %in% limma_dea$under_exp,
    #                                yes = 1, no = -1)
    
    # Indicate whether each gene is part of DESSeq2's under or over expressed set
    cna_genes$deseq_over = ifelse(test = cna_genes$ensembl_gene_id %in% deseq_dea$over_exp,
                                   yes = T, no = F)
    cna_genes$deseq_under = ifelse(test = cna_genes$ensembl_gene_id %in% deseq_dea$under_exp,
                                    yes = T, no = F)
    # table(cna_genes$deseq_over)
    # table(cna_genes$deseq_under)
    # Chi-square test of paradoxical significance
    # Ensure there are any genes that meet the cutoff before testing
    # limma_over_pos <- tryCatch(expr = chisq.test(x = cna_genes$limma_over, y = cna_genes$gain, correct = F),
    #                            error = function(x) list(p.value = 1))
    # deseq_over_pos <- tryCatch(expr = chisq.test(x = cna_genes$deseq_over,
    #                                              y = cna_genes$gain,
    #                                              correct = F),
    #                            error = function(x) list(p.value = 1))
    # cor(x = cna_genes$deseq_over, y = cna_genes$gain, method = "spearman")
    # cor(x = cna_genes$deseq_over, y = cna_genes$gain, method = "pearson")
    # 
    # cor(x = cna_genes$deseq_over, y = cna_genes$loss, method = "spearman")
    # cor(x = cna_genes$deseq_over, y = cna_genes$loss, method = "pearson")
    
    # limma_under_pos <- tryCatch(expr = chisq.test(x = cna_genes$limma_under, y = cna_genes$gain, correct = F),
    #                             error = function(x) list(p.value = 1))
    # deseq_under_pos <- tryCatch(expr = chisq.test(x = cna_genes$deseq_under,
    #                                               y = cna_genes$gain,
    #                                               correct = F),
    #                             error = function(x) list(p.value = 1))
    # 
    # # limma_over_neg <- tryCatch(expr = chisq.test(x = cna_genes$limma_over, y = cna_genes$loss, correct = F),
    # #                            error = function(x) list(p.value = 1))
    # deseq_over_neg <- tryCatch(expr = chisq.test(x = cna_genes$deseq_over, y = cna_genes$loss, correct = F),
    #                            error = function(x) list(p.value = 1))
    # # limma_under_neg <- tryCatch(expr = chisq.test(x = cna_genes$limma_under, y = cna_genes$loss, correct = F),
    # #                             error = function(x) list(p.value = 1))
    # deseq_under_neg <- tryCatch(expr = chisq.test(x = cna_genes$deseq_under, y = cna_genes$loss, correct = F),
    #                             error = function(x) list(p.value = 1))
    
    # Find and depict the ratio of paradoxical genes to all differentially expressed
    # limma_over_par <- cna_genes$Gene[cna_genes$limma_over & cna_genes$loss]
    deseq_over_par <- cna_genes$ensembl_gene_id[cna_genes$deseq_over & cna_genes$loss]
    
    # limma_under_par <- cna_genes$Gene[cna_genes$limma_under & cna_genes$gain]
    deseq_under_par <- cna_genes$ensembl_gene_id[cna_genes$deseq_under & cna_genes$gain]
    
    # Remove NAs
    # limma_over_par <- limma_over_par[!is.na(limma_over_par)]
    deseq_over_par <- deseq_over_par[!is.na(deseq_over_par)]
    # limma_under_par <- limma_under_par[!is.na(limma_under_par)]
    deseq_under_par <- deseq_under_par[!is.na(deseq_under_par)]
    
    # Save these paradoxical genes
    saveRDS(file = paste0("Paradoxical_Genes/Paradoxical_", proj, ".rds"),
            object = list(#limma_over_par = limma_over_par,
                          #limma_under_par = limma_under_par,
                          deseq_over_par = deseq_over_par,
                          deseq_under_par = deseq_under_par))
    
    # Show the paradoxical genes on the same plot as before
    # limma:
    # png(filename = paste0("CNA_Plots/limma_SegMean_Paradoxical_", proj, ".png"),
    #     res = 220, height = 1000, width = 2000)
    # par(pch = 20)
    # cna_genes$color[cna_genes$color == "black"] <- "white"
    # cna_genes$color[cna_genes$color == "green"] <- "black"
    # cna_genes$color[cna_genes$color == "red"] <- "black"
    # cna_genes$color[cna_genes$ensembl_gene_id %in% limma_over_par] <- "red"
    # cna_genes$color[cna_genes$ensembl_gene_id %in% limma_under_par] <- "red"
    # 
    # plot(cna_genes$Mean, col = cna_genes$color, xaxt = "n", ylab = "",
    #      xlab = "", ylim = c(-1, 1))
    # abline(h = upper_margin, col = "red", lwd = 1)
    # abline(h = log2(2/2), col = "red", lwd = 1)
    # abline(h = lower_margin, col = "red", lwd = 1)
    # text(x = nrow(cna_genes) / 2, y = lower_margin + 0.1,
    #      labels = paste0("over-exp p-value: ", signif(limma_over_neg$p.value, 2)))
    # text(x = nrow(cna_genes) / 2, y = upper_margin - 0.1,
    #      labels = paste0("under-exp p-value: ", signif(limma_under_pos$p.value, 2)))
    # 
    # title(main = paste0("limma determined paradoxical genes' segment means in\n", proj),
    #       ylab = "Segment Mean", mgp = c(2,0,0),
    #       xlab = paste0("Thresholds of ", upper_margin, " and ", lower_margin,
    #                     " were selected for gain and loss\nP-values using Chi-Squared Test for significance of paradoxical genes\namong differentially expressed genes"))
    # dev.off()
    
    # DESeq2:
    png(filename = paste0("CNA_Plots/DESeq2_SegMean_Paradoxical_", proj, ".png"),
        res = 220, height = 1000, width = 2000)
    par(pch = 20)
    # cna_genes$color[cna_genes$color == "black"] <- "white"
    # cna_genes$color[cna_genes$color == "green"] <- "black"
    NT_TP_segmeans$color[NT_TP_segmeans$color == "red"] <- "black"
    NT_TP_segmeans$color[NT_TP_segmeans$ensembl_gene_id %in% deseq_over_par] <- "red"
    NT_TP_segmeans$color[NT_TP_segmeans$ensembl_gene_id %in% deseq_under_par] <- "red"
    
    plot(x = NT_TP_segmeans$index, y = NT_TP_segmeans$NT_median, pch = 20,
         col = "blue", xaxt = "n", ylab = "", cex = 0.1,
         xlab = "", ylim = c(-.5, .5))
    points(x = NT_TP_segmeans$index, y = NT_TP_segmeans$TP_median, pch = 20,
           cex = 0.1, col = "black")
    points(x = NT_TP_segmeans[color == "red"]$index,
           y = NT_TP_segmeans[color == "red"]$TP_median,
           col = "red", cex = 0.5, pch = 20)
    for (i in 1:length(chrom_lines)) {
        abline(v = chrom_lines[i], col = "black", lwd = 0.1)
    }
    text(labels = chromosome_names,
         x = ((chrom_lines + chrom_lines[-1])/2)[-24],
         y = 0.55, srt = 45, xpd = T, adj = 0, cex = 0.6)
    abline(h = log2(2/2), col = "red", lwd = 0.5)
    
    title(main = paste0("DESeq2 determined paradoxical genes' medians of\nsegment means in ", proj),
          ylab = "Medians of segment mean", mgp = c(2,0,1),
          xlab = paste0("Margins of ", upper_margin, " and ", lower_margin,
                        " were selected for gain and loss\n", tumor_percentile,
                        " of genes in tumor samples have a margin with ",
                        normal_percentile, " of normal samples"))
    dev.off()
    
    # ==== Focus on the paradoxical genes' distribution ====
    cur_data <- assay(curFile, 1)
    cur_mrna <- get(load(file = paste0("SummarizedExperiments/mRNA/",
                                       proj, "_GeneExpression.rdata")))
    exprs <- as.data.table(assay(cur_mrna, 1), keep.rownames = T)
    mrna_NT_idx <- which(cur_mrna$shortLetterCode == "NT")
    mrna_NT_ids <- colnames(cur_mrna)[mrna_NT_idx]
    mrna_TP_idx = which(cur_mrna$shortLetterCode == "TP")
    mrna_TP_ids <- colnames(cur_mrna)[mrna_TP_idx]
    
    # Subset
    exprs <- exprs[, c(1, mrna_NT_idx + 1, mrna_TP_idx + 1), with = F]
    # Convert to long format    
    molten <- melt.data.table(exprs, id.vars = c("rn"))
    molten[, type := variable %in% mrna_NT_ids, by = "rn"]
    molten$value <- molten$value + 1
    molten$value <- log2(molten$value)
    
    dir.create(paste0("Paradoxical_Gene_Plots/", proj))
    deseq_over_par <- gene_dict[ensembl_gene_id %in% deseq_over_par]
    deseq_under_par <- gene_dict[ensembl_gene_id %in% deseq_under_par]
    # Plot paradoxical over expressed
    for (i in 1:nrow(deseq_over_par)) {
        cur_plot <- ggplot(data = molten[rn == deseq_over_par[i]$ensembl_gene_id]) +
            geom_density(mapping = aes(x = value, fill = type), alpha = 0.5) +
            ggtitle(paste0("Log2 gene counts of ",
                           deseq_over_par[i]$external_gene_name,
                           " (", deseq_over_par[i]$gene_biotype, ")"),
                    subtitle = paste0("Overexpressed in event of copy number loss",
                                      "\nProject: ", proj,
                                      "\nChromosome: ", deseq_over_par[i]$chromosome_name,
                                      "\nExon count: ", deseq_over_par[i]$exon_count)) +
            scale_fill_discrete(name = "Sample type", labels = c("Tumor", "Normal")) +
            xlab("log2(Counts + 1)") +
            ylab("Density")
        ggsave(filename = paste0("Paradoxical_Gene_Plots/", proj, "/", proj, "_over_par_", deseq_over_par[i], ".png"),
               plot = cur_plot, device = "png", dpi = 300)
    }
    # Plot paradoxical under expressed
    for (i in 1:nrow(deseq_under_par)) {
        cur_plot <- ggplot(data = molten[rn == deseq_under_par[i]$ensembl_gene_id]) +
            geom_density(mapping = aes(x = value, fill = type), alpha = 0.5) +
            ggtitle(paste0("Log2 gene counts of ",
                           deseq_under_par[i]$external_gene_name,
                           " (", deseq_under_par[i]$gene_biotype, ")"),
                    subtitle = paste0("Underexpressed in event of copy number gain",
                                      "\nProject: ", proj,
                                      "\nChromosome: ", deseq_under_par[i]$chromosome_name,
                                      "\nExon count: ",deseq_under_par[i]$exon_count)) +
            scale_fill_discrete(name = "Sample type", labels = c("Tumor", "Normal")) +
            xlab("log2(Counts + 1)") +
            ylab("Density")
        ggsave(filename = paste0("Paradoxical_Gene_Plots/", proj, "/", proj,
                                 "_under_par_", deseq_under_par[i], ".png"),
               plot = cur_plot, device = "png", dpi = 300)
    }
    
    
    # Find limma and DESeq2's overlap of paradoxical genes
    # over_overlap = sum(deseq_over_par %in% limma_over_par) / length(deseq_over_par)
    # under_overlap = sum(deseq_under_par %in% limma_under_par) / length(deseq_under_par)
    
    # Paradoxical to differentially expressed
    # sum_limma_par = sum(c(limma_over_par, limma_under_par) %in%
    #                         c(limma_dea$over_exp, limma_dea$under_exp))
    # sum_limma_dea = (length(limma_dea$over_exp) + length(limma_dea$under_exp))
    # 
    # sum_deseq_par = sum(c(deseq_over_par, deseq_under_par) %in% 
    #                         c(deseq_dea$over_exp, deseq_dea$under_exp))
    # sum_deseq_dea = (length(deseq_dea$over_exp) + length(deseq_dea$under_exp))
    
    # return(
    #     list(
    #         project = proj,
    #         over_overlap = over_overlap,
    #         under_overlap = under_overlap,
    #         sum_limma_dea = sum_limma_dea,
    #         sum_limma_par = sum_limma_par,
    #         sum_deseq_dea = sum_deseq_dea,
    #         sum_deseq_par = sum_deseq_par
    #     )
    # )
    
    # barplot(c(sum_limma_par, sum_limma_dea))
    # pie(x = c(sum_deseq_par, sum_deseq_dea), col = c("white", "black"))
    
    # # Find the mean and SD of segment means of limma and DESeq2 overexpressed genes
    # limma_over_segmeans <- cna_genes[Gene %in% limma_dea$over_exp, c("Gene", "Mean", "SD")]
    # deseq_over_segmeans <- cna_genes[Gene %in% deseq_dea$over_exp, c("Gene", "Mean", "SD")]
    
    # # Find the mean and SD of segment means of limma and DESeq2 underexpressed genes
    # limma_under_segmeans <- cna_genes[rn %in% limma_dea$under_exp, c("rn", "Mean", "SD")]
    # deseq_under_segmeans <- cna_genes[rn %in% deseq_dea$under_exp, c("rn", "Mean", "SD")]
    # 
    # # Find the sign of segnemt means
    # limma_over_segmeans$sign = sign(limma_over_segmeans$Mean)
    # limma_under_segmeans$sign = sign(limma_under_segmeans$Mean)
    # deseq_over_segmeans$sign = sign(deseq_over_segmeans$Mean)
    # deseq_under_segmeans$sign = sign(deseq_under_segmeans$Mean)
    # 
    # # Find paradoxical genes
    # limma_over_par = limma_over_segmeans[sign(limma_over_segmeans$Mean) == -1,]
    # deseq_over_par = deseq_over_segmeans[sign(deseq_over_segmeans$Mean) == -1,]
    # limma_over_par[, Mean := -1]
    # deseq_over_par[, Mean := -1]
    # 
    # limma_under_par = limma_under_segmeans[sign(limma_under_segmeans$Mean) == 1,]
    # deseq_under_par = deseq_under_segmeans[sign(deseq_under_segmeans$Mean) == 1,]
    # limma_under_par[, Mean := 1]
    # deseq_under_par[, Mean := 1]
    # table(limma_under_par$Mean, limma_over_par$Mean)
}

cna_analysis <- compiler::cmpfun(cna_analysis)
for (i in 1:length(full_cna_files)) {
    cat("Working on ", full_cna_files[i], "\n")
    cna_analysis(
        idx = i,
        upper_margin = 0.2,
        lower_margin = 0.2,
        p_thresh = 0.05,
        tumor_percentile = 0.4,
        normal_percentile = 0.9
    )
}

# mc.results <- mclapply(X = 1:length(full_cna_files), FUN = cna_analysis,
#                        lower_margin = -0.2, upper_margin = 0.2,
#                        mc.preschedule = T, mc.cores = detectCores())

# Subset projects with matching DEA files
valid_results <- sapply(X = 1:length(mc.results), FUN = function(x) is.list(mc.results[[x]]))
# Aggregate all results in one place
cna_results <- mc.results[valid_results]


cna_results <- as.data.table(cbind(
    sapply(cna_results, function(x) x$project),
    sapply(cna_results, function(x) x$over_overlap),
    sapply(cna_results, function(x) x$under_overlap),
    sapply(cna_results, function(x) x$sum_limma_dea),
    sapply(cna_results, function(x) x$sum_limma_par),
    sapply(cna_results, function(x) x$sum_deseq_dea),
    sapply(cna_results, function(x) x$sum_deseq_par)
    ))
colnames(cna_results) <- c(
    "project",
    "over_overlap",
    "under_overlap",
    "sum_limma_dea",
    "sum_limma_par",
    "sum_deseq_dea",
    "sum_deseq_par"
)

cna_results$limma_ratio <- as.integer(cna_results$sum_limma_par) /
    as.integer(cna_results$sum_limma_dea)
cna_results$deseq_ratio <- as.integer(cna_results$sum_deseq_par) /
    as.integer(cna_results$sum_deseq_dea)


# Show the ratio of the paradoxical genes per cancer type
paradoxical <- cna_results[, c("project", "limma_ratio", "deseq_ratio")]
# paradoxical <- cna_results[rn == "project" |
#                     rn == "sum_limma_dea" | rn == "sum_limma_par"]
                    # rn == "sum_deseq_dea" | rn == "sum_deseq_par"]
# colnames(paradoxical) <- unlist(paradoxical[1, ])
# paradoxical <- paradoxical[-1, ]
# colnames(paradoxical)[ncol(paradoxical)] <- "rn"
# Convert to long format
ggplot(data = paradoxical[, 1:2], mapping = aes(x = project, y = limma_ratio)) +
    geom_bar(stat = "identity") + 
    ggtitle(label = "Ratio of paradoxical to differentially expressed genes", 
            sub = "As determined by limma") + 
    xlab("Cancer Project") + ylab("Ratio") +
    ylim(c(0, 0.2)) +
    theme(axis.text.x = element_text(hjust = 0, angle = -45),
          plot.margin = unit(c(1,1,1,1), "cm"))
ggsave(filename = "CNA_Plots/overview_limma_ratios.png", device = "png", dpi = 320)

ggplot(data = paradoxical[, c(1, 3)], mapping = aes(x = project, y = deseq_ratio)) +
    geom_bar(stat = "identity") + 
    ggtitle(label = "Ratio of paradoxical to differentially expressed genes", 
            sub = "As determined by DESeq2") + 
    xlab("Cancer Project") + ylab("Ratio") +
    ylim(c(0, 0.2)) +
    theme(axis.text.x = element_text(hjust = 0, angle = -45),
          plot.margin = unit(c(1,1,1,1), "cm"))
ggsave(filename = "CNA_Plots/overview_deseq_ratios.png", device = "png", dpi = 320)


