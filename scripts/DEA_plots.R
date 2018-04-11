# DEA_plots.R
if (!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
}
if (!require(limma)) {
    install.packages("limma")
    library(limma)
}
if (!require(parallel)) {
    install.packages("parallel")
    library(parallel)
}
if (!require(compiler)) {
    install.packages("compiler")
    library(compiler)
}
if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}
if (!require(DESeq2)) {
    install.packages("DESeq2")
    library(DESeq2)
}

dir.create("SigGenes")

# ==== limma MA plots ====
# Note that the MA plot 
limma_files <- list.files("limma_Results/", full.names = T)
dir.create("MA_Plots")
dir.create("Volcano_Plots")

# PlotLimma_MA_Volcano <- function(idx, p_value = 0.001, min_lfc = 1) {
#     # Load DEA (list) file and extract project name
#     curDEA <- readRDS(limma_files[idx])
#     curInfo <- curDEA[[1]]
#     curDEA <- curDEA[[2]]
#     
#     proj <- gsub(pattern = ".*_(TCGA-.+)\\.rds", replacement = "\\1",
#                  x = limma_files[idx])
#     
#     # First filter by p-value, then by logFC, separately.
#     # The alternative way, i.e. ranking by fold change and then filtering by
#     # p-value, has the disadvantage of only ranking by a single parameter, and
#     # also fails to properly rank genes which have zero expression in group 1
#     # and a non-zero expression in group 2 - the fold change in this case would
#     # be theoretically infinite. In contrast, a p-value will give you a more
#     # tangible quantification of how significant such differential expression event is.
#     # The following sets the minimum logFC to 1, uses FDR and "separate" multiple 
#     # testing scheme.
#     # https://www.researchgate.net/post/FDR_or_log_fold_change_which_one_is_the_priority_for_selecting_the_DEGs
# 
#     # Find why the "separate" method must be used as opposed to "global"
#     deGenes <- decideTests(curDEA, method = "separate", lfc = min_lfc, 
#                            adjust.method = "fdr", p.value = p_value)
#     
#     # ==== MA Plot:
#     # Create a PNG device to save the MA plot
#     png(paste0("MA_Plots/limma_DEA_MA_Plot_", proj, ".png"),
#         width = 1500, height = 1000, res = 220)
#     par(mgp = c(2, 1, 0))
#     limma::plotMD(object = curDEA, status = deGenes[, 2], 
#            main = paste0("Averaged MA plot for ", proj),
#            hl.cex = 0.1, bg.cex = 0.1,
#            sub = paste0("# Normal",
#                         ": ", unname(curInfo[2]), " | ", "# Tumor",
#                         ": ", unname(curInfo[3]),
#                         " | max p-value: ", p_value,
#                         " | min logFC: ", min_lfc), legend = F)
#     abline(h = c(-min_lfc, min_lfc), col = "red")
#     legend(x = "topright", legend = c("Over Exp", "Under Exp"), pch = 20, pt.cex = 1,
#            col = c("green", "red"), cex = 0.8)
#     
#     dev.off()
#     
#     # ==== Volcano Plot:
#     deGenes <- as.data.table(deGenes, keep.rownames = T)
#     
#     # Get the over and under expressed genes
#     over <- deGenes[`mRNAdata$shortLetterCodeTP` == 1]$rn
#     under <- deGenes[`mRNAdata$shortLetterCodeTP` == -1]$rn
#     
#     # Save these genes for use later on
#     saveRDS(list(over_exp = over, under_exp = under),
#             file = paste0("SigGenes/limma_SigGenes_p_", p_value, "_lfc_",
#                           min_lfc, "_", proj, ".rds"))
#     
#     # Get current logFC and the adjusted p-values of all genes
#     curLFC <- topTable(fit = curDEA,
#                        number = Inf,
#                        adjust.method = "BH")[, c("logFC", "adj.P.Val"),]
#     # Convert to data.table
#     curLFC <- as.data.table(curLFC, keep.rownames = T)
#     # Convert adjusted p-values to -log base 10
#     curLFC$neg_base_ten_p <- -log(x = curLFC$adj.P.Val, base = 10)
#     # Delete the p-value column
#     curLFC[ , adj.P.Val := NULL]
#     
#     png(paste0("Volcano_Plots/limma_DEA_MA_Plot_", proj, ".png"),
#         width = 1500, height = 1000, res = 220)
#     par(mgp = c(2, 1, 0))
#     
#     plot(x = curLFC$logFC, y = curLFC$neg_base_ten_p, xlim = c(-20, 20),
#         pch = 16,
#         cex = 0.3,
#         main = paste0("Volcano plot of ", proj),
#         sub =  paste0("# Normal: ", nt_count,
#                       " | # Tumor: ", tp_count,
#                       " | min p-value: ", p_value,
#                       " | min log2FC: ", min_lfc),
#         #ylim = c(-10, 10),
#         xlab = "log2FC tumor vs normal",
#         ylab = "-log base 10 p-value")
# 
#     # Highlight significant genes
#     points(x = curLFC[rn %in% over]$logFC,
#            y = curLFC[rn %in% over]$neg_base_ten_p,
#            pch = 16, cex = 0.3, col = "green")
#     points(x = curLFC[rn %in% under]$logFC,
#            y = curLFC[rn %in% under]$neg_base_ten_p,
#            pch = 16, cex = 0.3, col = "red")
#     abline(v = c(-min_lfc, min_lfc), col = "red")
#     legend(x = "topright", legend = c("Over Exp", "Under Exp"), pch = 20, pt.cex = 1,
#            col = c("green", "red"), cex = .8)
#     
#     dev.off()
# }
# 
# PlotLimma_MA_Volcano <- compiler::cmpfun(PlotLimma_MA_Volcano)
# 
# mc.results <-
#     mclapply(
#         X = 1:length(limma_files),
#         FUN = PlotLimma_MA_Volcano,
#         p_value = 0.01,
#         min_lfc = 1,
#         mc.preschedule = T,
#         mc.cores = detectCores()
#     )

# ==== DESeq2 MA and Volcano plots ====
deseq_files <- list.files("DESeq_Results/", full.names = T)

# Plot MA and Volcano at the same time
plotDeseq_MA_Volcano <- function(idx, p_value = 0.001, min_lfc = 1) {
    curFile <- readRDS(deseq_files[idx])
    # Get DEA results
    res <- results(curFile, lfcThreshold = min_lfc, independentFiltering = T,
                    pAdjustMethod = "fdr", alpha = p_value)
    
    # Calculate the number of normal and tumor samples
    nt_count <- sum(curFile$shortLetterCode == "NT")
    tp_count <- sum(curFile$shortLetterCode == "TP")
    # Extract project name from meta data
    proj <- gsub(pattern = ".*(TCGA-\\w{1,4})\\s.*", replacement = "\\1",
                 metadata(curFile)$Title)
    
    # MA plot:
    # Creata PNG image
    png(filename = paste0("MA_Plots/DESeq2_DEA_MA_Plot_", proj, ".png"),
        res = 220, width = 2000, height = 2000)
    # Create an MA plot using DESeq2's own method, providing project info
    # and cutoffs used, while highlighting differentially expressed genes based 
    # on those cutoffs.
    DESeq2::plotMA(object = res, main = paste0("MA plot of ", proj),
           sub =  paste0("# Normal: ", nt_count,
                         " | # Tumor: ", tp_count,
                         " | min p-value: ", p_value,
                         " | min log2FC: ", min_lfc), ylim = c(-10, 10))
    abline(h = c(-min_lfc, min_lfc), col = "red")
    # Save MA plot
    dev.off()
    
    # TODO: Plot explanation:
    # Note that testing for |logFC| > 1 by TREAT is not the same as selecting
    # genes with |logFC| > 1. Genes will need to exceed this threshold by some
    # margin, depending on the data, before being declared statistically
    # significant. It would be better to interpret the threshold as
    # ‘the fold-change below which we are definitely not interested in the gene’
    # rather than ‘the fold-change above which we are interested in the gene’.
    
    # Volcano plot:
    # The NA values are a result of a threshold; a requirement for the number of 
    # counts high enough to have a good power in detecting differential expression
    volcano_data <-
        data.table(gene = rownames(res), logFC = res$log2FoldChange,
                   negLog10 = -log10(res$padj))
    sig_genes_sign <-
        sign(res[abs(res$log2FoldChange) > min_lfc &
                         (res$padj < p_value) %in% TRUE, "log2FoldChange"])
    signif_genes <- rownames((res[abs(res$log2FoldChange) > min_lfc &
                                      (res$padj < p_value) %in% TRUE, ]))
    
    # Save these genes for use later on
    under = signif_genes[sig_genes_sign == -1]
    over = signif_genes[sig_genes_sign == 1]
    
    saveRDS(list(under_exp = under, over_exp = over),
            file = paste0("SigGenes/DESeq2_SigGenes_p_", p_value, "_lfc_",
                          min_lfc, "_", proj, ".rds"))
    # Plot all genes
    png(filename = paste0("Volcano_Plots/Volcano_DESeq2_", proj, ".png"),
        res = 220, width = 2000, height = 2000)
    plot(x = volcano_data$logFC,
        y = volcano_data$negLog10,
        pch = 16,
        cex = 0.6,
        main = paste0("Volcano plot of ", proj),
        sub =  paste0("# Normal: ", nt_count,
                      " | # Tumor: ", tp_count,
                      " | min p-value: ", p_value,
                      " | min log2FC: ", min_lfc),
        #ylim = c(-10, 10),
        xlab = "log2FC tumor vs normal",
        ylab = "-log base 10 p-value")
    
    # Highlight significant genes
    points(x = volcano_data[gene %in% signif_genes]$logFC,
           y = volcano_data[gene %in% signif_genes]$negLog10,
           pch = 16, cex = 0.6, col = "red")
    # Save Volcano plot
    dev.off()
}

plotDeseq_MA_Volcano <- compiler::cmpfun(plotDeseq_MA_Volcano)
mc.results2 <-
    mclapply(
        X = 1:length(deseq_files),
        FUN = plotDeseq_MA_Volcano,
        p_value = 0.01,
        min_lfc = 1,
        mc.preschedule = T,
        mc.cores = detectCores()
    )


# Number of differentially expressed genes selected for each cancer type ====
sig_files <- list.files("SigGenes/", pattern = "lfc_1", full.names = T)
deseq_files <- list.files("DESeq_Results/", full.names = T)
parad_files <- list.files("Paradoxical_Genes/", full.names = T)
# counts <- data.table(
#     Project = character(),
#     OverExpressed = integer(),
#     UnderExpressed = integer(),
#     NormalSamples = integer(),
#     TumorSamples = integer()
# )
idx <- 1
cur_func <- function(idx) {
    cur_parad <- readRDS(parad_files[idx])
    proj_name <- gsub(pattern = ".*(TCGA-\\w{1,4})\\.rds", replacement = "\\1",
                      parad_files[idx])
    
    cur_mRNA <- readRDS(deseq_files[grepl(pattern = proj_name, x = deseq_files,
                                          fixed = T)])
    nt_count <- sum(cur_mRNA$shortLetterCode == "NT")
    tp_count <- sum(cur_mRNA$shortLetterCode == "TP")
    
    # if (sum(grepl(pattern = proj_name, x = parad_files, fixed = T)) == 0) {
    #     return ()
    # }
    cur_diff <- readRDS(sig_files[grepl(pattern = proj_name, x = sig_files,
                                           fixed = T)])
    
    cur_data <- data.table(Project = proj_name,
                           OverExpressed = length(cur_diff$over_exp),
                           UnderExpressed = length(cur_diff$under_exp),
                           NormalSamples = nt_count,
                           TumorSamples = tp_count,
                           ParadoxicallyOverExpressed = length(cur_parad$over_exp_parad),
                           ParadoxicallyUnderExpressed = length(cur_parad$under_exp_parad))
    return(cur_data)
}

mc.results <-
    mclapply(
        X = 1:length(parad_files),
        FUN = cur_func,
        mc.preschedule = T,
        mc.cores = 16
    )
all_counts <- rbindlist(mc.results)

# Save
fwrite(all_counts, "SigGene_Per_Tumor_Counts.csv")
