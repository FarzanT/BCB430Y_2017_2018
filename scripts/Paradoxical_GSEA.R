# Paradoxical_GSEA.R

if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}
if (!require(gage)) {
    BiocInstaller::biocLite("gage")
    library(gage)
}
if (!require(topGO)) {
    BiocInstaller::biocLite("topGO")
    library(topGO)
}
if (!require(Rgraphviz)) {
    BiocInstaller::biocLite("Rgraphviz")
    library(Rgraphviz)
}
if (!require(clusterProfiler)) {
    # Note that sudo apt-get install libudunits2-dev must be run in Linux
    devtools::install_github("GuangchuangYu/clusterProfiler")
    devtools::install_github("GuangchuangYu/DOSE")
    devtools::install_github("GuangchuangYu/enrichplot")
    library(clusterProfiler)
}
if (!require(clusterProfiler)) {
    BiocInstaller::biocLite("clusterProfiler")
    library(clusterProfiler)
}

if (!require(DESeq2)) {
    BiocInstaller::biocLite("DESeq2")
    library(DESeq2)
}
if (!require(org.Hs.eg.db)) {
    BiocInstaller::biocLite("org.Hs.eg.db")
    library(org.Hs.eg.db)
}
if (!require(fastmatch)) {
    install.packages("fastmatch")
    library(fastmatch)
}
if (!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
}



# ==== Run gene set enrichment analysis for paradoxically expressed genes ====
par_files <- list.files("Paradoxical_Genes/", full.names = T)
mrna_files <- list.files("SummarizedExperiments/mRNA/", full.names = T)
deseq_files <- list.files("DESeq_Results/", full.names = T)
dir.create("GSEA_files")
gene_dict <- fread("gene_dictionary.txt")
# Only keep rows where entrez ID is not NA
gene_dict <- gene_dict[!is.na(entrezgene)]

hsa <- search_kegg_organism('hsa', by='kegg_code')


idx = 1
gageAnalysis <- function(idx) {
    # Load paradoxical genes for each cancer project
    cur_parad <- readRDS(par_files[idx])
    # Find
    # Get the project name
    proj <- gsub(pattern = ".*(TCGA-\\w{1,4}).*",
                 replacement = "\\1",
                 x = par_files[idx])
    
    
    # # Convert to data.table
    # curLFC <- as.data.table(curLFC, keep.rownames = T)
    # # Get expression data
    # cur_assay <- assay(cur_mrna, 1)
    # # Calculate geometric mean of all genes per group (TP and NT only)
    # # Subset NT and TP
    # NT_cols <- colnames(cur_assay)[which(cur_mrna$shortLetterCode == "NT")]
    # TP_cols <- colnames(cur_assay)[which(cur_mrna$shortLetterCode == "TP")]
    # # Transpose for faster column wise operations
    # TP_assay <- t(cur_assay[, TP_cols])
    # NT_assay <- t(cur_assay[, NT_cols])
    # 
    # # Log, get column means and exponentiate to get geometric mean
    # TP_assay <- log(TP_assay)
    # TP_geom_mean <- exp(base::colMeans(TP_assay, na.rm = T))
    # NT_assay <- log(NT_assay)
    # NT_geom_mean <- exp(base::colMeans(NT_assay, na.rm = T))
    # TP_assay / NT_assay
    # cur_mrna <- cur_mrna[, cur_mrna$shortLetterCode %in% c("NT", "TP")]
    # normal_columns <- which(cur_mrna$shortLetterCode == "NT")
    # # Create a DESeq DataSet, design by tissue type
    # cur_des <- DESeq2::DESeqDataSet(se = cur_mrna, design = ~ shortLetterCode)
    # 
    # # We can remove the rows that have no or nearly no information about the
    # # amount of gene expression (across all samples)
    # cur_des <- cur_des[rowSums(counts(cur_des)) > 1, ]
    # # Get the normalized counts (first estimate size factors)
    # cur_des <- estimateSizeFactors(cur_des)
    # norm_counts <- counts(cur_des, normalized = T)
    # Subset for paradoxical genes and solid tissue samples
    # para_exp <- norm_counts[rownames(norm_counts) %in% c(cur_parad$under_exp_parad, cur_parad$over_exp_parad), ]
    all_parad <- c(cur_parad$under_exp_parad, cur_parad$over_exp_parad)
    
    # Convert to entrez gene IDs
    gene.df <- bitr(all_parad, fromType = "ENSEMBL", 
                    toType = c("ENTREZID", "SYMBOL"),
                    OrgDb = org.Hs.eg.db)
    
    # cur_assay[, NT_geom_mean := exp(mean(unlist(log(.SD[is.finite(log(unlist(.SD)))]), na.rm=T))), by = "rn", .SDcols = NT_cols]
    # cur_assay[, TP_geom_mean := exp(mean(unlist(log(.SD)), na.rm = T)), by = "rn", .SDcols = TP_cols]
    # cur_assay[, NT_geom_mean := exp(mean(unlist(log(.SD)), na.rm = T)), by = "rn", .SDcols = NT_cols]
    # 
    # all_entrez <- gene_dict[ensembl_gene_id %in% rownames(cur_mrna)][, c("ensembl_gene_id", "entrezgene")]
    # unique(all_entrez)
    
    kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                     organism     = 'hsa',
                     pvalueCutoff = 0.05)
    
    ego <- enrichGO(gene          = gene.df$ENTREZID,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.1,
                    readable      = TRUE)
        
    # kegg_result <- gage(exprs = para_exp, gsets = kegg.gs, ref = normal_columns, same.dir = F, compare = "unpaired", rank.test = T)
    # go_result <- gage(exprs = para_exp, gsets = go.gs, ref = normal_columns, same.dir = F, compare = "unpaired")
    
    
    # sel <- kegg_result$greater[, "p.val"] < 0.05 & !is.na(kegg_result$greater[, "p.val"])
    # path_ids <- rownames(kegg_result$greater)[sel]
    if (is.null(kk) | is.null(ego)) {
        return(proj)
    }
    # Save
    saveRDS(object = list(kegg = kk, go = ego),
            file = paste0("GSEA_files/GSEA_", proj, ".rds"))
}
gageAnalysis <- compiler::cmpfun(gageAnalysis)
gc()
mc.results <-
    mclapply(
        X = 1:length(par_files),
        FUN = gageAnalysis,
        mc.preschedule = T,
        mc.cores = 10
    )


# ==== 
dir.create("ClusterProfiler_Graphs")
dir.create("ClusterProfiler_Graphs/GO")
dir.create("ClusterProfiler_Graphs/KEGG")

gsea_files <- list.files("GSEA_files/", full.names = T)
f = gsea_files[2]
plot_func <- function(f) {
    # Get the project name
    proj <- gsub(pattern = ".*(TCGA-\\w{1,4}).*",
                 replacement = "\\1",
                 x = f)
    
    cur <- readRDS(f)
    cur_kegg <- cur$kegg
    cur_go <- cur$go
    
    # Find and load associated DESeq2 file
    cur_deseq <- readRDS(deseq_files[grepl(pattern = proj, x = deseq_files)])
    # Get current logFC and the adjusted p-values of all genes
    res <- results(cur_deseq, lfcThreshold = 0, independentFiltering = T,
                   pAdjustMethod = "BH", alpha = 0.05)
    cur_lfc <- data.table(ENSG = rownames(res), lfc = 2**res[, "log2FoldChange"])
    match_idx <- fastmatch::fmatch(x = cur_lfc$ENSG, table = gene_dict$ensembl_gene_id)
    # Change ENSG to ENTREZ
    cur_lfc$entrez <- gene_dict$entrezgene[match_idx]
    # Get the log fold change of paradoxical genes
    parad_lfc <- cur_lfc[entrez %in% cur_go@gene]$lfc
    names(parad_lfc) <- cur_lfc[entrez %in% cur_go@gene]$entrez
    # Find the median of gene expression values for that project
    # cur_mrna <- get(load(mrna_files[grepl(pattern = proj, x = mrna_files, fixed = T)]))
    # ==== GO plots ====
    # Create CNET plot
    go_cnet <- cnetplot(x = cur_go, foldChange = parad_lfc) +
        ggtitle(paste0("CNET plot of paradoxical genes in ", proj),
                subtitle = "GO based enrichment, biological process") 
    ggsave(
        filename = paste0("ClusterProfiler_Graphs/GO/", proj, "_cnet.png"),
        plot = go_cnet,
        device = "png",
        dpi = 320,
        width = 30,
        height = 25,
        units = "cm"
    )
    # Create emap plot
    go_emap <- emapplot(x = cur_go) +
        ggtitle(paste0("EMAP plot of paradoxical genes in ", proj),
                subtitle = "GO based enrichment, biological process")
    ggsave(
        filename = paste0("ClusterProfiler_Graphs/GO/", proj, "_emap.png"),
        plot = go_emap,
        device = "png",
        dpi = 320,
        width = 30,
        height = 25,
        units = "cm"
    )
    # ==== KEGG Plots ====
    go_cnet <- cnetplot(x = cur_go, foldChange = parad_lfc) +
        ggtitle(paste0("CNET plot of paradoxical genes in ", proj),
                subtitle = "KEGG based enrichment") 
    ggsave(
        filename = paste0("ClusterProfiler_Graphs/KEGG/", proj, "_cnet.png"),
        plot = go_cnet,
        device = "png",
        dpi = 320,
        width = 30,
        height = 25,
        units = "cm"
    )
    # Create emap plot
    go_emap <- emapplot(x = cur_go) +
        ggtitle(paste0("EMAP plot of paradoxical genes in ", proj),
                subtitle = "KEGG based enrichment")
    ggsave(
        filename = paste0("ClusterProfiler_Graphs/KEGG/", proj, "_emap.png"),
        plot = go_emap,
        device = "png",
        dpi = 320,
        width = 30,
        height = 25,
        units = "cm"
    )
    # ==== GO graph plot ====
    plotGOgraph(cur_go)
    dev.print(pdf, paste0("ClusterProfiler_Graphs/GO/", proj, "_GOgraph.pdf"))
}

plot_func <- compiler::cmpfun(plot_func)
mc.results <-
    mclapply(
        X = gsea_files,
        FUN = plot_func,
        mc.preschedule = T,
        mc.cores = 10
    )

