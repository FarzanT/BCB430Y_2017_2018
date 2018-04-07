# miRNA_CNA_correlation.R

# Test whether tumor deregulated miRNAs fall under abberated regions
if (!require(BSgenome.Hsapiens.UCSC.hg38)) {
    biocLite("BSgenome.Hsapiens.UCSC.hg38")
    library(BSgenome.Hsapiens.UCSC.hg38)
}

hg38 <- BSgenome.Hsapiens.UCSC.hg38

# Create chromosome name vector
chromosome_names <- c(paste0("chr", seq(1, 22)), "chrX")
# Get chromosome sizes from BSGenome's hg38
hg38 <- BSgenome::getSeq(x = hg38, names = chromosome_names)
chr_lengths <- seqlengths(hg38)
rm(hg38)

# Partition/bin the genome into 10 Kb segments
all_chroms <- list()
for (i in 1:23) {
    cur_seq <- seq(from = 1, to = chr_lengths[i] - 10000, by = 10000)
    all_chroms[[i]] <- data.table(seq(from = 1, to = length(cur_seq)),
                                  rep(x = i, times = length(cur_seq)),
                                  seq(from = 1, to = chr_lengths[i] - 10000, by = 10000),
                                  seq(from = 10001, to = chr_lengths[i], by = 10000))
}
bin_dt <- rbindlist(all_chroms)

colnames(bin_dt) <- c("SeqNum", "Chromosome", "Start", "End")
bin_dt[Chromosome == "chr23"]$Chromosome <- "chrX"
bin_gr <-
    makeGRangesFromDataFrame(
        df = bin_dt,
        keep.extra.columns = F,
        seqnames.field = "Chromosome",
        start.field = "Start",
        end.field = "End"
    )

dir.create("miRNA_CNA_Correlation")

# Load CNA data files, deregulated miRNA files and mirBase data
cna_files <- list.files("SummarizedExperiments/CNA/", full.names = T)
mirna_files <- list.files("Signif_miRNA/", full.names = T)
mirbase <- readRDS("Human_miRNA_GRanges.rds")
# Drop chrY and replace chrX with 23
seqlevels(mirbase) <- gsub(pattern = "chr", replacement = "", x = seqlevels(mirbase))
seqlevels(mirbase) <- gsub(pattern = "X", replacement = "23", x = seqlevels(mirbase))
seqnames(mirbase) <- gsub(pattern = "chr", replacement = "", x = seqnames(mirbase))
seqnames(mirbase) <- gsub(pattern = "X", replacement = "23", x = seqnames(mirbase))
mirbase <- mirbase[seqnames(mirbase) != "Y"]

idx <- 3
mirna_cna_func <- function(idx) {
    # Load the CNV file
    cnvMatrix <- as.data.table(get(load(cna_files[idx])))
    # TODO: REMOVE
    # cnvMatrix <- cnvMatrix[1:10000,]
    
    cnvMatrix[Chromosome == "X"]$Chromosome <- 23
    # Drop Num_Probes column
    cnvMatrix <- cnvMatrix[, !("Num_Probes")]
    
    # Discard invalid sequences with negative width
    cnvMatrix <- cnvMatrix[End > Start]
    
    # Get the project name
    proj <- gsub(pattern = ".*//(TCGA-\\w{1,4})_.*",
                 replacement = "\\1",
                 x = cna_files[idx])
    
    sample_names <- unique(cnvMatrix$Sample)
    sample_split <- str_split_fixed(sample_names, "-", n = 7)
    # Get tissue type
    sample_type <- as.integer(gsub(pattern = "(\\d\\d)\\w",
                                   replacement = "\\1",
                                   x = sample_split[, 4]
    ))
    # TP: Solid tumor 1
    # TM: Metastatic 6
    # NB: Blood derived normal 10
    # NT: Solid Tissue normal 11
    sample_type[sample_type == 1] <- "TP"
    # Type 3 is blood derived, we'll consider it as TP for this analysis
    sample_type[sample_type == 3] <- "TP"
    sample_type[sample_type == 6] <- "TM"
    sample_type[sample_type == 10] <- "NB"
    sample_type[sample_type == 11] <- "NT"
    
    # Create a 'dictionary'
    all_samples <- data.table(aliquot_id = sample_names, type = sample_type)
    
    # Find overlaps of each sample separately for less memory usage
    # sample <- all_samples$aliquot_id[1]
    find_ov <- function(sample) {
        cur_sample <- makeGRangesFromDataFrame(
            df = cnvMatrix[Sample == sample],
            keep.extra.columns = T,
            seqnames.field = "Chromosome",
            start.field = "Start",
            end.field = "End",
            ignore.strand = T
        )
        # Find the overlaps between the binned ranges
        overlap <- findOverlaps(query = cur_sample, subject = bin_gr, type = "any")
        final <- as.data.table(cur_sample[queryHits(overlap)]$Segment_Mean)
        
        final2 <- as.data.table(bin_gr[subjectHits(overlap), ])[, c(1, 2, 3)]
        final <- cbind(final2, final, rep(sample, nrow(final)))
        colnames(final) <- c("seqnames", "start", "end", "Segment_Mean", "sample")
        
        # Change shape to have samples as columns
        wide <- dcast.data.table(data = final, formula = seqnames + start + end ~ sample,
                                 value.var = "Segment_Mean", fun.aggregate = sum)
        data.table::setorder(wide, seqnames, start, end)
        
        return(wide)
    }
    
    find_ov <- compiler::cmpfun(find_ov)
    all_overlaps <- list()
    for (i in all_samples$aliquot_id) {
        all_overlaps[[i]] <- find_ov(i)
    }
    
    # Merge all columns
    all <- Reduce(function(...) merge(..., by = c("seqnames", "start", "end")), all_overlaps)
    rm(list = c("all_overlaps"))
    gc()
    
    # Calculate KS-test p-values, ensure to compare x: TP to y:NT
    wc_func <- function(row) {
        less <- wilcox.test(x = unlist(all[row, all_samples[type == "TP"]$aliquot_id, with = F]),
                            y = unlist(all[row, all_samples[type == "NT"]$aliquot_id, with = F]),
                            alternative = "less")$p.value
        greater <- wilcox.test(x = unlist(all[row, all_samples[type == "TP"]$aliquot_id, with = F]),
                               y = unlist(all[row, all_samples[type == "NT"]$aliquot_id, with = F]),
                               alternative = "greater")$p.value
        return(data.table(less = less, greater = greater))
    }
    # Compile and execute
    wc_func <- compiler::cmpfun(wc_func)
    system.time({
        wc_results <- mclapply(
            X = 1:nrow(all),
            FUN = wc_func,
            mc.preschedule = T,
            mc.cores = detectCores(),
            mc.cleanup = T
        )
        wc_results <- rbindlist(wc_results)
        
    })
    # Append to all segments
    all$wilcox_p_less <- wc_results$less
    all$wilcox_p_grtr <- wc_results$greater
    # Adjust p values
    all$wilcox_p_less_adj <- p.adjust(wc_results$less, method = "fdr")
    all$wilcox_p_grtr_adj <- p.adjust(wc_results$greater, method = "fdr")
    
    # Find the miRNAs in these regions
    all_less_gr <- tryCatch(expr = {
        makeGRangesFromDataFrame(df = all[wilcox_p_less_adj <= 0.1,
                                          c("seqnames", "start", "end"), with = F],
                                 seqnames.field = "seqnames")
    }, error = function(...) return(NULL))
    
    all_grtr_gr <- tryCatch(expr = {
        makeGRangesFromDataFrame(df = all[wilcox_p_grtr_adj < 0.1,
                                          c("seqnames", "start", "end")],
                                 seqnames.field = "seqnames")
    }, error = function(...) return(NULL))
    
    rm(all)
    
    # Find overlaps with gene dictionary
    loss <- NULL
    gain <- NULL
    if (!is.null(all_less_gr)) {
        all_less_ol <- findOverlaps(query = all_less_gr, subject = mirbase, type = "any")
        # Find genome hits, classify gain and loss
        loss <- names(mirbase[unique(subjectHits(all_less_ol)), ])
    }
    if (!is.null(all_grtr_gr)) {
        all_grtr_ol <- findOverlaps(query = all_grtr_gr, subject = mirbase, type = "any")
        gain <- names(mirbase[unique(subjectHits(all_grtr_ol)), ])
    }
    
    # Save
    saveRDS(list(gain = gain, loss = loss),
            paste0("miRNA_CNA_Correlation/", proj, "_mirna_cna_cor.rds"))
}
mirna_cna_func <- compiler::cmpfun(mirna_cna_func)
mirna_cna_func(1)
for (i in 4:length(cna_files)) {
    mirna_cna_func(i)
}
    # Add normal and tumor segment mean medians
    # all[, normal_med := median(unlist(.SD), na.rm = T),
    #                by = c("seqnames", "start", "end"), .SDcols = all_samples[type == "NT"]$aliquot_id]
    # all[, tumor_med := median(unlist(.SD), na.rm = T),
    #     by = c("seqnames", "start", "end"), .SDcols = all_samples[type == "TP"]$aliquot_id]
    
    # Find the genes in significant segments:}


# ==== Find which deregulated miRNAs are in CNA regions ====
mirna_cna_files <- list.files("miRNA_CNA_Correlation/", full.names = T)
mirna_difex_files <- list.files("miRNA_Differential_Expression/", full.names = T)

cur_func <- function(idx) {
    proj <- gsub(pattern = ".*(TCGA-\\w{1,4})\\.rds",
                 replacement = "\\1",
                 x = mirna_difex_files[idx])
    
    cur_mirna_difex <- readRDS(mirna_difex_files[idx])
    cur_results <-
        results(
            object = cur_mirna_difex,
            lfcThreshold = 2,
            alpha = 0.01,
            pAdjustMethod = "fdr",
            independentFiltering = T
        )
    temp_up <- rownames(cur_results)[cur_results$padj <= 0.01 & cur_results$log2FoldChange > 2]
    temp_down <- rownames(cur_results)[cur_results$padj <= 0.01 & cur_results$log2FoldChange < -2]
    up_mirna <- temp_up[!is.na(temp_up)]
    down_mirna <- temp_down[!is.na(temp_down)]
    
    # Load miRNA/CNA correlation file
    cur_mirna_cna <- readRDS(file = mirna_cna_files[grepl(pattern = proj,
                                             x = mirna_cna_files, fixed = T)])
    
    up_mirna %in% cur_mirna_cna$gain
    
    
}