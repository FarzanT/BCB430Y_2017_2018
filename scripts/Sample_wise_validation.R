# Sample_wise_validation.R

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

# List gene expression files, paradoxical genes, CNA data
all_exp <- list.files("SummarizedExperiments/mRNA/", full.names = T)
all_cna <- list.files("SummarizedExperiments/CNA", full.names = T)
all_parad <- list.files("Paradoxical_Genes/", pattern = "DESeq", full.names = T)
dir.create("Paradoxical_Plots")
# Load gene dictionary
gene_dict <- fread("gene_dictionary.txt")
# Discard chrY
gene_dict <- gene_dict[chromosome_name %in% c(1:22, "X")]
# Change chrX
gene_dict$chromosome_name[gene_dict$chromosome_name == "X"] <- "23"
gene_dict$chromosome_name <- as.integer(gene_dict$chromosome_name)

# Convert to GRanges
genome_gr <-
    makeGRangesFromDataFrame(
        df = gene_dict[, c("start_position",
                           "end_position",
                           "chromosome_name",
                           "ensembl_gene_id")],
        keep.extra.columns = T,
        ignore.strand = T,
        seqnames.field = "chromosome_name",
        start.field = "start_position",
        end.field = "end_position")

# Plot function
plot_exp_and_segmeans <- function(cur_parad,
                                  proj,
                                  mrna_match,
                                  cna_match,
                                  cur_mrna,
                                  cur_cna,
                                  tissue_type,
                                  norm_counts, parad_type, multi_sampled) {
    # ==== Plot gene expression data ====
    # Add tissue type
    sample_names <- unique(cur_mrna$barcode)
    sample_split <- str_split_fixed(sample_names, "-", n = 7)
    # Patient ID consists of the first three
    mrna_patients <- paste(sample_split[, 1], sample_split[, 2], sample_split[, 3], sep = "-")
    # Find the barcode of multi sampled patients for mRNA data
    multi_barcodes <- sample_split[mrna_patients %in% multi_sampled, ]
    
    # Only consider these samples
    valid_barcodes <- paste(
        multi_barcodes[, 1],
        multi_barcodes[, 2],
        multi_barcodes[, 3],
        multi_barcodes[, 4],
        multi_barcodes[, 5],
        multi_barcodes[, 6],
        multi_barcodes[, 7],
        sep = "-"
    )
    
    # Get tissue type
    gene_exp <- norm_counts[cur_parad$under_exp_parad, valid_barcodes]
    gene_exp <- as.data.table(gene_exp, keep.rownames = T)
    
    # Add tissue type and reshape for plotting
    molten <- melt.data.table(data = gene_exp, id.vars = "rn")
    molten <- merge(molten, tissue_type, by.x = "variable", by.y = "ID")
    colnames(molten)[2] <- "ensembl_gene_id"
    
    # TODO: LEARN THE FOLLOWING SYNTAX
    molten$gene_name <- gene_dict[molten, .SD, nomatch=0L, on="ensembl_gene_id",
                                  .SDcols=names(gene_dict)]$external_gene_name
    cur_plot <- ggplot(data = molten, mapping = aes(x = gene_name, y = value, col = type)) +
        geom_point() + 
        # Add mean bars
        stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean",
                     size= 0.3, geom = "crossbar", show.legend = F) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab("Paradoxical Genes") + ylab("DESeq2 Normalized Count") +
        ggtitle(paste0("Paradoxical Gene Expression Levels in ", proj),
                subtitle = paste0("Genes in tumor tissue are paradoxically ", parad_type, "-expressed"))
    
    ggsave(filename = paste0("Paradoxical_Plots/Gene_Exp_Plot_", proj, ".png"),
           plot = cur_plot, device = "png", dpi = 320, width = 20, units = "in")
    
    
    # ==== Plot CNA data ====
    # Discard invalid segments
    cur_cna <- cur_cna[cur_cna$End > cur_cna$Start, ]
    # Add tissue type
    sample_names <- unique(cur_cna$Sample)
    sample_split <- str_split_fixed(sample_names, "-", n = 7)
    
    # Get tissue type
    sample_type <- as.integer(gsub(pattern = "(\\d\\d)\\w",
                                   replacement = "\\1",
                                   x = sample_split[, 4]
    ))
    # Patient ID consists of the first three
    cna_patients <- paste(sample_split[, 1], sample_split[, 2], sample_split[, 3], sep = "-")
    # Find the barcode of multi sampled patients for mRNA data
    multi_barcodes <- sample_split[cna_patients %in% multi_sampled, ]
    
    # Only consider these samples
    valid_barcodes <- paste(
        multi_barcodes[, 1],
        multi_barcodes[, 2],
        multi_barcodes[, 3],
        multi_barcodes[, 4],
        multi_barcodes[, 5],
        multi_barcodes[, 6],
        multi_barcodes[, 7],
        sep = "-"
    )
    
    sample_type[sample_type == 1] <- "TP"
    sample_type[sample_type == 3] <- "TP"
    sample_type[sample_type == 6] <- "TM"
    sample_type[sample_type == 10] <- "NB"
    sample_type[sample_type == 11] <- "NT"
    # Only keep multi sampled data
    cur_cna <- cur_cna[cur_cna$Sample %in% valid_barcodes, ]
    
    # Create a 'dictionary'
    cna_tissues <- data.table(aliquot_id = sample_names, type = sample_type)
    # Convert CNA data to GRanges
    cna_gr <-
        makeGRangesFromDataFrame(
            df = cur_cna,
            keep.extra.columns = T,
            ignore.strand = T,
            seqnames.field = "Chromosome",
            start.field = "Start",
            end.field = "End"
        )
    # Find the positions of paradoxical genes in the gene dictionary
    parad_gr <- makeGRangesFromDataFrame(df = gene_dict[ensembl_gene_id %in% cur_parad$under_exp_parad,
                                                        c("chromosome_name", "start_position", "end_position", "external_gene_name")],
                                         keep.extra.columns = T, seqnames.field = "chromosome_name",
                                         start.field = "start_position", end.field = "end_position",
                                         ignore.strand = T)
    # Find overlaps with CNA data
    parad_ol <- findOverlaps(query = parad_gr, subject = cna_gr, type = "any")
    
    # Get segment mean values for all samples
    cna_hits <- cna_gr[subjectHits(parad_ol), c("Sample", "Segment_Mean")]
    gene_hits <- parad_gr[queryHits(parad_ol), ]$external_gene_name
    
    # Aggregate
    parad_cna <- data.table(
        sample = cna_hits$Sample,
        segment_mean = cna_hits$Segment_Mean,
        gene_name = gene_hits
    )
    
    # Merge
    parad_cna <- merge(parad_cna, cna_tissues, by.x = "sample", by.y = "aliquot_id")
    # Only keep NT and TP
    parad_cna <- parad_cna[type %in% c("TP", "NT")]
    
    cur_plot <- ggplot(data = parad_cna, mapping = aes(x = gene_name, y = segment_mean, col = type)) +
        geom_point(size = 0.5) + 
        # Add mean bars
        stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean",
                     size= 0.3, geom = "crossbar", show.legend = F) +
        ylab("Segment Mean") + xlab("Paradoxical Genes") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        ggtitle(paste0("Paradoxical Segment Means in ", proj),
                subtitle = paste0("Genes in tumor tissue are paradoxically ", parad_type, "-expressed")) 
    
    ggsave(filename = paste0("Paradoxical_Plots/CNA_Plot_", proj, ".png"),
           plot = cur_plot, device = "png", dpi = 320, width = 20, units = "in")
}



idx <- 2
# Find and depict samples with paradoxical gene expression
sample_analysis <- function(idx) {
    # Select a project with paradoxical genes
    cur_parad <- readRDS(all_parad[idx])
    # Get the project name
    proj <- gsub(pattern = ".*(TCGA-\\w{1,4}).*",
                 replacement = "\\1",
                 x = all_parad[idx])
    
    # Find and load associated CNA and mRNA files
    mrna_match <- grepl(pattern = proj, x = all_exp)
    cna_match <- grepl(pattern = proj, x = all_cna)
    
    cur_mrna <- get(load(all_exp[mrna_match]))
    cur_cna <- get(load(all_cna[cna_match]))
    
    # Subset mRNA data for normal and tumor solid tissues
    cur_mrna <- cur_mrna[, cur_mrna$shortLetterCode %in% c("NT", "TP")]
    # Get tissue type order
    tissue_type <- data.table(ID = cur_mrna$barcode, type = cur_mrna$shortLetterCode)
    
    # Create a DESeq DataSet, design by tissue type
    cur_des <- DESeq2::DESeqDataSet(se = cur_mrna, design = ~ shortLetterCode)
    
    # We can remove the rows that have no or nearly no information about the
    # amount of gene expression (across all samples)
    cur_des <- cur_des[rowSums(counts(cur_des)) > 1, ]
    # Get the normalized counts (first estimate size factors)
    cur_des <- estimateSizeFactors(cur_des)
    norm_counts <- counts(cur_des, normalized = T)
    
    # Find patients that have a normal and tumor tissue type
    patients <- data.table(cur_mrna$patient, cur_mrna$shortLetterCode)
    patients[, sum_types := length(unique(V2)), by = "V1"]
    multi_sampled <- patients[sum_types > 1]$V1
    
    if (length(multi_sampled) < 2) {
        return(paste0(proj, " does not have mutli-sampled patients"))
    }
    
    # Find all associated barcodes of these patients
    # cur_mrna$barcode[cur_mrna$patient %in% multi_sampled]
    
    # Find the expression values for the paradoxical genes:
    # Paradoxical gain
    if (length(cur_parad$over_exp_parad) > 0) {
        cur_parad$over_exp_parad
        plot_exp_and_segmeans(cur_parad,
                              proj,
                              mrna_match,
                              cna_match,
                              cur_mrna,
                              cur_cna,
                              tissue_type,
                              norm_counts, "over", multi_sampled)
    }
    if (length(cur_parad$under_exp_parad) > 0) {
        plot_exp_and_segmeans(cur_parad,
                              proj,
                              mrna_match,
                              cna_match,
                              cur_mrna,
                              cur_cna,
                              tissue_type,
                              norm_counts, "under", multi_sampled)
    }
}

plot_exp_and_segmeans <- compiler::cmpfun(plot_exp_and_segmeans)
sample_analysis <- compiler::cmpfun(sample_analysis)

mc.results <-
    mclapply(
        X = 1:length(all_parad),
        FUN = sample_analysis,
        mc.preschedule = T,
        mc.cores = detectCores()
    )

