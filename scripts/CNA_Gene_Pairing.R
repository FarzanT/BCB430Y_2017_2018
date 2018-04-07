# CNA_Gene_Pairing.R

# The GDC mRNA quantification analysis pipeline measures gene level expression
# in   HT-Seq raw read count, Fragments per Kilobase of transcript per Million
# mapped reads (FPKM), and FPKM-UQ (upper quartile normalization). These values
# are generated through this pipeline by first aligning reads to the GRCh38
# reference genome and then by quantifying the mapped reads. To facilitate
# harmonization across samples, all RNA-Seq reads are treated as unstranded
# during analyses.

# RNA-Seq expression level read counts are normalized using two related methods:
# FPKM and FPKM-UQ. Normalized values should be used only within the context of
# the entire gene set. Users are encouraged to normalize raw read count values
# if a subset of genes is investigated.

# FPKM:
# The Fragments per Kilobase of transcript per Million mapped reads (FPKM)
# calculation normalizes read count by dividing it by the gene length and the
# total number of reads mapped to protein-coding genes.

# ==== Package Load ====
if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}
if (!require(stringr)) {
    install.packages("stringr")
    library(stringr)
}
if (!require(fastmatch)) {
    install.packages("fastmatch")
    library(fastmatch)
}
if (!require(compiler)) {
    install.packages("compiler")
    library(compiler)
}

# BioConductor packages
# source("https://bioconductor.org/biocLite.R")
if (!require(biocInstaller)) {
    biocLite("BiocInstaller")
    library(BiocInstaller)
}
# Ensure that Curl is installed if running on Linux
if (!require(BSgenome.Hsapiens.UCSC.hg38)) {
    BiocInstaller::biocLite("BSgenome.Hsapiens.UCSC.hg38")
    library(BSgenome.Hsapiens.UCSC.hg38)
}
hg38 <- BSgenome.Hsapiens.UCSC.hg38
genomeInfo = SeqinfoForBSGenome(hg38)[c(paste0("chr", seq(1, 22)), "chrX", "chrY")]
rm(hg38)

if (!require(GenomicRanges)) {
    BiocInstaller::biocLite("GenomicRanges")
    library(GenomicRanges)
}
if (!require(biomaRt)) {
    BiocInstaller::biocLite("biomaRt")
    library(biomaRt)
}
if (!require(SummarizedExperiment)) {
    BiocInstaller::biocLite("SummarizedExperiment")
    library(SummarizedExperiment)
}
if (!require(DESeq2)) {
    BiocInstaller::biocLite("DESeq2")
    library(DESeq2)
}
if (!require(TCGAbiolinks)) {
    install.packages("TCGAbiolinks")
    library(TCGAbiolinks)
}
if (!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
}


# ==== Data Preparation ====
mrnaFiles = list.files("SummarizedExperiments/mRNA/", full.names = T)
cnaFiles = list.files("SummarizedExperiments/CNA/", full.names = T)
# mirnaFiles = list.files("SummarizedExperiments/miRNA/", full.names = T)
# dir.create("AnnotatedExperiments")
dir.create("CNA_Gene_Pairs")
dir.create("Plots_Gene_Expression_PCA")
dir.create("DESeq_Results")
dir.create("CNA_Gene_Pairs/Full/")
dir.create("CNA_Gene_Pairs/Partial/")

# CNA and gene pairings are created during the following script:
cnaExperiment <- function(idx) {
    # Append project information to metadata
    mRNAdata = get(load(mrnaFiles[idx]))
    
    # NOTE: If the expression samples do not have normal types, discard
    # TODO: Is there a way around this?
    # We can retrieve the tissue type data from the shortLetterCode attribute
    if (sum(mRNAdata$shortLetterCode == "NT") < 1) {
        return("This project does not have normal tissue gene expression data")
    }
    
    # Extract project name
    proj = strsplit(basename(mrnaFiles[idx]), "_GeneExpression")[[1]][1]
    # Add a title to meta data
    # metadata(mRNAdata) = list(data_release = metadata(mRNAdata)$data_release,
    #                           Title = paste0("Raw mRNA counts of project ",
    #                                          proj, " from the GDC."))
    
    # We can retrieve the gene ranges from the SummarizedExperiment object, making
    # it easy to append corresponding CNA information (via rowRanges())
    gRanges = rowRanges(mRNAdata)
    
    # Update seqinfo
    seqinfo(gRanges) = genomeInfo[seqlevels(gRanges)]
    
    # Only retain the samples that also have corresponding mRNA data:
    # Extract the project, tissue source site, participant ID, tissue type
    # and portion analyte type
    # These are the main identifiers for a unique participant (!)
    # https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode
    
    # # Available mRNA participants and tissues
    # mrnaBarcodes = mRNAdata$sample
    # 
    # # Split barcodes into sub-identifiers
    # mRNAparticipant = as.data.table(str_split(mRNAdata$sample, pattern = "-",
    #                                           simplify = T)[ , 1:4])
    # 
    # # Split tissue type
    # tissueSplit = str_split(mRNAparticipant$V4, "", simplify = T)
    # 
    # # Retain the unique identifiers mentioned above
    # mRNAparticipantSplit = sort(paste(mRNAparticipant$V1,
    #                                   mRNAparticipant$V2,
    #                                   mRNAparticipant$V3,
    #                                   paste0(tissueSplit[,1], tissueSplit[,2]),
    #                                   sep = "-"))
    
    # Load the corresponding CNA file, convert to GRanges object for overlapping
    # later on
    projCNA = as.data.table(get(load(cnaFiles[grep(pattern = proj, cnaFiles)])))
    
    # Categorize SegMean > 1 as "gain", < 1 as "loss" and everything else as
    # normal
    # TODO: Are the current cutoffs sensible?
    # projCNA$Segment_Category = character()
    # projCNA[projCNA$Segment_Mean > 1]$Segment_Category = "gain"
    # projCNA[projCNA$Segment_Mean < -1]$Segment_Category = "loss"
    # projCNA[is.na(projCNA$Segment_Category)]$Segment_Category = "norm"
    
    # # Available CNA participants and tissues
    # cnaBarcodes = sort(unique(projCNA$Sample))
    # # Do the same as above
    # CNAparticipant = as.data.table(str_split(cnaBarcodes, pattern = "-",
    #                                          simplify = T))
    # # TODO: If there are multiple vials/plates, only pick the first one!
    # # i.e. remove the participant duplicates
    # # Split the sample and vial column
    # CNAparticipant$V8 = str_split(CNAparticipant$V4, "\\D", simplify = T)[,1]
    # CNAparticipant = CNAparticipant[!duplicated(CNAparticipant[, c(1:3, 8)])]
    # 
    # tissueSplit = str_split(CNAparticipant$V4, "", simplify = T)
    # CNAparticipantSplit = sort(paste(CNAparticipant$V1,
    #                             CNAparticipant$V2,
    #                             CNAparticipant$V3,
    #                             paste0(tissueSplit[,1], tissueSplit[,2]),
    #                                    sep = "-"))
    # 
    # # Find CNA samples that correspond with our mRNA data
    # CNAmatch = CNAparticipantSplit[CNAparticipantSplit %in% mRNAparticipantSplit]
    # # Now that we have samples with both mRNA and CNA data, must match them back
    # # to their original full identifiers in TCGA
    # pattern = paste(CNAmatch, collapse = "|")
    # CNAparticipant = paste(
    #     CNAparticipant$V1, CNAparticipant$V2, CNAparticipant$V3,
    #     CNAparticipant$V4, CNAparticipant$V5, CNAparticipant$V6,
    #     CNAparticipant$V7,
    #     sep = "-")
    # 
    # fullCNAmatch = CNAparticipant[grep(pattern = pattern, x = CNAparticipant)]
    # 
    # # Subset CNA data by matching IDs
    # projCNA = projCNA[Sample %in% fullCNAmatch]
    
    # Check whether all segment ranges are valid
    # e.g. if end position < start position, the sample is invalid
    if (!all(projCNA$End >= projCNA$Start)) {
        # Find and remove the offending sample entirely
        invalids = unique(projCNA$Sample[projCNA$End < projCNA$Start])
        # Remove invalids from CNA data
        projCNA <- projCNA[!(Sample %in% invalids)]
    }
    # Prepend "chr" to the chromosome column of CNA data
    projCNA$Chromosome = paste0("chr", projCNA$Chromosome)
    # Create a GRanges for overlapping operations
    projCNA = makeGRangesFromDataFrame(df = projCNA, keep.extra.columns = T,
                                       start.field = "Start", end.field = "End",
                                       seqnames.field = "Chromosome")
    
    # Overlap with gene ranges, differentiate partial and full overlaps
    # (It is faster to find the overlaps across all samples first then subset/
    # match with the mRNA data rather than finding overlaps for each sample
    # at a time.)
    
    # Full overlaps
    fullOL = as.data.table(findOverlaps(query = gRanges, subject = projCNA,
                                        type = "within", ignore.strand = T))
    # Partial overlaps
    partialOL = as.data.table(findOverlaps(query = gRanges, subject = projCNA,
                                           type = "any", ignore.strand = T))
    
    # Find difference to be able to indicate which are partial _only_
    # Using the following trick: add a dummy column to each type of overlap,
    # and then merge. The ones that are differing will have an NA column
    partialOL$temp1_specific = T
    fullOL$temp2_specific = T
    dif = merge(partialOL, fullOL, all = T)
    rm(partialOL)
    dif = dif[is.na(dif$temp2_specific), c("queryHits", "subjectHits")]
    
    ## Create GRanges for full overlaps
    fOL = gRanges[fullOL$queryHits, 0]
    myMeta1 = as.data.table(elementMetadata(projCNA[fullOL$subjectHits, ]))
    # Only keep aliquot ID, Num_probes and segment mean
    myMeta1$`V` = "full"
    colnames(myMeta1) = c("aliquot_barcode", "num_probes", "segment_mean",
                          "overlap_type")
    elementMetadata(fOL) = cbind(elementMetadata(fOL), myMeta1)
    rm(list = c("fullOL", "myMeta1"))
    
    ## Create GRanges for partial overlaps
    pOL = gRanges[dif$queryHits, 0]
    myMeta2 = as.data.table(elementMetadata(projCNA[dif$subjectHits, ]))
    # Only keep aliquot ID, Num_probes and segment mean
    myMeta2$`V` = "partial"
    colnames(myMeta2) = c("aliquot_barcode", "num_probes", "segment_mean",
                          "overlap_type")
    elementMetadata(pOL) = cbind(elementMetadata(pOL), myMeta2)
    # Remove large 'dif' object to free memory
    rm(list = c("dif", "myMeta2"))
    
    # Concatenate the GRanges, remove unneeded files from the environment
    colnames(fOL)
    # Create a data.table with the overlaps
    full <- data.table(ensembl_gene_id = names(fOL),
                        aliquot_barcode = fOL$aliquot_barcode,
                        segment_mean = fOL$segment_mean)
    partial <- data.table(ensembl_gene_id = names(pOL),
                       aliquot_barcode = pOL$aliquot_barcode,
                       segment_mean = pOL$segment_mean)
    min(partial$segment_mean)
    min(full$segment_mean)
    
    # final <- append(fOL, pOL)
    # rm(list = c("fOL", "pOL"))
    
    # # Match corresponding participant and sample types in mRNA and CNA data
    # match = sort(mRNAparticipantSplit[mRNAparticipantSplit %in% CNAparticipantSplit])
    # pattern = paste(match, collapse = "|")
    # 
    # # Find matches in mRNA and CNA data
    # dict = data.table(CNAmatch = character(), mRNAmatch = character())
    # for (i in 1:length(match)) {
    #     # Only get the first matching element
    #     # TODO: Should other matching elements also be considered? i.e. for 
    #     # patients that have multiple samples
    #     cur_CNA = cnaBarcodes[grep(pattern = match[i],
    #                                           x = cnaBarcodes, ignore.case = T)[1]]
    #     cur_mRNA = mrnaBarcodes[grep(pattern = match[i],
    #                                             x = mrnaBarcodes, ignore.case = T)[1]]
    #     # Append to a 1-to-1 dictionary for mRNA and CNA matches
    #     dict <- rbindlist(list(dict, data.table(CNAmatch = cur_CNA,
    #                                        mRNAmatch = cur_mRNA)))
    # }
    # dict <- unique(dict)
    # rm(list = c("CNAmatch", "mRNAmatch"))
    
    # Note that some CNA segments overlap, resulting in multiple segment_mean
    # values recored for a gene in that sample
    # TODO: For now, take the largest segment mean for a gene that falls under
    # this condition
    # full_names <- names(final[final$overlap_type == "full"])
    # partial_names <- names(final[final$overlap_type == "partial"])
    
    # TODO: For now, consider the highest segment mean of a gene in a given
    # sample as the segment mean representative of that gene
    full[, max_seg_mean := max(segment_mean),
          by = c("ensembl_gene_id", "aliquot_barcode")]
    partial[, max_seg_mean := max(segment_mean),
         by = c("ensembl_gene_id", "aliquot_barcode")]
    
    # Subset
    full <- unique(full[, c(1, 2, 4)])
    partial <- unique(partial[, c(1, 2, 4)])
    
    # Convert to wide format, ready to convert to a SummarizedExperiment object
    full_wide <- dcast.data.table(data = full,
                               formula = ensembl_gene_id ~ aliquot_barcode,
                               value.var = "max_seg_mean")
    partial_wide <- dcast.data.table(data = partial,
                                     formula = ensembl_gene_id ~ aliquot_barcode,
                                     value.var = "max_seg_mean")
    
    rm(list = c("full", "partial"))
    
    # Merge assayResults to create another assay variable for the current
    # "SummarizedExperiment":
    # Convert to matrix and a SummarizedExperiment
    full_sum_expr <- SummarizedExperiment(as.matrix(full_wide[ , -1]))
    # Add feature names
    rownames(full_sum_expr) <- full_wide$ensembl_gene_id
    assayNames(full_sum_expr) <- "SNP 6.0 Copy NUmber Aberrations"
    metadata(full_sum_expr) <- list(title = paste0("Full gene overlap CNA data for project", proj))
    
    # Do the same for partial overlaps    
    partial_sum_expr <- SummarizedExperiment(as.matrix(partial_wide[ , -1]))
    rownames(partial_sum_expr) <- partial_wide$ensembl_gene_id
    assayNames(partial_sum_expr) <- "SNP 6.0 Copy NUmber Aberrations"
    metadata(partial_sum_expr) <- list(title = paste0("Partial gene overlap CNA data for project", proj))
    
    # Save each as an RDS file
    # saveRDS(mRNAdata, paste0("AnnotatedExperiments/", proj, "_GE_and_CNA.rds"))
    saveRDS(full_sum_expr, paste0("CNA_Gene_Pairs/Full/", proj,
                                  "_Full_Overlap_CNA_Gene_Pairs.rds"))
    saveRDS(partial_sum_expr, paste0("CNA_Gene_Pairs/Partial/",
                                     proj, "_Partial_Overlap_CNA_Gene_Pairs.rds"))
    
    return(paste0("Saved the results of project ", proj))
}


# Compile the function for slightly faster processing
cnaExperiment <- compiler::cmpfun(cnaExperiment)

# When memory becomes low, the Linux out of memory killer (oom killer)
# will start silently killing processes. It does not print anything to the
# console to let you know what it's doing, although the oom killer activities
# show up in the system log. In this situation the output of mclapply will
# appear to have been randomly contaminated with NULLS.

results <- parallel::mclapply(X = 1:length(cnaFiles), FUN = cnaExperiment,
                              mc.preschedule = T,
                              mc.cores = 2)


# [END]