# mRNA_Analysis.R

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


# A function to load a file to a temporary environment, to allow name assignment
load_obj <- function(f) {
    env <- new.env()
    nm <- load(f, env)[1]
    env[[nm]]
}
pcaPlot <- function(object, projectName, intgroup = NULL, ntop = 500) {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select,]))
    percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup,
                                                 drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = ":"))
    } else {
        colData(object)[[intgroup]]
    }
    d <- data.table(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = group,
        intgroup.df,
        name = colnames(object)
    )
    myPlot <- ggplot(data = d, aes(x = PC1, y = PC2, color = group)) +
        geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] *
                                                              100), "% variance")) + 
        ylab(paste0("PC2: ", round(percentVar[2] *
                                       100), "% variance")) + coord_fixed() +
        ggtitle(label = paste0("PCA plot of ", projectName),
                subtitle = paste0("Total number of samples used: ",
                                  ncol(object),
                                  "\nNumber of genes considered for PCA: ",
                                  ntop))+
        # Modify legend
        scale_color_discrete(name = "Tissue Type") +
        # Modify plot details
        theme(title = element_text(size = 10))
    return(myPlot)
}
pcaPlot <- compiler::cmpfun(pcaPlot)

# ==== Data Preparation ====
mrnaFiles = list.files("SummarizedExperiments/mRNA/", full.names = T)
cnaFiles = list.files("SummarizedExperiments/CNA/", full.names = T)
# mirnaFiles = list.files("SummarizedExperiments/miRNA/", full.names = T)
dir.create("AnnotatedExperiments")
dir.create("Plots Gene Expression PCA")
dir.create("DESeq Results")

# Two file types are updated/created during the following script:
#       CNA data tables
#       mRNA's annotated with corresponding CNA data

annotateExperiment <- function(idx) {
    # Append project information to metadata
    mRNAdata = load_obj(mrnaFiles[idx])
    
    # NOTE: If the expression samples do not have normal types, discard
    # TODO: Is there a way around this?
    # We can retrieve the tissue type data from the shortLetterCode attribute
    if (sum(mRNAdata$shortLetterCode == "NT") == 0) {
        stop("This project does not have a normal tissue gene expression ")
    }
    
    # Extract project name
    proj = strsplit(basename(mrnaFiles[idx]), "_GeneExpression")[[1]][1]
    # Add a title to meta data
    metadata(mRNAdata) = list(data_release = metadata(mRNAdata)$data_release,
                          Title = paste0("Raw mRNA counts of project ",
                                         proj, " from the GDC."))
    
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
    
    # Available mRNA participants and tissues
    mrnaBarcodes = mRNAdata$sample
    
    # Split barcodes into sub-identifiers
    mRNAparticipant = as.data.table(str_split(mRNAdata$sample, pattern = "-",
                                              simplify = T)[ , 1:4])
    
    # Split tissue type
    tissueSplit = str_split(mRNAparticipant$V4, "", simplify = T)
    
    # Retain the unique identifiers mentioned above
    mRNAparticipantSplit = sort(paste(mRNAparticipant$V1,
                                 mRNAparticipant$V2,
                                 mRNAparticipant$V3,
                                 paste0(tissueSplit[,1], tissueSplit[,2]),
                                 sep = "-"))
    
    # Use the following approach to normalize and run differential expression 
    # analysis
    # https://www.bioconductor.org/help/workflows/rnaseqGene/#differential-expression-analysis
    
    # Create a DESeq DataSet, indicate that the counts are dependent on the
    # tissue type
    # cur_des <- DESeq2::DESeqDataSet(se = mRNAdata, design = ~ shortLetterCode)
    
    # We can remove the rows that have no or nearly no information about the
    # amount of gene expression
    # cur_des <- cur_des[rowSums(counts(cur_des)) > 1, ]
    
    # We must properly adjust the variance of our data, before we can consider
    # the distances between the samples (for PCA and visualization)
    # Variance Stabilizing Transformation (VST) is faster and works better on
    # larger datasets compared to rlog: They correct for counts without inflating
    # variance
    # Sequencing depth correction is also automatically done
    # cur_vst <- vst(cur_des, blind = F)
    
    # Generate a PCA plot using the function defined above
    # curPCAplot <- pcaPlot(cur_vst, projectName = proj,
    #                       intgroup = "shortLetterCode", ntop = 5000)
    
    # Save the resulting PCA plot
    # ggsave(plot = curPCAplot,
    #        filename = paste0("Plots Gene Expression PCA/Plot_PCA_GeneExpression_",
    #                          proj, ".png"),
    #        device = "png", dpi = 220)
    
    # Run differential gene expression analysis
    # system.time({
    #     cur_difEx <- DESeq(object = cur_des, parallel = T)
    # })
    # saveRDS(cur_difEx, file = paste0("DESeq Results/", proj))
    
    
    # Load the corresponding CNA file, convert to GRanges object for overlapping
    # later on
    projCNA = as.data.table(load_obj(cnaFiles[grep(pattern = proj, cnaFiles)]))
    
    # Categorize SegMean > 1 as "gain", < 1 as "loss" and everything else as
    # normal
    # TODO: Are the current cutoffs sensible?
    projCNA$Segment_Category = character()
    projCNA[projCNA$Segment_Mean > 1]$Segment_Category = "gain"
    projCNA[projCNA$Segment_Mean < -1]$Segment_Category = "loss"
    projCNA[is.na(projCNA$Segment_Category)]$Segment_Category = "norm"
    
    # Available CNA participants and tissues
    cnaBarcodes = sort(unique(projCNA$Sample))
    # Do the same as above
    CNAparticipant = as.data.table(str_split(cnaBarcodes, pattern = "-",
                                             simplify = T))
    # TODO: If there are multiple vials/plates, only pick the first one!
    # i.e. remove the participant duplicates
    # Split the sample and vial column
    CNAparticipant$V8 = str_split(CNAparticipant$V4, "\\D", simplify = T)[,1]
    CNAparticipant = CNAparticipant[!duplicated(CNAparticipant[, c(1:3, 8)])]
    
    tissueSplit = str_split(CNAparticipant$V4, "", simplify = T)
    CNAparticipantSplit = sort(paste(CNAparticipant$V1,
                                CNAparticipant$V2,
                                CNAparticipant$V3,
                                paste0(tissueSplit[,1], tissueSplit[,2]),
                                       sep = "-"))
    
    # Find CNA samples that correspond with our mRNA data
    CNAmatch = CNAparticipantSplit[CNAparticipantSplit %in% mRNAparticipantSplit]
    # Now that we have samples with both mRNA and CNA data, must match them back
    # to their original full identifiers in TCGA
    pattern = paste(CNAmatch, collapse = "|")
    CNAparticipant = paste(
        CNAparticipant$V1, CNAparticipant$V2, CNAparticipant$V3,
        CNAparticipant$V4, CNAparticipant$V5, CNAparticipant$V6,
        CNAparticipant$V7,
        sep = "-")
    
    fullCNAmatch = CNAparticipant[grep(pattern = pattern, x = CNAparticipant)]
    
    # Subset CNA data by matching IDs
    projCNA = projCNA[Sample %in% fullCNAmatch]
    
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
    final = append(fOL, pOL)
    rm(list = c("fOL", "pOL"))
    
    # Match corresponding participant and sample types in mRNA and CNA data
    match = sort(mRNAparticipantSplit[mRNAparticipantSplit %in% CNAparticipantSplit])
    pattern = paste(match, collapse = "|")
    
    # Find matches in mRNA and CNA data
    dict = data.table(CNAmatch = character(), mRNAmatch = character())
    for (i in 1:length(match)) {
        # Only get the first matching element
        # TODO: Should other matching elements also be considered? i.e. for 
        # patients that have multiple samples
        cur_CNA = cnaBarcodes[grep(pattern = match[i],
                                              x = cnaBarcodes, ignore.case = T)[1]]
        cur_mRNA = mrnaBarcodes[grep(pattern = match[i],
                                                x = mrnaBarcodes, ignore.case = T)[1]]
        # Append to a 1-to-1 dictionary for mRNA and CNA matches
        dict <- rbindlist(list(dict, data.table(CNAmatch = cur_CNA,
                                           mRNAmatch = cur_mRNA)))
    }
    dict <- unique(dict)
    # rm(list = c("CNAmatch", "mRNAmatch"))
    
    # Note that some CNA segments overlap, resulting in multiple segment_mean
    # values recored for a gene in that sample
    # TODO: For now, take the largest segment mean for a gene that falls under
    # this condition
    largestMean = function(idx) {
        # Retrieve the sample
        curSample = final[final$aliquot_barcode == dict$CNAmatch[idx]]
        # Extract the max segment_mean of each gene
        curSample = data.table(ensembl_gene_id = names(curSample),
                               segment_mean = curSample$segment_mean)
        # Find the highest segment mean
        curSample[ , max_segment_mean := max(segment_mean),
                   by = "ensembl_gene_id"]
        # Discard segment_mean
        curSample = curSample[ , c(1,3)]
        curSample = unique(curSample)
        # Change column name to reflect the corresponding mRNA aliquot ID
        colnames(curSample) = c("ensembl_gene_id", dict$mRNAmatch[idx])
        
        return(curSample)
    }
    
    assayResults = mclapply(X = 1:nrow(dict),
                            FUN = largestMean, mc.preschedule = T,
                            mc.cores = detectCores())
    # rm(final)
    # Merge assayResults to create another assay variable for the current
    # "SummarizedExperiment":
    # Create place holder gene names that match the mRNA assay
    skeleton <- data.table(ensembl_gene_id = rowData(mRNAdata)$ensembl_gene_id)
    for (i in 1:length(assayResults)) {
        # Merge, keep all rows of skeleton for future merges (all.x)
        # Merge only if the assay result appears in mRNA data columns
        if (!is.null(assayResults[[i]]) && colnames(assayResults[[i]])[2] %in% colnames(mRNAdata)) {
            skeleton <- merge(skeleton, assayResults[[i]], all.x = T)
        }
    }
    # rm(assayResults)
    # Append the unmatched aliquot_id's to the skeleton
    unmatched = colnames(mRNAdata)[!(colnames(mRNAdata) %in% colnames(skeleton))]
    if (length(unmatched) > 0) {
        skeleton[ , unmatched] <- NA
    }
    sum(colnames(skeleton) %in% colnames(mRNAdata))
    ncol(mRNAdata); ncol(skeleton)
    colnames(mRNAdata)[1:10]; colnames(skeleton)[1:10]
    # Order columns to match SummarizedExperiment, and add unmatched samples at
    # the end
    setcolorder(x = skeleton, neworder = c("ensembl_gene_id", colnames(mRNAdata)))
    
    # Convert to matrix
    assay2 = as.matrix(skeleton[ , -1])
    rownames(assay2) = skeleton$ensembl_gene_id
    # sum(colnames(assay2) %in% colnames(mRNAdata)); ncol(mRNAdata)
    
    rm(skeleton)
    
    # Append the CNA segment information to the SummarizedExperiment
    # Test:
    # all(dimnames(assay(mRNAdata))[[1]] == dimnames(assay2)[[1]])
    # all(dimnames(assay(mRNAdata))[[2]] == dimnames(assay2)[[2]])
    
    assay(mRNAdata, 2) = assay2
    assayNames(mRNAdata) = c(assayNames(mRNAdata)[1], "CNA - DNAcopy")
    rm(assay2)
    
    # Reorder assay columns
    # assay
    # Save as RDS

    saveRDS(mRNAdata, paste0("AnnotatedExperiments/", proj, "_GE_and_CNA.rds"))
}

# results <- parallel::mclapply(X = 1:length(cnaFiles), FUN = annotateExperiment,
                              # mc.preschedule = T, mc.cores = detectCores(logical = F))

# Compile the function for slightly faster processing
annotateExperiment <- compiler::cmpfun(annotateExperiment)

for (i in 1:length(cnaFiles)) {
    annotateExperiment(i)
}


# [END]