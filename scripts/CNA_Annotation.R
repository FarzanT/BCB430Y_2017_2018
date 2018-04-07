# CNA_Annotation.R

# ==== Preparation ====
# CRAN packages
if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}
if (!require(parallel)) {
    install.packages("parallel")
    library(parallel)
}

# BioConductor packages
source("https://bioconductor.org/biocLite.R")
if (!require(biocInstaller)) {
    biocLite("BiocInstaller")
    library(BiocInstaller)
}

if (!require(GenomicRanges)) {
    biocLite("GenomicRanges")
    library(GenomicRanges)
}
if (!require(biomaRt)) {
    biocLite("biomaRt")
    library(biomaRt)
}
if (!require(BSgenome.Hsapiens.UCSC.hg38)) {
    biocLite("BSgenome.Hsapiens.UCSC.hg38")
    library(BSgenome.Hsapiens.UCSC.hg38)
}
hg38 = BSgenome.Hsapiens.UCSC.hg38

# Use the biomaRt package to retrieve any overlapping genes given start
# and positions:
myMart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                 dataset = "hsapiens_gene_ensembl")

dir.create("CNA_Annotated")
# List all non-log files
allFiles <-  list.files("GDC_Data/", recursive = T, full.names = T, pattern = "(seg.txt)$")
length(allFiles)

annotateSample <- function(idx) {
    curSample = fread(allFiles[idx])
    # First, create a GenomicRanges object for the sample data
    curSampleGR = makeGRangesFromDataFrame(df = curSample, keep.extra.columns = T,
                                         ignore.strand = T,
                                         seqnames.field = "Chromosome",
                                         start.field = "Start", 
                                         end.field = "End")
    # Add sequence information
    seqinfo(curSampleGR) = Seqinfo(seqnames = as.character(seqnames(curSampleGR)@values),
                                 genome = rep("hg38", length(unique(curSample$Chromosome))),
                                 isCircular = rep(F, length(unique(curSample$Chromosome))))
    
    # Query:
    queryResult = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol",
                                       "gene_biotype", "strand", "chromosome_name",
                                       "start_position", "end_position"),
                        filters = c("chromosome_name", "start", "end"),
                        values = list(as.character(seqnames(curSampleGR)),
                                      start(curSampleGR), end(curSampleGR)),
                        mart = myMart, checkFilters = T)
    
    # Create a GRanges object from this result, so that it can be overlapped
    # with the sample data:
    result = as.data.table(queryResult)
    result$strand = ifelse(test = result$strand == 1, yes = "+", no = "-")
    resultGR = makeGRangesFromDataFrame(df = result, keep.extra.columns = T,
                                        ignore.strand = F,
                                        seqnames.field = "chromosome_name",
                                        start.field = "start_position",
                                        end.field = "end_position",
                                        strand.field = "strand")
    # Overlap:
    anyOL = findOverlaps(query = resultGR, subject = curSampleGR, type = "any")
    fullOL = findOverlaps(query = resultGR, subject = curSampleGR,
                          type = "within")
    
    # Find difference
    temp1 = as.data.table(anyOL)
    temp1$temp1_specific = T
    temp2 = as.data.table(fullOL)
    temp2$temp2_specific = T
    dif = merge(temp1, temp2, all = T)
    dif = dif[is.na(dif$temp2_specific), c("queryHits", "subjectHits")]
    
    # Create GRanges for full overlaps with annotations
    temp1 = resultGR[queryHits(fullOL), ]
    myMeta1 = as.data.table(curSampleGR[subjectHits(fullOL), ])
    myMeta1 = cbind(myMeta1[ , !c("seqnames", "strand"), with = F],
                    rep("full", nrow(myMeta1)))
    colnames(myMeta1) = c("segment_start", "segment_end", "segment_width",
                          "sample_id", "num_probes", "segment_mean",
                          "overlap_type")
    elementMetadata(temp1) = cbind(elementMetadata(temp1), myMeta1)
    
    # Create GRanges for partial overlaps with annotations
    temp2 = resultGR[dif$queryHits, ]
    myMeta2 = as.data.table(curSampleGR[dif$subjectHits, ])
    myMeta2 = cbind(myMeta2[ , !c("seqnames", "strand"), with = F],
                    rep("partial", nrow(myMeta2)))
    colnames(myMeta2) = c("segment_start", "segment_end", "segment_width",
                          "sample_id", "num_probes", "segment_mean",
                          "overlap_type")
    elementMetadata(temp2) = cbind(elementMetadata(temp2), myMeta2)
    
    # Concatenate the GRanges
    gList = GRangesList(temp1, temp2)
    final = unlist(x = gList, recursive = TRUE, use.names = TRUE)
    
    # Append sequence information
    chrs = (seqlengths(hg38)[1:length(seqinfo(final))])
    # Discard chr in the names
    names(chrs) = substring(text = names(chrs), first = 4)
    finalSeqInfo = Seqinfo(seqnames = seqlevels(final),
                           isCircular = rep(F, length(seqinfo(final))),
                           seqlengths = chrs,
                           genome = rep("hg38", length(seqinfo(final))))
    seqinfo(final) = finalSeqInfo
    
    # Save with original file name
    id = paste(strsplit(x = allFiles[idx], "/")[[1]][c(3,4)], collapse = "_")
    path = paste0("CNA_Annotated/", id, ".rds")
    saveRDS(final, file = path)
}

# brcaFiles = list.files("GDC_Data/", recursive = T, full.names = T, pattern = "seg.txt")
# neckFiles = list.files("GDC_Data/TCGA-Head_and_Neck/", recursive = T, full.names = T)
# kidneyFiles = list.files("GDC_Data/TCGA-Kidney/", recursive = T, full.names = T)
# liverFiles = list.files("GDC_Data/TCGA-Liver/", recursive = T, full.names = T)
# lungFiles = list.files("GDC_Data/TCGA-Lung/", recursive = T, full.names = T)
# prostateFiles = list.files("GDC_Data/TCGA-Prostate/", recursive = T, full.names = T)
# skinFiles = list.files("GDC_Data/TCGA-Skin/", recursive = T, full.names = T)
# thyroidFiles = list.files("GDC_Data/TCGA-Thyroid/", recursive = T, full.names = T)

annotationResult = mclapply(X = 1:length(allFiles), FUN = annotateSample,
                            mc.preschedule = T, mc.cores = detectCores())

# [END]