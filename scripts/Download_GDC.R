# Download_GDC.R

# ==== Preperation ====
if (!require(devtools)) {
    install.packages("devtools")
    library(devtools)
}
if (!require(SummarizedExperiment)) {
    install.packages("SummarizedExperiment")
    library(SummarizedExperiment)
}
if (!require(parallel)) {
    install.packages("parallel")
    library(parallel)
}
if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}
# Note that on Linux, libmariadb-client-lgpl-dev must be installed
# This allows for configuration of RMySQL
if (!require(TCGAbiolinks)) {
    devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
    library(TCGAbiolinks)
}

getGDCInfo()
# Data Release 10
# Move all data to a single directory later on
dir.create("SummarizedExperiments")
dir.create("SummarizedExperiments/miRNA")
dir.create("SummarizedExperiments/CNA")
dir.create("SummarizedExperiments/mRNA")

projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA', projects, perl=TRUE)]

# ==== Download miRNA data ====
plyr::alply(sort(projects), 1, function(proj) {
    tryCatch(
        query <- GDCquery(
            project = proj,
            data.category = "Transcriptome Profiling",
            data.type = "miRNA Expression Quantification"
        )
    )
    query <- GDCquery(project = proj,
                      data.category = "Transcriptome Profiling",
                      data.type = "miRNA Expression Quantification")
    GDCdownload(query, files.per.chunk = 20)
    GDCprepare(
        query,
        save = T,
        save.filename = paste0("SummarizedExperiments/miRNA/", proj,
                               "_microRNAexpression.rdata")
    )
})

# ==== Download CNA data ====
plyr::alply(sort(projects), 1, function(proj) {
    tryCatch(
        query <- GDCquery(
            project = proj,
            data.category = "Copy Number Variation",
            data.type = "Copy Number Segment",
            platform = "Affymetrix SNP 6.0"
        )
    )
    query <- GDCquery(
        project = proj,
        data.category = "Copy Number Variation",
        data.type = "Copy Number Segment",
        platform = "Affymetrix SNP 6.0"
    )
    GDCdownload(query, files.per.chunk = 20)
    GDCprepare(query,
               save = T,
               save.filename = paste0("SummarizedExperiments/CNA/", proj,
                                      "_CopyNumber.rdata"))
})

# ==== Download gene expression data ====
plyr::alply(projects, 1, function(proj) {
    # tryCatch(
    #     query <- GDCquery(
    #         project = proj,
    #         data.category = "Transcriptome Profiling",
    #         data.type = "Gene Expression Quantification",
    #         workflow.type = "HTSeq - Counts"
    #     )
    # )
    query <- GDCquery(
        project = proj,
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "HTSeq - Counts"
    )
    GDCdownload(query, files.per.chunk = 20)
    GDCprepare(
        query,
        save = T,
        save.filename = paste0("SummarizedExperiments/mRNA/", proj, "_GeneExpression.rdata")
    )
})


# We can use the rename function to move the file to a new path (note that the
# original files are deleted, use file.move() if this is not desired)
# oldPaths = list.files(pattern = ".*microrna*",
#                          full.names = T, ignore.case = T)
# newPaths = paste0("SummarizedExperiments/miRNA/", basename(oldPaths))
# file.rename(from = oldPaths, to = newPaths)
# 
# oldPaths = list.files(pattern = ".*copynumber*",
#                          full.names = T, ignore.case = T)
# newPaths = paste0("SummarizedExperiments/CNA/", basename(oldPaths))
# file.rename(from = oldPaths, to = newPaths)
# 
# oldPaths = list.files(pattern = ".*geneexpression.*",
#                          full.names = T, ignore.case = T)
# newPaths = paste0("SummarizedExperiments/mRNA/", basename(oldPaths))
# file.rename(from = oldPaths, to = newPaths)

# [END]