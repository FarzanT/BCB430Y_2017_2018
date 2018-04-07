# DEA.R
# The following script runs Differential Expression Analysis on the gene
# expression counts obtained from TCGA

# ==== Package load ====
if (!require(SummarizedExperiment)) {
    BiocInstaller::biocLite("SummarizedExperiment")
    library(SummarizedExperiment)
}
if (!require(DESeq2)) {
    BiocInstaller::biocLite("DESeq2")
    library(DESeq2)
}
if (!require(limma)) {
    BiocInstaller::biocLite("limma")
    library(limma)
}
if (!require(BiocParallel)) {
    BiocInstaller::biocLite("BiocParallel")
    library(BiocParallel)
}
if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}
if (!require(compiler)) {
    install.packages("compiler")
    library(compiler)
}
if (!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
}

# ==== Function input preparation ====
SumExpFiles <- list.files("SummarizedExperiments/mRNA/", full.names = T)

# A function to create a PCA plot from gene expression data
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

# ==== DEA using DESeq2 ====
# The following function prepares the gene expression file, saves a PCA plot and
# then runs DEA
dir.create("DESeq_Results/")
RunDifExp <- function(idx) {
    # Load a SummarizedExperiment file containing the gene expression values
    mRNAdata <- get(load(SumExpFiles[idx]))
    
    # Extract project name
    proj <- strsplit(basename(SumExpFiles[idx]), "_GeneExpression")[[1]][1]
    
    # Add a title to meta data
    # metadata(mRNAdata) <- list(data_release = metadata(mRNAdata)$data_release,
    #                           Title = paste0("Raw mRNA counts of project ",
    #                                          proj, " from the GDC."))
    
    # NOTE: If the expression samples do not have normal types, discard
    # We can retrieve the tissue type data from the shortLetterCode attribute
    if (sum(mRNAdata$shortLetterCode == "NT") == 0) {
        cat("\nThis project does not have a normal tissue gene expression\n")
        return()
    }
    
    # Only compare between normal and solid tumor tissues (remove others)
    mRNAdata <- mRNAdata[, mRNAdata$shortLetterCode %in% c("NT", "TP")]
    
    # Use the following approach to normalize and run differential expression 
    # analysis
    # https://www.bioconductor.org/help/workflows/rnaseqGene/#differential-expression-analysis
    
    # Create a DESeq DataSet, indicate that the counts are dependent on the
    # tissue type
    cur_des <- DESeq2::DESeqDataSet(se = mRNAdata, design = ~ shortLetterCode)
    
    # We can remove the rows that have no or nearly no information about the
    # amount of gene expression (across all samples)
    cur_des <- cur_des[rowSums(counts(cur_des)) > 1, ]
    
    # We must properly adjust the variance of our data, before we can consider
    # the distances between the samples (for PCA and visualization)
    # Variance Stabilizing Transformation (VST) is faster and works better on
    # larger datasets compared to rlog: They correct for counts without inflating
    # variance
    # Sequencing depth correction is also automatically done
    cat("\nVariance Scaling for the PCA plot of", proj)
    cur_vst <- vst(cur_des, blind = F)
    
    # Generate a PCA plot using the function defined above
    cat("\nGenerating the PCA plot of", proj)
    curPCAplot <- pcaPlot(cur_vst, projectName = proj,
                          intgroup = "shortLetterCode", ntop = 20000)
    
    # Save the resulting PCA plot
    ggsave(plot = curPCAplot,
           filename = paste0("Plots Gene Expression PCA/Plot_PCA_GeneExpression_",
                             proj, ".png"),
           device = "png", dpi = 220)
    
    # Clear environment
    rm(list = c("cur_vst", "curPCAplot"))
    
    # Run differential gene expression analysis
    cat("\nRunning differential gene expression analysis for", proj)
    
    if (file.size(SumExpFiles[idx]) > 200e6) {
        # Tune the following based on the ammount of RAM available
        myBPPARAM <-
            BiocParallel::MulticoreParam(workers = detectCores(logical = F) - 2)
    } else {
        myBPPARAM <- BiocParallel::MulticoreParam(workers = detectCores())
    }
    
    cur_difEx <- DESeq(object = cur_des, parallel = T, BPPARAM = myBPPARAM)
    
    # Save results
    cat("\nSaving results for", proj)
    saveRDS(cur_difEx, file = paste0("DESeq_Results/DESeq_Results_", proj, ".rds"))
    
    # Stop the cluster and run a garbage collection
    bpstop(myBPPARAM)
    gc()
}

RunDifExp <- compiler::cmpfun(RunDifExp)

# Run differential expression analysis for the selected projects
# TCGA-HNSC failed
for (i in 1:length(SumExpFiles)) {
    # Only process files less than 200 MB for now
    if (file.size(SumExpFiles[i]) > 200e6) {
        RunDifExp(2)
    }
}

# ==== DEA using limma ====
dir.create("limma_Results")
LimmaDEA <- function(idx) {
    # Same procedure as the above function
    mRNAdata <- get(load(SumExpFiles[idx]))
    proj <- strsplit(basename(SumExpFiles[idx]), "_GeneExpression")[[1]][1]
    
    # Append project name to object in case the file name is changed
    if (!"Project" %in% names(metadata(mRNAdata))) {
        metadata(mRNAdata)[["Project"]] <- proj
    }
    
    # Save file (in rdata format)
    save(mRNAdata, file = SumExpFiles[idx])
    
    if (sum(mRNAdata$shortLetterCode == "NT") == 0) {
        cat("\nThis project does not have a normal tissue gene expression:", proj)
        return(NULL)
    }
    mRNAdata <- mRNAdata[, mRNAdata$shortLetterCode %in% c("NT", "TP")]
    
    # Convert the shortLetterCode to factors (limma requirement)
    mRNAdata$shortLetterCode <- factor(mRNAdata$shortLetterCode)
    
    # Create a design matrix
    curDesign <- model.matrix(~ mRNAdata$shortLetterCode)

    # Extract the counts from the SummarizedExperiment
    countData <- assay(mRNAdata, 1)
    
    # Remove genes with little or no expression values
    countData <- countData[rowSums(countData) > 1, ]
    
    # Transform RNA-Seq data for linear modelling
    # i.e. by assigning weights so that mean-variance relationship (e.g. higher
    # variance in higher count data) is offset
    cat("\nTransforming counts to CPM for project", proj)
    transformedCounts <- limma::voom(counts = countData,
                                     design = curDesign)
    
    # Fit the transformed object
    myLMFit <- lmFit(object = transformedCounts, design = curDesign)
    
    # Use empirical Bayes to calculate a moderated t-statistic...
    cat("\nRunning empirical Bayes for project", proj)
    eBayesFit <- eBayes(fit = myLMFit)
    
    # Append project related data (as a named vector)
    prop <- c(proj, sum(mRNAdata$shortLetterCode == "NT"),
              sum(mRNAdata$shortLetterCode == "TP"))
    names(prop) <- c("ProjectName", "NormalTissueCount", "TumorTissueCount")
    
    resultList <- list(prop, eBayesFit)
    cat("\nSaving raw results for project", proj)
    saveRDS(resultList, paste0("limma_Results/limma_Results_", proj, ".rds"))
}

# Compile
LimmaDEA <- compiler::cmpfun(LimmaDEA)

results <-
    mclapply(
        X = 1:length(SumExpFiles),
        FUN = LimmaDEA,
        mc.preschedule = T,
        mc.cores = detectCores(logical = F)
    )
