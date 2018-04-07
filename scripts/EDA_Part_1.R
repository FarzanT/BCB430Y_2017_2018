# EDA_Part_1.R
# ==== Package load =====
if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}
if (!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
}
if (!require(SummarizedExperiment)) {
    install.packages("SummarizedExperiment")
    library(SummarizedExperiment)
}

# ==== Data preparation ====
# List annotated files
annotated_experiments <- list.files("AnnotatedExperiments/", full.names = T)
load_obj <- function(f) {
    env <- new.env()
    nm <- load(f, env)[1]
    env[[nm]]
}


# ==== Plot by chromosome for each tumor type ====
# List all available files
all_files <- list.files("SummarizedExperiments/CNA/", full.names = T)
dir.create("Plot SegMean vs Chromosome Position")
dir.create("Plot Gain Loss per Tumor")

SegvChrom <- function(i) {
    curFile <- load_obj(f = all_files[i])
    # Extract the current tumor type from the file name
    curTumor <- gsub(pattern = "TCGA-(.*)_CopyNumber.rdata",
                     replacement = "\\1",
                     x = basename(all_files[i]))
    
    ### For each chromosome, plot a dot at "Start" on X-axis and "Segment_Mean"
    # on the Y-axis:
    # Extract the required columns, convert to long format
    curDT <- as.data.table(curFile[, c(1, 2, 3, 6)])
    # Shrink X axis to lie between 0 and 1
    curDT$Start <- curDT$Start / max(curDT$Start)
    curPlot <- ggplot(data = curDT,
                      mapping = aes(x = Start, y = Segment_Mean)) +
        geom_point(size = 0.1) +
        # Draw a red horizontal line on y = 0
        geom_hline(yintercept = 0, col = "red") + 
        facet_wrap(facets = ~ Chromosome,
                   scales = "free_x", ncol = 4) +
        xlab("Relative Chromosome Position") + ylab("Log2 Ratio Segment Mean") +
        ggtitle(label = paste0("Segment Mean vs Chromosome Position in ", curTumor),
                subtitle = paste0("Total number of samples: ",
                                  length(unique(curFile$Sample)))) +
        theme(axis.text.x = element_blank(),
              strip.text.x = element_text(size = 6, margin = margin(1, 1, 1, 1)))
    # Save plot to current tumor's folder
    ggsave(filename = paste0("Plot SegMean vs Chromosome Position/",
                             "/Plot_SegMean_vs_ChrPos_",
                             curTumor, ".png"),
           plot = curPlot, device = "png", dpi = 220)
    
    ### Create another plot for gain/loss ratio for each tumor type, across 
    # different chromosomes (the current cutoffs are -1 and 1)
    curDT[, Gain := sum(Segment_Mean > 1), by = Chromosome]
    curDT[, Loss := sum(Segment_Mean < -1), by = Chromosome]
    curDT[, Normal := sum(Segment_Mean < 1 &
                                  Segment_Mean > -1), by = Chromosome]
    
    curDT <- unique(curDT[, c(2, 3, 5, 6, 7)])
    # Find fraction of gain, loss and normal
    curDT$GainRatio = curDT$Gain / (curDT$Gain + curDT$Loss + curDT$Normal)
    curDT$LossRatio = curDT$Loss / (curDT$Gain + curDT$Loss + curDT$Normal)
    curDT$NormalRatio = curDT$Normal / (curDT$Gain + curDT$Loss + curDT$Normal)

    pie_data <- curDT[, c("Chromosome", "GainRatio", "LossRatio", "NormalRatio")]
    
    pie_molten <- melt.data.table(pie_data, id.vars = "Chromosome",
                                  variable.name = "Type", value.name = "Ratio")
    pie_molten <- unique(pie_molten)
    
    # Plot a pie chart chromosome-wise
    piePlot <- ggplot(pie_molten, mapping = aes(x = "", y = Ratio, fill = Type)) +
        # Bar plot
        geom_bar(stat = "identity", position = "stack") +
        # Polar coordinates convert bar charts to pie charts
        coord_polar(theta = "y", start = 0) +
        facet_wrap(facets = ~ Chromosome, ncol = 6) +
        # Axis labels
        xlab("") + ylab("") + 
        # Title
        ggtitle(paste0("Relative ratio of copy number gains vs losses\nper chromosome in ", curTumor),
                subtitle = paste0("Total number of samples: ",
                                  length(unique(curFile$Sample)))) +
        # Modify legend labels
        scale_fill_discrete(name = "Segment change",
                            labels = c("Gain", "Loss", "Normal")) +
        # Plot settings
        theme(title = element_text(size = 9),
              axis.text.y = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              strip.text.x = element_text(size = 6, margin = margin(1, 1, 1, 1)),
              legend.position = "bottom")
    ggsave(filename = paste0("Plot Gain Loss per Tumor/PieChart_GainLoss_Ratio_",
                             curTumor, ".png"),
           plot = , piePlot, device = "png", dpi = 220)
}

# Compile for better speed
SegvChrom <- compiler::cmpfun(SegvChrom)
results <- mclapply(X = 1:length(all_files), FUN = SegvChrom,
                    mc.preschedule = T, mc.cores = detectCores(logical = F))


# ==== Depict segment mean precedence across all genes in each tumor type ====
# Create a directory for the resulting plots
dir.create("Plots Segment Mean Density")

SegPlot <- function(i) {
    # Read an annotated tumor file
    curFile <- readRDS(annotated_experiments[i])
    # Extract the current tumor type from the file name
    curTumor <- gsub(pattern = "_GE_and_CNA.rds", replacement = "",
                     x = basename(annotated_experiments[i]))
    # Extract the segment means from the second "assay" in the summarized experiment
    max_seg_means <- as.data.table(assay(curFile, 2), keep.rownames = T)
    # The first column has the row (feature) names
    
    # Convert assay data to long format
    seg_molten <- melt.data.table(data = max_seg_means, id.vars = c("rn"),
                                variable.name = "SampleID", value.name = "SegMean",
                                na.rm = T)
    colnames(seg_molten)[1] <- "GeneID"
    
    # Plot a distribution of the number of genes vs the max segment mean
    graph <- ggplot(data = seg_molten) +
        geom_density(mapping = aes(x = SegMean)) +
        ggtitle(paste0("Density Plot of Segment Mean in ", curTumor),
                subtitle = paste0("Number of samples: ",
                                  length(unique(seg_molten$SampleID)))) +
        xlab("Segment Means") + ylab("Density") +
        scale_x_continuous(limits = c(-2.5, 2.5))
    ggsave(plot = graph,
           filename = paste0("Plots Segment Mean Density/Density_", curTumor, ".png"),
           device = "png", dpi = 220)
}
SegPlot <- compiler::cmpfun(SegPlot)
results <- parallel::mclapply(X = 1:length(annotated_experiments), FUN = SegPlot,
                              mc.preschedule = T,
                              mc.cores = detectCores(logical = F))

# ==== Depict differential gene expression across different tumor types ====

