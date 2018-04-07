# mirDIP_data_modification.R
require(data.table)
require(biomaRt)

# Load data from mirDIP, keep
mirdip <- fread("Fa_Ta_mirDIP_v_3_1520371014/mirDIP-All-Data-Version3.0.txt",
                showProgress = T, nrows = 10)
mirdip <- mirdip[V6 %in% "Top_Third"]
mirdip <- mirdip[, c(1, 2)]
mirdip <- unique(mirdip)
colnames(mirdip) <- c("GENENAME", "miRNA")
mirdip$miRNA <- gsub(pattern = "hsa-mir-(\\d{1,}\\w*)-*.*", replacement = "hsa-miR-\\1", x = mirdip$miRNA)
gene_dict <- fread("gene_dictionary.txt")

# Add ENSG IDs to mirDIP data
myMart = useMart(biomart = "ensembl",
                 dataset = "hsapiens_gene_ensembl")
cur_ensg <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "external_gene_name",
    values = unique(mirdip$GENENAME),
    mart = myMart
)
gene_dict <- gene_dict[!is.na(ensembl_gene_id)]
idx <- fastmatch::fmatch(x = mirdip$GENENAME, table = cur_ensg$external_gene_name)
mirdip$ENSG <- cur_ensg$ensembl_gene_id[idx]
fwrite(mirdip, "Fa_Ta_mirDIP_v_3_1520371014/Unique_mirdip.txt")