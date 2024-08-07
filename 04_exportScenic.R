##### Export matrix for SCENIC
library(Seurat)
# This object was generated from script 03_cluster_markers.Rmd
object <- readRDS("git_obj.RDS")
scenic <- object[["RNA"]]@counts

# Filter to retain only genes that are expressed more than minimal counts per gene threshold
minCountsPerGene <- 3*0.01*ncol(scenic)
sum_genes <- apply(scenic, 1, sum)
keep <- names(sum_genes)[which(sum_genes > minCountsPerGene)]
# Filter to retain only genes that are expressed in more than the minimal number of cells
minSamples <- ncol(scenic)*0.01
zeros <- apply(scenic, 1, function(x) sum(x>0))
keep2 <- names(zeros)[which(zeros > minSamples)]
keep <- keep[keep %in% keep2]
scenic <- as.matrix(scenic[keep,])

## Genes and cell info for AnnData
genes <- data.frame(row.names = keep, Gene = keep)
cells <- object@meta.data
######-------------
# write.csv(genes, "exported_matrices/geneInfo.csv", quote = FALSE)
# write.csv(cells, "exported_matrices/cellInfo.csv", quote = FALSE)
# write.csv(scenic, paste0("exported_matrices/forScenic.csv"), quote = FALSE)








