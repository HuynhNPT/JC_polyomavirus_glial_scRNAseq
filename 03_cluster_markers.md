---
title: "Find cluster markers and cell type markers"
date: "`r format(Sys.Date())`"
output:
  rmarkdown::html_document:
    highlight: pygments
    toc: true
    toc_depth: 3
    fig_width: 5
    toc_float: true
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding[utf8]{inputenc}
abstract: "Script requres 01 and 02 to be run first. bluehive r/4.1.1/b1"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, eval = TRUE, message = FALSE, warning = FALSE)
```

## PART 0: Load libraries

```{r}
.libPaths("/gpfs/fs1/sfw2/r/4.1.1/b1/lib64/R/library")
library(Seurat)
library(SeuratDisk)
library(cowplot)
library(ggplot2)
library(rhdf5)
main_pal <- c("grey", viridis::viridis_pal(option = "B")(300))
source("_helper.R")
```

## PART 1: Data 

Note that the imported object was produced from `02_scVI_integration`

```{r}
Convert("git_integrated_object_postscVI.h5ad", "git_postScanpy.h5seurat", assay = "RNA", overwrite = TRUE)
# Updates in October 2023 gives error with old codes without misc and meta.data specified to FALSE
obj <- LoadH5Seurat("git_postScanpy.h5seurat", assay = "RNA")
file.remove("git_postScanpy.h5seurat")
```

## PART 2: Cluster markers 

```{r}
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- SetIdent(obj, value = "leiden")
```

```{r eval = FALSE, echo = TRUE, include=TRUE}
cluster_markers <- FindAllMarkers(obj, assay = "RNA", test.use = "MAST", only.pos = TRUE)
write.table(cluster_markers, "tables/clusterMarkers/leidenMarkers.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
```

Leiden clusters were annotated for cell types based on the following markers: <br>
0: APC: PDGFRA, MDK, PTN <br>
1: Astrocyte: GFAP, AQP4, APOE <br>
2: cAPC:  <br>
3: rxAstrocyte: APOE, S100A10 <br>

```{r}
leiden_clusters <- c("0", "1", "2", "3")
cellType <- c("APC", "Astrocyte", "cAPC", "reactive_Astrocyte")
obj$cellType <- plyr::mapvalues(obj$leiden, leiden_clusters, cellType)
obj$cellType <- factor(obj$cellType, levels = c("cAPC", "APC", "Astrocyte", "reactive_Astrocyte"))
```

```{r, fig.width=15, fig.height=8}
goi <- c("MKI67", "TOP2A", "PDGFRA", "PTN", "GFAP", "AQP4", "APOE", "S100A10")
p1 <- DimPlot(obj, group.by = "cellType", label = TRUE)
p_list <- FeaturePlot(obj, goi, order = TRUE, combine = FALSE)
p_list <- lapply(p_list, function(x) x + scale_color_gradientn(colors = main_pal))
plot_grid(p1, p_list[[1]], p_list[[2]], p_list[[3]], p_list[[4]], p_list[[5]], p_list[[6]], p_list[[7]], p_list[[8]], ncol = 3)
```

```{r, fig.width=15, fig.height=8}
VlnPlot(obj, goi, group.by = "cellType", pt.size = 0, ncol = 4)
```

```{r}
obj <- SetIdent(obj, value = "cellType")
```

## PART 3: Add more meta data information to this object 

### Cell cycle 

```{r}
# Read in cc markers
s.genes = cc.genes$s.genes
s.genes[s.genes == "MLF1IP"] <- "CENPU"
g2m.genes = cc.genes$g2m.genes
g2m.genes[g2m.genes == "FAM64A"] <- "PIMREG"
g2m.genes[g2m.genes == "HN1"] <- "JPT1"

obj<- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
```

### Infection

```{r}
# Cells that are infected and cells that are not infected 
Jvgp <- obj[["RNA"]]@counts[grep("Jvgp", row.names(obj[["RNA"]]@counts)),]
Jvgp <- apply(Jvgp, 2, sum)
# Call cells that have a total sum of all detected viral RNAs less than 15 - uninfected 
# Why did we decide on 15? Histogram 
ggplot(as.data.frame(Jvgp), aes(x=Jvgp)) + geom_histogram(bins = 50) + xlim(0, 50)
Jvgp <- ifelse(Jvgp < 15, "Uninfected", "Infected")
obj$Infection <- Jvgp
```

## Save objects

```{r}
saveRDS(obj, "git_obj.RDS")
```

## Software

```{r echo=FALSE}

sessionInfo()
```