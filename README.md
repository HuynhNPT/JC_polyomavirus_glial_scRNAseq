# JC polyomavirus propagation in glial cells - scRNAseq analysis

## **Abstract**
Progressive multifocal leukoencephalopathy (PML) is a demyelinating infection of the immunosuppressed brain, mediated by the gliotropic polyomavirus JCV. JCV replicates in human glial progenitor cells and astrocytes, which undergo viral T antigen-triggered mitosis, enabling viral replication. We asked if JCV spread might therefore be accelerated by glial proliferation. We report here that dividing human astrocytes support JCV propagation to a substantially greater degree than do mitotically quiescent cells. Accordingly, bulk and single cell RNA-seq reveal that JCV-infected glia exhibit cycle-linked disruption ofboth DNA damage response and transcriptional regulatory pathways. In vivo, JCV infection of humanized glial chimeras is greatly accentuated by cuprizone-induced demyelination and its associated mobilization of GPCs. Importantly, such infection triggers the death of uninfected as well as infected glia in vivo, reflecting significant bystander death. Together, these results suggest that JCV propagation in PML may be potentiated by glial cell division. As a result, the accentuated glial proliferation attending disease-associateddemyelination may provide an especially favorable environment for JCV propagation, thus potentiating oligodendrocytic bystander death and further accelerating demyelination in susceptible hosts.

## **Scripts**
Follow the scripts below in order to generate results of this manuscript. 
1. [Import data]( https://rawcdn.githack.com/HuynhNPT/JC_polyomavirus_glial_scRNAseq/main/01_Import_data.html)</br>

2. [Perform integration]( https://rawcdn.githack.com/HuynhNPT/JC_polyomavirus_glial_scRNAseq/main/02_scVI_integration.html)</br>

3. [Find differentially expressed markers for each cell cluster]( https://rawcdn.githack.com/HuynhNPT/JC_polyomavirus_glial_scRNAseq/main/03_cluster_markers.html)</br>

4. Run pySCENIC to derive the gene regulatory network with scripts `04_exportSCENIC.R`, [05]( https://rawcdn.githack.com/HuynhNPT/JC_polyomavirus_glial_scRNAseq/main/05_saveModules.html), [06]( https://rawcdn.githack.com/HuynhNPT/JC_polyomavirus_glial_scRNAseq/main/06_AUCCell.html)</br>

5. Find [differential activity]( https://rawcdn.githack.com/HuynhNPT/JC_polyomavirus_glial_scRNAseq/main/07_Differential_Activity.html) and [differential expression]( https://rawcdn.githack.com/HuynhNPT/JC_polyomavirus_glial_scRNAseq/main/08_Differential_Expression.html) </br>

6. Deconstruct the [cell-cell communication network]( https://rawcdn.githack.com/HuynhNPT/JC_polyomavirus_glial_scRNAseq/main/09_CellChat.html)</br>



## **Data**
Processed `filtered_feature_bc_matrix` can be downloaded from GEO - Accession number: GSE144626
