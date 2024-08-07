library(CellChat)
# Set options to allow 1Gb of RAM
options(future.globals.maxSize = 1000*1024^2)
# Register parallel
future::plan("multiprocess", workers = 4)
# Set ligand-receptor database
CellChatDB <- CellChatDB.human

##### Merged cell chats with infected and uinfected meta tagged to cell types
# This object is created from PART 2 of script 09_CellChat
cellChat
print(">>>>>>>>>>>>>>>>>>>>")
cat("\n")
token <- "merge"
cat(token)
cat("\n")
print(cellChat)
cellChat@DB <- CellChatDB
# Subset the expression data of signaling genes for saving computational cost
# This step is necessary even if using the whole database
# In this example, we are using all genes
cellChat <- subsetData(cellChat, features = NULL)
# Identification of overexpressed genes and interactions
cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)

# # Compute the communication probability and infer cellular communication network
cellChat <- computeCommunProb(cellChat) # Default uses trimean method
# # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellChat <- filterCommunication(cellChat, min.cells = 5)
# cell-Cell communication on a signlaing pathway level (i.e. bigger picture)
cellChat <- computeCommunProbPathway(cellChat)
# Aggregated ccc network
cellChat <- aggregateNet(cellChat)
# Save work and clean up
saveRDS(cellChat, "cellChat_merge.RDS")