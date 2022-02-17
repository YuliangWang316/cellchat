library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

counts<-read.table("D:/jyc/AllSample.NormalizedCounts.txt",sep = "\t",header = TRUE,row.names = 1)
#counts_1<-counts[!duplicated(counts$Gene_symbol),]
#rownames(counts_1)<-counts_1$Gene_symbol
#counts_1<-counts_1[,-1]
#counts_1<-as.matrix(counts_1)
metadata<-read.table("f:/cellphoneDB/metadata.txt",sep = "\t",header = TRUE,row.names = 1)
counts_1<-counts[,rownames(metadata)]
counts_1<-as.matrix(counts_1)
#rownames(counts_1)<-tolower(rownames(counts_1))
#library(Hmisc)
#rownames(counts_1)<-capitalize(rownames(counts_1))

cellchat <- createCellChat(object = counts_1, meta = metadata,group.by = "cell_type")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
#future::plan("multiprocess", workers = 20)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
