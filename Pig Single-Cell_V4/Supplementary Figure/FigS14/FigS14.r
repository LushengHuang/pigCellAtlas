

rm(list=ls())
library("ArchR")
addArchRThreads(threads=10)
library("org.Ss.eg.db")
library("BSgenome")
library("rtracklayer")
library("GenomicRanges")
library("Biostrings")
library("GenomicFeatures")
library("Cairo")
library("ggrastr")
library("GenomeInfoDb")
library("ensembldb")
library("biomaRt")
library("Seurat")
library("SummarizedExperiment")
library("BSgenome.Sscrofa.ensembl.susScr11")

# ---------------------Step 1 create annotation for the genome and genes 
geneAnnotation = readRDS("geneAnnotation_new.rds")
genomeAnnotation = readRDS("genomeAnnotation.rds")

# -----------------    Step 2 import data
# ArchR uses "scanTabix" to read fragment files and "scanBam" to read BAM files. 
cellRangerDirs <- c(
  "data/HS4M64a1W_fragments.tsv.gz",
  "data/HS4M63a1R_fragments.tsv.gz",
  "data/HS4M57a1R_fragments.tsv.gz", 
  "data/HS4M04a1W_fragments.tsv.gz",
  "data/HS4M07b1W_fragments.tsv.gz",
  "data/HS4M08c1L_fragments.tsv.gz",
  "data/HS4M08e1L_fragments.tsv.gz",
  "data/HS4M08h1L_fragments.tsv.gz",
  "data/HS4M08j1L_fragments.tsv.gz",
  "data/HS4M08k1L_fragments.tsv.gz"
)

sampleNames <- c("HS4M64a1W","HS4M63a1R","HS4M57a1R","HS4M04a1W","HS4M07b1W",
                 "HS4M08c1L","HS4M08e1L","HS4M08h1L","HS4M08j1L","HS4M08k1L")
addArchRChrPrefix(chrPrefix = FALSE)

proj <- readRDS("./proj/Save-ArchR-Project.rds")
proj <- filterDoublets(ArchRProj = proj)
#sample_names <- unique(proj$Sample)
geneNames <- data.frame(getGeneAnnotation(proj)$genes)$symbol

ec_clusters <- list("C8","C12","C11","C10","C2","C5","C4","C14","C1","C2")
names(ec_clusters) <- sampleNames

#-------------------------------------------------------------------------------
cellInfo <- NULL

for(i in 1:length(sampleNames)){

target_sample <- sampleNames[i]
cells_to_keep <- getCellNames(proj)[proj$Sample == target_sample]

if("proj_subset" %in% ls()){
    rm("proj_subset")
}
proj_subset <- proj[cells_to_keep,]

# ---------------Dimensionality Reduction 
proj_subset <- addIterativeLSI(
    ArchRProj = proj_subset,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 2,
    clusterParams = list(
        resolution = c(0.2),
        sampleCells = 10000,
        n.start = 10 
    ),
    varFeatures = 20000,
    dimsToUse = 1:30
)

#----------------------- performed using Seurat::FindClusters() function
proj_subset <- addClusters(
    input=proj_subset,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 1,
    force = TRUE)

tmp <- data.frame(cellnames=getCellNames(proj_subset),clusters=proj_subset$Clusters,sample=target_sample)
cellInfo <- rbind(cellInfo,tmp)

proj_subset <- addUMAP(
    ArchRProj = proj_subset, 
    reducedDims="IterativeLSI",
    name="UMAP", 
    nNeighbors=30, 
    minDist = 0.5, 
    metric="cosine"
    )


pdf(paste0("./figures/alignedUMAP_harmony.",target_sample,".pdf"))
p1 <- plotEmbedding(ArchRProj = proj_subset, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_subset, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1,p2, type="h")
dev.off()

markerGenes <- c("ENSSSCG00000061957","PECAM1","CDH5","CLDN5","TEK","RGCC","LYVE1",
                 "CCL21","SPARCL1","SLC39A10","SLC16A4")
proj_subset <- addImputeWeights(proj_subset)

pdf(paste0("./figures/UMAP_markerGenes.",target_sample,".pdf"))
p <- plotEmbedding(
    ArchRProj=proj_subset,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_subset)
  )
print(p)
dev.off()

#pdf(paste0("./figures/TrackPlot_markerGenes_CD36.",target_sample,".pdf"))
for(gene in markerGenes){
p <- plotBrowserTrack(
    ArchRProj=proj_subset,
    groupBy = "Clusters",
    geneSymbol = gene,
    upstream = 100000,
    downstream = 100000
  )
plotPDF(plotList=p, name = paste0("TrackPlot_markerGenes_",gene,"_",target_sample,".pdf"), ArchRProj=proj_subset, addDOC=FALSE, width=7.5, height=5)
}
save(cellInfo,file="cellInfo.Rdata")
}


EC2keep <- NULL
for(i in sampleNames){
    subdat <- subset(cellInfo,sample==i)
    subdat2 <- subset(subdat,clusters==ec_clusters[[i]])
    ECs <- subdat2[,"cellnames"]
    EC2keep <- c(EC2keep,ECs)
}


proj_EC <- proj[EC2keep,]
saveArchRProject(ArchRProj = proj_EC, outputDirectory = "ArchRProj_EC", load = FALSE)

proj_EC <- readRDS("./ArchRProj_EC/Save-ArchR-Project.rds")

# ----------------------- add group info（case/control）
proj_EC$Group <- getCellColData(proj_EC, "Sample")$Sample
proj_EC$Group <- ifelse(proj_EC$Group %in% c("HS4M64a1W", "HS4M63a1R", "HS4M57a1R", "HS4M04a1W", "HS4M07b1W"),"Peripheral", "CNS")

# -----------------------（pseudo-bulk replicates）
proj_EC <- addGroupCoverages(ArchRProj = proj_EC, groupBy = "Group")

# ------------------------Call Peaks (using MACS2)
pathToMacs2 <- '/opt/conda/envs/r4/bin/macs2'
totalGenomeSize <- sum(width(genomeAnnotation[[3]]))
pdf("./figures/tmp.pdf")
proj_EC <- addReproduciblePeakSet(ArchRProj = proj_EC,groupBy = "Group",pathToMacs2 = pathToMacs2,genomeSize = totalGenomeSize)
dev.off()

system("rm ./figures/tmp.pdf")

# ------------------ Add Peak Matrix
proj_EC <- addPeakMatrix(proj_EC)

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj_EC,
  useMatrix = "PeakMatrix",
  groupBy = "Group",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)")
)


markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
peripheral <- as.data.frame(markerList$Peripheral)

saveArchRProject(ArchRProj = proj_EC, outputDirectory = "ArchRProj_EC", load = FALSE)

proj_EC$Sample <- factor(proj_EC$Sample, levels = sampleNames)
proj_EC$Condition <- ifelse(proj_EC$Sample %in% sampleNames[1:5], "Peripheral", "CNS")



# 添加一个新字段 Sample_ordered，顺序为 factor
prefix <- c(rep("01_",396),rep("02_",85),rep("03_",228),rep("04_",105),rep("05_",464),
            rep("06_",44),rep("07_",134),rep("08_",53),rep("09_",42),rep("10_",96))

proj_EC$Sample_ordered <- paste0(prefix,proj_EC$Sample)


pdf("./figures/tmp.pdf")
p <- plotBrowserTrack(
    ArchRProj=proj_EC,
    groupBy = "Sample_ordered",
    geneSymbol = "ENSSSCG00000061957",
    upstream = 150000,
    downstream = 50000
  )
plotPDF(plotList=p, name = paste0("TrackPlot_markerGenes_CD36_EC.pdf"), 
    ArchRProj=proj_EC, addDOC=FALSE, width=7.5, height=5)




