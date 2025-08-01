library(ArchR)

sampleid = c("sampleid1","sampleid2")

for (i in 1:length(sampleid)){
proj=loadArchRProject(path = "./Save-proj1_filterDoublets_2.5_VAR15000", force = FALSE, showLogo = TRUE)
proj=addImputeWeights(proj)

markersPeaks <- getMarkerFeatures(
  ArchRProj = tmp, 
  useMatrix = "PeakMatrix", 
  groupBy = "celltype1",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerGenes <- c("IRX6")
## 绘制在一个pdf里
for (j in 1:length(markerGenes)){
p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "celltype1", 
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  geneSymbol = markerGenes[j],
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", returnGR = TRUE),
  upstream = 9000,
  downstream = 2000,
  facetbaseSize = 15,
  baseSize =15,
)
plotPDF(p, name = paste0(sampleid[i],'_',markerGenes[j], "_BrowserTrack_Plot.pdf"), width = 7, height = 10, ArchRProj = Merged.proj2, addDOC = FALSE)
}