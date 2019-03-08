## ---- echo=TRUE, eval=TRUE, results=FALSE, message=FALSE-----------------
library(Seurat)
library(MTGOsc)

## ---- echo=TRUE, eval=TRUE-----------------------------------------------
#this is a (simplified) Seurat object
print(bladder)

#this come from applying Seurat function FindAllMarkers() to the bladder dataset
head(markers)

#this is a simple gene -> pathway dictionary
head(mouse.pathways)

## ---- echo=TRUE, eval=TRUE-----------------------------------------------
table(bladder@ident)

## ---- echo=TRUE, eval=TRUE-----------------------------------------------
cells.selected = bladder@ident == 'Basal epithelial cell(Bladder)'
genes.selected = subset(markers, cluster == 'Basal epithelial cell(Bladder)')$gene

my.bladder = bladder
my.bladder@data = my.bladder@data[genes.selected, cells.selected]

## ---- echo=TRUE, eval=TRUE-----------------------------------------------
root = tempdir() #change this to your preferred local path
dir.create(root, recursive = TRUE, showWarnings = FALSE)

## ---- echo=TRUE, eval=TRUE-----------------------------------------------
#building a genes-pathways dictionary
dict = write.dictionary(genes=mouse.pathways$gene, terms = mouse.pathways$pathway, outfolder = root)

#computign gene coexpression (default function is 'cor')
coexp = write.coexpressionMatrix(geneExpression = my.bladder, outfolder = root)

#thinning coexpression network via scale criterion
edges = write.edges(coexpression = coexp, outfolder = root, keep.weights = FALSE, fun = scale_free_threshold)

#writing a parameter file, useful for MGTO
write.paramFile(outfolder = root)

#actual call to MTGO
call.MTGO(outfolder = root, verbose = TRUE)

#building and saving representation of resulting network
network.collapsed = export.network.modules(infolder = root, collapse.modules = TRUE) 
network.full = export.network.modules(infolder = root, collapse.modules = FALSE) 

## ----eval=TRUE, include=TRUE, message=FALSE------------------------------
# load libraries for gene enrichment on Reactome (those are on Bioconductor, not CRAN)
library(ReactomePA)      
library(clusterProfiler)
library(org.Mm.eg.db)

#a support function to take care of gene upper/lower case convention
firstup = function(x) {
  x = tolower(x)
  substr(x, 1, 1) = toupper(substr(x, 1, 1))
  return(x)
}

#the list of all genes involved in the cluster
genes = unique(c(as.character(edges$gene1), as.character(edges$gene2)))

#correct casing of gene names
genes = firstup(genes)

#translating gene names to ENTREZID via org.Mm.eg.db database
genes = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

#the actual enrichment
enriched = enrichPathway(gene=genes$ENTREZID, pvalueCutoff=0.05, readable=T, organism = "mouse", pAdjustMethod = "BH")

## ----eval=FALSE, include=FALSE-------------------------------------------
#  # In this example we search the literature for the Basal Epithelial cell type and the pathway terms
#  # retrieved by both MTGO-SC (pathway labelling gene modules) and ReactomePA (pathways enriched for the whole basal epithelial gene network
#  
#  library(RISmed)
#  
#  # set the search parameters
#  file_MTGO<-"Modules_Best_QGO.txt" #path
#  file_enriched<-"Basal epithelial cell" #passare la variabile
#  tissue<-"Bladder"
#  celltype<-"Basal Epithelial Cell"
#  
#  modules<-read.table(file_MTGO,sep="\t",header=T)
#  terms_MTGO<-as.character(modules[,4])
#  terms_enriched<-setdiff(as.character(read.table(file_enriched,header=T)[,1]),NA)
#  
#  terms<-union(terms_MTGO,terms_enriched)
#  
#  
#  # Perform all the PUBMED searchings
#  
#  pubmed_search_both<-c()
#  pubmed_search_term<-c()
#  for(k in terms)
#  {
#    print(k)
#    res <- EUtilsSummary(paste(tissue,celltype,k), type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2019, retmax=500)
#    pubmed_search_both[[k]]<-QueryCount(res)
#    res <- EUtilsSummary(paste(tissue,k), type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2019, retmax=500)
#    pubmed_search_term[[k]]<-QueryCount(res)
#    Sys.sleep(0.5)  # a connection error might occur without this pause
#  }
#  
#  pubmed_search_universe<-QueryCount(EUtilsSummary(tissue, type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2019, retmax=500))
#  pubmed_search_celltype<-QueryCount(EUtilsSummary(paste(tissue,celltype), type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2019, retmax=500))
#  
#  
#  # Compute p.values
#  
#  pvals<-c()
#  
#  M<-pubmed_search_celltype
#  N<-pubmed_search_universe
#  
#  for(i in terms)
#  {
#      x<-pubmed_search_both[[i]]
#      k<-pubmed_search_term[[i]]
#      pvals[[i]]<-signif(phyper(x-1,M,N-M,k,lower.tail=F),2)
#  }
#  
#  PUBMEDsearch_results<-cbind(terms,pvals,pubmed_search_both,pubmed_search_term,pubmed_search_celltype,pubmed_search_universe)

