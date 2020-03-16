#######################################################################################################################

###### GO term enrichment using topGO ######
##### useful link: http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html #####
##### looped based on: https://github.com/PBR/scripts/blob/master/topGO.R #####

### README: topGO takes GO hierarchy into account when calculating enrichment, leading to fewer false positives. ###

#######################################################################################################################


### set working directory
setwd("~/Desktop/nvir_final/enrichment/") #adjust to subdirectory with DGE tables of interest (e.g. skin)

### install the latest Bioconductor core packages and topGO
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  + install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("topGO")

### load packages
library("topGO")
library("BiocGenerics")
library("AnnotationDbi")
library("S4Vectors")
library("data.table")

### create list of files to use in loop
list_files=dir(pattern='.*txt')
#list_files[1]

### run topGO and generate output files: 1) significantly enriched terms for each GO category, 2) list of genes related to those sig enriched terms
for (i in 1:length(list_files))
{
  ## Required #1: gene-to-GO mapping (tab-delimited with first column gene name, second column list of GO terms separated by comma and space) ##
  geneIDs2GO <- readMappings(file="../Newt_GOterms_only_final.txt") #remember to delete header if present
  L=lapply(geneIDs2GO, function(x) {gsub('\"',"",x)}) #grep to remove '"'
  geneID2GO=L
  names(geneID2GO) <- gsub('\"',"",names(geneID2GO)) #grep to remove '"'
  
  ## Required #2: list of all genes in genome/transcriptome with GO annotation. Use names from the gene-to-GO mapping file. ##
  geneUniverse <- names(geneID2GO)
  
  ## Required #3: list of genes of interest, optionally with a gene-wise score (e.g. DGE significance p-value from edgeR) ##
  tableOfInterest=fread(list_files[[i]]) #read in all files as list of files
  genesOfInterest=tableOfInterest$gene #only interested in gene column
  genesOfInterest
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest)) #create object telling topGO which genes in the gene universe are your genes of interest
  names(geneList) <- geneUniverse
  geneList
  
  ### Create topGOdata objects for each GO category
  myGOdata_BP <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=5)
  print(myGOdata_BP)
  myGOdata_CC <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=5)
  print(myGOdata_CC)
  myGOdata_MF <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=5)
  print(myGOdata_MF)
  
  ### Perform analysis
  # Biological process
  resultFisher <- runTest(myGOdata_BP, algorithm="weight01", statistic="fisher") #use weighted Fisher
  allRes <- GenTable(myGOdata_BP, weight01Fisher = resultFisher, topNodes = length(resultFisher@score)) 
  #print(allRes)
  allRes_2=allRes[as.numeric(allRes$weight01Fisher)<0.05,] #subset only significant enrichment
  print(allRes_2)
  write.table(allRes_2,file=paste(sep="",strsplit(list_files[[i]], "\\.")[[1]][1],"_BP_enrichment.txt"), quote=F, sep="\t")
  myterms=allRes_2[,1]
  mygenes=genesInTerm(myGOdata_BP, myterms)
  print(allRes_2)
  print (myterms)
  print(mygenes)
  x=data.frame(GOs = rep(names(mygenes), lapply(mygenes, length)), Genes = unlist(mygenes), stringsAsFactors=FALSE) #generate table of genes linked to sig enriched terms
  write.table(x, paste(sep="",strsplit(list_files[[i]], "\\.")[[1]][1],"genes_corresponding_overrepresent_BP_GO.txt"), append=T, row.names=F, quote=F, sep="\t")
  rm(resultFisher)

  #Cellular component
  resultFisher <- runTest(myGOdata_CC, algorithm="weight01", statistic="fisher") #use weighted Fisher
  allRes <- GenTable(myGOdata_CC, weight01Fisher = resultFisher, topNodes = length(resultFisher@score))
  #print(allRes)
  allRes_2=allRes[as.numeric(allRes$weight01Fisher)<0.05,] #subset only significant enrichment
  #print(allRes_2)
  write.table(allRes_2,file=paste(sep="",strsplit(list_files[[i]], "\\.")[[1]][1],"_CC_enrichment.txt"), quote=F, sep="\t")
  myterms=allRes_2[,1]
  mygenes=genesInTerm(myGOdata_CC, myterms)
  print(allRes_2)
  print (myterms)
  print(mygenes)
  x=data.frame(GOs = rep(names(mygenes), lapply(mygenes, length)), Genes = unlist(mygenes), stringsAsFactors=FALSE) 
  write.table(x, paste(sep="",strsplit(list_files[[i]], "\\.")[[1]][1],"genes_corresponding_overrepresent_CC_GO.txt"), append=T, row.names=F, quote=F, sep="\t")
  rm(resultFisher)
  
  #Molecular function
  resultFisher <- runTest(myGOdata_MF, algorithm="weight01", statistic="fisher") #use weighted Fisher
  allRes <- GenTable(myGOdata_MF, weight01Fisher = resultFisher, topNodes = length(resultFisher@score))
  #print(allRes)
  allRes_2=allRes[as.numeric(allRes$weight01Fisher)<0.05,] #subset only significant enrichment
  #print(allRes_2)
  write.table(allRes_2,file=paste(sep="",strsplit(list_files[[i]], "\\.")[[1]][1],"_MF_enrichment.txt"), quote=F, sep="\t")
  myterms=allRes_2[,1]
  mygenes=genesInTerm(myGOdata_MF, myterms)
  print(allRes_2)
  print (myterms)
  print(mygenes)
  x=data.frame(GOs = rep(names(mygenes), lapply(mygenes, length)), Genes = unlist(mygenes), stringsAsFactors=FALSE) #generate table of genes linked to sig enriched terms
  write.table(x, paste(sep="",strsplit(list_files[[i]], "\\.")[[1]][1],"genes_corresponding_overrepresent_MF_GO.txt"), append=T, row.names=F, quote=F, sep="\t")
  rm(resultFisher)
  
}
  