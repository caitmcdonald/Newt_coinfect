### Add UniprotKB annotation data to annotated Nvir genes ###

setwd("~/nvir/nvir_euro_reference/annotations/")

library(dplyr)

### Read in annotations file created from only swissprot data from reference (see "nvir_add_uniprot_annot.txt")
annotations <- read.delim("newt_swissprot.txt", header=T, stringsAsFactors = F)
  #n_occur <- data.frame(table(annotations$gene)) #gives number of occurences of each $gene
  #how_many <- n_occur[n_occur$Freq >1,] #gives # of $gene with more than 1 occurrence
  #annotations[annotations$gene %in% n_occur$Var1[n_occur$Freq>1],] #lists the annotations with more than 1 $gene occurrence
str(annotations)

# the annotations file has multiple copies of genes with different uniprot descriptions (isoforms hit different uniprot ids). need to collapse annotations file to only unique genes (result is genes with multiple uniprot ids separated by ', ')
collapsed_annot <- aggregate(annotations[,2:4], list(annotations[,1]), function(x) paste0(unique(x)))
str(collapsed_annot)
names(collapsed_annot) <- c("gene","sp_gene","sp_accession","uniprot_id")
head(collapsed_annot)

# convert lists to character vectors
collapsed_annot$gene <- vapply(collapsed_annot$gene, paste, collapse=", ", character(1L))
collapsed_annot$sp_gene <- vapply(collapsed_annot$sp_gene, paste, collapse=", ", character(1L))
collapsed_annot$sp_accession <- vapply(collapsed_annot$sp_accession, paste, collapse=", ", character(1L))
collapsed_annot$uniprot_id <- vapply(collapsed_annot$uniprot_id, paste, collapse=", ", character(1L))
str(collapsed_annot)

### Read in metadata from UniprotKB database search
uni_kb <- read.delim("./uniprot_kb/uniprot_kb_all.txt", header=T, stringsAsFactors = F, sep="\t")
head(uni_kb)
str(uni_kb)

# collapse annotations so only unique uniprot_ids
collapsed_uni_kb <- aggregate(uni_kb[,2:12], list(uni_kb[,1]), function(x) paste0(unique(x)))
names(collapsed_uni_kb) <- c("uniprot_id","Entry","Entry_name","GO_ID","GO_bio_process","GO_cellular","GO_molecular","Gene_names","Protein_names","Annotation","Interacts_with","Tissue_specificity")
str(collapsed_uni_kb)

collapsed_uni_kb$uniprot_id <- vapply(collapsed_uni_kb$uniprot_id, paste, collapse=", ", character(1L))
collapsed_uni_kb$Entry <- vapply(collapsed_uni_kb$Entry, paste, collapse=", ", character(1L))
collapsed_uni_kb$Entry_name <- vapply(collapsed_uni_kb$Entry_name, paste, collapse=", ", character(1L))
collapsed_uni_kb$GO_ID <- vapply(collapsed_uni_kb$GO_ID, paste, collapse=", ", character(1L))
collapsed_uni_kb$GO_bio_process <- vapply(collapsed_uni_kb$GO_bio_process, paste, collapse=", ", character(1L))
collapsed_uni_kb$GO_cellular <- vapply(collapsed_uni_kb$GO_cellular, paste, collapse=", ", character(1L))
collapsed_uni_kb$GO_molecular <- vapply(collapsed_uni_kb$GO_molecular, paste, collapse=", ", character(1L))
collapsed_uni_kb$Gene_names <- vapply(collapsed_uni_kb$Gene_names, paste, collapse=", ", character(1L))
collapsed_uni_kb$Protein_names <- vapply(collapsed_uni_kb$Protein_names, paste, collapse=", ", character(1L))
collapsed_uni_kb$Annotation <- vapply(collapsed_uni_kb$Annotation, paste, collapse=", ", character(1L))
collapsed_uni_kb$Interacts_with <- vapply(collapsed_uni_kb$Interacts_with, paste, collapse=", ", character(1L))
collapsed_uni_kb$Tissue_specificity <- vapply(collapsed_uni_kb$Tissue_specificity, paste, collapse=", ", character(1L))

str(collapsed_uni_kb)


### Merge annotation and uniprot annot dataframes
merged <- left_join(collapsed_annot,collapsed_uni_kb)
str(merged)
names(merged)
head(merged)

# Write table of all UniprotKB annotation data to append to tximport counts matrix for DGE in edgeR
write.table(merged, "~/nvir/nvir_tximport/Newt_gene_swissprot_all_annot.txt", quote=F, sep="\t", row.names=F)

# Write table of GO IDs for use in topGO
merged2 <- merged[,c(1,7)]
head(merged2)
merged2_format <- sapply(merged2, gsub, pattern=";", replacement=",")
merged2_format <- gsub(";",",", merged2$GO_ID)
head(merged2_format)
write.table(merged2_format, "~/nvir/nvir_enrichment/Newt_GOterms_only_final.txt", quote=F, sep="\t", row.names = F)
