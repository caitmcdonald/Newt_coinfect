###

# Generate dataframes with all edgeR output + all topGO output
# Project: Nvir coinfection
#README: generates tables for each pairwise edgeR comparison that have topGO output appended (GO accessions, GO category, GO term). Use these to search for immune keywords as in: nvir_all_immune_terms_search.sh ####

###


# get packages
library(dplyr)
library(tidyr)
library(data.table)

# set working directory (change for each tissue)
setwd("~/nvir_final/enrichment/spleen/immune_search/")

### STEPS 1 and 2: done in nvir_immune_gene_search.sh


### STEP 3: collapse enriched GO term files by gene ID. Steps 1 & 2 generate *with_GOname.txt files, which list all significantly enriched GO terms, the corresponding gene_IDs, GO_categories, and GO terms. However, each gene_ID may be enriched for multiple GO terms, so I need to collapse the files by gene_ID ###

# load overrepresented GO terms as list
enrichment_files <- dir(pattern='*GOname.txt')
enrichment <- lapply(enrichment_files, read.delim, header=F)
names(enrichment) <-  c(tools::file_path_sans_ext(basename(enrichment_files)))

# collapse overrepresented GO term files by gene_ID
enriched_collapse <- list()
for (i in 1:length(enrichment))
{
  enriched_collapse[[i]] <- aggregate(enrichment[[i]][,c(1,3:4)], list(enrichment[[i]][,2]), function(x) paste0(x))
  names(enriched_collapse[[i]]) <- c("gene_ID","GO_accession","GO_category","GO_term")
  enriched_collapse[[i]]$gene_ID <- vapply(enriched_collapse[[i]]$gene_ID, paste, collapse=", ", character(1L))
  enriched_collapse[[i]]$GO_accession <- vapply(enriched_collapse[[i]]$GO_accession, paste, collapse=", ", character(1L))
  enriched_collapse[[i]]$GO_category <- vapply(enriched_collapse[[i]]$GO_category, paste, collapse=", ", character(1L))
  enriched_collapse[[i]]$GO_term <- vapply(enriched_collapse[[i]]$GO_term, paste, collapse=", ", character(1L))
}
names(enriched_collapse) <- c(tools::file_path_sans_ext(basename(enrichment_files)))

### STEP 4: add GO information back to DGE output. For my DGE data to be searchable, I need to have the gene name (which is output by edgeR), plus all the GO terms. So, need to append all the GO info to the original DGE dataframes, generating new, complete dfs with all information ###

# load down files as list 
dge_down_files <- dir(pattern='*down.txt')
dge_down <- lapply(dge_down_files, read.delim, header=F)
names(dge_down) <-  c(tools::file_path_sans_ext(basename(dge_down_files)))

# convert lists to character vectors
dge_down_names <- list()
for (i in 1:length(dge_down))
{
  dge_down_names[[i]] <- dge_down[[i]]
  names(dge_down_names[[i]]) <- c("gene_ID","uniprot_ID","uniprot_name","uniprot_GO","Eff_length","logFC","logCPM","LR","PValue","FDR")
  dge_down_names[[i]]$gene_ID <- vapply(dge_down_names[[i]]$gene_ID, paste, collapse=", ", character(1L))
  dge_down_names[[i]]$uniprot_ID <- vapply(dge_down_names[[i]]$uniprot_ID, paste, collapse=", ", character(1L))
  dge_down_names[[i]]$uniprot_name <- vapply(dge_down_names[[i]]$uniprot_name, paste, collapse=", ", character(1L))
  dge_down_names[[i]]$uniprot_GO <- vapply(dge_down_names[[i]]$uniprot_GO, paste, collapse=", ", character(1L))
}
names(dge_down_names) <- c(tools::file_path_sans_ext(basename(dge_down_files)))

# load up files
dge_up_files <- dir(pattern='*up.txt')
dge_up <- lapply(dge_up_files, read.delim, header=F)
names(dge_up) <-  c(tools::file_path_sans_ext(basename(dge_up_files)))

# convert lists to character vectors
dge_up_names <- list()
for (i in 1:length(dge_up))
{
  dge_up_names[[i]] <- dge_up[[i]]
  names(dge_up_names[[i]]) <- c("gene_ID","uniprot_ID","uniprot_name","uniprot_GO","Eff_length","logFC","logCPM","LR","PValue","FDR")
  dge_up_names[[i]]$gene_ID <- vapply(dge_up_names[[i]]$gene_ID, paste, collapse=", ", character(1L))
  dge_up_names[[i]]$uniprot_ID <- vapply(dge_up_names[[i]]$uniprot_ID, paste, collapse=", ", character(1L))
  dge_up_names[[i]]$uniprot_name <- vapply(dge_up_names[[i]]$uniprot_name, paste, collapse=", ", character(1L))
  dge_up_names[[i]]$uniprot_GO <- vapply(dge_up_names[[i]]$uniprot_GO, paste, collapse=", ", character(1L))
}
names(dge_up_names) <- c(tools::file_path_sans_ext(basename(dge_up_files)))

# merge down and up lists back and order by name so that it matches enriched_collapse list
dge_merged <- c(dge_down_names, dge_up_names)
dge_merged_order <- dge_merged[order(names(dge_merged))]

names(enriched_collapse)
names(dge_merged_order)

# join GO information from enriched_collapse list to dge_merged_order list
dge_with_enr <- list()
for (i in 1:length(dge_merged_order))
{
  dge_with_enr[[i]] <- left_join(dge_merged_order[[i]], enriched_collapse[[i]], by="gene_ID", all=TRUE)
}
names(dge_with_enr) <- names(dge_merged_order)

# write output files: all DGE information, plus all information about GO terms to which the particular gene_IDs were related (according to topGO)

# spleen
lapply(1:length(dge_with_enr),
       function(i) write.table(dge_with_enr[[i]],
                               file=paste0("spl_",names(dge_with_enr[i]),"_enr.txt"),
                               row.names=F, quote=F, sep="\t"))

