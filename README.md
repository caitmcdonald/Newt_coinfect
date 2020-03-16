# Chytridiomycosis (*Bsal-Bd*) coinfection in Eastern newts (*Notophthalmus viridescens*)

## This repository contains files related to RNA-sequencing and differential gene expression analyses for McDonald et al. 2020

#### Overview:
*Bsal* is an invasive fungal disease of salamanders that is currently spreading to naïve populations. North American salamanders are naïve to *Bsal*, but may have high *Bd* prevalence. To investigate host immune responses under coinfected vs. single infection scenarios, we performed an infection trial and harvested immune-related tissues for differential gene expression analyses.

#### Infection trial:
Eastern newts were divided into four groups:
1. __*Bd-infected*:__  natural *Bd* infection
1. __*Bsal-infected*:__  *Bd*-negative + experimental *Bsal* inoculation
1. __Coinfected:__  natural *Bd* infection + experimental *Bsal* inoculation
1. __Controls:__  *Bd*-negative + *Bsal*-negative

Skin, spleen, and liver tissue was harvested upon trial completion for RNA-seq.

#### Sequencing:
- Illumina NextSeq 500, SE 75-bp
- SRA:
- Sample metadata: [nvir_seq_metadata.txt](data/nvir_seq_metadata.txt)

#### Analysis:
__1. Read trimming (Trimmomatic):__
- Standard trimming according to suggestions in [MacManes et al. 2014](https://www.frontiersin.org/articles/10.3389/fgene.2014.00013/full)
- Script: [nvir_trimmomatic.sh](nvir_trimmomatic.sh)

__2. Transcriptome and annotation data:__
- Note: transcriptome provided by authors of Looso et al. 2013 via associated website (http://newtomics.mpi-bn.mpg.de/index.php) is NOT most current
- I had to email the authors to obtain the most recent transcriptome and annotation data. I was sent the following files between Jul 6 and Aug 1, 2018: "newt_transcriptome_trinity_clean.fa", "newt_annotation.2.0.tsv"
- From here, I used the author-provided Swissprot annotations to mine the UniprotKB database for all associated annotation information (e.g. GO terms, gene names, protein names)
- Procedure/scripts: [nvir_add_uniprot_annot.txt](nvir_add_uniprot_annot.txt), [nvir_add_annotation.R](nvir_add_annotation.R)
- Output files: [Newt_GOterms_only_final.txt](Newt_GOterms_only_final.txt), [Newt_gene_swissprot_all_annot.txt](Newt_gene_swissprot_all_annot.txt)

__3. Read quantification (Salmon):__
- Used newt transcriptome sent from Looso et al.: "newt_transcriptome_trinity_clean.fa"
- Estimated SE fragment length and SD based on Fragment Analyzer results
- Script: [nvir_salmon.sh](nvir_salmon.sh)

__4. Differential expression (tximport and edgeR):__
- Followed recommendations in tximport and edgeR documentation
- Script: [nvir_tximport_edger.R](nvir_tximport_edger.R)
- Compiled output file: [FileS1_all_edger_results.txt](FileS1_all_edger_results.txt)

__5. Functional enrichment (topGO):__
- Followed recommendations in topGO documentation
- Script: [nvir_topgo_looped.R](nvir_topgo_looped.R)
- Compiled output file: [FileS2_all_enrichment_results.txt](FileS2_all_enrichment_results.txt)

__6. Immune gene search:__
- Searched significantly differentially expressed genes for those pertinent to immune response
- Step 1: add GO term names back to DGE results
    - Scripts: [nvir_immune_gene_search.sh](nvir_immune_gene_search.sh), [nvir_immune_gene_search.R](nvir_immune_gene_search.R)
- Step 2: search all DGE results files (with GO names added) for custom immune search terms
    - Script: [nvir_all_immune_terms_search.sh](nvir_all_immune_terms_search.sh)
    - Compiled output file: [FileS3_immune_genes.txt](FileS3_immune_genes.txt)
