###

# Salmon alignment-free transcript quantification #
# Following documentation: https://salmon.readthedocs.io/en/latest/salmon.html#using-salmon 
# Project: Nvir coinfection

###


## Index transcriptome
/programs/salmon-0.11.3/bin/salmon index -t newt_transcriptome_trinity_clean.fa -i transcripts_index --type quasi -k 31

## Quantification 
/programs/salmon-0.11.3/bin/salmon quant -i transcripts_index --libType U -r bdlivrep1.fastq --fldMean 320 --fldSD 40 -o bd_liver_rep1

# fldMean: for single-end reads, salmon cannot calculate the fragment length distribution, so it must be set by the user (based on Bioanalyzer readout)
# fldSD: same as above


## Or, can also run using Trinity util align_and_estimate_abundance.pl instead

## prep reference and perform alignment-free quantification
#nohup /programs/trinityrnaseq-Trinity-v2.8.4/util/align_and_estimate_abundance.pl --transcripts newt_transcriptome_trinity_clean.fa --seqType fq --samples_file trinity_samples_file.txt  --est_method salmon --output_dir Nvir_salmon_all --fragment_length 320 --fragment_std 40 --thread 23 --prep_reference >& quant.logfile &

### compile library counts to matrix ###
#nohup /programs/trinityrnaseq-Trinity-v2.8.4/util/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map new_mapping_file.txt --name_sample_by_basedir bd_liver_rep1/quant.sf bd_liver_rep2/quant.sf bd_liver_rep3/quant.sf bd_liver_rep4/quant.sf bd_skin_rep1/quant.sf bd_skin_rep2/quant.sf bd_skin_rep3/quant.sf bd_skin_rep4/quant.sf bd_spleen_rep1/quant.sf bd_spleen_rep2/quant.sf bd_spleen_rep3/quant.sf bd_spleen_rep4/quant.sf both_liver_rep1/quant.sf both_liver_rep2/quant.sf both_liver_rep3/quant.sf both_liver_rep4/quant.sf both_skin_rep1/quant.sf both_skin_rep2/quant.sf both_skin_rep3/quant.sf both_skin_rep4/quant.sf both_spleen_rep1/quant.sf both_spleen_rep2/quant.sf both_spleen_rep3/quant.sf both_spleen_rep4/quant.sf bsal_liver_rep1/quant.sf bsal_liver_rep2/quant.sf bsal_liver_rep3/quant.sf bsal_liver_rep4/quant.sf bsal_skin_rep1/quant.sf bsal_skin_rep2/quant.sf bsal_skin_rep3/quant.sf bsal_skin_rep4/quant.sf bsal_spleen_rep1/quant.sf bsal_spleen_rep2/quant.sf bsal_spleen_rep3/quant.sf bsal_spleen_rep4/quant.sf control_liver_rep1/quant.sf control_liver_rep2/quant.sf control_liver_rep3/quant.sf control_liver_rep4/quant.sf control_skin_rep1/quant.sf control_skin_rep2/quant.sf control_skin_rep3/quant.sf control_skin_rep4/quant.sf control_spleen_rep1/quant.sf control_spleen_rep2/quant.sf control_spleen_rep3/quant.sf control_spleen_rep4/quant.sf >& matrix.logfile &