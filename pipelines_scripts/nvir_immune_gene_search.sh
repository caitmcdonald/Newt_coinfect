###

# Append GO IDs and GO names to DGE results files for immune gene search
# Project: Nvir coinfection

###


# Before filtering: strip header from all files in directory
for f in *.txt; do
    tail -n +2 "$f" > "${f}".tmp && mv "${f}".tmp "$f"
    echo "Processing $f"
done

# Before filtering: cat all the overrepresented files (BP, CC, MF) into a single file
# a: add a column with filename to each overrepresented file
for f in *GO.txt; do sed -i '' "s/$/    $f/" $f; done

# b. cat the BP, CC, and MF files into one, swap column order to match edgeR file
cat spl_bd_bsal_downgenes*.txt | awk '{print $2"\t"$1"\t"$3}' > spl_bd_bsal_down_gene_and_GO.txt

cat spl_bd_bsal_upgenes*.txt | awk '{print $2"\t"$1"\t"$3}' > spl_bd_bsal_up_gene_and_GO.txt

cat spl_bd_ctl_downgenes*.txt | awk '{print $2"\t"$1"\t"$3}' > spl_bd_ctl_down_gene_and_GO.txt

cat spl_bd_ctl_upgenes*.txt | awk '{print $2"\t"$1"\t"$3}' > spl_bd_ctl_up_gene_and_GO.txt

cat spl_both_bd_downgenes*.txt | awk '{print $2"\t"$1"\t"$3}' > spl_both_bd_down_gene_and_GO.txt

cat spl_both_bd_upgenes*.txt | awk '{print $2"\t"$1"\t"$3}' > spl_both_bd_up_gene_and_GO.txt

cat spl_both_bsal_downgenes*.txt | awk '{print $2"\t"$1"\t"$3}' > spl_both_bsal_down_gene_and_GO.txt

cat spl_both_bsal_upgenes*.txt | awk '{print $2"\t"$1"\t"$3}' > spl_both_bsal_up_gene_and_GO.txt

cat spl_both_ctl_downgenes*.txt | awk '{print $2"\t"$1"\t"$3}' > spl_both_ctl_down_gene_and_GO.txt

cat spl_both_ctl_upgenes*.txt | awk '{print $2"\t"$1"\t"$3}' > spl_both_ctl_up_gene_and_GO.txt

cat spl_both_vs_single_downgenes*.txt | awk '{print $2"\t"$1"\t"$3}' > spl_both_vs_single_down_gene_and_GO.txt

cat spl_both_vs_single_upgenes*.txt | awk '{print $2"\t"$1"\t"$3}' > spl_both_vs_single_up_gene_and_GO.txt

cat spl_bsal_ctl_downgenes*.txt | awk '{print $2"\t"$1"\t"$3}' > spl_bsal_ctl_down_gene_and_GO.txt

cat spl_bsal_ctl_upgenes*.txt | awk '{print $2"\t"$1"\t"$3}' > spl_bsal_ctl_up_gene_and_GO.txt

cat spl_inf_ctl_downgenes*.txt | awk '{print $2"\t"$1"\t"$3}' > spl_inf_ctl_down_gene_and_GO.txt

cat spl_inf_ctl_upgenes*.txt | awk '{print $2"\t"$1"\t"$3}' > spl_inf_ctl_up_gene_and_GO.txt


# 1. Filter topGO overrepresented terms by my edgeR significant DGE results (match gene id in latter to gene id in former, filter former)
awk -F\t 'NR==FNR{c[$1]++;next};c[$1]' spl_bd_bsal_down.txt spl_bd_bsal_down_gene_and_GO.txt | awk '{print $2"\t"$1"\t"$3}'> spl_bd_bsal_down_gene_and_GO.filtered.txt

awk -F\t 'NR==FNR{c[$1]++;next};c[$1]' spl_bd_bsal_up.txt spl_bd_bsal_up_gene_and_GO.txt | awk '{print $2"\t"$1"\t"$3}'> spl_bd_bsal_up_gene_and_GO.filtered.txt

awk -F\t 'NR==FNR{c[$1]++;next};c[$1]' spl_bd_ctl_down.txt spl_bd_ctl_down_gene_and_GO.txt | awk '{print $2"\t"$1"\t"$3}'> spl_bd_ctl_down_gene_and_GO.filtered.txt

awk -F\t 'NR==FNR{c[$1]++;next};c[$1]' spl_bd_ctl_up.txt spl_bd_ctl_up_gene_and_GO.txt | awk '{print $2"\t"$1"\t"$3}'> spl_bd_ctl_up_gene_and_GO.filtered.txt

awk -F\t 'NR==FNR{c[$1]++;next};c[$1]' spl_both_bd_down.txt spl_both_bd_down_gene_and_GO.txt | awk '{print $2"\t"$1"\t"$3}'> spl_both_bd_down_gene_and_GO.filtered.txt

awk -F\t 'NR==FNR{c[$1]++;next};c[$1]' spl_both_bd_up.txt spl_both_bd_up_gene_and_GO.txt | awk '{print $2"\t"$1"\t"$3}'> spl_both_bd_up_gene_and_GO.filtered.txt

awk -F\t 'NR==FNR{c[$1]++;next};c[$1]' spl_both_bsal_down.txt spl_both_bsal_down_gene_and_GO.txt | awk '{print $2"\t"$1"\t"$3}'> spl_both_bsal_down_gene_and_GO.filtered.txt

awk -F\t 'NR==FNR{c[$1]++;next};c[$1]' spl_both_bsal_up.txt spl_both_bsal_up_gene_and_GO.txt | awk '{print $2"\t"$1"\t"$3}'> spl_both_bsal_up_gene_and_GO.filtered.txt

awk -F\t 'NR==FNR{c[$1]++;next};c[$1]' spl_both_ctl_down.txt spl_both_ctl_down_gene_and_GO.txt | awk '{print $2"\t"$1"\t"$3}'> spl_both_ctl_down_gene_and_GO.filtered.txt

awk -F\t 'NR==FNR{c[$1]++;next};c[$1]' spl_both_ctl_up.txt spl_both_ctl_up_gene_and_GO.txt | awk '{print $2"\t"$1"\t"$3}'> spl_both_ctl_up_gene_and_GO.filtered.txt

awk -F\t 'NR==FNR{c[$1]++;next};c[$1]' spl_both_vs_single_down.txt spl_both_vs_single_down_gene_and_GO.txt | awk '{print $2"\t"$1"\t"$3}'> spl_both_vs_single_down_gene_and_GO.filtered.txt

awk -F\t 'NR==FNR{c[$1]++;next};c[$1]' spl_both_vs_single_up.txt spl_both_vs_single_up_gene_and_GO.txt | awk '{print $2"\t"$1"\t"$3}'> spl_both_vs_single_up_gene_and_GO.filtered.txt

awk -F\t 'NR==FNR{c[$1]++;next};c[$1]' spl_bsal_ctl_down.txt spl_bsal_ctl_down_gene_and_GO.txt | awk '{print $2"\t"$1"\t"$3}'> spl_bsal_ctl_down_gene_and_GO.filtered.txt

awk -F\t 'NR==FNR{c[$1]++;next};c[$1]' spl_bsal_ctl_up.txt spl_bsal_ctl_up_gene_and_GO.txt | awk '{print $2"\t"$1"\t"$3}'> spl_bsal_ctl_up_gene_and_GO.filtered.txt

awk -F\t 'NR==FNR{c[$1]++;next};c[$1]' spl_inf_ctl_down.txt spl_inf_ctl_down_gene_and_GO.txt | awk '{print $2"\t"$1"\t"$3}'> spl_inf_ctl_down_gene_and_GO.filtered.txt

awk -F\t 'NR==FNR{c[$1]++;next};c[$1]' spl_inf_ctl_up.txt spl_inf_ctl_up_gene_and_GO.txt | awk '{print $2"\t"$1"\t"$3}'> spl_inf_ctl_up_gene_and_GO.filtered.txt


# 2. Add GO term matching GO ID
awk -F\t 'NR==FNR {a[$1]=$2; next} {$(NF+1)=a[$1]} 1' goID_and_names.txt spl_bd_bsal_down_gene_and_GO.filtered.txt > spl_bd_bsal_down_with_GOname.txt #problem: separates by space, not tab, when try to separate by tab, replaces every space in string of GO name column

awk -F\t 'NR==FNR {a[$1]=$2; next} {$(NF+1)=a[$1]} 1' goID_and_names.txt spl_bd_bsal_up_gene_and_GO.filtered.txt > spl_bd_bsal_up_with_GOname.txt

awk -F\t 'NR==FNR {a[$1]=$2; next} {$(NF+1)=a[$1]} 1' goID_and_names.txt spl_bd_ctl_down_gene_and_GO.filtered.txt > spl_bd_ctl_down_with_GOname.txt

awk -F\t 'NR==FNR {a[$1]=$2; next} {$(NF+1)=a[$1]} 1' goID_and_names.txt spl_bd_ctl_up_gene_and_GO.filtered.txt > spl_bd_ctl_up_with_GOname.txt

awk -F\t 'NR==FNR {a[$1]=$2; next} {$(NF+1)=a[$1]} 1' goID_and_names.txt spl_bsal_ctl_down_gene_and_GO.filtered.txt > spl_bsal_ctl_down_with_GOname.txt

awk -F\t 'NR==FNR {a[$1]=$2; next} {$(NF+1)=a[$1]} 1' goID_and_names.txt spl_bsal_ctl_up_gene_and_GO.filtered.txt > spl_bsal_ctl_up_with_GOname.txt

awk -F\t 'NR==FNR {a[$1]=$2; next} {$(NF+1)=a[$1]} 1' goID_and_names.txt spl_both_bd_down_gene_and_GO.filtered.txt > spl_both_bd_down_with_GOname.txt

awk -F\t 'NR==FNR {a[$1]=$2; next} {$(NF+1)=a[$1]} 1' goID_and_names.txt spl_both_bd_up_gene_and_GO.filtered.txt > spl_both_bd_up_with_GOname.txt

awk -F\t 'NR==FNR {a[$1]=$2; next} {$(NF+1)=a[$1]} 1' goID_and_names.txt spl_both_bsal_down_gene_and_GO.filtered.txt > spl_both_bsal_down_with_GOname.txt

awk -F\t 'NR==FNR {a[$1]=$2; next} {$(NF+1)=a[$1]} 1' goID_and_names.txt spl_both_bsal_up_gene_and_GO.filtered.txt > spl_both_bsal_up_with_GOname.txt

awk -F\t 'NR==FNR {a[$1]=$2; next} {$(NF+1)=a[$1]} 1' goID_and_names.txt spl_both_ctl_down_gene_and_GO.filtered.txt > spl_both_ctl_down_with_GOname.txt

awk -F\t 'NR==FNR {a[$1]=$2; next} {$(NF+1)=a[$1]} 1' goID_and_names.txt spl_both_ctl_up_gene_and_GO.filtered.txt > spl_both_ctl_up_with_GOname.txt

awk -F\t 'NR==FNR {a[$1]=$2; next} {$(NF+1)=a[$1]} 1' goID_and_names.txt spl_both_vs_single_down_gene_and_GO.filtered.txt > spl_both_vs_single_down_with_GOname.txt

awk -F\t 'NR==FNR {a[$1]=$2; next} {$(NF+1)=a[$1]} 1' goID_and_names.txt spl_both_vs_single_up_gene_and_GO.filtered.txt > spl_both_vs_single_up_with_GOname.txt

awk -F\t 'NR==FNR {a[$1]=$2; next} {$(NF+1)=a[$1]} 1' goID_and_names.txt spl_inf_ctl_down_gene_and_GO.filtered.txt > spl_inf_ctl_down_with_GOname.txt

awk -F\t 'NR==FNR {a[$1]=$2; next} {$(NF+1)=a[$1]} 1' goID_and_names.txt spl_inf_ctl_up_gene_and_GO.filtered.txt > spl_inf_ctl_up_with_GOname.txt

### change the filename column manually to only BP, CC, MF in text editor

### change all spaces to tabs manually in text editor


### STEP 3: collapse enriched GO term files by gene ID. Steps 1 & 2 generate *with_GOname.txt files, which list all significantly enriched GO terms, the corresponding gene_IDs, GO_categories, and GO terms. However, each gene_ID may be enriched for multiple GO terms, so I need to collapse the files by gene_ID ###
	#done in R: see nvir_immune_gene_search.R
    #output: geneID    GO1, GO2    GOcat1, GOcat2    GOname1, GOname2
    
    
### STEP 4: add GO information back to DGE output. For my DGE data to be searchable, I need to have the gene name (which is output by edgeR), plus all the GO terms. So, need to append all the GO info to the original DGE data, generating new, complete dfs with all information ###
	#done in R: see nvir_immune_gene_search.R
	#output: gene_id, uniprot_id, uniprot_name, DGE results ... GO1, GO2   GOcat1, GOcat2    GOterm1, GOterm2
	
	
### STEP 5: filter compiled DGE results by immune search terms, generate new file for each immune term. ###
	# use grep to search all compiled DGE files for immune terms, add column with search term, and write to new file: grep -i "antigen" *enr.txt |  awk '{print $0, "\tantigen"}' > antigen_enriched.txt
	# see: nvir_all_immune_terms_search.sh
	# cat all *enriched.txt files to generate a single file with all immune-related, significantly DGE genes per each tissue: cat *enriched.txt > spleen_all_immune_enriched.txt
