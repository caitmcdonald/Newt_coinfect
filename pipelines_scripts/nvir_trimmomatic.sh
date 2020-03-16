####Trimmomatic parameters for Nvir Bsal-Bd####
##based on optimal trimming parameters described in MacManes 2014: https://doi.org/10.3389/fgene.2014.00013

mkdir trimmed
for file in *.fastq.gz
do
    java -jar /programs/trimmomatic/trimmomatic-0.36.jar SE -phred33 -threads 20 "$file" "trimmed/${file%.fastq.gz}.trimmed.fastq.gz" ILLUMINACLIP:/programs/trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:25
done
