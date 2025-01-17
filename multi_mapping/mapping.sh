#!/bin/bash
#SBATCH --job-name=mapping
#SBATCH --error=mapping_err
#SBATCH --output=mapping_log
#SBATCH --time=50:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=150GB
#SBATCH --partition=all
source ~/anaconda3/bin/activate
conda activate samtools_env
start_path="/home/miham/vavilon_save"
index="${start_path}/t2t_2_0_ncbi/t2t_index/t2t"

dna_path="${start_path}/reads/SRR3633289.fastq"
file_name=$(basename "$dna_path") 
dna="${file_name%.*}" # Removing the file extension

rna_path="${start_path}/reads/SRR3633288.fastq"
file_name=$(basename "$rna_path") 
rna="${file_name%.*}" # Removing the file extension
ss="${start_path}/t2t_2_0_ncbi/splice_site"
#для пересечений 
repeats="${start_path}/repeat_t2t_v2/individual_repeats/repeatmasker_chr_ids.bed"
genes="${start_path}/t2t_2_0_ncbi/GCF_009914755.1_T2T-CHM13v2.0_genomic_short.bed"
regulatory="${start_path}/t2t_2_0_ncbi/GCF_009914755.1_T2T-CHM13v2.0_genomic_regulatory.bed"

hisat="/home/miham/gpfs/tools/hisat2/hisat2"
hi_out="${start_path}/mm1s_t2t"
log="${start_path}/mm1s_t2t/mapping_log.txt"


echo "Script starts" 1> ${log} 

echo "Mapping ${dna_path} to the genome" 1>> ${log}
${hisat} -x ${index} -U ${dna_path} -S ${hi_out}/${dna}.sam --no-softclip -k 20 --no-spliced-alignment 2>> ${log}
echo "Converting sam to bam and sorting" 1>> ${log} 
samtools sort ${hi_out}/${dna}.sam -o ${hi_out}/${dna}.bam 2>> ${log}

echo "Indexing ${dna}" 1>> ${log} 
samtools index ${hi_out}/${dna}.bam 2>> ${log}
echo "Mapping info" 1>> ${log}
samtools flagstat ${dna}.bam 1>> ${dna}_info.txt  2>>${log}

echo "Obtaining only mapped reads" 1>> ${log} 
samtools view -@ 3 -h -b -F 4 ${hi_out}/${dna}.bam > ${hi_out}/${dna}_mapped.bam 2>>${log}
echo "Obtaining only uniquely mapped reads" 1>> ${log} 
samtools view -@ 3 -h ${hi_out}/${dna}_mapped.bam | grep -P "^@|NH:i:1\b" | samtools view -Sb - > ${hi_out}/${dna}_unique.bam 2>>${log}
echo "Obtaining only multiply mapped reads" 1>> ${log} 
samtools view -@ 3 -h ${hi_out}/${dna}_mapped.bam | grep -w -v "NH:i:1" | samtools view -Sb - > ${hi_out}/${dna}_multi.bam 2>>${log}

echo "Converting Bam to Bed" 1>> ${log} 
bedtools bamtobed -i ${hi_out}/${dna}_unique.bam 1> ${hi_out}/coressponding_beds/${dna}_unique.bed 2>>${log}
bedtools bamtobed -i ${hi_out}/${dna}_multi.bam 1> ${hi_out}/coressponding_beds/${dna}_multi.bed 2>>${log}

#awk " $3 = $2+1 1" SRR3633288_multi.bed  > trash.txt
#cat trash.txt > SRR3633288_multi.bed


echo "Intersecting with repeats" 1>> ${log} 
bedtools intersect -a ${hi_out}/coressponding_beds/${dna}_unique.bed -b ${repeats} -wa -wb -loj 1> trash.txt 2>>${log}
awk -F '\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $8 "\t" $9 "\t" $10}' < trash.txt 1> ${hi_out}/repeats_intersected/${dna}_unique_repeats.bed 2>>${log}

bedtools intersect -a ${hi_out}/coressponding_beds/${dna}_multi.bed -b ${repeats} -wa -wb -loj 1> trash.txt 2>>${log}
awk -F '\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $8 "\t" $9 "\t" $10}' < trash.txt 1> ${hi_out}/repeats_intersected/${dna}_multi_repeats.bed 2>>${log}


echo "Intersecting with genes" 1>> ${log}
bedtools intersect -a ${hi_out}/coressponding_beds/${dna}_unique.bed -b ${genes} -wa -wb -loj 1> trash.txt 2>>${log}
awk -F '\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $8 "\t" $9 "\t" $10 "\t" $11}' < trash.txt 1>${hi_out}/genes_intersected/${dna}_unique_genes.bed 2>>${log}

bedtools intersect -a ${hi_out}/coressponding_beds/${dna}_multi.bed -b ${genes} -wa -wb -loj 1> trash.txt 2>>${log}
awk -F '\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $8 "\t" $9 "\t" $10 "\t" $11}' < trash.txt 1> ${hi_out}/genes_intersected/${dna}_multi_genes.bed 2>>${log}


echo "Intersecting with regulatory " 1>> ${log} 
bedtools intersect -a ${hi_out}/coressponding_beds/${dna}_unique.bed -b ${regulatory} -wa -wb -loj 1> trash.txt 2>>${log}
awk -F '\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $8 "\t" $9 "\t" $11}' < trash.txt 1>${hi_out}/genes_intersected/${dna}_unique_regulatory.bed 2>>${log}

bedtools intersect -a ${hi_out}/coressponding_beds/${dna}_multi.bed -b ${regulatory} -wa -wb -loj 1> trash.txt 2>>${log}
awk -F '\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $8 "\t" $9 "\t" $11}' < trash.txt 1>${hi_out}/genes_intersected/${dna}_multi_regulatory.bed 2>>${log}


echo "Mapping ${rna_path} to the genome" 1>> ${log}
${hisat} -x ${index} -U ${rna_path} -S ${hi_out}/${rna}.sam --no-softclip -k 20 --no-spliced-alignment --known-splicesite-infile ${ss} 2>> ${log}
echo "Converting sam to bam and sorting" 1>> ${log} 
samtools sort -@ 3 ${hi_out}/${rna}.sam -o ${hi_out}/${rna}.bam -T ${hi_out}/temps 2>> ${log}

echo "Indexing ${rna}" 1>> ${log} 
samtools index ${hi_out}/${rna}.bam 2>> ${log}
echo "Mapping info" 1>> ${log}
samtools flagstat ${rna}.bam 1>> ${rna}_info.txt  2>>${log}

echo "Obtaining only mapped reads" 1>> ${log} 
samtools view -@ 5 -h -b -F 4 ${hi_out}/${rna}.bam > ${hi_out}/${rna}_mapped_other_strand.bam 2>>${log}

echo "You need to change strand for rna data" 1>> ${log} 
#samtools view -h ${hi_out}/${rna}_mapped_other_strand.bam | awk -F '\t' 'BEGIN {OFS = FS} {if ($1 !~ /^@/){$2 = ($2 == 0) ? 16 : 0} print}' - | samtools view -Sb - 1> ${hi_out}/${rna}_mapped.bam  2>> ${log} 
conda deactivate
conda activate pysam_env
python3 change_strand.py ${hi_out}/${rna}_mapped_other_strand.bam ${hi_out}/${rna}_mapped.bam
#rm ${hi_out}/${rna}_mapped_other_strand.bam 2>> ${log} 
conda activate samtools_env

echo "Obtaining only uniquely mapped reads" 1>> ${log} 
samtools view -@ 5 -h ${hi_out}/${rna}_mapped.bam | grep -P "^@|NH:i:1\b" | samtools view -Sb - > ${hi_out}/${rna}_unique.bam 2>>${log}

echo "Obtaining only multiply mapped reads" 1>> ${log} 
samtools view -@ 5 -h ${hi_out}/${rna}_mapped.bam | grep -w -v "NH:i:1" | samtools view -Sb - > ${hi_out}/${rna}_multi.bam 2>>${log}

echo "Converting Bam to Bed" 1>> ${log} 
bedtools bamtobed -i ${hi_out}/${rna}_unique.bam 1> ${hi_out}/coressponding_beds/${rna}_unique.bed 2>>${log}
bedtools bamtobed -i ${hi_out}/${rna}_multi.bam 1> ${hi_out}/coressponding_beds/${rna}_multi.bed 2>>${log}


echo "Intersecting with genes" 1>> ${log}
bedtools intersect -a ${hi_out}/coressponding_beds/${rna}_unique.bed -b ${genes} -wa -wb -loj 1> trash.txt 2>>${log}
awk -F '\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $8 "\t" $9 "\t" $10 "\t" $11}' < trash.txt 1>${hi_out}/genes_intersected/${rna}_unique_genes.bed 2>>${log}

bedtools intersect -a ${hi_out}/coressponding_beds/${rna}_multi.bed -b ${genes} -wa -wb -loj 1> trash.txt 2>>${log}
awk -F '\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $8 "\t" $9 "\t" $10 "\t" $11}' < trash.txt 1> ${hi_out}/genes_intersected/${rna}_multi_genes.bed 2>>${log}

echo "Script finished" 1>> ${log} 

