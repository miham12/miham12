#!/bin/bash
#SBATCH --job-name=raw_assembly
#SBATCH --error=raw_assembly.err
#SBATCH --output=raw_assembly.log
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=100GB

source ~/.bashrc

bb="/home/miham/gpfs/tools/bbmap/bbduk.sh"
hisat="/home/miham/gpfs/tools/hisat2/hisat2"
stringtie="/home/miham/gpfs/tools/stringtie-2.2.1.Linux_x86_64/stringtie"

annotation="/home/miham/gpfs/tribolium_castaneum/genome/GCA_000002335.3_Tcas5.2_genomic_no_genes.gtf"
index="/home/miham/gpfs/tribolium_castaneum/genome/castaneum_index/castaneum"
log="/home/miham/gpfs/tribolium_castaneum/reads/raw_assembly.txt"
stages=("egg" "larvae" "pupae" "adult")

echo "Скрипт начался" 1> ${log}

for stage in "${stages[@]}"
do
    echo "Жук на стадии ${stage}" 1>> ${log}
    reads="/home/miham/gpfs/tribolium_castaneum/reads/${stage}"
    trim_out="/home/miham/gpfs/tribolium_castaneum/reads/${stage}_trimmed"
    hi_out="/home/miham/gpfs/tribolium_castaneum/mapped/${stage}"
    stringtie_out="/home/miham/gpfs/tribolium_castaneum/assembly/${stage}"
    
    mkdir ${trim_out}
    mkdir ${hi_out}
    mkdir ${stringtie_out}
    for f in "$reads"/*.fastq
    do
        no_path=$(basename "$f") # Getting the filename with extension
        file_name="${no_path%.*}" # Removing the file extension
        
        echo "Триммирую ${f}" 1>> ${log}
        $bb -Xmx1g in="$f" out=$trim_out/${no_path} ftl=15 qtrim=rl trimq=10 minlen=20
        echo "Картирую $trim_out/${no_path}" 1>> ${log}
        ${hisat} -x ${index} -U $trim_out/${no_path} -S ${hi_out}/${file_name}.sam 2>> ${log}
        echo "Конвертирую сам в бам и сортирую" 1>> ${log} 
        samtools sort ${hi_out}/${file_name}.sam -o ${hi_out}/${file_name}.bam 2>> ${log}
        echo "Индексирую" 1>> ${log} 
        samtools index ${hi_out}/${file_name}.bam 2>> ${log}
        echo "Сохраняем инфу о выравнивании" 1>> ${log}
        samtools flagstat ${hi_out}/${file_name}.bam 1>> ${hi_out}/${file_name}_info.txt  2>>${log}
        echo "Сборка транскриптома" 1>> ${log}
        ${stringtie} -o ${stringtie_out}/${file_name}.gtf -G ${annotation} ${hi_out}/${file_name}.bam 2>>${log}
    done
done
echo "Скрипт окончен" 1>> ${log}