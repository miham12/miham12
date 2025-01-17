#!/bin/bash
#SBATCH --job-name=download_sra
#SBATCH --error=download_err
#SBATCH --output=download_log
#SBATCH --time=900
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=40GB
index_list=("SRR2812572" "SRR2812573" "SRR2812574" "SRR2812891" "SRR2812959" "SRR2813002" "SRR2811815" "SRR2811826" "SRR2811827")

# Цикл по каждому индексу
for index in "${index_list[@]}"
do
    # Используйте sra-toolkit для скачивания чтений в формате fastq
    ~/gpfs/tools/sratoolkit.3.0.10-ubuntu64/bin/prefetch $index
    ~/gpfs/tools/sratoolkit.3.0.10-ubuntu64/bin/fasterq-dump $index
    rm -r $index
done
