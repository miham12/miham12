#!/bin/bash
#SBATCH --job-name=indexing
#SBATCH --error=indexing_err
#SBATCH --output=indexing_log
#SBATCH --time=900
#SBATCH --nodes=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=10GB
source ~/.bashrc
FASTA_File="/home/miham/gpfs/tribolium_castaneum/genome/GCA_000002335.3_Tcas5.2_genomic.fna"
GFF="/home/miham/gpfs/tribolium_castaneum/genome/GCA_000002335.3_Tcas5.2_genomic.fna"
SPLICE="splice_site"
EXON="exon"
BASE_NAME="castaneum"
/home/miham/gpfs/tools/hisat2/hisat2_extract_splice_sites.py ${GFF} > ${SPLICE}
/home/miham/gpfs/tools/hisat2/hisat2_extract_exons.py ${GFF} > ${EXON}
/home/miham/gpfs/tools/hisat2/hisat2-build -p 18 --exon ${EXON} --ss ${SPLICE} ${FASTA_File} ${BASE_NAME}
echo 'Script finished' >> index_log
