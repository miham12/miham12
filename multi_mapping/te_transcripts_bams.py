import pysam
import pandas as pd
import subprocess
import multiprocessing
import logging
import itertools
logging.basicConfig(filename='te_transcripts_bams.log', level=logging.INFO, format='%(asctime)s - %(message)s')
  
te="/home/miham/gpfs/mm1s_t2t/te_transcripts/TEtranscripts-master/bin/TEtranscripts"
repeats="/home/miham/gpfs/repeat_t2t_v2/te_annotation_2024/T2T-CHM13v2_rmsk_TE.gtf"
genes="/home/miham/gpfs/t2t_2_0_ncbi/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf"
file = '/home/miham/gpfs/mm1s_t2t_voting/voting_result_final/rna_mappings_no_mismatch.4-voted_SRR3633288.bed' # Замените на путь к вашему файлу с айдишниками

rna_voted = pd.read_csv(file,sep='\t',usecols=[3,10,11])
rna_voted = rna_voted.rename(columns={'id': 'read_id'})
rna_voted['gene_name']=rna_voted['gene_name']+'_'+rna_voted['gene_type']
del rna_voted['gene_type']

counts = rna_voted['gene_name'].value_counts()
gene_ids_to_keep = counts[counts > 500].index
rna_voted = rna_voted[rna_voted['gene_name'].isin(gene_ids_to_keep)]

logging.info(f"Dict prepared")

gene_read_dict = rna_voted.groupby('gene_name')['read_id'].agg(set).to_dict()
chunk_size = 4000
gene_read_chunks = [dict(itertools.islice(gene_read_dict.items(), i, i + chunk_size)) for i in range(0, len(gene_read_dict), chunk_size)]
logging.info(f"Dict was devided into chunks")

input_bam_file = "SRR3633289_voted.bam"

for chunk_index, chunk in enumerate(gene_read_chunks):
    input_bam = pysam.AlignmentFile(input_bam_file, "rb")
    logging.info(f"Starting the chunk {chunk_index}")
    output_bam_files={}
    for gene_id in chunk.keys():
        output_bam_files[gene_id] = pysam.AlignmentFile(f"./genes_bams_500/{gene_id}.bam", "wb", template=input_bam)

    for read in input_bam:
        read_id = read.query_name
        for gene_id, read_ids in chunk.items():
            if read_id in read_ids:
                output_bam_files[gene_id].write(read)
                break  # Exit loop once the read is written to the appropriate file
                
    for output_bam in output_bam_files.values():
        output_bam.close()
    input_bam.close()
logging.info(f"All bam files are ready")