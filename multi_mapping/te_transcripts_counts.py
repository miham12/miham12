import pandas as pd 
import re 
import os
import time
import logging
import subprocess
import multiprocessing
from multiprocessing import Pool
#start = time.time()
logging.basicConfig(filename='te_transcripts_counts.log', level=logging.INFO, format='%(asctime)s - %(message)s')
logging.info(f"Script for TEtranscripts starts")
def give_time(start):
    end = time.time() - start
    print(end)
    
def process_gene(gene_id):
    logging.info(f"Starting process for {gene_id}")
    try:
        subprocess.run([te, "-b", f"./genes_bams_500/{gene_id}.bam",'--outdir','./genes_te_transcripts_again/','--project',f'{gene_id}', "--GTF", genes, "--TE", repeats,'--mode',"multi"])
        logging.info(f"Process for {gene_id} completed successfully")
    except Exception as e:
        logging.error(f"Error in process for {gene_id}: {str(e)}")
    
 #   os.unlink(f'./genes_te_transcripts/{gene_id}_DESeq2.R')
    
        
        
te="/home/miham/gpfs/mm1s_t2t/te_transcripts/TEtranscripts-master/bin/TEcount"
repeats="/home/miham/gpfs/repeat_t2t_v2/te_annotation_2024/T2T-CHM13v2_rmsk_TE.ind"
genes="/home/miham/gpfs/t2t_2_0_ncbi/GCF_009914755.1_T2T-CHM13v2.0_genomic.ind"

dir_path = './genes_bams_500/'
gene_ids = [f.rstrip('.bam') for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]

num_processes = 20 
with Pool(processes=num_processes) as pool:
    pool.map(process_gene, gene_ids)
  #  for gene_id in gene_ids:
  #      pool.map(process_gene, gene_ids)
        
logging.info(f"Закончили создавать таблицы с каунтами, го мерджить")
dfs = []
dfs.append(pd.read_csv(f'./genes_te_transcripts_again/{gene_ids[0]}.cntTable', sep='\t',skiprows=1, usecols=[0, 1],names=['gene/TE',f'{gene_ids[0]}']))
           
for gene_id in gene_ids[1:]:
    df = pd.read_csv(f'./genes_te_transcripts_again/{gene_id}.cntTable', sep='\t', usecols=[1],names=[f'{gene_id}'],skiprows=1)
    dfs.append(df)
logging.info(f"Смерджили")


merged_df = pd.concat(dfs, axis=1)
merged_df.to_csv('te_transcripts_all.csv',sep='\t',header=True,index=None)
logging.info(f"Скрипт окончен ура ура ")