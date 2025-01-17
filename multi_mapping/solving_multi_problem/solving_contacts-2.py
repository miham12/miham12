import pandas as pd
from intervaltree import Interval, IntervalTree
import os
import gc
import argparse
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import logging
import shutil
from prepare_chr_probs import *


def cleanup_directory(directory):
    logging.info(f"Deleting directory: {directory}")
    shutil.rmtree(directory)
    logging.info(f"Deleted directory: {directory}")
    
def preprocess_contacts(filename):
    logging.info(f"Preprocessing contacts from file: {filename}")
    processed_chunks = []
    for chunk in pd.read_csv(filename, sep='\t', header=0, chunksize=30000000):
        logging.info(f"Reading a chunk of size 50000000 from {filename}")
        chunk = chunk[chunk['dna_chr'] != '.']
        processed_chunks.append(chunk)
        
    contacts = pd.concat(processed_chunks, ignore_index=True)
    logging.info(f"Filtered . contacts")
    
    contacts['dna_start'] = contacts['dna_start'].astype(float).astype(int)
    contacts['dna_end'] = contacts['dna_end'].astype(float).astype(int)
    read_id_counts = contacts['read_id'].value_counts()
    unique_read_ids = read_id_counts[read_id_counts == 1].index
    multiple_occurrences_read_ids = read_id_counts[read_id_counts > 1].index
    
    logging.info(f"Divided unique and multi contacts ")
    unique_contacts = contacts[contacts['read_id'].isin(unique_read_ids)]
    multi_contacts = contacts[contacts['read_id'].isin(multiple_occurrences_read_ids)]
    logging.info(f"Preprocessed contacts from file: {filename}")
    
    return unique_contacts, multi_contacts

def build_interval_trees(unique_df):
    interval_trees = {}
    for chr_key, group in unique_df.groupby('dna_chr'):
        interval_trees[chr_key] = IntervalTree(
            Interval(row['dna_start'], row['dna_end'], row['read_id']) for _, row in group.iterrows()
        )
    return interval_trees

def calculate_weights(unique_df, multi_df, size, interval_trees):
    def count_intersections(row):
        chr_key = row['dna_chr']
        if chr_key in interval_trees:
            start = row['dna_start']
            end = row['dna_end']
            intersected_intervals = interval_trees[chr_key][start-size:end+size]
            read_center = ((row['dna_end'] + row['dna_start']) // 2)
            weight_sum = 1000 * sum(1 / (abs(read_center - ((interval.end + interval.begin) // 2))+1) for interval in intersected_intervals)
            return weight_sum
        else:
            return 0

    weights_series = multi_df.apply(lambda row: count_intersections(row), axis=1)
    
    return weights_series

def process_gene(gene_name, gene_unique_contacts, gene_multi_contacts, chr_contacts_probs, size, gene_results_dir):
    try:
        logging.info(f"Processing gene: {gene_name}")
        anything_saved=True
        if gene_unique_contacts.empty:
            logging.info(f"No unique-mapping reads for gene: {gene_name}. Writing Nothing.")
            return
        if gene_multi_contacts.empty:
            logging.info(f"No multi-mapping reads available for gene: {gene_name}. Writing only unique contacts.")
            output_filename = os.path.join(gene_results_dir, f'{gene_name}_results.csv')
            with open(output_filename, 'w') as f:
                gene_unique_contacts[f'intersections_count_{size}'] = float('inf')
                gene_unique_contacts.to_csv(f, index=False, sep='\t')
            return
        interval_trees = build_interval_trees(gene_unique_contacts)
        
        weights_series = calculate_weights(gene_unique_contacts, gene_multi_contacts, size, interval_trees)

        gene_multi_contacts[f'intersections_count_{size}'] = weights_series
        gene_unique_contacts[f'intersections_count_{size}'] = float('inf')
    
        gene_multi_contacts = gene_multi_contacts[gene_multi_contacts[f'intersections_count_{size}'] != 0]
        
        if len(gene_multi_contacts) > 0:
            gene_multi_contacts.loc[:, f'intersections_count_{size}'] = gene_multi_contacts.apply(lambda row: row[f'intersections_count_{size}'] * chr_contacts_probs[(row['dna_chr'], row['rna_chr'])],axis=1)
            gene_multi_contacts = gene_multi_contacts.loc[gene_multi_contacts.groupby('read_id')[f'intersections_count_{size}'].idxmax()]
        else:
            logging.info(f"No multi-mapping saved for gene: {gene_name}. Writing only unique contacts.")
            anything_saved=False
        output_filename = os.path.join(gene_results_dir, f'{gene_name}_results.csv')
    
        with open(output_filename, 'w') as f:
            gene_unique_contacts.to_csv(f, index=False, sep='\t')
            if anything_saved:
                gene_multi_contacts.to_csv(f, index=False, header=False, sep='\t')
        
        del gene_unique_contacts, gene_multi_contacts, weights_series, interval_trees
        logging.info(f"Processed gene: {gene_name}")
        gc.collect()
        
    except Exception as e:
        logging.error(f"Error processing gene {gene_name}: {e}")

def process_genes(filename, divisors, size, max_workers):
    unique_contacts, multi_contacts = preprocess_contacts(filename)

    logging.info(f"Counting probabilities of contacts based on chrs")
    chr_contacts_probs = chr_contacts_normalizing(count_chr_contacts(unique_contacts), divisors)
    chr_contacts_probs = chr_contacts_probs.stack().to_dict()
    logging.info(f"Counted probabilities of contacts based on chrs")
    
    base_filename = os.path.splitext(os.path.basename(filename))[0]
    gene_results_dir = f'gene_results_{base_filename}_{size}'
    
    logging.info(f"Creating directory {gene_results_dir}")
    os.makedirs(gene_results_dir, exist_ok=True)
    logging.info(f"Created directory {gene_results_dir}")

    gene_names = unique_contacts['gene_name'].unique()
    with mp.Pool(processes=max_workers) as pool:
        for gene_name in gene_names:
            gene_unique_contacts = unique_contacts[unique_contacts['gene_name'] == gene_name]
            gene_multi_contacts = multi_contacts[multi_contacts['gene_name'] == gene_name]
            pool.apply_async(process_gene, (gene_name, gene_unique_contacts, gene_multi_contacts, chr_contacts_probs, size, gene_results_dir))
        logging.info(f"Processed all genes")
        pool.close()
        pool.join() 

def combine_results(output_file, gene_results_dir):
    logging.info(f"Combining results")
    files = [f for f in os.listdir(gene_results_dir) if f.endswith('_results.csv')]
    combined_df = pd.concat([pd.read_csv(os.path.join(gene_results_dir, f), sep='\t') for f in files], ignore_index=True)
    combined_df.to_csv(output_file, sep='\t', header=True, index=False)
    logging.info(f"Script finished")
    
def main():
    parser = argparse.ArgumentParser(description="Process and analyze contact data.")
    parser.add_argument('input_file', type=str, help='Path to the input contact file')
    parser.add_argument('divisors_path', type=str, help='Divisors for normalizing the num of contacts')
    parser.add_argument('-s', '--size', type=int, default=10000, help='Range for reads consideration (default: 10000)')
    parser.add_argument('-p', '--processes', type=int, default=1, help='Number of parallel processes to use (default: 1)')
    args = parser.parse_args()

    input_file = args.input_file
    size = args.size
    max_workers = args.processes
    divisors = args.divisors_path
    divisors = pd.read_csv(divisors, sep='\t', index_col=0)

    base_filename = os.path.splitext(os.path.basename(input_file))[0]
    log_filename = f'solving_contacts_{base_filename}_{size}.log'

    logging.basicConfig(
        filename=log_filename,
        filemode='w',  
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s:%(message)s'
    )
    
    process_genes(input_file, divisors, size, max_workers)
    
    output_file = input_file.replace('pre', str(size)) 
    gene_results_dir = f'gene_results_{base_filename}_{size}'
    combine_results(output_file, gene_results_dir)
    cleanup_directory(gene_results_dir)
if __name__ == '__main__':
    main()

