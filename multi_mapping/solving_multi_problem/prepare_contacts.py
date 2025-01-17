import pandas as pd
import argparse
import logging
import os 
logging.basicConfig(filename='prepare_contacts.log', level=logging.INFO, format='%(asctime)s - %(message)s')
'''
input : voted rna file,
'''

def main(rna_path, dna_mapped_path, min_num_contacts):
    logging.info(f"Reading RNA voted file ")
    rna = pd.read_csv(rna_path, sep='\t', header=0, usecols=[0, 1, 2, 3, 4, 10, 11])
    logging.info(f"Leaving high contacting RNA")
    gene_counts = rna['gene_name'].value_counts()
    valid_genes = gene_counts[gene_counts >= min_num_contacts].index
    rna = rna[rna['gene_name'].isin(valid_genes)]
    rna.rename(columns={'id': 'read_id'}, inplace=True)
    del gene_counts
    del valid_genes
    logging.info(f"Readng DNA mappings")
    dna_mapped = pd.read_csv(dna_mapped_path, sep='\t', usecols=[0, 1, 2, 3, 5], names=['dna_chr', 'dna_start', 'dna_end', 'read_id', 'dna_strand'])
    logging.info(f"Merging RNA with DNA ")
    contacts = rna.merge(dna_mapped, on='read_id', how='left')
    logging.info(f"Filling NAs")
    contacts = contacts.fillna('.')

    output_dir = os.path.dirname(rna_path)
    output_filename = os.path.basename(rna_path).replace('.bed', '_pre_solved_contacts.bed')
 #   output_path = os.path.join(output_dir, output_filename)
    output_path = output_filename
    logging.info(f"Saving result")
    contacts.to_csv(output_path, header=True, index=None, sep='\t')
    logging.info(f"Script finished")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process RNA and DNA mapping data.')
    parser.add_argument('rna_path', type=str, help='Path to the RNA file')
    parser.add_argument('dna_mapped_path', type=str, help='Path to the DNA mapped file')
    parser.add_argument('min_num_contacts', type=int, help='Minimum number of contacts')

    args = parser.parse_args()
    main(args.rna_path, args.dna_mapped_path, args.min_num_contacts)