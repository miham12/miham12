from pyfaidx import Fasta
import pandas as pd
import argparse
import logging
import os 
def obtain_chr_CR_lengths(bed_path):
    '''
    input: path to bed file (converted from gtf file)
    output: length of coding regions of chrs
    '''
    bed=pd.read_csv(bed_path,sep='\t',header=None)
    bed=bed[bed[7]=='gene']
    bed['length']= bed[2]-bed[1]
    bed=bed.groupby(by=[0])['length'].sum()
    bed = bed.reset_index()
    bed.columns = ['rna_chr', 'size']
    return bed

def obtain_chr_lengths(genome_path):
    '''
    input: path to fasts file 
    output: chr lengths
    '''
    fasta = Fasta(genome_path)
    chromosome_sizes = {chrom: len(fasta[chrom]) for chrom in fasta.keys()}
    df = pd.DataFrame(list(chromosome_sizes.items()), columns=['dna_chr', 'size'])
    return df

def multiply_lengths(chr_lengths, chr_CR_lengths):
    '''
    input: lengths of chrs, CR lengths of chrs
    output: multiplyes lengths matrix 
    '''
    result_df = pd.DataFrame(index=chr_lengths['dna_chr'], columns=chr_CR_lengths['rna_chr'])
    
    for i, row_L in chr_lengths.iterrows():
        for j, row_CR in chr_CR_lengths.iterrows():
            result_df.at[row_L['dna_chr'], row_CR['rna_chr']] = row_L['size'] * row_CR['size']

    return result_df

    
def main(bed_path, genome_path):
    chr_CR_lengths=obtain_chr_CR_lengths(bed_path)
    chr_lengths=obtain_chr_lengths(genome_path)
    result_df=multiply_lengths(chr_lengths,chr_CR_lengths)
    result_df.to_csv("divisors.tab",sep='\t',header=True,index=True)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process RNA and DNA mapping data.')
    parser.add_argument('bed_path', type=str, help='Path to bed file converted from gtf')
    parser.add_argument('genome_path', type=str, help='Path to genome file')

    args = parser.parse_args()
    main(args.bed_path, args.genome_path)