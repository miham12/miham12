import pandas as pd
def count_chr_contacts(unique_contacts):
    '''
    input: prepared pre-solved unique contacts 
    output: matrix of number of contacts among chromosomes
    '''
    сombinations_count = unique_contacts.groupby(['dna_chr', 'rna_chr']).size().reset_index(name='count')
    pivot_table = сombinations_count.pivot_table(index='dna_chr', columns='rna_chr', values='count', fill_value=0)
    return pivot_table

def chr_contacts_normalizing(df, divisors):
    '''
    input: matrix of number of contacts among chromosomes, normalization matrix 
    output: probability of contact rna_chr with dna_chr
    '''
    df_safe = df.copy()
    df = df * 100000000000
    df = df / divisors
    column_sums = df.sum()
    df_normalized = df.div(column_sums)
    return df_normalized