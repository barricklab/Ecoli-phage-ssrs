#!/usr/bin/env python3

from Bio import SeqIO
from random import shuffle
import argparse
import csv
import random
import copy

#Returns a list of rows with columns length and count
max_ssr_length = 20

def count_SSRs(s, minimum_length):

    ssr_count_list = {}
    ssr_count_list['AT'] = [0] * max_ssr_length
    ssr_count_list['GC'] = [0] * max_ssr_length

    return_rows = []
    i = 0
    while i <  len(s):

        #biopython positions are 0-indexed
        homopolymer_start_position = i
        homopolymer_length = 1
        homopolymer_base = s[i]
        #print(homopolymer_start_position, homopolymer_base)
        i += 1
        while i < len(s) and s[i] == homopolymer_base:
            homopolymer_length += 1
            i += 1
        if homopolymer_base != 'N':
            if homopolymer_base in ['G', 'C']:
                ssr_count_list['GC'][homopolymer_length-1] += 1 
            if homopolymer_base in ['A', 'T']:
                ssr_count_list['AT'][homopolymer_length-1] += 1 

    return ssr_count_list

#Assigns a dictionary with keys as amino acids and items as lists of codons
def add_to_codon_list_per_amino_acid(aa_seq_string, nt_seq_string, codon_list_per_amino_acid):
    if codon_list_per_amino_acid is None:
        codon_list_per_amino_acid = {}

    remaining_nt_seq_string = nt_seq_string
    for this_aa in aa_seq_string:
        #print(this_aa)
        #print(remaining_nt_seq_string)
        this_codon=str(remaining_nt_seq_string[:3])
        #print (this_codon) 
        remaining_nt_seq_string = remaining_nt_seq_string[3:]
        if not this_aa in codon_list_per_amino_acid:
            codon_list_per_amino_acid[this_aa] = []
        codon_list_per_amino_acid[this_aa] += [this_codon]

def shuffle_codon_list_per_amino_acid(codon_list_per_amino_acid):
    #Randomize within each codon
    for this_aa in codon_list_per_amino_acid.keys():
        #print(codon_list_per_amino_acid[this_aa])
        random.shuffle(codon_list_per_amino_acid[this_aa])
        #print(codon_list_per_amino_acid[this_aa])


#Returns the codon list with the used codons removed
def sample_aa_seq_from_codon_list_per_amino_acid(aa_seq_string, codon_list_per_amino_acid):
    # Reassemble the new nucleotide sequence
    resampled_nt_seq = ''
    for this_aa in aa_seq_string:
        #print(this_aa)
        resampled_nt_seq += codon_list_per_amino_acid[this_aa].pop()
        #print(resampled_nt_seq)

    return (resampled_nt_seq)


def main():
    parser = argparse.ArgumentParser(
                        prog='randomize_genome_within_genes.py',
                        description='Shuffle codons within all genes in genome keeping amino acid sequence constant and examine homopolymer repeats',
                        epilog='')

    parser.add_argument("-i", "--input", type=str, required = True, help="Input nucleotide FASTA with one gene sequence")
    parser.add_argument("-o", "--output-actual", type=str, default = "genome_actual_within_gene.csv", help="Output actual SSR CSV file path")
    parser.add_argument("-p", "--output-resampled", type=str, default = "genome_resampled_within_gene.csv", help="Output resampled SSR CSV file path")
    parser.add_argument("-g", "--output-gc-content", type=str, default = "genome_gc_content.csv", help="Output resampled SSR CSV file path")
    parser.add_argument("-l", "--minimum-length", default=6, type=int, help="Minimum length of homopolymer repeat to report")
    parser.add_argument("-r", "--resamplings", type=int, default = 1000, help="Number of resamplings")
    parser.add_argument("-s", "--random-seed", type=int, default = None, help="Random seed. (0=system)")
    parser.add_argument("-m", "--randomization-method", type=str, default = "gene", help="'gene' or 'genome' to determine level of shuffling")

    args = parser.parse_args()
    input_filename = args.input
    minimum_length = args.minimum_length
    output_actual_filename = args.output_actual
    output_resampled_filename = args.output_resampled
    output_gc_content_filename = args.output_gc_content
    num_replicates = args.resamplings
    randomization_method = args.randomization_method
    if not randomization_method in ['gene', 'genome']:
        print("Invalid --randomization-method setting: " + randomization_method)

    if not args.random_seed is None:
        random.seed(args.random_seed)

    output_rows = []

    # Read in nucleotide sequences of all genes
    nt_seq_strings = []
    aa_seq_strings = []
    locus_tags = []
    seq_ids = []
    overall_seq_ids = []
    gc_base_count = 0
    at_base_count = 0
    other_base_count = 0
    total_base_count = 0
    for record in SeqIO.parse(input_filename, "genbank"):
        seq_id = record.id
        overall_seq_ids.append(seq_id)
        for feature in record.features:
            if feature.type == 'CDS':

                if 'protein_id' in feature.qualifiers.keys():

                    locus_tag = feature.qualifiers['protein_id'][0]

                    # extract amino acid and nucleotide sequences
                    nt_seq = feature.extract(record)
                    aa_seq = nt_seq.translate()

                    # Use strings
                    nt_seq_string = nt_seq.seq.upper()
                    aa_seq_string = aa_seq.seq.upper()

                    # Remove stop codon? - not really necessary
                    #if aa_seq_string[-1] == '*':
                    #    aa_seq_string = aa_seq_string[:-1]
                    #    nt_seq_string = nt_seq_string[:-3]

                    nt_seq_strings.append(nt_seq_string)
                    aa_seq_strings.append(aa_seq_string)
                    locus_tags.append(locus_tag)
                    seq_ids.append(seq_id)

                    total_base_count = total_base_count + len(nt_seq_string)
                    for this_nt in nt_seq:
                        if this_nt in ['A', 'T']:
                            at_base_count = at_base_count + 1
                        elif this_nt in ['G', 'C']:
                            gc_base_count = gc_base_count + 1
                        else:
                            other_base_count = other_base_count + 1

    # Record GC content (lumping all seq_ids in genome)
    gc_percent = 100 * (gc_base_count / (gc_base_count + at_base_count))
    gc_content_file = open(output_gc_content_filename, mode='w', newline='')
    gc_content_writer = csv.DictWriter(gc_content_file, fieldnames=['seq_id', 'gc_percent', "total_bases", "other_bases"])
    gc_content_writer.writeheader() 
    for seq_id in overall_seq_ids:
        gc_content_writer.writerows([{'seq_id' : seq_id, 'gc_percent' : gc_percent, 'total_bases' : total_base_count, 'other_bases' : other_base_count}]) 

    # Count the actual totals per gene
    actual_file = open(output_actual_filename, mode='w', newline='')
    actual_writer = csv.DictWriter(actual_file, fieldnames=['seq_id', 'locus_tag', 'length', 'base_pair', 'count'])
    actual_writer.writeheader()
    for i in range(len(aa_seq_strings)):

        aa_seq_string = aa_seq_strings[i]
        nt_seq_string = nt_seq_strings[i]
        locus_tag = locus_tags[i]
        seq_id = seq_ids[i]

        ssr_counts = count_SSRs(nt_seq_string, minimum_length)
        for base_pair in ['GC', 'AT']:
            for i in range(minimum_length-1, max_ssr_length):
                count = ssr_counts[base_pair][i]
                if count == 0:
                    continue
                row = {}
                row['seq_id'] = seq_id
                row['locus_tag'] = locus_tag
                row['length'] = i+1
                row['base_pair'] = base_pair
                row['count'] = count
                actual_writer.writerow(row)
        

    actual_file.close()



    resampled_file = open(output_resampled_filename, mode='w', newline='')
    writer = csv.DictWriter(resampled_file, fieldnames=['replicate', 'seq_id', 'locus_tag', 'length', 'base_pair', 'count'])
    writer.writeheader()

    # Count the resampled totals per gene, writing as we go


    if randomization_method == "gene":

        # Set up master codon lists per gene. We will deep copy shuffle and pop from these
        codon_list_per_amino_acid_list = []
        for i in range(len(aa_seq_strings)):
            codon_list_per_amino_acid = {}
            aa_seq_string = aa_seq_strings[i]
            nt_seq_string = nt_seq_strings[i]
            add_to_codon_list_per_amino_acid(aa_seq_string, nt_seq_string, codon_list_per_amino_acid)
            codon_list_per_amino_acid_list.append(codon_list_per_amino_acid)

        for r in range(num_replicates):
            for i in range(len(aa_seq_strings)):

                aa_seq_string = aa_seq_strings[i]
                locus_tag = locus_tags[i]
                seq_id = seq_ids[i]

                codon_list_per_amino_acid = copy.deepcopy(codon_list_per_amino_acid_list[i])
                shuffle_codon_list_per_amino_acid(codon_list_per_amino_acid)
                resampled_nt_seq_string = sample_aa_seq_from_codon_list_per_amino_acid(aa_seq_string, codon_list_per_amino_acid)

                ssr_counts = count_SSRs(resampled_nt_seq_string, minimum_length)
                for base_pair in ['GC', 'AT']:
                    for i in range(minimum_length-1, max_ssr_length):
                        count = ssr_counts[base_pair][i]
                        if count == 0:
                            continue
                        row = {}
                        row['replicate'] = r+1
                        row['seq_id'] = seq_id
                        row['locus_tag'] = locus_tag
                        row['length'] = i+1
                        row['base_pair'] = base_pair
                        row['count'] = count
                        writer.writerow(row)
               
        
            if (r+1) % 10==0:
                print("Processed replicates: " + str(r+1))

    elif randomization_method == "genome":

        # Set up master codon list per genome. We will deep copy shuffle and pop from this one list
        master_codon_list_per_amino_acid = {}
        for i in range(len(aa_seq_strings)):
            aa_seq_string = aa_seq_strings[i]
            nt_seq_string = nt_seq_strings[i]
            add_to_codon_list_per_amino_acid(aa_seq_string, nt_seq_string, master_codon_list_per_amino_acid)

        for r in range(num_replicates):
            
            codon_list_per_amino_acid = copy.deepcopy(master_codon_list_per_amino_acid)
            shuffle_codon_list_per_amino_acid(codon_list_per_amino_acid)

            for i in range(len(aa_seq_strings)):

                aa_seq_string = aa_seq_strings[i]
                locus_tag = locus_tags[i]
                seq_id = seq_ids[i]

                resampled_nt_seq_string = sample_aa_seq_from_codon_list_per_amino_acid(aa_seq_string, codon_list_per_amino_acid)

                ssr_counts = count_SSRs(resampled_nt_seq_string, minimum_length)
                for base_pair in ['GC', 'AT']:
                    for i in range(minimum_length-1, max_ssr_length):
                        count = ssr_counts[base_pair][i]
                        if count == 0:
                            continue
                        row = {}
                        row['replicate'] = r+1
                        row['seq_id'] = seq_id
                        row['locus_tag'] = locus_tag
                        row['length'] = i+1
                        row['base_pair'] = base_pair
                        row['count'] = count
                        writer.writerow(row)
               
        
            if (r+1) % 10==0:
                print("Processed replicates: " + str(r+1))




if __name__ == '__main__':
    main()
