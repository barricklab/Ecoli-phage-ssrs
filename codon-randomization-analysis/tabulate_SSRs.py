#!/usr/bin/env python3

import argparse
import pandas
import glob
import os.path
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SimpleLocation
from Bio.SeqFeature import ExactPosition

parser = argparse.ArgumentParser(description="Survey homopolymer CL candidates across genomes")
parser.add_argument("-i", "--input", help="Input folder of Genbank files", action="store", required=True)
parser.add_argument("-o", "--output", help="Output CSV file", action="store", default="output.csv")
parser.add_argument("-l", "--minimum-length", help="Minimum length of homopolymer repeat to report", action="store", default=8, type=int)

args = parser.parse_args()

full_list = []
genbank_input_files=glob.glob(os.path.join(args.input, "*.gbk"))
for g in genbank_input_files:
    p = re.compile("/(.+).gbk")
    genome_id = p.search(g).group(1)
    print(genome_id)
    this_list = []
    for record in SeqIO.parse(g, "genbank"):
        seq_id = record.id
        s = str(record.seq)
        current_base = s[0]
        i = 0
        while i <  len(s):

            #biopython positions are 0-indexed
            homopolymer_start_position = i
            homopolymer_length = 1
            homopolymer_base = s[i]
            #print(homopolymer_start_position, homopolymer_base)
            i = i + 1
            while i < len(s) and s[i] == homopolymer_base:
                homopolymer_length = homopolymer_length + 1
                i = i + 1
            if homopolymer_base != 'N' and homopolymer_length >= args.minimum_length:

                homopolymer_end_position = i - 1
                # Now assign genes that it's within or upstream
                overlapping_gene = None
                promoter_gene = None
                gene_position = None
                original_homopolymer_base = homopolymer_base
                
                for feature in record.features:
                    if feature.type == 'CDS':
                        locus_tag = "undefined"
                        if 'protein_id' in feature.qualifiers.keys():
                            locus_tag = feature.qualifiers['protein_id'][0]
                        current_gene_length = feature.location.end - feature.location.start + 1
                        if feature.location.strand == +1:
                            upstream_region = SimpleLocation(ExactPosition(feature.location.start-201), ExactPosition(feature.location.start-1))
                            current_gene_position = homopolymer_start_position - feature.location.start
                            current_homopolymer_base = original_homopolymer_base
                        elif feature.location.strand == -1:
                            upstream_region = SimpleLocation(ExactPosition(feature.location.end+1), ExactPosition(feature.location.end+201))
                            current_gene_position = feature.location.end - (homopolymer_start_position + homopolymer_length)
                            #base should be on gene strand
                            current_homopolymer_base = str(Seq(original_homopolymer_base).reverse_complement())
                            
                        if homopolymer_start_position in feature:
                            if overlapping_gene == None or gene_position < gene_position:
                                promoter_gene = None
                                overlapping_gene = locus_tag
                                gene_position = current_gene_position
                                gene_length = current_gene_length
                                homopolymer_base = current_homopolymer_base

                        elif homopolymer_end_position in feature:
                            if overlapping_gene == None or gene_position < gene_position:
                                promoter_gene = None
                                overlapping_gene = locus_tag
                                gene_position = current_gene_position
                                gene_length = current_gene_length
                                homopolymer_base = current_homopolymer_base

                        if overlapping_gene == None:
                            if homopolymer_start_position in upstream_region:
                                if promoter_gene == None or gene_position > gene_position:
                                    promoter_gene = locus_tag
                                    gene_position = current_gene_position
                                    gene_length = current_gene_length
                                    homopolymer_base = current_homopolymer_base

                            elif homopolymer_end_position in upstream_region:
                                if promoter_gene == None or gene_position > gene_position:
                                    promoter_gene = locus_tag
                                    gene_position = current_gene_position
                                    gene_length = current_gene_length
                                    homopolymer_base = current_homopolymer_base



                location_type = ""
                if promoter_gene != None:
                    location_type = "promoter"
                if overlapping_gene != None:
                    location_type = "gene"

                if overlapping_gene == None:
                    overlapping_gene = ""

                if promoter_gene == None:
                    promoter_gene = ""

                gene = overlapping_gene
                if gene == "":
                    gene = promoter_gene

                if gene_position == None:
                    gene_position = ""
                    gene_length = ""

                base_pair = "NA"
                if homopolymer_base in ['G', 'C']:
                    base_pair = 'GC'
                elif homopolymer_base in ['A', 'T']:
                    base_pair = 'AT'

                new_data = { 
                    'genome_id' : genome_id,
                    'seq_id' : record.id,
                    'position' : homopolymer_start_position+1, #switch bach to 1-indexed
                    'base' : homopolymer_base,
                    'base_pair' : base_pair,
                    'length' : homopolymer_length,
                    'gene' : gene,
                    'gene_position' : gene_position,
                    'gene_length' : gene_length,
                    'location_type' : location_type
                }
                #print(new_data)
                new_df = pandas.DataFrame(new_data, index=[0])
                #print(new_df)
                this_list.append(new_df)

    if len(this_list) > 0:
        this_df = pandas.concat(this_list)
        full_list.append(this_df)


full_df = pandas.concat(full_list, ignore_index=True)
full_df.to_csv(args.output, index=False)
