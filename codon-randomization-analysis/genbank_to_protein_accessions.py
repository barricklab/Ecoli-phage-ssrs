#!/usr/bin/env python3

import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import csv
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True)
parser.add_argument("-o", "--output", type=str, required=True)
args = parser.parse_args()

all_entries = []
for filename in os.listdir(args.input):
    if not filename.endswith(".gbk"):
        continue
    full_path = os.path.join(args.input, filename)
    if os.path.isfile(full_path):
        print(full_path)
        for record in SeqIO.parse(full_path, "genbank"):
            for feature in record.features:
                if feature.type == 'CDS':
                    if 'protein_id' in feature.qualifiers.keys():
                        #print(feature)
                        protein_id = feature.qualifiers['protein_id'][0]
                        all_entries.append({'seq_id' : record.id, 'protein_id' : protein_id})
                        # print(cds)

with open(args.output, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=all_entries[0].keys())
    writer.writeheader()
    writer.writerows(all_entries)