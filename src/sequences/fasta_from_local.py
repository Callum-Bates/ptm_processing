# create fasta file here based on the chosen df
import os
from urllib.request import urlretrieve
import gzip
import matplotlib.pyplot as plt
import numpy as np
from pyteomics import fasta, parser, mass, achrom, electrochem, auxiliary
import pandas as pd
import re 
from pathlib import Path
import argparse

'''
This script makes FASTA file from list of UniProt IDs
Note - this is for organisms for which the proteome fasta file is already downloaded. 
If it is not, then use other function which gets sequence by querying uniprot
'''


# create fasta file here based on the chosen df

#out_fasta = "/home/callum/sno_data/processed/P03372.fasta"
#input_file = '/home/callum/sno_data/processed/pcys/pcys_org.csv'

DEFAULT_PROTEOME = "/home/callum/human_proteome/proteome.gz"

#wanted = df[df["organism"] == "Homo sapiens"]

def get_ids(csv_path, id_col="id"):
    df = pd.read_csv(csv_path)
    return set(df[id_col].dropna().astype(str).str.strip())


# example of df path
#x = get_id_list("df", input_file, lst = None)
# print(x)
#example of list
# x = get_id_list("lst", lst = ["P03372"])

# print(x)
#wanted = "P03372"

def accession(description):
    # 'sp|P03372|ESR1_HUMAN Estrogen receptor ...' -> 'P03372'
    return description.split('|')[1]

def open_fasta(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)




def fasta_from_local(wanted, proteome):
	records = []

	with open_fasta(proteome) as f:
		for description, sequence in fasta.FASTA(f):
			if accession(description) in wanted:
				records.append((description, sequence))
	found = {accession(d) for d, _ in records}
	return records, wanted - found

def output_paths(csv_path):
    """yang/yang_processed.csv -> yang/yang.fasta, yang/yang_missing.txt"""
    csv_path = Path(csv_path)
    name = csv_path.stem.removesuffix("_processed")
    return csv_path.with_name(f"{name}.fasta"), csv_path.with_name(f"{name}_missing.txt")
 
 
def main(csv_path, proteome):
    out_fasta, out_missing = output_paths(csv_path)
    wanted = get_ids(csv_path)
    records, missing = fasta_from_local(wanted, proteome)
 
    fasta.write(records, str(out_fasta), file_mode="w")
    if missing:
        out_missing.write_text("\n".join(sorted(missing)) + "\n")
 
    print(f"{len(records)} / {len(wanted)} found -> {out_fasta}")
    if missing:
        print(f"{len(missing)} not in proteome -> {out_missing}")
 
 
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", help="processed CSV with an 'id' column")
    ap.add_argument("--proteome", default=DEFAULT_PROTEOME)
    args = ap.parse_args()
    main(args.csv, args.proteome)

 
