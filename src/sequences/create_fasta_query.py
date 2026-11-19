import requests
import pandas as pd
from tqdm import tqdm
from pathlib import Path


'''
This script takes a list of IDs and returns a fasta file associated to a specific dataset
Note - this is for organisms that don't have a proteome fasta saved locally 
'''

UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{}.fasta"

def fetch_uniprot_sequences(ids, batch = 50):
    missing = []
    id_dict = {}
    headers = {}
    remapped = {}
    for i in tqdm(range(0, len(ids), batch)):
        
        chunk = ids[i : i + batch]
        #print(chunk)
        
        for uniprot_id in chunk:
            url = UNIPROT_FASTA.format(uniprot_id)
            response = requests.get(url, timeout=30)
            
            #print(response.text)

            if response.status_code == 200 and response.text.strip():                
                fasta = response.text
                #print(fasta)
                for line in fasta.splitlines():
                    if line.startswith(">"):
                        parts = line.split("|")
                        current_id = parts[1]

                        id_dict[current_id] = ""

                        headers[current_id] = line
                        if current_id != uniprot_id:
                            remapped[uniprot_id] = current_id

                    else:
                        # Append sequence to the current protein ID
                        if line.strip():
                            id_dict[current_id] += line.strip()
                            
                        # f.write(f"{line}\n")

            else:
                print(f"Failed to retrieve {uniprot_id}")
                missing.append(uniprot_id)
    return id_dict, headers, missing, remapped


def create_fasta(sequences, headers, output_file):
    with open(output_file, "w") as f:
        for uniprot_id, seq in sequences.items():
            f.write(f"{headers[uniprot_id]}\n{seq}\n")



if __name__ == "__main__":
    import argparse
    from fasta_from_local import get_ids, output_paths

    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    args = ap.parse_args()

    out_fasta, out_missing = output_paths(args.csv)
    sequences, headers, missing, remapped = fetch_uniprot_sequences(sorted(get_ids(args.csv)))
    create_fasta(sequences, headers, out_fasta)

    if missing:
        out_missing.write_text("\n".join(missing) + "\n")
    if remapped:
        pd.DataFrame(remapped.items(), columns=["source_id", "id"]).to_csv(
            out_fasta.with_name(out_fasta.stem + "_id_map.csv"), index=False)
    print(f"{len(sequences)} fetched, {len(missing)} missing, {len(remapped)} remapped -> {out_fasta}")