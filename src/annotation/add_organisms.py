import argparse
import io
import os
import time
import requests
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("input_file")
parser.add_argument("--dir_name", default=".")
parser.add_argument("--out_file", default=".")

args = parser.parse_args()


BASE = "https://rest.uniprot.org/idmapping"
BATCH_SIZE = 500
MASTER_PATH = "/home/callum/sno_data/processed/mapping/uniprot_organisms.csv"  # absolute path, adjust as needed
PROCESSED_PATH = "/home/callum/sno_data/processed/"






def get_ids(file):
    
    df = pd.read_csv(file)
    ids = df["id"].dropna().astype(str).unique().tolist()
    
    return ids


def fetch_organisms(ids):
    """Query UniProt for a list of IDs.

    Returns (mapping, failed): a dict of {input_id: scientific_name}
    and a list of IDs that UniProt could not map.
    """
    mapping, failed = {}, []

    for start in range(0, len(ids), BATCH_SIZE):
        batch = ids[start:start + BATCH_SIZE]
        batch = [i.strip() for i in batch if i and i.strip() and i.lower() != "nan"]

        # Submit job
        r = requests.post(
            f"{BASE}/run",
            data={"from": "UniProtKB_AC-ID", "to": "UniProtKB", "ids": ",".join(batch)},
        )
        r.raise_for_status()
        job_id = r.json()["jobId"]

        # Poll until finished
        while True:
            s = requests.get(f"{BASE}/status/{job_id}")
            s.raise_for_status()
            job_status = s.json().get("jobStatus")

            if job_status in ("NEW", "RUNNING"):
                time.sleep(1)
                continue
            if job_status is not None:  # e.g. "ERROR"
                raise RuntimeError(f"Job {job_id} ended with status {job_status}")
            break  # no jobStatus key means results are ready

        # Fetch accession + organism as JSON
        res = requests.get(
            f"{BASE}/uniprotkb/results/stream/{job_id}",
            params={"format": "json", "fields": "accession,organism_name"},
        )
        res.raise_for_status()
        data = res.json()


        for entry in data.get("results", []):
            org = entry["to"].get("organism", {}).get("scientificName")
            if org is None:
                failed.append(entry["from"])
                continue
            mapping[entry["from"]] = org
            
            
        # for entry in data.get("results", []):
        #     mapping[entry["from"]] = entry["to"]["organism"]["scientificName"]
        # failed.extend(data.get("failedIds", []))

    return mapping, failed





def add_organisms(df, id_col="id", out_col="organism", master_path=MASTER_PATH):
    """Add an organism column to df, querying UniProt only for IDs not already
    in the master CSV. Updates the master CSV with any newly mapped IDs."""
    # Load existing master mapping if there is one
    if os.path.exists(master_path):
        master = pd.read_csv(master_path)
        print(master.head())
    else:
        master = pd.DataFrame(columns=["id", "organism"])

    known = dict(zip(master["id"], master["organism"]))
    

    # Only query IDs we haven't mapped before
    ids = df[id_col].dropna().astype(str).unique().tolist()
    new_ids = [i for i in ids if i not in known]
    print(f"{len(known)} already mapped, {len(new_ids)} to query")

    if new_ids:
        new_mapping, failed = fetch_organisms(new_ids)
        known.update(new_mapping)

        if failed:
            print(f"{len(failed)} IDs failed to map, e.g. {failed[:5]}")

        pd.DataFrame(known.items(), columns=["id", "organism"]).to_csv(
            master_path, index=False
        )

    df[out_col] = df[id_col].map(known)
    return df


def save_new_df(df, dir_name, file_name, dir_path = PROCESSED_PATH):
    df.to_csv(os.path.join(dir_path,dir_name,file_name), index = False)


df = pd.read_csv(os.path.join(PROCESSED_PATH, args.dir_name, args.input_file))
new_df = add_organisms(df)
save_new_df(new_df, args.dir_name, args.out_file)


