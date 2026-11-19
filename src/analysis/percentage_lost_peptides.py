"""Lost-cysteine summary for a chosen peptide length window.
 
Set the paths and N_MIN / N_MAX below, then run:  python lost_cys_summary.py
 
Reports: total Cys, Cys in lost (short / long) peptides, how many of those are modified
in your data (conflicts), and the final lost count with conflicts removed.
Assumes trypsin, 0 missed cleavages, and a CSV with columns: id, position (1-based).
"""
import gzip
 
import pandas as pd
from pyteomics import fasta, parser
 
# ---------------- settings ----------------
input_file = "/home/callum/sno_data/processed/pcys/pcys_processed.csv"
proteome = "/home/callum/human_proteome/proteome.gz"
out_fasta = "/home/callum/sno_data/processed/pcys.fasta"
N_MIN = 5      # peptides shorter than this are lost
N_MAX = 50     # peptides longer than this are lost
RESIDUE = "C"
# ------------------------------------------
 
 
def accession(description):
    # 'sp|P03372|ESR1_HUMAN Estrogen receptor ...' -> 'P03372'
    return description.split("|")[1]
 
 
def load_sites(path):
    df = pd.read_csv(path)
    if "organism" in df.columns:
        df = df[df["organism"] == "Homo sapiens"]
    df = df.dropna(subset=["id", "position"]).copy()
    df["id"] = df["id"].astype(str).str.strip()
    df["position"] = df["position"].astype(int)
    return df.drop_duplicates(["id", "position"])
 
 
def build_fasta(wanted, proteome_path, output_path):
    found, subset = set(), []
    with gzip.open(proteome_path, mode="rt") as f:
        for description, sequence in fasta.FASTA(f):
            acc = accession(description)
            if acc in wanted:
                subset.append((description, sequence))
                found.add(acc)
    fasta.write(subset, output_path)
    print(f"{len(found)} proteins found, {len(wanted - found)} not in proteome")
 
 
def digest(seqs, mod_set, residue):
    """One row per Cys: protein, position, length of the tryptic peptide it sits in, is it modified."""
    records = []
    for acc, seq in seqs.items():
        for start, pep in parser.icleave(seq, "trypsin", missed_cleavages=0):
            for i, aa in enumerate(pep):
                if aa == residue:
                    pos = start + i + 1                      # start is 0-based
                    records.append((acc, pos, len(pep), (acc, pos) in mod_set))
    return pd.DataFrame(records, columns=["id", "pos", "pep_len", "is_mod"])
 
 
def lost_summary(cys_df, n_min, n_max):
    short = cys_df[cys_df["pep_len"] < n_min]
    long_ = cys_df[cys_df["pep_len"] > n_max]
    res = {
        "total_cys": len(cys_df),
        "short_lost": len(short),
        "long_lost": len(long_),
        "total_lost": len(short) + len(long_),
        "short_conflicts": int(short["is_mod"].sum()),
        "long_conflicts": int(long_["is_mod"].sum()),
    }
    res["total_conflicts"] = res["short_conflicts"] + res["long_conflicts"]
    res["final_lost"] = res["total_lost"] - res["total_conflicts"]
    res["final_lost_pct_of_cys"] = round(100 * res["final_lost"] / res["total_cys"], 2)
    return res
 
 
if __name__ == "__main__":
    # 1. sites and FASTA
    df = load_sites(input_file)
    build_fasta(set(df["id"]), proteome, out_fasta)
    with open(out_fasta) as f:
        seqs = {accession(d): s for d, s in fasta.FASTA(f)}
 
    # 2. keep only sites that really are a Cys at the stated position
    sites = df[df["id"].isin(seqs)]
    sites = sites[[0 < p <= len(seqs[i]) and seqs[i][p - 1] == RESIDUE
                   for i, p in zip(sites["id"], sites["position"])]]
    mod_set = set(zip(sites["id"], sites["position"]))
    print(f"{len(sites)} modified sites are a {RESIDUE} at the stated position "
          f"({len(df[df['id'].isin(seqs)]) - len(sites)} dropped)")
 
    # 3. digest and summarise
    cys_df = digest(seqs, mod_set, RESIDUE)
    print(f"\nWindow: {N_MIN}-{N_MAX} residues")
    for k, v in lost_summary(cys_df, N_MIN, N_MAX).items():
        print(f"{k}: {v}")


print("number of cysteine sites in pcysmod:", len(df))
print("number of proteins in pcysmod:", df["id"].nunique())

