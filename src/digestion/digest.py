"""Check how many modified Cys sites fall on peptides trypsin would make undetectable.

Usage:
	python tryptic_check.py --dataset dn --input /path/to/sites.csv --min_len 7 --max_len 50 --missed 1

Input CSV needs columns: id (UniProt accession), position (1-based residue number).
An optional 'organism' column is used to keep only Homo sapiens rows.
"""
import argparse
import gzip
import json
from dataclasses import dataclass, asdict
from pathlib import Path

import pandas as pd
from pyteomics import fasta, parser


@dataclass
class Config:
	dataset: str
	input_csv: str
	proteome: str = "/home/callum/human_proteome/proteome.gz"
	out_dir: str = "/home/callum/sno_data/processed/cleavage"
	enzyme: str = "trypsin"
	min_len: int = 7           # peptides shorter than this are assumed undetectable
	max_len: int = 50          # peptides longer than this are assumed undetectable
	missed_cleavages: int = 0
	residue: str = "C"


def load_sites(cfg):
	df = pd.read_csv(cfg.input_csv)
	if "organism" in df.columns:
		df = df[df["organism"] == "Homo sapiens"]
	df = df.dropna(subset=["id", "position"]).copy()
	df["id"] = df["id"].astype(str).str.strip()
	df["position"] = df["position"].astype(int)
	return df.drop_duplicates(subset=["id", "position"])


def load_sequences(cfg, wanted):
	seqs = {}
	with gzip.open(cfg.proteome, mode="rt") as f:
		for description, sequence in fasta.FASTA(f):
			acc = description.split("|")[1]
			if acc in wanted:
				seqs[acc] = sequence
	return seqs


def lost_positions(seq, cfg):
	"""1-based positions of residues not covered by any peptide in the length window."""
	all_pos = [i + 1 for i, aa in enumerate(seq) if aa == cfg.residue]
	covered = set()
	for start, pep in parser.icleave(seq, cfg.enzyme, missed_cleavages=cfg.missed_cleavages):
		if cfg.min_len <= len(pep) <= cfg.max_len:
			for i, aa in enumerate(pep):
				if aa == cfg.residue:
					covered.add(start + i + 1)   # start is 0-based
	return [p for p in all_pos if p not in covered], len(all_pos)


def run(cfg):
	out = Path(cfg.out_dir) / cfg.dataset
	out.mkdir(parents=True, exist_ok=True)

	df = load_sites(cfg)
	seqs = load_sequences(cfg, set(df["id"]))
	missing = sorted(set(df["id"]) - set(seqs))

	# write the subset FASTA for the record
	fasta.write(list(seqs.items()), str(out / "subset.fasta"))

	# per protein: which residues are lost
	lost = {}
	total_res = 0
	for acc, seq in seqs.items():
		lost[acc], n = lost_positions(seq, cfg)
		total_res += n
	n_lost = sum(len(v) for v in lost.values())

	# per site checks
	df = df[df["id"].isin(seqs)].copy()
	df["is_residue"] = [seqs[i][p - 1] == cfg.residue if 0 < p <= len(seqs[i]) else False
						for i, p in zip(df["id"], df["position"])]
	df["in_lost_peptide"] = [p in lost[i] for i, p in zip(df["id"], df["position"])]

	bad = df[~df["is_residue"]]
	conflicts = df[df["in_lost_peptide"]]

	df.to_csv(out / "sites_checked.csv", index=False)
	conflicts.to_csv(out / "conflicts.csv", index=False)
	bad.to_csv(out / "bad_positions.csv", index=False)
	
	mod_ok = df[df["is_residue"]]
	n_mod = len(mod_ok)
	n_mod_lost = int(mod_ok["in_lost_peptide"].sum())
	n_neg = total_res - n_mod
	n_neg_lost = n_lost - n_mod_lost

	summary = {
		"proteins_requested": len(set(df["id"])) + len(missing),
		"proteins_found": len(seqs),
		"ids_not_in_proteome": len(missing),
		"residues_total": total_res,
		"residues_lost": n_lost,
		"fraction_lost": round(n_lost / total_res, 4) if total_res else None,
		"modified_sites": len(df),
		"sites_not_the_stated_residue": len(bad),
		"modified_sites_in_lost_peptides": len(conflicts),
		"config": asdict(cfg),
		"modified_sites_valid": n_mod,
		"negatives_total": n_neg,
		"negatives_lost": n_neg_lost,
		"fraction_negatives_lost": round(n_neg_lost / n_neg, 4) if n_neg else None,
		
	}
	with open(out / "summary.json", "w") as f:
		json.dump(summary, f, indent=2)

	for k, v in summary.items():
		if k != "config":
			print(f"{k}: {v}")
	print(f"Results written to {out}")


if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("--dataset", required=True)
	ap.add_argument("--input", required=True)
	ap.add_argument("--proteome", default=Config.proteome)
	ap.add_argument("--out_dir", default=Config.out_dir)
	ap.add_argument("--enzyme", default=Config.enzyme)
	ap.add_argument("--min_len", type=int, default=Config.min_len)
	ap.add_argument("--max_len", type=int, default=Config.max_len)
	ap.add_argument("--missed", type=int, default=Config.missed_cleavages)
	ap.add_argument("--residue", default=Config.residue)
	a = ap.parse_args()

	run(Config(dataset=a.dataset, input_csv=a.input, proteome=a.proteome, out_dir=a.out_dir,
			enzyme=a.enzyme, min_len=a.min_len, max_len=a.max_len,
			missed_cleavages=a.missed, residue=a.residue))