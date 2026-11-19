#!/usr/bin/env python3
"""Run CD-HIT on protein sequences and return clustered"""
import argparse
import subprocess
from pathlib import Path


def run_cdhit(fasta: Path, output: Path, identity: float = 0.4,
              word_size: int = 2, desc_len: int = 0) -> None:
    if not fasta.is_file():
        raise FileNotFoundError(fasta)
    output.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "cd-hit",
        "-i", str(fasta),
        "-o", str(output),
        "-c", str(identity),
        "-n", str(word_size),
        "-d", str(desc_len),
    ]
    subprocess.run(cmd, check=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-i", "--input", type=Path, required=True, help="input FASTA")
    p.add_argument("-o", "--output", type=Path, required=True, help="output path")
    p.add_argument("-c", "--identity", type=float, default=0.4, help="sequence identity threshold")
    p.add_argument("-n", "--word_size", type=int, default=2)
    p.add_argument("-d", "--desc_len", type=int, default=0)
    args = p.parse_args()

    run_cdhit(args.input, args.output, args.identity, args.word_size, args.desc_len)


if __name__ == "__main__":
    main()