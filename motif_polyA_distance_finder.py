#!/usr/bin/env python3

"""
Compute distances between a poly(A) site and the first occurrence of a motif
within fixed-length sequence windows, and visualize the distribution as a histogram.

This script assumes input sequences are polyA-centered windows of length 101 bp,
with the poly(A) site fixed at index 50. For each sequence, it searches for the
first occurrence of a user-specified motif using fuzzy matching (substitutions
only, no indels). The distance is calculated as:

  start_of_motif - 50

Negative values indicate upstream motifs; positive values indicate downstream
motifs.

A histogram of motif distances is saved as a PNG.

Workflow:
1) Read sequences from FASTA (ignoring headers).
2) Search each sequence for the first fuzzy motif match.
3) Compute distances relative to the polyA center.
4) Plot and save a histogram of distances.

Inputs:
  -f / --fasta   FASTA file of polyA-centered windows (101 bp assumed)
  --motif        Motif to search (default: TGTTTGTT)
  --subs         Maximum substitutions allowed (default: 0)

Output:
  polyA_motif_hist.png
    Histogram of motif distances relative to the polyA site.

Dependencies:
  regex, matplotlib

Author: Tania Iranpour
"""

import argparse
import regex
import matplotlib.pyplot as plt

polyA_index = 50  # fixed center index for 101-bp windows

def read_fasta(path):
    """Return a list of sequences (skip '>' lines)."""
    seqs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue  # skip headers and empty lines
            seqs.append(line.upper())
    return seqs


def main():
    ap = argparse.ArgumentParser(
        description="Measure distance from polyA (index 50) to start of a fuzzy motif."
    )
    ap.add_argument("-f", "--fasta", required=True, help="Input FASTA file.")
    ap.add_argument("--motif", default="TGTTTGTT", help="Motif to search (default: TGTTTGTT).")
    ap.add_argument("--subs", type=int, default=0, help="Max substitutions allowed (default: 0).")
    args = ap.parse_args()

    # fuzzy regex: allow up to args.subs substitutions (no indels)
    pat = regex.compile(f"({args.motif})" + f"{{s<={args.subs}}}", flags=regex.IGNORECASE)

    distances = []
    for seq in read_fasta(args.fasta):
        m = pat.search(seq)  # FIRST match only
        if not m:
            continue
        start = m.start()
        dist = start - polyA_index
        distances.append(dist)

    if not distances:
        print("No motif hits found; no histogram will be created.")
        return

    # histogram
    png_path = "polyA_motif_hist.png"
    plt.figure(figsize=(8, 4.5))
    plt.hist(distances, bins=41, edgecolor="black")
    plt.xlabel("Distance (bp) from polyA (start_of_motif - 50)")
    plt.ylabel("Number of sequences")
    plt.title(f"Motif {args.motif} — distances to polyA@50")
    plt.tight_layout()
    plt.savefig(png_path, dpi=200)
    plt.close()
    print(f"Done. Matches: {len(distances)}")
    print(f"Wrote: {png_path}")

if __name__ == "__main__":
    main()
