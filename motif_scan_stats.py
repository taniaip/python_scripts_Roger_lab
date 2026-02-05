#!/usr/bin/env python3

"""
motif_scan_stats.py

Scan nucleotide FASTA sequences for a user-defined motif allowing a specified
number of substitutions, and report summary statistics on motif occurrence.

This script searches each sequence for overlapping motif matches using fuzzy
regex (substitutions only), then computes:

- Total number of sequences
- Number and percentage of sequences containing ≥1 hit
- Total motif hits across all sequences
- Mean and standard deviation of motif start and end positions (1-based)

Inputs:
  --fasta      Input FASTA file
  --motif      Motif sequence (e.g., TGTTTGTT)
  --pmismatch  Maximum number of allowed substitutions

Outputs:
  Printed summary statistics to stdout.

Intended for exploratory motif prevalence and positional analysis in genomic
or transcriptomic datasets.

Dependencies:
  Biopython, regex

Author: Tania Iranpour
"""

import argparse
from Bio import SeqIO
import regex
import statistics


def parse_args():
    p = argparse.ArgumentParser(
        description="Scan FASTA for a motif allowing up to N substitutions; report averages, std, and prevalence."
    )
    p.add_argument("-f", "--fasta", dest="fasta_file", required=True, help="Input FASTA file.")
    p.add_argument("--motif", dest="motif_seq", required=True, help="Motif to search for (e.g., TGTTTGTT).")
    p.add_argument("--pmismatch", dest="pmismatch", type=int, required=True, help="Max substitutions allowed.")
    return p.parse_args()


def get_motif_info(fasta_file: str, motif_seq: str, pmismatch: int):
    motif_seq = motif_seq.upper()
    motif_len = len(motif_seq)

    # substitutions only, allow overlaps
    pattern = regex.compile(fr"(?b)(?=({motif_seq}){{s<={pmismatch}}})")

    start_positions = []
    end_positions = []
    total_seqs = 0
    total_hits = 0
    seqs_with_hit = 0

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        total_seqs += 1

        hits_this_seq = 0
        for m in pattern.finditer(seq):
            start_1 = m.start() + 1
            end_1 = start_1 + motif_len - 1

            start_positions.append(start_1)
            end_positions.append(end_1)

            total_hits += 1
            hits_this_seq += 1

        if hits_this_seq > 0:
            seqs_with_hit += 1

    if total_hits == 0:
        return {
            "total_sequences": total_seqs,
            "sequences_with_hit": 0,
            "hitcount": 0,
            "avg_start": None,
            "avg_end": None,
            "std_start": None,
            "std_end": None,
            "percent_sequences_with_hit": 0.0,
        }

    avg_start = statistics.mean(start_positions)
    avg_end = statistics.mean(end_positions)

    # sample std requires at least 2 values
    std_start = statistics.stdev(start_positions) if total_hits > 1 else 0.0
    std_end = statistics.stdev(end_positions) if total_hits > 1 else 0.0

    return {
        "total_sequences": total_seqs,
        "sequences_with_hit": seqs_with_hit,
        "hitcount": total_hits,
        "avg_start": avg_start,
        "avg_end": avg_end,
        "std_start": std_start,
        "std_end": std_end,
        "percent_sequences_with_hit": (seqs_with_hit / total_seqs) * 100.0 if total_seqs else 0.0,
    }


def main():
    args = parse_args()
    stats = get_motif_info(args.fasta_file, args.motif_seq, args.pmismatch)

    print(f"Total_sequences: {stats['total_sequences']}")
    print(f"Sequences with hit (>=1): {stats['sequences_with_hit']}")
    print(f"Hitcount: {stats['hitcount']}")
    print(f"Percent_sequences_with_hit: {stats['percent_sequences_with_hit']:.2f}%")

    if stats["hitcount"] == 0:
        print("No motif hits found.")
        return

    print(f"Average_start (1-based): {stats['avg_start']:.3f} ± {stats['std_start']:.3f}")
    print(f"Average_end (1-based): {stats['avg_end']:.3f} ± {stats['std_end']:.3f}")


if __name__ == "__main__":
    main()


# $ python3 ../motif_pos_ave_finder.py -f filtered_ST7B_all_sequence.fa --motif TGTTTGTT --pmismatch 1
# Total_sequences: 5103
# Sequences with hit (>=1): 4386
# Hitcount: 4814
# Percent_sequences_with_hit: 85.95%
# Average_start (1-based): 54.753 ± 8.881
# Average_end (1-based): 61.753 ± 8.881


# $ python3 ../motif_pos_ave_finder.py -f filtered_ST2_ALL_sequence.fa --motif TGTTTGTT --pmismatch 1
# Total_sequences: 5596
# Sequences with hit (>=1): 4647
# Hitcount: 4733
# Percent_sequences_with_hit: 83.04%
# Average_start (1-based): 54.776 ± 7.561
# Average_end (1-based): 61.776 ± 7.561


# $ python3 ../motif_pos_ave_finder.py -f blastoise_filtered_sequence.fa --motif TTGTTTTGTT --pmismatch 1 
# Total_sequences: 2301
# Sequences with hit (>=1): 984
# Hitcount: 1000
# Percent_sequences_with_hit: 42.76%
# Average_start (1-based): 55.323 ± 11.173
# Average_end (1-based): 64.323 ± 11.173



# $ python3 ../motif_pos_ave_finder.py -f Rus-s_filtered_sequence.fa --motif TTGTTTTGTT --pmismatch 1 
# Total_sequences: 4170
# Sequences with hit (>=1): 1870
# Hitcount: 1919
# Percent_sequences_with_hit: 44.84%
# Average_start (1-based): 54.999 ± 14.516
# Average_end (1-based): 63.999 ± 14.516


# $ python3 ../motif_pos_ave_finder.py -f Hermann_filtered_sequence.fa --motif TTGTTTTGTT --pmismatch 1 
# Total_sequences: 1353
# Sequences with hit (>=1): 753
# Hitcount: 782
# Percent_sequences_with_hit: 55.65%
# Average_start (1-based): 56.455 ± 12.426
# Average_end (1-based): 65.455 ± 12.426
