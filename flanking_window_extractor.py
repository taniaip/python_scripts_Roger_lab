#!/usr/bin/env python3

"""
Extract fixed-length sequence windows centered on poly(A) site coordinates from
a genomic FASTA file and write them to a FASTA output.

This script reads poly(A) site positions from a GFF3 file, retrieves upstream and
downstream flanking sequence from the corresponding reference genome, and outputs
each window as an individual FASTA entry. Strand is respected: sites on the minus
strand are reverse-complemented so all sequences are reported in 5′→3′ orientation.

Each output window has length (2 * flank + 1), centered on the poly(A) coordinate.

Inputs:
  --gff3    GFF3 file containing poly(A) site coordinates
  --fasta   Genomic FASTA file
  --flank   Number of bases upstream and downstream (default: 50)

Output:
  sequences.fasta
    FASTA file containing extracted poly(A)-centered windows.

Intended for motif discovery and positional analysis of polyadenylation signals.

Dependencies:
  Biopython

Author: Tania Iranpour
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract polyA-centered sequence windows and write to FASTA."
    )
    parser.add_argument(
        "-g", "--gff3", dest="gff3_file", required=True,
        help="Input GFF3 file with polyA site coordinates."
    )
    parser.add_argument(
        "-f", "--fasta", dest="fasta_file", required=True,
        help="Input genomic FASTA file."
    )
    parser.add_argument(
        "--flank", type=int, default=50,
        help="Number of bases upstream and downstream of polyA site (default: 50)."
    )
    return parser.parse_args()


def extract_polyA_sequences(gff3_file, fasta_file, flank):
    sequences = []

    fasta_seqs = {
        record.id: str(record.seq)
        for record in SeqIO.parse(fasta_file, "fasta")
    }

    with open(gff3_file) as gf:
        for line in gf:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 7:
                continue

            seqid = parts[0]
            strand = parts[6]

            try:
                start = int(parts[3]) - 1  # convert to 0-based
            except ValueError:
                continue

            if seqid not in fasta_seqs:
                continue

            full_seq = fasta_seqs[seqid]

            region = full_seq[start - flank : start + flank + 1].upper()
            if len(region) != 2 * flank + 1:
                continue

            if strand == "-":
                region = str(Seq(region).reverse_complement())

            sequences.append(region)

    return sequences


def write_fasta(sequences, outfile="sequences.fasta"):
    with open(outfile, "w") as out:
        for i, seq in enumerate(sequences, start=1):
            out.write(f">sequence_{i}\n")
            out.write(f"{seq}\n")


def main():
    args = parse_args()
    sequences = extract_polyA_sequences(
        gff3_file=args.gff3_file,
        fasta_file=args.fasta_file,
        flank=args.flank
    )
    write_fasta(sequences)
    print(f"Wrote {len(sequences)} sequences to sequences.fasta")


if __name__ == "__main__":
    main()

