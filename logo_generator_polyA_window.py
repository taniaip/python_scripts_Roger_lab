#!/usr/bin/env python

"""
Generate a GC/background-normalized sequence logo from poly(A)-centered genomic
windows.

This script:
1) Computes genome-wide nucleotide background frequencies (A/C/G/T) from an input
   genomic FASTA (ignoring non-ACGT characters).
2) Extracts fixed-length windows centered on poly(A) site coordinates from a GFF3
   file (default flank=50, producing windows of length 2*flank+1). Minus-strand
   sites are reverse-complemented so all windows are in 5′→3′ orientation.
3) Builds a position-wise probability matrix from the extracted windows and
   converts it to an information-content matrix (bits) using the genome-wide
   background frequencies.
4) Saves the resulting sequence logo as a PNG and also writes the extracted
   windows to sequences.fasta for downstream use.

Inputs:
  --gff3   GFF3 file containing poly(A) site coordinates
  --fasta  Genomic FASTA file
  -o       Output image basename (PNG extension added)
  -t       Title text for the logo figure

Outputs:
  <out>.png          Sequence logo (information content, bits)
  sequences.fasta    Extracted poly(A)-centered windows

Dependencies:
  Biopython, pandas, logomaker, matplotlib

>>>>>> The "U" version is the first script and the 2nd one is the "T" version of teh script.

Author: Tania Iranpour
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import logomaker
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description="Generate sequence logos (PNG).")
    parser.add_argument("-g", "--gff3", dest="gff3_file", type=str, required=True,
                        help="Input GFF3 file with polyA site info.")
    parser.add_argument("-f", "--fasta", dest="fasta_file", type=str, required=True,
                        help="Input FASTA file of genomic sequences.")
    parser.add_argument("-o", dest="out", required=True,
                        help="Logo picture result (basename, no extension).")
    parser.add_argument("-t", dest="title", required=True,
                        help="Title for the output figure.")
    parser.add_argument("--flank", type=int, default=25,
                        help="Number of bases upstream/downstream (default: 50).")
    parser.add_argument("--rna", action="store_true",
                        help="Convert extracted windows to RNA (T->U) before logo.")
    return parser.parse_args()

def extract_polyA_sequences(gff3_file, fasta_file, flank=50, to_rna=False):
    sequences = []
    fasta_seqs = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

    with open(gff3_file, "r") as gf:
        for line in gf:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            seqid = parts[0]
            try:
                start0 = int(parts[3]) - 1  # 0-based
            except ValueError:
                continue

            strand = parts[6]
            if seqid not in fasta_seqs:
                continue

            full_seq = fasta_seqs[seqid]

            left = start0 - flank
            right = start0 + flank + 1
            if left < 0 or right > len(full_seq):
                continue

            region_seq = full_seq[left:right].upper()
            if len(region_seq) != 2 * flank + 1:
                continue

            # Reverse-complement on DNA first
            if strand == "-":
                region_seq = str(Seq(region_seq).reverse_complement())

            # Optionally convert DNA -> RNA for downstream logo work
            if to_rna:
                region_seq = region_seq.replace("T", "U")

            sequences.append(region_seq)

    # Write sequences to FASTA (match what you're actually using: DNA or RNA)
    with open("sequences.fasta", "w") as out:
        for i, seq in enumerate(sequences, start=1):
            out.write(f">sequence_{i}\n{seq}\n")

    return sequences

def compute_genome_background(fasta_file):
    """
    Compute genome-wide nucleotide background frequencies (A, C, G, T),
    treating U as T if it appears (so it still works on RNA-like inputs).
    Ignores non-ACGTU characters.
    """
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        for base in seq:
            if base == "U":
                base = "T"  # treat U as T
            if base in counts:
                counts[base] += 1

    total = sum(counts.values())
    if total == 0:
        raise ValueError("No valid A/C/G/T/U bases found in FASTA for background computation.")

    background = {b: counts[b] / total for b in counts}
    return background

def save_logo_logomaker(sequences, output_file, background, title):
    if not sequences:
        print(f"No sequences found for {output_file}, skipping logo.")
        return

    # If sequences contain U, build a U-alphabet logo; otherwise use T-alphabet.
    has_u = any("U" in s for s in sequences)

    if has_u:
        valid_bases = ["A", "C", "G", "U"]

        # Counts directly from RNA sequences (with U)
        counts_df = logomaker.alignment_to_matrix(sequences, to_type="counts")

        # Ensure all A/C/G/U columns exist
        for b in valid_bases:
            if b not in counts_df.columns:
                counts_df[b] = 0
        counts_df = counts_df[valid_bases]

        # Convert counts -> probabilities
        prob_df = counts_df.div(counts_df.sum(axis=1), axis=0)

        # Background: use genome's T frequency as U frequency
        bg_u = {
            "A": background.get("A", 0.25),
            "C": background.get("C", 0.25),
            "G": background.get("G", 0.25),
            "U": background.get("T", 0.25),
        }
        bg_list = [bg_u[b] for b in valid_bases]

    else:
        valid_bases = ["A", "C", "G", "T"]

        counts_df = logomaker.alignment_to_matrix(sequences, to_type="counts")
        for b in valid_bases:
            if b not in counts_df.columns:
                counts_df[b] = 0
        counts_df = counts_df[valid_bases]

        prob_df = counts_df.div(counts_df.sum(axis=1), axis=0)
        bg_list = [background[b] for b in valid_bases]

    info_df = logomaker.transform_matrix(
        prob_df,
        from_type="probability",
        to_type="information",
        background=bg_list
    )

    plt.figure(figsize=(len(info_df) / 3.0, 3))
    logo = logomaker.Logo(info_df, color_scheme="classic")
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left", "bottom"], visible=True)
    logo.ax.set_ylabel("bits")
    logo.ax.set_ylim([0, 2.0])
    plt.title(f"Sequence Logo for {title}")
    plt.tight_layout()
    plt.savefig(output_file + ".png", dpi=200)
    plt.close()

def main():
    args = parse_args()

    background = compute_genome_background(args.fasta_file)

    seqs = extract_polyA_sequences(
        gff3_file=args.gff3_file,
        fasta_file=args.fasta_file,
        flank=args.flank,
        to_rna=args.rna
    )

    save_logo_logomaker(seqs, args.out, background, args.title)

if __name__ == "__main__":
    main()


# import argparse
# from Bio import SeqIO
# from Bio.Seq import Seq
# import pandas as pd
# import logomaker
# import matplotlib.pyplot as plt

# def parse_args():
#     parser = argparse.ArgumentParser(description="Generate sequence logos (PNG) and write 'else' GFF3 lines.")
#     parser.add_argument("-g", "--gff3", dest="gff3_file", type=str, required=True,
#                         help="Input GFF3 file with polyA site info.")
#     parser.add_argument("-f", "--fasta", dest="fasta_file", type=str, required=True,
#                         help="Input FASTA file of genomic sequences.")
#     parser.add_argument("-o", dest="out", required=True,
#                         help="Logo picture result (basename, no extension).")
#     parser.add_argument("-t", dest="title", required=True,
#                         help="title for the output figure.")
#     return parser.parse_args()

# def extract_polyA_sequences(gff3_file, fasta_file, flank50):
#     sequences = []

#     fasta_seqs = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

#     with open(gff3_file, "r") as gf:
#         for line in gf:
#             line = line.strip()
#             if not line or line.startswith("#"):
#                 continue

#             parts = line.split("\t")
#             if len(parts) < 7:
#                 continue

#             seqid = parts[0]
#             try:
#                 start = int(parts[3]) - 1
#             except ValueError:
#                 continue

#             strand = parts[6]

#             if seqid not in fasta_seqs:
#                 continue

#             full_seq = fasta_seqs[seqid]

#             region_seq = full_seq[start - flank: start + flank + 1].upper()
#             # Expect 2*flank + 1 bp window
#             if len(region_seq) != 2 * flank + 1:
#                 continue

#             if strand == "-":
#                 region_seq = str(Seq(region_seq).reverse_complement())

#             if region_seq:
#                 sequences.append(region_seq)

#     # Write sequences to FASTA
#     with open("sequences.fasta", "w") as out:
#         for i, seq in enumerate(sequences, start=1):
#             out.write(f">sequence_{i}\n")
#             out.write(f"{seq}\n")

#     return sequences

# def compute_genome_background(fasta_file):
#     """
#     Compute genome-wide nucleotide background frequencies (A, C, G, T)
#     from the input FASTA file. Ignores any non-ACGT characters.
#     """
#     counts = {"A": 0, "C": 0, "G": 0, "T": 0}

#     for record in SeqIO.parse(fasta_file, "fasta"):
#         seq = str(record.seq).upper()
#         for base in seq:
#             if base in counts:
#                 counts[base] += 1

#     total = sum(counts.values())

#     background = {b: counts[b] / total for b in counts}
#     gc = (counts["G"] + counts["C"]) / float(total)

#     return background

# def save_logo_logomaker(sequences, output_file, background, title):
#     if not sequences:
#         print(f"No sequences found for {output_file}, skipping logo.")
#         return

#     # First get a counts matrix from the alignment
#     counts_df = logomaker.alignment_to_matrix(sequences, to_type="counts")

#     # Convert counts to position-wise probabilities
#     prob_df = counts_df.div(counts_df.sum(axis=1), axis=0)

#     # 3) Keep only A, C, G, T columns (drop N or anything else)
#     valid_bases = ["A", "C", "G", "T"]
#     prob_df = prob_df[[b for b in prob_df.columns if b in valid_bases]]

#     # 4) Background list in EXACT same order as prob_df columns
#     bg_list = [background[b] for b in prob_df.columns]

#     info_df = logomaker.transform_matrix(
#         prob_df,
#         from_type="probability",
#         to_type="information",
#         background=bg_list
#     )

#     # Set up the logo figure
#     plt.figure(figsize=(len(info_df) / 3.0, 3))
#     logo = logomaker.Logo(info_df, color_scheme="classic")
#     logo.style_spines(visible=False)
#     logo.style_spines(spines=["left", "bottom"], visible=True)
#     logo.ax.set_ylabel("bits")
#     # Max information with 4 bases is 2 bits, so 0-2 is fine
#     logo.ax.set_ylim([0, 2.0])
#     plt.title(f"Sequence Logo for {title}")
#     plt.tight_layout()

#     # Save as PNG
#     plt.savefig(output_file + ".png", dpi=200)
#     plt.close()

# def main():
#     args = parse_args()

#     # 1) Compute genome-wide background from the FASTA
#     background = compute_genome_background(args.fasta_file)

#     # 2) Extract polyA-centered windows
#     seqs = extract_polyA_sequences(
#         gff3_file=args.gff3_file,
#         fasta_file=args.fasta_file,
#         flank=50
#     )

#     # 3) Make GC-normalized logo
#     save_logo_logomaker(seqs, args.out, background, args.title)

# if __name__ == "__main__":
#     main()




### Fully functionla but doesnt account for background frequencies:

# #!/usr/bin/env python

# import argparse
# from Bio import SeqIO
# from Bio.Seq import Seq
# import pandas as pd
# import logomaker
# import matplotlib.pyplot as plt
# import io

# def parse_args():
#     parser = argparse.ArgumentParser(description="Generate sequence logos (PNG) and write 'else' GFF3 lines.")
#     parser.add_argument("-g", "--gff3", dest="gff3_file", type=str, required=True,
#                         help="Input GFF3 file with polyA site info.")
#     parser.add_argument("-f", "--fasta", dest="fasta_file", type=str, required=True,
#                         help="Input FASTA file of genomic sequences.")
#     parser.add_argument("-o",  dest="out",  required=True, help="Logo picture result")
#     return parser.parse_args()

# def extract_polyA_sequences(gff3_file, fasta_file, flank=50):
#     sequences  = []

#     fasta_seqs = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

#     with open(gff3_file, "r") as gf:
#         for line in gf:
#             line = line.strip()
#             parts = line.split("\t")

#             seqid = parts[0]
#             start = int(parts[3]) - 1
#             strand = parts[6]
#             if seqid not in fasta_seqs:
#                 continue

#             full_seq = fasta_seqs[seqid]

#             region_seq = full_seq[start - flank: start + flank + 1].upper()
#             if len(region_seq) != 101:
#                 continue
#             if strand == "-":
#                 region_seq = str(Seq(region_seq).reverse_complement())
#             if region_seq:
#                 # region_seq = f"{strand}{region_seq}"
#                 sequences.append(region_seq)    
    # # Write sequences to FASTA
    # with open("sequences.fasta", "w") as out:
    #     for i, seq in enumerate(sequences, start=1):
    #         out.write(f">sequence_{i}\n")
    #         out.write(f"{seq}\n")
            
    # return sequences

# def save_logo_logomaker(sequences, output_file):
#     if not sequences:
#         print(f"No sequences found for {output_file}, skipping logo.")
#         return
#     # Convert the alignment to an information content matrix (bits)
#     info_df = logomaker.alignment_to_matrix(sequences, to_type='information')

#     # Set up the logo figure
#     plt.figure(figsize=(len(info_df)/3, 3))
#     logo = logomaker.Logo(info_df, color_scheme='classic')
#     logo.style_spines(visible=False)
#     logo.style_spines(spines=['left', 'bottom'], visible=True)
#     logo.ax.set_ylabel('bits')
#     logo.ax.set_ylim([0, 2.0])  
#     plt.title(f"Sequence Logo: {output_file}")
#     plt.tight_layout()
#     # Save as PNG
#     plt.savefig(output_file + ".png", dpi=200)
#     plt.close()

# def main():
#     args = parse_args()
#     seq = extract_polyA_sequences(
#         gff3_file =args.gff3_file,
#         fasta_file=args.fasta_file,
#         flank=50
#     )

#     save_logo_logomaker(seq, args.out)

# if __name__ == "__main__":
#     main()
